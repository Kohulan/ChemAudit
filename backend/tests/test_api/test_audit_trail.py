"""
Tests for Validation Audit Trail

Tests append-only audit logging, batch event logging, paginated
history with filtering, and verification that no UPDATE/DELETE exists.
Uses SQLite in-memory for test isolation.
"""

import pytest
from sqlalchemy import select

from app.db.models.audit import ValidationAuditEntry
from app.services.audit.service import log_batch_event, log_validation_event


@pytest.fixture
async def db_session():
    """Create an in-memory SQLite async session for testing."""
    from sqlalchemy.ext.asyncio import async_sessionmaker, create_async_engine
    from sqlalchemy.pool import StaticPool

    from app.db import Base

    engine = create_async_engine(
        "sqlite+aiosqlite:///:memory:",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.create_all)

    session_factory = async_sessionmaker(engine, expire_on_commit=False)
    async with session_factory() as session:
        yield session

    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.drop_all)
    await engine.dispose()


class TestLogValidationEvent:
    """Tests for single validation audit logging."""

    @pytest.mark.asyncio
    async def test_log_single_validation(self, db_session):
        """Single validation event is logged with correct fields."""
        entry = await log_validation_event(
            db_session,
            smiles="CCO",
            inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            outcome="pass",
            score=95,
            source="single",
        )
        assert entry.id is not None
        assert entry.smiles == "CCO"
        assert entry.outcome == "pass"
        assert entry.score == 95
        assert entry.source == "single"

    @pytest.mark.asyncio
    async def test_log_validation_warn(self, db_session):
        """Warn outcome is stored correctly."""
        entry = await log_validation_event(
            db_session,
            smiles="c1ccccc1",
            inchikey=None,
            outcome="warn",
            score=55,
            source="single",
        )
        assert entry.outcome == "warn"
        assert entry.score == 55

    @pytest.mark.asyncio
    async def test_log_validation_fail(self, db_session):
        """Fail outcome is stored correctly."""
        entry = await log_validation_event(
            db_session,
            smiles="invalid",
            inchikey=None,
            outcome="fail",
            score=10,
            source="single",
        )
        assert entry.outcome == "fail"


class TestLogBatchEvent:
    """Tests for batch completion audit logging."""

    @pytest.mark.asyncio
    async def test_log_batch_event_pass(self, db_session):
        """Batch with no failures gets outcome='pass'."""
        entry = await log_batch_event(
            db_session,
            job_id="batch-123",
            molecule_count=100,
            pass_count=100,
            fail_count=0,
        )
        assert entry.outcome == "pass"
        assert entry.source == "batch"
        assert entry.job_id == "batch-123"
        assert entry.molecule_count == 100

    @pytest.mark.asyncio
    async def test_log_batch_event_warn(self, db_session):
        """Batch with <50% failures gets outcome='warn'."""
        entry = await log_batch_event(
            db_session,
            job_id="batch-456",
            molecule_count=100,
            pass_count=70,
            fail_count=30,
        )
        assert entry.outcome == "warn"

    @pytest.mark.asyncio
    async def test_log_batch_event_fail(self, db_session):
        """Batch with >=50% failures gets outcome='fail'."""
        entry = await log_batch_event(
            db_session,
            job_id="batch-789",
            molecule_count=100,
            pass_count=30,
            fail_count=70,
        )
        assert entry.outcome == "fail"


class TestHistoryPagination:
    """Tests for paginated history queries."""

    @pytest.mark.asyncio
    async def test_pagination(self, db_session):
        """Inserting 10 entries and querying page_size=3 returns correct subset."""
        from sqlalchemy import func

        for i in range(10):
            await log_validation_event(
                db_session,
                smiles=f"C{'C' * i}O",
                inchikey=None,
                outcome="pass" if i % 2 == 0 else "fail",
                score=50 + i * 5,
                source="single",
            )

        # Count total
        count_result = await db_session.execute(
            select(func.count(ValidationAuditEntry.id))
        )
        total = count_result.scalar()
        assert total == 10

        # Page 1 of size 3
        result = await db_session.execute(
            select(ValidationAuditEntry)
            .order_by(ValidationAuditEntry.created_at.desc())
            .offset(0)
            .limit(3)
        )
        page = result.scalars().all()
        assert len(page) == 3

    @pytest.mark.asyncio
    async def test_filter_by_outcome(self, db_session):
        """Filtering by outcome returns correct subset."""
        await log_validation_event(
            db_session, smiles="CCO", inchikey=None,
            outcome="pass", score=90, source="single",
        )
        await log_validation_event(
            db_session, smiles="C=O", inchikey=None,
            outcome="fail", score=20, source="single",
        )
        await log_validation_event(
            db_session, smiles="CC", inchikey=None,
            outcome="pass", score=85, source="single",
        )

        result = await db_session.execute(
            select(ValidationAuditEntry).where(
                ValidationAuditEntry.outcome == "fail"
            )
        )
        entries = result.scalars().all()
        assert len(entries) == 1
        assert entries[0].smiles == "C=O"

    @pytest.mark.asyncio
    async def test_smiles_search(self, db_session):
        """SMILES substring search returns matching entries."""
        await log_validation_event(
            db_session, smiles="CCO", inchikey=None,
            outcome="pass", score=90, source="single",
        )
        await log_validation_event(
            db_session, smiles="c1ccccc1", inchikey=None,
            outcome="pass", score=80, source="single",
        )
        await log_validation_event(
            db_session, smiles="CCCO", inchikey=None,
            outcome="pass", score=85, source="single",
        )

        result = await db_session.execute(
            select(ValidationAuditEntry).where(
                ValidationAuditEntry.smiles.contains("CCO")
            )
        )
        entries = result.scalars().all()
        assert len(entries) == 2  # CCO and CCCO both contain "CCO"

    @pytest.mark.asyncio
    async def test_audit_trail_is_append_only(self, db_session):
        """Verify that the history router has no UPDATE or DELETE endpoints."""
        from app.api.routes.history import router

        for route in router.routes:
            methods = getattr(route, "methods", set())
            assert "PUT" not in methods, "History should have no PUT endpoint"
            assert "PATCH" not in methods, "History should have no PATCH endpoint"
            assert "DELETE" not in methods, "History should have no DELETE endpoint"
