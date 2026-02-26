"""
Tests for Session-Based Data Isolation

Tests that bookmarks and audit entries are scoped to sessions,
invisible across sessions, and properly purged.
Uses SQLite in-memory for test isolation.
"""

import uuid
from datetime import datetime, timedelta

import pytest
from sqlalchemy import delete, func, select
from sqlalchemy.ext.asyncio import async_sessionmaker, create_async_engine
from sqlalchemy.pool import StaticPool

from app.core.session import SESSION_MAX_AGE, create_session_id
from app.db import Base
from app.db.models.audit import ValidationAuditEntry
from app.db.models.bookmark import Bookmark


@pytest.fixture
async def db_session():
    """Create an in-memory SQLite async session for testing."""
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


class TestBookmarkSessionScoping:
    """Bookmarks are isolated between sessions."""

    @pytest.mark.asyncio
    async def test_session_a_invisible_to_session_b(self, db_session):
        """Session A's bookmarks are not visible to session B."""
        # Session A creates bookmarks
        bm1 = Bookmark(smiles="CCO", name="Ethanol", session_id="session-aaa")
        bm2 = Bookmark(smiles="CC(=O)O", name="Acetic acid", session_id="session-aaa")
        # Session B creates a bookmark
        bm3 = Bookmark(smiles="c1ccccc1", name="Benzene", session_id="session-bbb")
        db_session.add_all([bm1, bm2, bm3])
        await db_session.commit()

        # Query as session A
        result = await db_session.execute(
            select(Bookmark).where(Bookmark.session_id == "session-aaa")
        )
        session_a_bookmarks = result.scalars().all()
        assert len(session_a_bookmarks) == 2
        assert {bm.name for bm in session_a_bookmarks} == {"Ethanol", "Acetic acid"}

        # Query as session B
        result = await db_session.execute(
            select(Bookmark).where(Bookmark.session_id == "session-bbb")
        )
        session_b_bookmarks = result.scalars().all()
        assert len(session_b_bookmarks) == 1
        assert session_b_bookmarks[0].name == "Benzene"

    @pytest.mark.asyncio
    async def test_api_key_scoped_bookmark(self, db_session):
        """API key-scoped bookmarks use api_key_hash, not session_id."""
        bm = Bookmark(
            smiles="CCO",
            name="API Ethanol",
            api_key_hash="hash-abc123",
            session_id=None,
        )
        db_session.add(bm)
        await db_session.commit()

        # Visible by api_key_hash
        result = await db_session.execute(
            select(Bookmark).where(Bookmark.api_key_hash == "hash-abc123")
        )
        assert len(result.scalars().all()) == 1

        # Not visible by session
        result = await db_session.execute(
            select(Bookmark).where(Bookmark.session_id == "session-aaa")
        )
        assert len(result.scalars().all()) == 0

    @pytest.mark.asyncio
    async def test_no_session_no_bookmarks(self, db_session):
        """Bookmarks with session_id are not visible without a session filter."""
        bm = Bookmark(smiles="CCO", name="Ethanol", session_id="session-aaa")
        db_session.add(bm)
        await db_session.commit()

        # Query without session or api_key filter returns none
        result = await db_session.execute(
            select(Bookmark).where(
                Bookmark.session_id.is_(None),
                Bookmark.api_key_hash.is_(None),
            )
        )
        assert len(result.scalars().all()) == 0


class TestAuditSessionScoping:
    """Audit entries are isolated between sessions."""

    @pytest.mark.asyncio
    async def test_audit_entries_scoped_to_session(self, db_session):
        """Audit entries for session A are invisible to session B."""
        entry_a = ValidationAuditEntry(
            smiles="CCO", outcome="pass", score=95, source="single",
            session_id="session-aaa",
        )
        entry_b = ValidationAuditEntry(
            smiles="c1ccccc1", outcome="warn", score=60, source="single",
            session_id="session-bbb",
        )
        db_session.add_all([entry_a, entry_b])
        await db_session.commit()

        # Session A sees only its entry
        result = await db_session.execute(
            select(ValidationAuditEntry).where(
                ValidationAuditEntry.session_id == "session-aaa"
            )
        )
        entries = result.scalars().all()
        assert len(entries) == 1
        assert entries[0].smiles == "CCO"

    @pytest.mark.asyncio
    async def test_no_session_returns_empty(self, db_session):
        """Without session/API key, no audit entries are visible."""
        entry = ValidationAuditEntry(
            smiles="CCO", outcome="pass", score=95, source="single",
            session_id="session-aaa",
        )
        db_session.add(entry)
        await db_session.commit()

        result = await db_session.execute(
            select(ValidationAuditEntry).where(
                ValidationAuditEntry.session_id.is_(None),
                ValidationAuditEntry.api_key_hash.is_(None),
            )
        )
        assert len(result.scalars().all()) == 0

    @pytest.mark.asyncio
    async def test_api_key_scoped_audit(self, db_session):
        """API key audit entries are scoped by api_key_hash."""
        entry = ValidationAuditEntry(
            smiles="CCO", outcome="pass", score=90, source="single",
            api_key_hash="hash-key-1", session_id=None,
        )
        db_session.add(entry)
        await db_session.commit()

        result = await db_session.execute(
            select(ValidationAuditEntry).where(
                ValidationAuditEntry.api_key_hash == "hash-key-1"
            )
        )
        assert len(result.scalars().all()) == 1

        result = await db_session.execute(
            select(ValidationAuditEntry).where(
                ValidationAuditEntry.api_key_hash == "hash-key-other"
            )
        )
        assert len(result.scalars().all()) == 0


class TestPurgeLogic:
    """Tests for the purge data operation."""

    @pytest.mark.asyncio
    async def test_purge_deletes_session_bookmarks(self, db_session):
        """Purge removes all bookmarks for a session."""
        # Create bookmarks for two sessions
        db_session.add_all([
            Bookmark(smiles="CCO", name="A1", session_id="session-purge"),
            Bookmark(smiles="CC", name="A2", session_id="session-purge"),
            Bookmark(smiles="CCC", name="B1", session_id="session-keep"),
        ])
        await db_session.commit()

        # Purge session-purge
        await db_session.execute(
            delete(Bookmark).where(Bookmark.session_id == "session-purge")
        )
        await db_session.commit()

        # Verify session-purge is gone
        result = await db_session.execute(
            select(func.count()).select_from(Bookmark).where(
                Bookmark.session_id == "session-purge"
            )
        )
        assert result.scalar() == 0

        # Verify session-keep is untouched
        result = await db_session.execute(
            select(func.count()).select_from(Bookmark).where(
                Bookmark.session_id == "session-keep"
            )
        )
        assert result.scalar() == 1

    @pytest.mark.asyncio
    async def test_purge_deletes_session_audit_entries(self, db_session):
        """Purge removes all audit entries for a session."""
        db_session.add_all([
            ValidationAuditEntry(
                smiles="CCO", outcome="pass", score=90, source="single",
                session_id="session-purge",
            ),
            ValidationAuditEntry(
                smiles="CC", outcome="pass", score=85, source="single",
                session_id="session-keep",
            ),
        ])
        await db_session.commit()

        await db_session.execute(
            delete(ValidationAuditEntry).where(
                ValidationAuditEntry.session_id == "session-purge"
            )
        )
        await db_session.commit()

        result = await db_session.execute(
            select(func.count()).select_from(ValidationAuditEntry).where(
                ValidationAuditEntry.session_id == "session-purge"
            )
        )
        assert result.scalar() == 0

        result = await db_session.execute(
            select(func.count()).select_from(ValidationAuditEntry).where(
                ValidationAuditEntry.session_id == "session-keep"
            )
        )
        assert result.scalar() == 1

    @pytest.mark.asyncio
    async def test_purge_does_not_affect_api_key_data(self, db_session):
        """Purging a session doesn't touch API-key-scoped data."""
        db_session.add_all([
            Bookmark(smiles="CCO", name="Session", session_id="session-purge"),
            Bookmark(
                smiles="CC", name="API Key",
                api_key_hash="hash-abc", session_id=None,
            ),
        ])
        await db_session.commit()

        await db_session.execute(
            delete(Bookmark).where(Bookmark.session_id == "session-purge")
        )
        await db_session.commit()

        result = await db_session.execute(
            select(func.count()).select_from(Bookmark).where(
                Bookmark.api_key_hash == "hash-abc"
            )
        )
        assert result.scalar() == 1


class TestCleanupTask:
    """Tests for the expired session cleanup logic."""

    @pytest.mark.asyncio
    async def test_cleanup_removes_old_session_data(self, db_session):
        """Rows older than 30 days with session_id are purged."""
        now = datetime.utcnow()
        cutoff = now - timedelta(days=30)

        # Old session-scoped bookmark
        old_bm = Bookmark(smiles="CCO", name="Old", session_id="session-old")
        db_session.add(old_bm)
        await db_session.commit()

        # Manually set created_at to 31 days ago (simulating age)
        old_bm.created_at = now - timedelta(days=31)
        await db_session.commit()

        # Recent session-scoped bookmark
        recent_bm = Bookmark(smiles="CC", name="Recent", session_id="session-new")
        db_session.add(recent_bm)
        await db_session.commit()

        # API-key bookmark (old but should survive)
        api_bm = Bookmark(smiles="CCC", name="API", api_key_hash="hash-api")
        db_session.add(api_bm)
        await db_session.commit()
        api_bm.created_at = now - timedelta(days=31)
        await db_session.commit()

        # Run cleanup logic (same as cleanup.py)
        await db_session.execute(
            delete(Bookmark).where(
                Bookmark.session_id.isnot(None),
                Bookmark.api_key_hash.is_(None),
                Bookmark.created_at < cutoff,
            )
        )
        await db_session.commit()

        # Old session bookmark: deleted
        result = await db_session.execute(
            select(func.count()).select_from(Bookmark).where(
                Bookmark.session_id == "session-old"
            )
        )
        assert result.scalar() == 0

        # Recent session bookmark: kept
        result = await db_session.execute(
            select(func.count()).select_from(Bookmark).where(
                Bookmark.session_id == "session-new"
            )
        )
        assert result.scalar() == 1

        # API-key bookmark: kept (even though old)
        result = await db_session.execute(
            select(func.count()).select_from(Bookmark).where(
                Bookmark.api_key_hash == "hash-api"
            )
        )
        assert result.scalar() == 1

    @pytest.mark.asyncio
    async def test_cleanup_removes_old_audit_entries(self, db_session):
        """Old session-scoped audit entries are purged."""
        now = datetime.utcnow()
        cutoff = now - timedelta(days=30)

        old_entry = ValidationAuditEntry(
            smiles="CCO", outcome="pass", score=90, source="single",
            session_id="session-old",
        )
        db_session.add(old_entry)
        await db_session.commit()
        old_entry.created_at = now - timedelta(days=31)
        await db_session.commit()

        recent_entry = ValidationAuditEntry(
            smiles="CC", outcome="pass", score=85, source="single",
            session_id="session-new",
        )
        db_session.add(recent_entry)
        await db_session.commit()

        await db_session.execute(
            delete(ValidationAuditEntry).where(
                ValidationAuditEntry.session_id.isnot(None),
                ValidationAuditEntry.api_key_hash.is_(None),
                ValidationAuditEntry.created_at < cutoff,
            )
        )
        await db_session.commit()

        result = await db_session.execute(
            select(func.count()).select_from(ValidationAuditEntry)
        )
        assert result.scalar() == 1  # Only recent entry survives


class TestSessionModule:
    """Tests for the session.py module functions."""

    def test_create_session_id_is_uuid(self):
        """create_session_id returns a valid UUID string."""
        sid = create_session_id()
        parsed = uuid.UUID(sid)
        assert str(parsed) == sid

    def test_create_session_id_unique(self):
        """Each call produces a unique session ID."""
        ids = {create_session_id() for _ in range(100)}
        assert len(ids) == 100

    def test_session_max_age_is_30_days(self):
        """SESSION_MAX_AGE is exactly 30 days in seconds."""
        assert SESSION_MAX_AGE == 30 * 24 * 3600
