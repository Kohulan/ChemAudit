"""
Tests for Permalink API

Tests batch permalink creation, resolution, expiry enforcement,
and single molecule stateless permalinks.
Uses SQLite in-memory for test isolation.
"""

import json
from datetime import datetime, timedelta, timezone

import pytest
from sqlalchemy import select

from app.db.models.permalink import BatchPermalink


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


class TestBatchPermalinks:
    """Tests for batch permalink create/resolve."""

    @pytest.mark.asyncio
    async def test_create_batch_permalink(self, db_session):
        """Creating a permalink stores short_id and job_id."""
        permalink = BatchPermalink(
            short_id="abc12345",
            job_id="test-job-id-123",
            snapshot_data=json.dumps({"results": []}),
        )
        db_session.add(permalink)
        await db_session.commit()
        await db_session.refresh(permalink)

        assert permalink.id is not None
        assert permalink.short_id == "abc12345"
        assert permalink.job_id == "test-job-id-123"

    @pytest.mark.asyncio
    async def test_resolve_permalink(self, db_session):
        """Resolving a permalink returns job_id and snapshot."""
        permalink = BatchPermalink(
            short_id="resolve01",
            job_id="job-to-resolve",
            snapshot_data=json.dumps({"count": 42}),
        )
        db_session.add(permalink)
        await db_session.commit()

        # Resolve
        result = await db_session.execute(
            select(BatchPermalink).where(BatchPermalink.short_id == "resolve01")
        )
        found = result.scalar_one_or_none()
        assert found is not None
        assert found.job_id == "job-to-resolve"
        assert json.loads(found.snapshot_data) == {"count": 42}

    @pytest.mark.asyncio
    async def test_resolve_nonexistent(self, db_session):
        """Resolving a non-existent short_id returns None."""
        result = await db_session.execute(
            select(BatchPermalink).where(BatchPermalink.short_id == "nope1234")
        )
        assert result.scalar_one_or_none() is None

    @pytest.mark.asyncio
    async def test_resolve_expired_permalink(self, db_session):
        """Expired permalink is detected via expires_at check."""
        past = datetime.now(timezone.utc) - timedelta(hours=1)
        permalink = BatchPermalink(
            short_id="expired1",
            job_id="old-job",
            expires_at=past,
        )
        db_session.add(permalink)
        await db_session.commit()

        # Resolve and check expiry
        result = await db_session.execute(
            select(BatchPermalink).where(BatchPermalink.short_id == "expired1")
        )
        found = result.scalar_one_or_none()
        assert found is not None

        # Expiry check (as done in the route)
        expires = found.expires_at
        if expires.tzinfo is None:
            expires = expires.replace(tzinfo=timezone.utc)
        assert expires < datetime.now(timezone.utc)

    @pytest.mark.asyncio
    async def test_unique_short_id(self, db_session):
        """Short IDs must be unique."""
        p1 = BatchPermalink(short_id="unique01", job_id="job-1")
        db_session.add(p1)
        await db_session.commit()

        p2 = BatchPermalink(short_id="unique01", job_id="job-2")
        db_session.add(p2)
        with pytest.raises(Exception):
            await db_session.commit()


class TestSingleMoleculePermalink:
    """Tests for stateless single molecule permalinks."""

    def test_single_molecule_permalink_encoding(self):
        """SMILES is URL-encoded in the permalink."""
        from urllib.parse import quote

        smiles = "CC(=O)Oc1ccccc1C(O)=O"
        encoded = quote(smiles, safe="")
        url = f"http://localhost:3002/?smiles={encoded}"

        assert "CC" in url
        assert "%28" in url  # ( is encoded
        assert "%3D" in url  # = is encoded
