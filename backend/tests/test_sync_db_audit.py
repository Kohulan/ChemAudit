"""Synchronous Celery DB path: audit writes and session purge use a sync engine.

Verifies the dual-engine design (3.5) without asyncio.run() bridging. Uses an
in-memory SQLite engine so it runs in any environment.
"""

from datetime import datetime, timedelta

import pytest
from sqlalchemy import create_engine, select
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from app.db import Base, _sync_db_url
from app.db.models.audit import ValidationAuditEntry
from app.db.models.bookmark import Bookmark
from app.services.audit.service import _batch_outcome, log_batch_event_sync


@pytest.fixture
def sync_db(monkeypatch):
    engine = create_engine(
        "sqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    Base.metadata.create_all(engine)
    session_local = sessionmaker(engine, expire_on_commit=False)
    # Point app.db.sync_session at the in-memory DB for code under test.
    monkeypatch.setattr("app.db.sync_session", session_local)
    yield session_local
    engine.dispose()


def test_sync_db_url_maps_async_drivers():
    assert _sync_db_url("postgresql+asyncpg://u:p@h/db") == "postgresql+psycopg://u:p@h/db"
    assert _sync_db_url("sqlite+aiosqlite:///x.db") == "sqlite+pysqlite:///x.db"


@pytest.mark.parametrize(
    "molecule_count,fail_count,expected",
    [(10, 0, "pass"), (10, 2, "warn"), (10, 8, "fail"), (2, 1, "fail")],
)
def test_batch_outcome(molecule_count, fail_count, expected):
    assert _batch_outcome(molecule_count, fail_count) == expected


def test_log_batch_event_sync_inserts_row(sync_db):
    with sync_db() as db:
        entry = log_batch_event_sync(
            db, job_id="job-x", molecule_count=10, pass_count=8, fail_count=2
        )
        assert entry.id is not None
        assert entry.source == "batch"
        assert entry.outcome == "warn"

    with sync_db() as db:
        rows = (
            db.execute(
                select(ValidationAuditEntry).where(ValidationAuditEntry.job_id == "job-x")
            )
            .scalars()
            .all()
        )
    assert len(rows) == 1
    assert rows[0].molecule_count == 10


def test_tasks_log_batch_audit_uses_sync_session(sync_db):
    """The Celery audit helper writes via the sync session (no asyncio.run)."""
    from app.services.batch import tasks

    tasks._log_batch_audit(job_id="job-y", molecule_count=4, pass_count=4, fail_count=0)

    with sync_db() as db:
        rows = (
            db.execute(
                select(ValidationAuditEntry).where(ValidationAuditEntry.job_id == "job-y")
            )
            .scalars()
            .all()
        )
    assert len(rows) == 1
    assert rows[0].outcome == "pass"


def test_purge_expired_sessions_sync(sync_db):
    """The purge task deletes only old, session-scoped rows via the sync engine."""
    from app.services.session.cleanup import purge_expired_sessions

    old = datetime.utcnow() - timedelta(days=60)
    recent = datetime.utcnow()

    with sync_db() as db:
        # Old anonymous audit row — should be purged.
        db.add(
            ValidationAuditEntry(
                smiles="", outcome="pass", source="batch", session_id="s1",
                api_key_hash=None, created_at=old,
            )
        )
        # Recent anonymous audit row — should survive.
        db.add(
            ValidationAuditEntry(
                smiles="", outcome="pass", source="batch", session_id="s1",
                api_key_hash=None, created_at=recent,
            )
        )
        # Old API-key row — must be preserved (not session-scoped).
        db.add(
            ValidationAuditEntry(
                smiles="", outcome="pass", source="batch", session_id=None,
                api_key_hash="hash", created_at=old,
            )
        )
        # Old anonymous bookmark — should be purged.
        db.add(Bookmark(smiles="CCO", session_id="s1", api_key_hash=None, created_at=old))
        db.commit()

    result = purge_expired_sessions()

    assert result["audit_entries_deleted"] == 1
    assert result["bookmarks_deleted"] == 1
    with sync_db() as db:
        remaining = db.execute(select(ValidationAuditEntry)).scalars().all()
    assert len(remaining) == 2  # recent anonymous + old api-key
