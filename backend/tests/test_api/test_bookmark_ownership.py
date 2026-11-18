"""
HTTP-level tests that bookmark by-id routes only reach the caller's own rows.

Regression tests for Claude Security findings F1, F3, F5, F6, F8: these routes used to
look bookmarks up by primary key only and relied on Postgres RLS, which is not in force
on SQLite or on a fresh-DB deployment.
"""

from unittest.mock import patch

import pytest
from httpx import ASGITransport, AsyncClient
from sqlalchemy import select

from app.core.security import generate_csrf_token
from app.db import get_db
from app.db.models.bookmark import Bookmark
from app.main import app


@pytest.fixture
async def db_session():
    """In-memory SQLite session shared by the app and the test."""
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
    factory = async_sessionmaker(engine, expire_on_commit=False)
    async with factory() as session:
        yield session
    await engine.dispose()


@pytest.fixture
async def api(db_session):
    """App client wired to the SQLite session, with the Postgres-only RLS call stubbed."""

    async def _get_db():
        yield db_session

    async def _no_rls(*_args, **_kwargs):
        return None

    app.dependency_overrides[get_db] = _get_db
    with patch("app.api.routes.bookmarks.set_rls_context", _no_rls):
        async with AsyncClient(transport=ASGITransport(app=app), base_url="http://test") as c:
            yield c
    app.dependency_overrides.pop(get_db, None)


async def _add(db, session_id, smiles="CCO"):
    bm = Bookmark(smiles=smiles, name=smiles, session_id=session_id)
    db.add(bm)
    await db.commit()
    await db.refresh(bm)
    return bm.id


def _as(session_id):
    """Cookie and CSRF header for a browser session."""
    return {
        "cookies": {"chemaudit_sid": session_id},
        "headers": {"X-CSRF-Token": generate_csrf_token()},
    }


async def _exists(db, bookmark_id):
    db.expire_all()
    res = await db.execute(select(Bookmark).where(Bookmark.id == bookmark_id))
    return res.scalar_one_or_none()


async def test_get_other_sessions_bookmark_404(api, db_session):
    bid = await _add(db_session, "A")
    assert (await api.get(f"/api/v1/bookmarks/{bid}", **_as("B"))).status_code == 404
    assert (await api.get(f"/api/v1/bookmarks/{bid}")).status_code == 404
    assert (await api.get(f"/api/v1/bookmarks/{bid}", **_as("A"))).status_code == 200


async def test_update_other_sessions_bookmark_404(api, db_session):
    bid = await _add(db_session, "A")
    r = await api.put(f"/api/v1/bookmarks/{bid}", json={"name": "pwned"}, **_as("B"))
    assert r.status_code == 404
    assert (await _exists(db_session, bid)).name == "CCO"


async def test_delete_other_sessions_bookmark_404(api, db_session):
    bid = await _add(db_session, "A")
    assert (await api.delete(f"/api/v1/bookmarks/{bid}", **_as("B"))).status_code == 404
    assert await _exists(db_session, bid) is not None


async def test_bulk_delete_only_own(api, db_session):
    mine = await _add(db_session, "A")
    theirs = await _add(db_session, "B")
    r = await api.delete(f"/api/v1/bookmarks/bulk?ids={mine}&ids={theirs}", **_as("A"))
    assert r.status_code == 204
    assert await _exists(db_session, mine) is None
    assert await _exists(db_session, theirs) is not None


async def test_bulk_delete_rejects_over_1000_ids(api, db_session):
    ids = "&".join(f"ids={i}" for i in range(1001))
    assert (await api.delete(f"/api/v1/bookmarks/bulk?{ids}", **_as("A"))).status_code == 422


async def test_batch_submit_only_own(api, db_session):
    mine = await _add(db_session, "A")
    theirs = await _add(db_session, "B", smiles="CCN")
    with (
        patch("app.api.routes.bookmarks.process_batch_job") as run,
        patch("app.api.routes.bookmarks.store_job_owner") as own,
    ):
        r = await api.post(
            "/api/v1/bookmarks/batch-submit", json={"bookmark_ids": [theirs]}, **_as("A")
        )
        assert r.status_code == 404
        r = await api.post(
            "/api/v1/bookmarks/batch-submit", json={"bookmark_ids": [mine, theirs]}, **_as("A")
        )
    assert r.status_code == 200
    assert r.json()["molecule_count"] == 1
    assert [m["smiles"] for m in run.call_args.args[1]] == ["CCO"]
    own.assert_called_once_with(r.json()["job_id"], "A")


async def test_batch_submit_rejects_over_1000_ids(api):
    r = await api.post(
        "/api/v1/bookmarks/batch-submit", json={"bookmark_ids": list(range(1001))}, **_as("A")
    )
    assert r.status_code == 422
