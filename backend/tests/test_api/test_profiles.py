"""
Tests for Scoring Profile CRUD API

Tests the ProfileService and profile API endpoints with 8 preset templates.
Uses SQLite in-memory for test isolation.
"""

import json

import pytest

from app.services.profiles.service import PRESET_PROFILES, ProfileService

A = ("session-a", None)
B = ("session-b", None)
ANON = (None, None)


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


@pytest.fixture
def service():
    """Create a ProfileService instance."""
    return ProfileService()


class TestPresetSeeding:
    """Tests for preset profile seeding."""

    @pytest.mark.asyncio
    async def test_seed_presets_creates_8(self, db_session, service):
        """Verify 8 presets are created after seeding."""
        created = await service.seed_presets(db_session)
        assert created == 8

        profiles = await service.list_all(db_session, owner=A)
        presets = [p for p in profiles if p.is_preset]
        assert len(presets) == 8

    @pytest.mark.asyncio
    async def test_seed_presets_idempotent(self, db_session, service):
        """Seeding twice does not create duplicates."""
        first = await service.seed_presets(db_session)
        second = await service.seed_presets(db_session)
        assert first == 8
        assert second == 0

        profiles = await service.list_all(db_session, owner=A)
        presets = [p for p in profiles if p.is_preset]
        assert len(presets) == 8

    @pytest.mark.asyncio
    async def test_preset_names_match(self, db_session, service):
        """Verify preset names match the constant definitions."""
        await service.seed_presets(db_session)
        profiles = await service.list_all(db_session, owner=A)
        preset_names = {p.name for p in profiles if p.is_preset}
        expected_names = {p["name"] for p in PRESET_PROFILES}
        assert preset_names == expected_names


class TestProfileCRUD:
    """Tests for user profile CRUD operations."""

    @pytest.mark.asyncio
    async def test_create_user_profile(self, db_session, service):
        """Creating a user profile sets is_preset=False."""
        profile = await service.create(
            db_session,
            name="My Custom Profile",
            description="Test description",
            thresholds={"mw": {"min": 100, "max": 600}},
            weights={"mw": 1.0},
            owner=A,
        )
        assert profile.name == "My Custom Profile"
        assert profile.is_preset is False
        assert profile.is_active is True
        assert json.loads(profile.thresholds) == {"mw": {"min": 100, "max": 600}}

    @pytest.mark.asyncio
    async def test_get_profile(self, db_session, service):
        """Get profile by ID."""
        profile = await service.create(
            db_session, name="Test", description=None, thresholds={}, weights={}, owner=A
        )
        fetched = await service.get(db_session, profile.id, owner=A)
        assert fetched is not None
        assert fetched.name == "Test"

    @pytest.mark.asyncio
    async def test_get_nonexistent_profile(self, db_session, service):
        """Get returns None for non-existent ID."""
        result = await service.get(db_session, 99999, owner=A)
        assert result is None

    @pytest.mark.asyncio
    async def test_update_user_profile(self, db_session, service):
        """Update changes the specified fields."""
        profile = await service.create(
            db_session, name="Original", description=None, thresholds={}, weights={}, owner=A
        )
        updated = await service.update(
            db_session, profile.id, {"name": "Updated", "thresholds": {"logp": {"max": 5}}}, owner=A
        )
        assert updated is not None
        assert updated.name == "Updated"
        assert json.loads(updated.thresholds) == {"logp": {"max": 5}}

    @pytest.mark.asyncio
    async def test_cannot_update_preset(self, db_session, service):
        """Updating a preset raises ValueError."""
        await service.seed_presets(db_session)
        profiles = await service.list_all(db_session, owner=A)
        preset = next(p for p in profiles if p.is_preset)

        with pytest.raises(ValueError, match="Cannot modify a preset"):
            await service.update(db_session, preset.id, {"name": "Hacked"}, owner=A)

    @pytest.mark.asyncio
    async def test_delete_user_profile(self, db_session, service):
        """Deleting soft-deletes by setting is_active=False."""
        profile = await service.create(
            db_session, name="ToDelete", description=None, thresholds={}, weights={}, owner=A
        )
        deleted = await service.delete(db_session, profile.id, owner=A)
        assert deleted is True

        # Should not be visible anymore
        fetched = await service.get(db_session, profile.id, owner=A)
        assert fetched is None

    @pytest.mark.asyncio
    async def test_cannot_delete_preset(self, db_session, service):
        """Deleting a preset raises ValueError."""
        await service.seed_presets(db_session)
        profiles = await service.list_all(db_session, owner=A)
        preset = next(p for p in profiles if p.is_preset)

        with pytest.raises(ValueError, match="Cannot delete a preset"):
            await service.delete(db_session, preset.id, owner=A)

    @pytest.mark.asyncio
    async def test_duplicate_preset(self, db_session, service):
        """Duplicating a preset creates a non-preset copy."""
        await service.seed_presets(db_session)
        profiles = await service.list_all(db_session, owner=A)
        preset = next(p for p in profiles if p.is_preset)

        copy = await service.duplicate(db_session, preset.id, "My Custom Lipinski", owner=A)
        assert copy is not None
        assert copy.name == "My Custom Lipinski"
        assert copy.is_preset is False
        assert copy.thresholds == preset.thresholds
        assert copy.weights == preset.weights

    @pytest.mark.asyncio
    async def test_export_import_profile(self, db_session, service):
        """Export and re-import produces identical thresholds/weights."""
        profile = await service.create(
            db_session,
            name="Export Test",
            description="For export",
            thresholds={"mw": {"min": 200, "max": 500}},
            weights={"mw": 1.5, "logp": 0.8},
            owner=A,
        )
        exported = await service.export_json(db_session, profile.id, owner=A)
        assert exported is not None

        imported = await service.import_json(db_session, exported, owner=A)
        assert imported.name == "Export Test"
        assert json.loads(imported.thresholds) == {"mw": {"min": 200, "max": 500}}
        assert json.loads(imported.weights) == {"mw": 1.5, "logp": 0.8}


class TestProfileOwnership:
    """Custom profiles are visible and mutable only by their creator (F7, F16)."""

    async def _make(self, db, service, owner=A):
        return await service.create(
            db, name="mine", description=None, thresholds={}, weights={}, owner=owner
        )

    async def test_other_owner_cannot_get_update_delete(self, db_session, service):
        p = await self._make(db_session, service)
        assert await service.get(db_session, p.id, owner=B) is None
        assert await service.update(db_session, p.id, {"name": "x"}, owner=B) is None
        assert await service.delete(db_session, p.id, owner=B) is False
        assert await service.export_json(db_session, p.id, owner=B) is None
        assert await service.duplicate(db_session, p.id, "copy", owner=B) is None
        assert (await service.get(db_session, p.id, owner=A)).name == "mine"
        assert (await service.update(db_session, p.id, {"name": "x"}, owner=A)).name == "x"
        assert await service.delete(db_session, p.id, owner=A) is True

    async def test_list_all_shows_presets_and_own_only(self, db_session, service):
        await service.seed_presets(db_session)
        await self._make(db_session, service)
        assert len(await service.list_all(db_session, owner=B)) == len(PRESET_PROFILES)
        assert len(await service.list_all(db_session, owner=A)) == len(PRESET_PROFILES) + 1

    async def test_anonymous_sees_presets_only(self, db_session, service):
        await service.seed_presets(db_session)
        p = await self._make(db_session, service)
        assert len(await service.list_all(db_session, owner=ANON)) == len(PRESET_PROFILES)
        assert await service.get(db_session, p.id, owner=ANON) is None
        preset = (await service.list_all(db_session, owner=ANON))[0]
        assert await service.get(db_session, preset.id, owner=ANON) is not None

    async def test_api_key_owner_isolated(self, db_session, service):
        p = await self._make(db_session, service, owner=(None, "h1"))
        assert await service.get(db_session, p.id, owner=(None, "h2")) is None
        assert await service.get(db_session, p.id, owner=(None, "h1")) is not None

    async def test_legacy_unowned_profile_hidden(self, db_session, service):
        from app.db.models.profile import ScoringProfile

        legacy = ScoringProfile(name="old", thresholds="{}", weights="{}", is_preset=False)
        db_session.add(legacy)
        await db_session.commit()
        assert await service.get(db_session, legacy.id, owner=A) is None

    async def test_new_rows_store_owner(self, db_session, service):
        p = await self._make(db_session, service, owner=(None, "h1"))
        copy = await service.duplicate(db_session, p.id, "copy", owner=(None, "h1"))
        imported = await service.import_json(db_session, {"name": "imp"}, owner=A)
        assert (p.api_key_hash, copy.api_key_hash) == ("h1", "h1")
        assert imported.session_id == "session-a"


async def test_purge_my_data_deletes_own_profiles(db_session, service):
    """DELETE /me/data also removes the caller's custom profiles, not other users'."""
    from httpx import ASGITransport, AsyncClient

    from app.core.security import generate_csrf_token
    from app.db import get_db
    from app.db.models.profile import ScoringProfile
    from app.main import app

    await service.create(db_session, "a", None, {}, {}, owner=A)
    await service.create(db_session, "b", None, {}, {}, owner=B)

    async def _get_db():
        yield db_session

    app.dependency_overrides[get_db] = _get_db
    try:
        async with AsyncClient(transport=ASGITransport(app=app), base_url="http://t") as c:
            r = await c.delete(
                "/api/v1/me/data",
                cookies={"chemaudit_sid": "session-a"},
                headers={"X-CSRF-Token": generate_csrf_token()},
            )
    finally:
        app.dependency_overrides.pop(get_db, None)
    assert r.status_code == 200
    assert r.json()["deleted"]["profiles"] == 1
    db_session.expire_all()
    from sqlalchemy import select

    names = (await db_session.execute(select(ScoringProfile.name))).scalars().all()
    assert names == ["b"]
