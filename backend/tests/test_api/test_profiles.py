"""
Tests for Scoring Profile CRUD API

Tests the ProfileService and profile API endpoints with 8 preset templates.
Uses SQLite in-memory for test isolation.
"""

import json

import pytest

from app.services.profiles.service import PRESET_PROFILES, ProfileService


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

        profiles = await service.list_all(db_session)
        presets = [p for p in profiles if p.is_preset]
        assert len(presets) == 8

    @pytest.mark.asyncio
    async def test_seed_presets_idempotent(self, db_session, service):
        """Seeding twice does not create duplicates."""
        first = await service.seed_presets(db_session)
        second = await service.seed_presets(db_session)
        assert first == 8
        assert second == 0

        profiles = await service.list_all(db_session)
        presets = [p for p in profiles if p.is_preset]
        assert len(presets) == 8

    @pytest.mark.asyncio
    async def test_preset_names_match(self, db_session, service):
        """Verify preset names match the constant definitions."""
        await service.seed_presets(db_session)
        profiles = await service.list_all(db_session)
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
        )
        assert profile.name == "My Custom Profile"
        assert profile.is_preset is False
        assert profile.is_active is True
        assert json.loads(profile.thresholds) == {"mw": {"min": 100, "max": 600}}

    @pytest.mark.asyncio
    async def test_get_profile(self, db_session, service):
        """Get profile by ID."""
        profile = await service.create(
            db_session, name="Test", description=None, thresholds={}, weights={}
        )
        fetched = await service.get(db_session, profile.id)
        assert fetched is not None
        assert fetched.name == "Test"

    @pytest.mark.asyncio
    async def test_get_nonexistent_profile(self, db_session, service):
        """Get returns None for non-existent ID."""
        result = await service.get(db_session, 99999)
        assert result is None

    @pytest.mark.asyncio
    async def test_update_user_profile(self, db_session, service):
        """Update changes the specified fields."""
        profile = await service.create(
            db_session, name="Original", description=None, thresholds={}, weights={}
        )
        updated = await service.update(
            db_session, profile.id, {"name": "Updated", "thresholds": {"logp": {"max": 5}}}
        )
        assert updated is not None
        assert updated.name == "Updated"
        assert json.loads(updated.thresholds) == {"logp": {"max": 5}}

    @pytest.mark.asyncio
    async def test_cannot_update_preset(self, db_session, service):
        """Updating a preset raises ValueError."""
        await service.seed_presets(db_session)
        profiles = await service.list_all(db_session)
        preset = next(p for p in profiles if p.is_preset)

        with pytest.raises(ValueError, match="Cannot modify a preset"):
            await service.update(db_session, preset.id, {"name": "Hacked"})

    @pytest.mark.asyncio
    async def test_delete_user_profile(self, db_session, service):
        """Deleting soft-deletes by setting is_active=False."""
        profile = await service.create(
            db_session, name="ToDelete", description=None, thresholds={}, weights={}
        )
        deleted = await service.delete(db_session, profile.id)
        assert deleted is True

        # Should not be visible anymore
        fetched = await service.get(db_session, profile.id)
        assert fetched is None

    @pytest.mark.asyncio
    async def test_cannot_delete_preset(self, db_session, service):
        """Deleting a preset raises ValueError."""
        await service.seed_presets(db_session)
        profiles = await service.list_all(db_session)
        preset = next(p for p in profiles if p.is_preset)

        with pytest.raises(ValueError, match="Cannot delete a preset"):
            await service.delete(db_session, preset.id)

    @pytest.mark.asyncio
    async def test_duplicate_preset(self, db_session, service):
        """Duplicating a preset creates a non-preset copy."""
        await service.seed_presets(db_session)
        profiles = await service.list_all(db_session)
        preset = next(p for p in profiles if p.is_preset)

        copy = await service.duplicate(db_session, preset.id, "My Custom Lipinski")
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
        )
        exported = await service.export_json(db_session, profile.id)
        assert exported is not None

        imported = await service.import_json(db_session, exported)
        assert imported.name == "Export Test"
        assert json.loads(imported.thresholds) == {"mw": {"min": 200, "max": 500}}
        assert json.loads(imported.weights) == {"mw": 1.5, "logp": 0.8}
