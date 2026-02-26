"""
Tests for Bookmark CRUD API

Tests bookmark creation with InChIKey auto-computation, tag filtering,
and bookmark-to-batch workflow.
Uses SQLite in-memory for test isolation.
"""

import pytest

from app.api.routes.bookmarks import _compute_inchikey, _str_to_tags, _tags_to_str
from app.db.models.bookmark import Bookmark


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


class TestTagSerialization:
    """Tests for tag serialization helpers."""

    def test_tags_to_str(self):
        """Tags list serializes to comma-separated string."""
        assert _tags_to_str(["drug", "lead"]) == "drug,lead"

    def test_tags_to_str_none(self):
        """None input returns None."""
        assert _tags_to_str(None) is None

    def test_tags_to_str_empty(self):
        """Empty list returns empty string."""
        assert _tags_to_str([]) == ""

    def test_str_to_tags(self):
        """Comma-separated string deserializes to tag list."""
        assert _str_to_tags("drug,lead") == ["drug", "lead"]

    def test_str_to_tags_none(self):
        """None input returns empty list."""
        assert _str_to_tags(None) == []

    def test_str_to_tags_empty(self):
        """Empty string returns empty list."""
        assert _str_to_tags("") == []


class TestInChIKeyComputation:
    """Tests for InChIKey auto-computation."""

    def test_valid_smiles_ethanol(self):
        """Valid SMILES computes InChIKey."""
        key = _compute_inchikey("CCO")
        assert key is not None
        assert len(key) == 27  # Standard InChIKey length

    def test_valid_smiles_aspirin(self):
        """Aspirin SMILES computes InChIKey."""
        key = _compute_inchikey("CC(=O)Oc1ccccc1C(O)=O")
        assert key is not None

    def test_invalid_smiles(self):
        """Invalid SMILES returns None."""
        key = _compute_inchikey("not_a_smiles_string!!!")
        assert key is None

    def test_empty_smiles(self):
        """Empty SMILES returns None."""
        key = _compute_inchikey("")
        assert key is None


class TestBookmarkCRUD:
    """Tests for bookmark CRUD operations via ORM."""

    @pytest.mark.asyncio
    async def test_create_bookmark(self, db_session):
        """Creating a bookmark persists all fields."""
        bm = Bookmark(
            smiles="CCO",
            name="Ethanol",
            inchikey=_compute_inchikey("CCO"),
            tags="drug,solvent",
            notes="Test molecule",
            source="manual",
        )
        db_session.add(bm)
        await db_session.commit()
        await db_session.refresh(bm)

        assert bm.id is not None
        assert bm.smiles == "CCO"
        assert bm.name == "Ethanol"
        assert bm.inchikey is not None

    @pytest.mark.asyncio
    async def test_create_bookmark_invalid_smiles(self, db_session):
        """Bookmark with invalid SMILES still persists (inchikey=None)."""
        bm = Bookmark(
            smiles="invalid!!!",
            name="Bad molecule",
            inchikey=None,
            tags=None,
        )
        db_session.add(bm)
        await db_session.commit()
        await db_session.refresh(bm)

        assert bm.id is not None
        assert bm.inchikey is None

    @pytest.mark.asyncio
    async def test_update_bookmark(self, db_session):
        """Updating bookmark changes specified fields."""
        bm = Bookmark(smiles="CCO", name="Old Name", tags="old")
        db_session.add(bm)
        await db_session.commit()

        bm.name = "New Name"
        bm.tags = "new,updated"
        await db_session.commit()
        await db_session.refresh(bm)

        assert bm.name == "New Name"
        assert _str_to_tags(bm.tags) == ["new", "updated"]

    @pytest.mark.asyncio
    async def test_delete_bookmark(self, db_session):
        """Deleting a bookmark removes it from the database."""
        from sqlalchemy import select

        bm = Bookmark(smiles="CCO", name="To Delete")
        db_session.add(bm)
        await db_session.commit()
        bm_id = bm.id

        await db_session.delete(bm)
        await db_session.commit()

        result = await db_session.execute(
            select(Bookmark).where(Bookmark.id == bm_id)
        )
        assert result.scalar_one_or_none() is None

    @pytest.mark.asyncio
    async def test_list_bookmarks_paginated(self, db_session):
        """Pagination returns correct subset."""
        from sqlalchemy import func, select

        for i in range(5):
            bm = Bookmark(smiles=f"C{'C' * i}O", name=f"Mol {i}")
            db_session.add(bm)
        await db_session.commit()

        # Count total
        count_result = await db_session.execute(select(func.count(Bookmark.id)))
        total = count_result.scalar()
        assert total == 5

        # Paginate: page 1, size 2
        result = await db_session.execute(
            select(Bookmark).order_by(Bookmark.id).offset(0).limit(2)
        )
        page1 = result.scalars().all()
        assert len(page1) == 2

    @pytest.mark.asyncio
    async def test_filter_by_tags(self, db_session):
        """Tag filtering returns only matching bookmarks."""
        from sqlalchemy import select

        bm1 = Bookmark(smiles="CCO", name="Ethanol", tags="drug,solvent")
        bm2 = Bookmark(smiles="CC(=O)O", name="Acetic acid", tags="reagent")
        bm3 = Bookmark(smiles="c1ccccc1", name="Benzene", tags="solvent")
        db_session.add_all([bm1, bm2, bm3])
        await db_session.commit()

        # Filter by "solvent"
        result = await db_session.execute(
            select(Bookmark).where(Bookmark.tags.contains("solvent"))
        )
        filtered = result.scalars().all()
        assert len(filtered) == 2
        names = {bm.name for bm in filtered}
        assert names == {"Ethanol", "Benzene"}
