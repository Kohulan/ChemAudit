"""
Bookmark API Routes

CRUD endpoints for molecule bookmarks with tag filtering and batch-submit workflow.
"""

import uuid
from typing import List, Optional

from fastapi import APIRouter, Depends, HTTPException, Query
from rdkit import Chem
from rdkit.Chem import inchi as rdkit_inchi
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from app.db import get_db
from app.db.models.bookmark import Bookmark
from app.schemas.bookmarks import (
    BookmarkBatchSubmit,
    BookmarkCreate,
    BookmarkResponse,
    BookmarkUpdate,
)
from app.services.batch.tasks import process_batch_job

router = APIRouter()


def _tags_to_str(tags: Optional[List[str]]) -> Optional[str]:
    """Serialize list of tags to comma-separated string for DB storage."""
    if tags is None:
        return None
    return ",".join(t.strip() for t in tags if t.strip())


def _str_to_tags(tags_str: Optional[str]) -> List[str]:
    """Deserialize comma-separated tag string from DB to list."""
    if not tags_str:
        return []
    return [t.strip() for t in tags_str.split(",") if t.strip()]


def _compute_inchikey(smiles: str) -> Optional[str]:
    """Compute InChIKey from SMILES using RDKit. Returns None if invalid."""
    try:
        if not smiles or not smiles.strip():
            return None
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol_inchi = rdkit_inchi.MolToInchi(mol)
        if not mol_inchi:
            return None
        key = rdkit_inchi.MolToInchiKey(mol)
        return key if key else None
    except Exception:
        return None


def _bookmark_to_response(bm: Bookmark) -> dict:
    """Convert a Bookmark ORM instance to a response dict."""
    return {
        "id": bm.id,
        "smiles": bm.smiles,
        "name": bm.name,
        "inchikey": bm.inchikey,
        "tags": _str_to_tags(bm.tags),
        "notes": bm.notes,
        "source": bm.source,
        "job_id": bm.job_id,
        "created_at": bm.created_at,
    }


@router.get("/bookmarks", response_model=list[BookmarkResponse])
async def list_bookmarks(
    page: int = Query(default=1, ge=1),
    page_size: int = Query(default=50, ge=1, le=200),
    tags: Optional[str] = Query(None, description="Comma-separated tag filter"),
    source: Optional[str] = Query(None, description="Source filter"),
    search: Optional[str] = Query(None, description="SMILES substring search"),
    db: AsyncSession = Depends(get_db),
):
    """List bookmarks with optional filters and pagination."""
    query = select(Bookmark)

    # Apply filters
    if tags:
        tag_list = [t.strip() for t in tags.split(",") if t.strip()]
        for tag in tag_list:
            query = query.where(Bookmark.tags.contains(tag))
    if source:
        query = query.where(Bookmark.source == source)
    if search:
        query = query.where(Bookmark.smiles.contains(search))

    # Order and paginate
    query = query.order_by(Bookmark.created_at.desc())
    offset = (page - 1) * page_size
    query = query.offset(offset).limit(page_size)

    result = await db.execute(query)
    bookmarks = result.scalars().all()
    return [_bookmark_to_response(bm) for bm in bookmarks]


@router.get("/bookmarks/{bookmark_id}", response_model=BookmarkResponse)
async def get_bookmark(bookmark_id: int, db: AsyncSession = Depends(get_db)):
    """Get a single bookmark by ID."""
    result = await db.execute(select(Bookmark).where(Bookmark.id == bookmark_id))
    bm = result.scalar_one_or_none()
    if bm is None:
        raise HTTPException(status_code=404, detail="Bookmark not found")
    return _bookmark_to_response(bm)


@router.post("/bookmarks", response_model=BookmarkResponse, status_code=201)
async def create_bookmark(body: BookmarkCreate, db: AsyncSession = Depends(get_db)):
    """Create a new bookmark. Auto-computes InChIKey from SMILES."""
    inchikey = _compute_inchikey(body.smiles)

    bm = Bookmark(
        smiles=body.smiles,
        name=body.name,
        inchikey=inchikey,
        tags=_tags_to_str(body.tags),
        notes=body.notes,
        source=body.source,
        job_id=body.job_id,
    )
    db.add(bm)
    await db.commit()
    await db.refresh(bm)
    return _bookmark_to_response(bm)


@router.put("/bookmarks/{bookmark_id}", response_model=BookmarkResponse)
async def update_bookmark(
    bookmark_id: int,
    body: BookmarkUpdate,
    db: AsyncSession = Depends(get_db),
):
    """Update bookmark metadata (name, tags, notes)."""
    result = await db.execute(select(Bookmark).where(Bookmark.id == bookmark_id))
    bm = result.scalar_one_or_none()
    if bm is None:
        raise HTTPException(status_code=404, detail="Bookmark not found")

    if body.name is not None:
        bm.name = body.name
    if body.tags is not None:
        bm.tags = _tags_to_str(body.tags)
    if body.notes is not None:
        bm.notes = body.notes

    await db.commit()
    await db.refresh(bm)
    return _bookmark_to_response(bm)


@router.delete("/bookmarks/{bookmark_id}", status_code=204)
async def delete_bookmark(bookmark_id: int, db: AsyncSession = Depends(get_db)):
    """Delete a bookmark (hard delete)."""
    result = await db.execute(select(Bookmark).where(Bookmark.id == bookmark_id))
    bm = result.scalar_one_or_none()
    if bm is None:
        raise HTTPException(status_code=404, detail="Bookmark not found")

    await db.delete(bm)
    await db.commit()


@router.post("/bookmarks/batch-submit")
async def bookmark_batch_submit(
    body: BookmarkBatchSubmit, db: AsyncSession = Depends(get_db)
):
    """
    Submit bookmarked molecules as a new batch job.

    Looks up SMILES for each bookmark ID and starts batch processing.
    Returns the new job_id.
    """
    # Fetch bookmarks
    result = await db.execute(
        select(Bookmark).where(Bookmark.id.in_(body.bookmark_ids))
    )
    bookmarks = result.scalars().all()

    if not bookmarks:
        raise HTTPException(
            status_code=404, detail="No bookmarks found for the given IDs"
        )

    # Build molecule dicts for batch processing
    mol_dicts = [
        {
            "smiles": bm.smiles,
            "name": bm.name or f"bookmark_{bm.id}",
            "index": i,
            "properties": {},
            "parse_error": None,
        }
        for i, bm in enumerate(bookmarks)
    ]

    job_id = str(uuid.uuid4())
    process_batch_job(job_id, mol_dicts, safety_options={})

    return {
        "job_id": job_id,
        "molecule_count": len(mol_dicts),
        "message": f"Batch job created from {len(mol_dicts)} bookmarks",
    }


@router.delete("/bookmarks/bulk", status_code=204)
async def bulk_delete_bookmarks(
    ids: List[int] = Query(..., description="Bookmark IDs to delete"),
    db: AsyncSession = Depends(get_db),
):
    """Bulk delete bookmarks by IDs."""
    result = await db.execute(select(Bookmark).where(Bookmark.id.in_(ids)))
    bookmarks = result.scalars().all()

    for bm in bookmarks:
        await db.delete(bm)
    await db.commit()
