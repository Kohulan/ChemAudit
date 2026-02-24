"""
Bookmark Schemas

Pydantic schemas for molecule bookmarks.
"""

from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel, Field


class BookmarkCreate(BaseModel):
    """Schema for creating a bookmark."""

    smiles: str = Field(..., min_length=1, description="SMILES string")
    name: Optional[str] = Field(None, max_length=200, description="Molecule name")
    tags: Optional[List[str]] = Field(None, description="List of tag strings")
    notes: Optional[str] = Field(None, description="User notes")
    source: Optional[str] = Field(None, max_length=50, description="Source context")
    job_id: Optional[str] = Field(None, max_length=50, description="Batch job ID")


class BookmarkUpdate(BaseModel):
    """Schema for updating a bookmark (all fields optional)."""

    name: Optional[str] = Field(None, max_length=200)
    tags: Optional[List[str]] = None
    notes: Optional[str] = None


class BookmarkResponse(BaseModel):
    """Schema for bookmark response."""

    id: int
    smiles: str
    name: Optional[str] = None
    inchikey: Optional[str] = None
    tags: List[str] = Field(default_factory=list)
    notes: Optional[str] = None
    source: Optional[str] = None
    job_id: Optional[str] = None
    created_at: datetime

    model_config = {"from_attributes": True}


class BookmarkBatchSubmit(BaseModel):
    """Request body for submitting bookmarks as a new batch job."""

    bookmark_ids: List[int] = Field(
        ..., min_length=1, description="IDs of bookmarks to submit"
    )
