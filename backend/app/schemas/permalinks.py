"""
Permalink Schemas

Pydantic schemas for shareable result permalinks.
"""

from datetime import datetime
from typing import Any, Dict, Optional

from pydantic import BaseModel, Field


class PermalinkCreateRequest(BaseModel):
    """Request body for creating a batch permalink."""

    job_id: str = Field(..., description="Batch job ID")
    snapshot_data: Optional[Dict[str, Any]] = Field(
        None, description="Optional snapshot of results at creation time"
    )
    settings: Optional[Dict[str, Any]] = Field(
        None, description="Optional settings used for this job"
    )


class PermalinkResponse(BaseModel):
    """Response for a created/resolved permalink."""

    short_id: str
    job_id: str
    url: str
    created_at: datetime
    expires_at: Optional[datetime] = None


class PermalinkResolveResponse(BaseModel):
    """Response when resolving a permalink."""

    job_id: str
    snapshot_data: Optional[Dict[str, Any]] = None
    settings: Optional[Dict[str, Any]] = None
