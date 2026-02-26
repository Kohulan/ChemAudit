"""
Permalink API Routes

Endpoints for creating and resolving shareable permalinks for results.
Batch permalinks are stored in the database with optional expiry.
Single molecule permalinks are stateless (URL-encoded).
"""

import json
import secrets
from datetime import datetime, timezone
from typing import Optional
from urllib.parse import quote

from fastapi import APIRouter, Depends, HTTPException, Query, Request
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from app.core.config import settings
from app.core.ownership import verify_job_access
from app.core.rate_limit import get_rate_limit_key, limiter
from app.db import get_db
from app.db.models.permalink import BatchPermalink
from app.schemas.permalinks import (
    PermalinkCreateRequest,
    PermalinkResolveResponse,
    PermalinkResponse,
)

router = APIRouter()


@router.post("/permalinks", response_model=PermalinkResponse, status_code=201)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def create_permalink(
    request: Request,
    body: PermalinkCreateRequest,
    db: AsyncSession = Depends(get_db),
):
    """
    Create a batch permalink with a unique short ID.

    Stores the job_id and optional snapshot in the database.
    Returns a shareable URL.
    """
    verify_job_access(request, body.job_id)

    short_id = secrets.token_urlsafe(8)

    permalink = BatchPermalink(
        short_id=short_id,
        job_id=body.job_id,
        snapshot_data=json.dumps(body.snapshot_data) if body.snapshot_data else None,
        settings=json.dumps(body.settings) if body.settings else None,
    )
    db.add(permalink)
    await db.commit()
    await db.refresh(permalink)

    return PermalinkResponse(
        short_id=short_id,
        job_id=body.job_id,
        url=f"{settings.BASE_URL}/report/{short_id}",
        created_at=permalink.created_at,
        expires_at=permalink.expires_at,
    )


@router.get("/report/{short_id}", response_model=PermalinkResolveResponse)
async def resolve_permalink(
    short_id: str,
    db: AsyncSession = Depends(get_db),
):
    """
    Resolve a batch permalink by short ID.

    Returns the job_id and snapshot data for the frontend to load results.
    Returns 404 if not found, 410 Gone if expired.
    """
    result = await db.execute(
        select(BatchPermalink).where(BatchPermalink.short_id == short_id)
    )
    permalink = result.scalar_one_or_none()

    if permalink is None:
        raise HTTPException(status_code=404, detail="Permalink not found")

    # Check expiry
    if permalink.expires_at is not None:
        now = datetime.now(timezone.utc)
        expires = permalink.expires_at
        # Handle timezone-naive datetimes from SQLite
        if expires.tzinfo is None:
            from datetime import timezone as tz

            expires = expires.replace(tzinfo=tz.utc)
        if expires < now:
            raise HTTPException(
                status_code=410, detail="This report link has expired"
            )

    return PermalinkResolveResponse(
        job_id=permalink.job_id,
        snapshot_data=json.loads(permalink.snapshot_data) if permalink.snapshot_data else None,
        settings=json.loads(permalink.settings) if permalink.settings else None,
    )


@router.get("/permalinks/single")
async def single_molecule_permalink(
    smiles: str = Query(..., description="SMILES string"),
    settings_json: Optional[str] = Query(
        None, alias="settings", description="JSON settings string"
    ),
):
    """
    Generate a stateless single molecule permalink.

    No database storage needed. URL-encodes SMILES and optional settings.
    """
    encoded_smiles = quote(smiles, safe="")
    url = f"{settings.BASE_URL}/?smiles={encoded_smiles}"
    if settings_json:
        encoded_settings = quote(settings_json, safe="")
        url += f"&settings={encoded_settings}"

    return {"url": url, "smiles": smiles}
