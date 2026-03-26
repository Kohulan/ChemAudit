"""
Permalink API Routes

Endpoints for creating and resolving shareable permalinks for results.
Batch permalinks are stored in the database with optional expiry.
Single molecule permalinks are stateless (URL-encoded).
"""

import json
import logging
import secrets
from datetime import datetime, timezone
from typing import Optional
from urllib.parse import quote

from fastapi import APIRouter, Depends, HTTPException, Path, Query, Request
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

logger = logging.getLogger(__name__)

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

    # Auto-snapshot results from Redis so shared links survive past Redis TTL.
    # If caller already provided snapshot_data, use that instead.
    snapshot_data = body.snapshot_data
    if not snapshot_data:
        try:
            from dataclasses import asdict

            from app.services.batch.result_aggregator import result_storage

            all_results = result_storage.get_all_results(body.job_id)
            stats = result_storage.get_statistics(body.job_id)
            if all_results:
                snapshot_data = {
                    "results": all_results,
                    "statistics": asdict(stats) if stats else None,
                    "total_results": len(all_results),
                }
        except Exception:
            logger.warning(
                "Failed to auto-snapshot results for job %s", body.job_id, exc_info=True
            )

    permalink = BatchPermalink(
        short_id=short_id,
        job_id=body.job_id,
        snapshot_data=json.dumps(snapshot_data) if snapshot_data else None,
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
    short_id: str = Path(..., max_length=16, pattern=r"^[A-Za-z0-9_-]+$"),
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
            expires = expires.replace(tzinfo=timezone.utc)
        if expires < now:
            raise HTTPException(
                status_code=410, detail="This report link has expired"
            )

    snapshot = None
    if permalink.snapshot_data:
        try:
            snapshot = json.loads(permalink.snapshot_data)
        except (json.JSONDecodeError, TypeError):
            logger.warning("Corrupt snapshot_data for permalink %s", short_id)

    link_settings = None
    if permalink.settings:
        try:
            link_settings = json.loads(permalink.settings)
        except (json.JSONDecodeError, TypeError):
            logger.warning("Corrupt settings for permalink %s", short_id)

    return PermalinkResolveResponse(
        job_id=permalink.job_id,
        snapshot_data=snapshot,
        settings=link_settings,
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
