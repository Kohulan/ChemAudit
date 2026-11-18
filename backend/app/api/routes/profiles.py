"""
Scoring Profile API Routes

CRUD endpoints for custom scoring profiles with 8 immutable presets.
"""

import logging
from typing import Optional

from fastapi import APIRouter, Depends, HTTPException, Request, Response
from sqlalchemy.ext.asyncio import AsyncSession

from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.core.session import ensure_session_cookie, get_data_scope
from app.db import get_db
from app.schemas.profiles import (
    DuplicateRequest,
    ScoringProfileCreate,
    ScoringProfileExport,
    ScoringProfileResponse,
    ScoringProfileUpdate,
)
from app.services.profiles.service import Owner, ProfileService, _profile_to_response_dict

logger = logging.getLogger(__name__)

router = APIRouter()
profile_service = ProfileService()


async def _owner_for_write(request: Request, response: Response) -> Owner:
    """Caller's scope for a new profile; starts a session cookie when there is none."""
    owner = await get_data_scope(request)
    if owner == (None, None):
        owner = (ensure_session_cookie(response, request), None)
    return owner


@router.get("/profiles", response_model=list[ScoringProfileResponse])
@limiter.limit("60/minute", key_func=get_rate_limit_key)
async def list_profiles(
    request: Request,
    db: AsyncSession = Depends(get_db),
    api_key: Optional[str] = Depends(get_api_key),
):
    """List all active scoring profiles (user + presets)."""
    profiles = await profile_service.list_all(db, owner=await get_data_scope(request))
    return [_profile_to_response_dict(p) for p in profiles]


@router.get("/profiles/{profile_id}", response_model=ScoringProfileResponse)
@limiter.limit("60/minute", key_func=get_rate_limit_key)
async def get_profile(
    request: Request,
    profile_id: int,
    db: AsyncSession = Depends(get_db),
    api_key: Optional[str] = Depends(get_api_key),
):
    """Get a single scoring profile by ID."""
    profile = await profile_service.get(db, profile_id, owner=await get_data_scope(request))
    if profile is None:
        raise HTTPException(status_code=404, detail="Profile not found")
    return _profile_to_response_dict(profile)


@router.post("/profiles", response_model=ScoringProfileResponse, status_code=201)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def create_profile(
    request: Request,
    response: Response,
    body: ScoringProfileCreate,
    db: AsyncSession = Depends(get_db),
    api_key: Optional[str] = Depends(get_api_key),
):
    """Create a new user scoring profile."""
    profile = await profile_service.create(
        db,
        name=body.name,
        description=body.description,
        thresholds=body.thresholds,
        weights=body.weights,
        owner=await _owner_for_write(request, response),
    )
    return _profile_to_response_dict(profile)


@router.put("/profiles/{profile_id}", response_model=ScoringProfileResponse)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def update_profile(
    request: Request,
    profile_id: int,
    body: ScoringProfileUpdate,
    db: AsyncSession = Depends(get_db),
    api_key: Optional[str] = Depends(get_api_key),
):
    """Update a user scoring profile. Cannot update presets (returns 400)."""
    try:
        updates = body.model_dump(exclude_none=True)
        profile = await profile_service.update(
            db, profile_id, updates, owner=await get_data_scope(request)
        )
    except ValueError as e:
        logger.debug("Profile update rejected: %s", e)
        raise HTTPException(status_code=400, detail="Cannot modify preset profiles")

    if profile is None:
        raise HTTPException(status_code=404, detail="Profile not found")
    return _profile_to_response_dict(profile)


@router.delete("/profiles/{profile_id}", status_code=204)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def delete_profile(
    request: Request,
    profile_id: int,
    db: AsyncSession = Depends(get_db),
    api_key: Optional[str] = Depends(get_api_key),
):
    """Soft-delete a user scoring profile. Cannot delete presets (returns 400)."""
    try:
        deleted = await profile_service.delete(db, profile_id, owner=await get_data_scope(request))
    except ValueError as e:
        logger.debug("Profile delete rejected: %s", e)
        raise HTTPException(status_code=400, detail="Cannot delete preset profiles")

    if not deleted:
        raise HTTPException(status_code=404, detail="Profile not found")


@router.post(
    "/profiles/{profile_id}/duplicate",
    response_model=ScoringProfileResponse,
    status_code=201,
)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def duplicate_profile(
    request: Request,
    response: Response,
    profile_id: int,
    body: DuplicateRequest,
    db: AsyncSession = Depends(get_db),
    api_key: Optional[str] = Depends(get_api_key),
):
    """Duplicate any profile (including presets) as a new user profile."""
    profile = await profile_service.duplicate(
        db, profile_id, body.name, owner=await _owner_for_write(request, response)
    )
    if profile is None:
        raise HTTPException(status_code=404, detail="Source profile not found")
    return _profile_to_response_dict(profile)


@router.get("/profiles/{profile_id}/export", response_model=ScoringProfileExport)
@limiter.limit("60/minute", key_func=get_rate_limit_key)
async def export_profile(
    request: Request,
    profile_id: int,
    db: AsyncSession = Depends(get_db),
    api_key: Optional[str] = Depends(get_api_key),
):
    """Export a scoring profile as JSON for file sharing."""
    data = await profile_service.export_json(db, profile_id, owner=await get_data_scope(request))
    if data is None:
        raise HTTPException(status_code=404, detail="Profile not found")
    return data


@router.post(
    "/profiles/import",
    response_model=ScoringProfileResponse,
    status_code=201,
)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def import_profile(
    request: Request,
    response: Response,
    body: ScoringProfileExport,
    db: AsyncSession = Depends(get_db),
    api_key: Optional[str] = Depends(get_api_key),
):
    """Import a scoring profile from JSON."""
    profile = await profile_service.import_json(
        db, body.model_dump(), owner=await _owner_for_write(request, response)
    )
    return _profile_to_response_dict(profile)
