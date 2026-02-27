"""
Scoring Profile API Routes

CRUD endpoints for custom scoring profiles with 8 immutable presets.
"""

from typing import Optional

from fastapi import APIRouter, Depends, HTTPException, Request
from sqlalchemy.ext.asyncio import AsyncSession

from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.db import get_db
from app.schemas.profiles import (
    DuplicateRequest,
    ScoringProfileCreate,
    ScoringProfileExport,
    ScoringProfileResponse,
    ScoringProfileUpdate,
)
from app.services.profiles.service import ProfileService, _profile_to_response_dict

router = APIRouter()
profile_service = ProfileService()


@router.get("/profiles", response_model=list[ScoringProfileResponse])
@limiter.limit("60/minute", key_func=get_rate_limit_key)
async def list_profiles(
    request: Request,
    db: AsyncSession = Depends(get_db),
    api_key: Optional[str] = Depends(get_api_key),
):
    """List all active scoring profiles (user + presets)."""
    profiles = await profile_service.list_all(db)
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
    profile = await profile_service.get(db, profile_id)
    if profile is None:
        raise HTTPException(status_code=404, detail="Profile not found")
    return _profile_to_response_dict(profile)


@router.post("/profiles", response_model=ScoringProfileResponse, status_code=201)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def create_profile(
    request: Request,
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
        profile = await profile_service.update(db, profile_id, updates)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

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
        deleted = await profile_service.delete(db, profile_id)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

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
    profile_id: int,
    body: DuplicateRequest,
    db: AsyncSession = Depends(get_db),
    api_key: Optional[str] = Depends(get_api_key),
):
    """Duplicate any profile (including presets) as a new user profile."""
    profile = await profile_service.duplicate(db, profile_id, body.name)
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
    data = await profile_service.export_json(db, profile_id)
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
    body: ScoringProfileExport,
    db: AsyncSession = Depends(get_db),
    api_key: Optional[str] = Depends(get_api_key),
):
    """Import a scoring profile from JSON."""
    profile = await profile_service.import_json(db, body.model_dump())
    return _profile_to_response_dict(profile)
