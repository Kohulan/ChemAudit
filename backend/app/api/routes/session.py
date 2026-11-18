"""
Session Data Management Routes

GDPR-compliant data purge endpoint.
"""

from fastapi import APIRouter, Depends, HTTPException, Request, Response
from sqlalchemy import and_, delete, func, select
from sqlalchemy.ext.asyncio import AsyncSession

from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.core.session import SESSION_COOKIE, get_data_scope
from app.db import get_db
from app.db.models.audit import ValidationAuditEntry
from app.db.models.bookmark import Bookmark
from app.db.models.profile import ScoringProfile

router = APIRouter()


@router.delete("/me/data")
@limiter.limit("5/minute", key_func=get_rate_limit_key)
async def purge_my_data(
    request: Request,
    response: Response,
    db: AsyncSession = Depends(get_db),
    api_key: str | None = Depends(get_api_key),
):
    """
    Delete all data associated with the current session or API key.

    GDPR Article 17 — Right to Erasure.
    """
    session_id, api_key_hash = await get_data_scope(request)

    if not session_id and not api_key_hash:
        raise HTTPException(status_code=404, detail="No active session or API key")

    # Build scope filter
    bm_filter = (
        Bookmark.api_key_hash == api_key_hash
        if api_key_hash
        else Bookmark.session_id == session_id
    )
    audit_filter = (
        ValidationAuditEntry.api_key_hash == api_key_hash
        if api_key_hash
        else ValidationAuditEntry.session_id == session_id
    )

    # Count before delete
    bm_count = (
        await db.execute(select(func.count()).select_from(Bookmark).where(bm_filter))
    ).scalar() or 0
    audit_count = (
        await db.execute(
            select(func.count()).select_from(ValidationAuditEntry).where(audit_filter)
        )
    ).scalar() or 0

    profile_filter = and_(
        ScoringProfile.is_preset.is_(False),
        ScoringProfile.api_key_hash == api_key_hash
        if api_key_hash
        else ScoringProfile.session_id == session_id,
    )
    profile_count = (
        await db.execute(select(func.count()).select_from(ScoringProfile).where(profile_filter))
    ).scalar() or 0

    # Delete
    await db.execute(delete(ScoringProfile).where(profile_filter))
    await db.execute(delete(Bookmark).where(bm_filter))
    await db.execute(delete(ValidationAuditEntry).where(audit_filter))
    await db.commit()

    # Clear session cookie
    if session_id:
        response.delete_cookie(SESSION_COOKIE)

    return {
        "status": "purged",
        "deleted": {"bookmarks": bm_count, "history": audit_count, "profiles": profile_count},
    }
