"""
Batch Job Ownership

Session-based ownership for batch jobs. Stores the session ID of the job creator
in Redis and verifies access on subsequent requests. Returns 404 (not 403) to
prevent job ID enumeration.
"""

import logging
from typing import Optional

from fastapi import HTTPException, Request

from app.core.config import settings
from app.core.session import get_session_id

logger = logging.getLogger(__name__)

_OWNER_KEY_PREFIX = "batch:owner:"


def _get_redis():
    """Get synchronous Redis client."""
    import redis
    return redis.from_url(settings.REDIS_URL, decode_responses=True)


def store_job_owner(job_id: str, session_id: Optional[str]) -> None:
    """Store the session that created a batch job."""
    if not session_id:
        return
    try:
        r = _get_redis()
        r.set(
            f"{_OWNER_KEY_PREFIX}{job_id}",
            session_id,
            ex=settings.BATCH_RESULT_TTL,
        )
    except Exception:
        logger.warning("Failed to store job owner for %s", job_id)


def verify_job_access(request: Request, job_id: str) -> None:
    """
    Verify the requester owns the batch job.

    Skips check for API-key authenticated requests (they use a different
    auth model). Returns 404 (not 403) to prevent job ID enumeration.
    Degrades gracefully if Redis is unavailable or no session cookie exists.
    """
    # API-key users bypass session ownership
    if request.headers.get("X-API-Key"):
        return

    session_id = get_session_id(request)
    if not session_id:
        # No session cookie — degrade gracefully (e.g. first visit, curl)
        return

    try:
        r = _get_redis()
        owner = r.get(f"{_OWNER_KEY_PREFIX}{job_id}")
    except Exception:
        # Redis down — degrade gracefully
        return

    if owner is None:
        # No owner recorded (legacy job or Redis expired) — allow access
        return

    if owner != session_id:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
