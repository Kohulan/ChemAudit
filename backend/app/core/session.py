"""
Session Management for Data Isolation

Provides per-browser session scoping via HttpOnly cookie.
API key users are scoped by api_key_hash instead.
"""

import uuid
from typing import Optional

from fastapi import Request, Response
from sqlalchemy import text
from sqlalchemy.ext.asyncio import AsyncSession

from app.core.config import settings
from app.core.security import hash_api_key_for_lookup

SESSION_COOKIE = "chemaudit_sid"
SESSION_MAX_AGE = 30 * 24 * 3600  # 30 days


def get_session_id(request: Request) -> Optional[str]:
    """Read session ID from cookie. Returns None if no cookie."""
    return request.cookies.get(SESSION_COOKIE)


def create_session_id() -> str:
    """Generate a new cryptographically random session ID."""
    return str(uuid.uuid4())


def ensure_session_cookie(response: Response, session_id: str) -> None:
    """Set session cookie on response if not already present."""
    is_dev = (
        "localhost" in settings.CORS_ORIGINS_STR
        or "127.0.0.1" in settings.CORS_ORIGINS_STR
    )
    response.set_cookie(
        key=SESSION_COOKIE,
        value=session_id,
        max_age=SESSION_MAX_AGE,
        httponly=True,
        secure=not is_dev,
        samesite="lax",
    )


async def get_data_scope(request: Request) -> tuple[Optional[str], Optional[str]]:
    """
    Determine data scope from request.

    Returns:
        (session_id, api_key_hash) — API key takes precedence.
    """
    api_key = request.headers.get("X-API-Key")
    if api_key:
        api_key_hash = hash_api_key_for_lookup(api_key)
        return None, api_key_hash

    session_id = get_session_id(request)
    return session_id, None


async def set_rls_context(
    db: AsyncSession, session_id: Optional[str], api_key_hash: Optional[str]
) -> None:
    """Set PostgreSQL session variables for RLS policy evaluation.

    Uses set_config() instead of SET LOCAL because SET doesn't support
    bind parameters in PostgreSQL (asyncpg translates :param to $1).
    The third argument (true) makes it transaction-local, same as SET LOCAL.
    """
    await db.execute(
        text("SELECT set_config('app.session_id', :sid, true)"),
        {"sid": session_id or ""},
    )
    await db.execute(
        text("SELECT set_config('app.api_key_hash', :akh, true)"),
        {"akh": api_key_hash or ""},
    )
