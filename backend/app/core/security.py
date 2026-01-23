"""
API Key authentication and validation.

Provides optional API key authentication with usage tracking.
"""
from fastapi import Depends, HTTPException, status, Security, Request
from fastapi.security import APIKeyHeader
import secrets
import hashlib
from datetime import datetime, timezone
from typing import Optional
import redis.asyncio as redis

from app.core.config import settings

api_key_header = APIKeyHeader(name="X-API-Key", auto_error=False)


async def get_redis_client():
    """Get async Redis client."""
    return redis.from_url(settings.REDIS_URL, decode_responses=True)


def hash_api_key(key: str) -> str:
    """
    Hash API key for storage (one-way).

    Uses SHA256 for secure, one-way hashing.
    """
    return hashlib.sha256(key.encode()).hexdigest()


async def validate_api_key(api_key: str) -> Optional[dict]:
    """
    Validate API key against Redis storage.

    Args:
        api_key: The API key to validate

    Returns:
        Dictionary with key metadata if valid, None if invalid
    """
    client = await get_redis_client()
    try:
        key_hash = hash_api_key(api_key)
        key_data = await client.hgetall(f"apikey:{key_hash}")
        return dict(key_data) if key_data else None
    finally:
        await client.aclose()


async def get_api_key(
    request: Request,
    api_key: Optional[str] = Security(api_key_header)
) -> Optional[str]:
    """
    Validate API key if provided. Returns None for anonymous access.

    This dependency allows both anonymous and authenticated access.
    Invalid API keys result in 401 Unauthorized.

    Args:
        request: FastAPI request object
        api_key: API key from X-API-Key header

    Returns:
        API key if valid, None if not provided

    Raises:
        HTTPException: 401 if API key is invalid
    """
    if api_key is None:
        return None  # Anonymous access allowed

    key_data = await validate_api_key(api_key)
    if key_data is None:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid API key",
            headers={"WWW-Authenticate": "ApiKey"}
        )

    # Update usage stats (async, don't wait)
    import asyncio
    asyncio.create_task(_update_usage_stats(api_key))

    return api_key


async def _update_usage_stats(api_key: str):
    """
    Update API key usage statistics.

    Increments request count and updates last_used timestamp.
    """
    client = await get_redis_client()
    try:
        key_hash = hash_api_key(api_key)
        await client.hset(
            f"apikey:{key_hash}",
            "last_used",
            datetime.now(timezone.utc).isoformat()
        )
        await client.hincrby(f"apikey:{key_hash}", "request_count", 1)
    finally:
        await client.aclose()
