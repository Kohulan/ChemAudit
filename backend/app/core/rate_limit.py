"""
Rate limiting configuration using SlowAPI with Redis storage.

Provides different rate limits for anonymous users vs API key users.
"""
from fastapi import Request
from fastapi.responses import JSONResponse
from slowapi import Limiter
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded

from app.core.config import settings

# Rate limit tiers
RATE_LIMITS = {
    "anonymous": "10/minute",
    "api_key": "300/minute",
}


def get_rate_limit_key(request: Request) -> str:
    """
    Get rate limit key - API key if present, else IP address.

    This allows API key users to get a higher rate limit tier.
    """
    api_key = request.headers.get("X-API-Key")
    if api_key:
        return f"apikey:{api_key}"
    return get_remote_address(request)


# Configure limiter with Redis storage
# IMPORTANT: Use separate Redis DB (db=1) to avoid conflicts with Celery
limiter = Limiter(
    key_func=get_rate_limit_key,
    storage_uri=f"{settings.REDIS_URL}/1",  # Use db=1 for rate limiting
    default_limits=["10/minute"],
    enabled=settings.RATE_LIMIT_ENABLED,
)


def rate_limit_exceeded_handler(request: Request, exc: RateLimitExceeded):
    """
    Custom handler for rate limit exceeded.

    Returns 429 Too Many Requests with Retry-After header.
    """
    return JSONResponse(
        status_code=429,
        content={
            "error": "rate_limit_exceeded",
            "message": f"Rate limit exceeded: {exc.detail}",
            "retry_after": getattr(exc, 'retry_after', 60),
        },
        headers={
            "Retry-After": str(getattr(exc, 'retry_after', 60)),
            "X-RateLimit-Limit": (
                request.state.view_rate_limit
                if hasattr(request.state, 'view_rate_limit')
                else "unknown"
            ),
        }
    )
