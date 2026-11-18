"""
Rate limiting configuration using SlowAPI with Redis storage.

Provides:
- Different rate limits for anonymous users vs API key users
- IP banning for repeated rate limit violations
"""

import ipaddress
import logging
from datetime import datetime, timezone

import redis
from fastapi import HTTPException, Request
from fastapi.responses import JSONResponse
from limits import parse as parse_limit
from slowapi import Limiter
from slowapi.errors import RateLimitExceeded

from app.core.config import settings

logger = logging.getLogger(__name__)

# Redis key prefixes for banning
BAN_KEY_PREFIX = "ratelimit:ban:"
VIOLATION_KEY_PREFIX = "ratelimit:violations:"


_SYNC_REDIS: "redis.Redis | None" = None


def client_ip(request: Request) -> str:
    """
    Return the address of the real client.

    In production nginx sets X-Real-IP to the connecting address (overwriting any
    value the client sent). The header is trusted only when the direct peer is in
    settings.TRUSTED_PROXY_CIDRS; otherwise the peer address is used, so clients
    cannot spoof their address. Without this every user behind nginx shares the
    nginx container's address for rate limits and IP bans.
    """
    peer = request.client.host if request.client else "unknown"
    cidrs = [c.strip() for c in settings.TRUSTED_PROXY_CIDRS.split(",") if c.strip()]
    if not cidrs:
        return peer
    try:
        peer_addr = ipaddress.ip_address(peer)
        if not any(peer_addr in ipaddress.ip_network(c, strict=False) for c in cidrs):
            return peer
        return str(ipaddress.ip_address(request.headers.get("X-Real-IP", "").strip()))
    except ValueError:
        return peer


def get_sync_redis_client():
    """Return a process-wide pooled synchronous Redis client (lazy singleton)."""
    global _SYNC_REDIS
    if _SYNC_REDIS is None:
        _SYNC_REDIS = redis.from_url(
            settings.REDIS_URL, decode_responses=True, max_connections=20
        )
    return _SYNC_REDIS


def is_ip_banned(ip_address: str) -> bool:
    """
    Check if an IP address is currently banned.

    Args:
        ip_address: The IP address to check

    Returns:
        True if banned, False otherwise
    """
    try:
        client = get_sync_redis_client()
        ban_key = f"{BAN_KEY_PREFIX}{ip_address}"
        return client.exists(ban_key) > 0
    except Exception as e:
        logger.error(f"Error checking IP ban status: {e}")
        return False  # Fail open for availability


def get_ban_info(ip_address: str) -> dict | None:
    """
    Get ban information for an IP address.

    Args:
        ip_address: The IP address to check

    Returns:
        Dictionary with ban info, or None if not banned
    """
    try:
        client = get_sync_redis_client()
        ban_key = f"{BAN_KEY_PREFIX}{ip_address}"
        ban_data = client.hgetall(ban_key)
        if ban_data:
            ttl = client.ttl(ban_key)
            return {
                "banned_at": ban_data.get("banned_at"),
                "reason": ban_data.get("reason"),
                "violations": int(ban_data.get("violations", 0)),
                "remaining_seconds": max(0, ttl),
            }
        return None
    except Exception as e:
        logger.error(f"Error getting ban info: {e}")
        return None


def record_rate_limit_violation(ip_address: str) -> int:
    """
    Record a rate limit violation and potentially ban the IP.

    Uses progressive banning: after threshold violations within the window,
    the IP is banned for the configured duration.

    Args:
        ip_address: The IP address that violated the rate limit

    Returns:
        Current violation count for this IP
    """
    try:
        client = get_sync_redis_client()
        violation_key = f"{VIOLATION_KEY_PREFIX}{ip_address}"

        # Increment violation count
        violations = client.incr(violation_key)

        # Set expiry on first violation (sliding window)
        if violations == 1:
            # Violations expire after 2x the ban duration
            client.expire(
                violation_key, settings.RATE_LIMIT_BAN_DURATION_MINUTES * 60 * 2
            )

        # Check if threshold reached. Never ban loopback: in-process callers (the
        # MCP bridge, health checks) use it, and banning it would lock them all out.
        if violations >= settings.RATE_LIMIT_BAN_THRESHOLD and not _is_loopback(ip_address):
            _ban_ip(ip_address, violations)

        return violations

    except Exception as e:
        logger.error(f"Error recording rate limit violation: {e}")
        return 0


def _is_loopback(ip_address: str) -> bool:
    """True when the address parses as a loopback address."""
    try:
        return ipaddress.ip_address(ip_address).is_loopback
    except ValueError:
        return False


def _ban_ip(ip_address: str, violations: int):
    """
    Ban an IP address for repeated rate limit violations.

    Args:
        ip_address: The IP to ban
        violations: Number of violations that triggered the ban
    """
    try:
        client = get_sync_redis_client()
        ban_key = f"{BAN_KEY_PREFIX}{ip_address}"

        # Store ban information
        client.hset(
            ban_key,
            mapping={
                "banned_at": datetime.now(timezone.utc).isoformat(),
                "reason": "rate_limit_violations",
                "violations": str(violations),
            },
        )

        # Set ban expiry
        client.expire(ban_key, settings.RATE_LIMIT_BAN_DURATION_MINUTES * 60)

        # Clear violation counter
        violation_key = f"{VIOLATION_KEY_PREFIX}{ip_address}"
        client.delete(violation_key)

        logger.warning(
            f"IP {ip_address} banned for {settings.RATE_LIMIT_BAN_DURATION_MINUTES} minutes "
            f"after {violations} rate limit violations"
        )

    except Exception as e:
        logger.error(f"Error banning IP: {e}")


def get_rate_limit_key(request: Request) -> str:
    """
    Get rate limit key - API key if present, else IP address.

    This allows API key users to get a higher rate limit tier.

    Returns:
        Rate limit key string based on API key or IP address.
    """
    # Only keys already validated in this worker get their own bucket. An unknown
    # header value is throttled by client IP, so random keys cannot mint fresh
    # buckets, and nothing here runs the slow key hash.
    from app.core.security import _VALIDATED_KEY_HASHES

    key_hash = _VALIDATED_KEY_HASHES.get(request.headers.get("X-API-Key") or "")
    if key_hash:
        return f"apikey:{key_hash}"
    return client_ip(request)


def check_ip_ban_middleware(request: Request):
    """
    Check if the requesting IP is banned.

    This should be called early in request processing.

    Args:
        request: The FastAPI request

    Raises:
        HTTPException: 403 if IP is banned
    """
    ip_address = client_ip(request)
    if is_ip_banned(ip_address):
        ban_info = get_ban_info(ip_address)
        remaining = ban_info.get("remaining_seconds", 0) if ban_info else 0

        raise HTTPException(
            status_code=403,
            detail={
                "error": "ip_banned",
                "message": "Your IP has been temporarily banned due to repeated rate limit violations",
                "retry_after": remaining,
            },
            headers={"Retry-After": str(remaining)},
        )


# Configure limiter with Redis storage
# IMPORTANT: Use separate Redis DB (db=1) to avoid conflicts with Celery
limiter = Limiter(
    key_func=get_rate_limit_key,
    storage_uri=f"{settings.REDIS_URL}/1",  # Use db=1 for rate limiting
    default_limits=["10/minute"],
    enabled=settings.RATE_LIMIT_ENABLED,
)


# Attempts with an API key this worker has not validated, per client IP. Checked in
# get_api_key before the slow key hash: route limits run only after dependencies,
# so without this a stream of bad keys is never throttled.
_UNVALIDATED_KEY_LIMIT = parse_limit("10/minute")


def allow_unvalidated_key_attempt(request: Request) -> bool:
    """
    Count one unvalidated-key attempt for the client IP; False once over the limit.

    Records a rate-limit violation (towards an IP ban) when the limit is exceeded.
    Fails open when rate limiting is off or its storage is unavailable, like the
    rest of the limiter.
    """
    if not limiter.enabled:
        return True
    ip = client_ip(request)
    try:
        allowed = limiter.limiter.hit(_UNVALIDATED_KEY_LIMIT, "unvalidated-api-key", ip)
    except Exception as e:
        logger.error(f"Unvalidated-key limit check failed: {e}")
        return True
    if not allowed:
        record_rate_limit_violation(ip)
    return allowed


def rate_limit_exceeded_handler(
    request: Request, exc: RateLimitExceeded
) -> JSONResponse:
    """
    Custom handler for rate limit exceeded.

    Records the violation for potential IP banning.
    Returns 429 Too Many Requests with Retry-After header.
    """
    # Record the violation
    ip_address = client_ip(request)
    violations = record_rate_limit_violation(ip_address)

    # Get rate limit info
    rate_limit = "unknown"
    if hasattr(request.state, "view_rate_limit"):
        limit_value = request.state.view_rate_limit
        if isinstance(limit_value, tuple):
            rate_limit = str(limit_value[0])
        else:
            rate_limit = str(limit_value)

    # Check if now banned
    ban_info = get_ban_info(ip_address)
    if ban_info:
        return JSONResponse(
            status_code=403,
            content={
                "error": "ip_banned",
                "message": "Your IP has been temporarily banned due to repeated rate limit violations",
                "retry_after": ban_info.get("remaining_seconds", 60),
            },
            headers={
                "Retry-After": str(ban_info.get("remaining_seconds", 60)),
            },
        )

    return JSONResponse(
        status_code=429,
        content={
            "error": "rate_limit_exceeded",
            "message": f"Rate limit exceeded: {exc.detail}",
            "retry_after": getattr(exc, "retry_after", 60),
            "violations": violations,
            "ban_threshold": settings.RATE_LIMIT_BAN_THRESHOLD,
        },
        headers={
            "Retry-After": str(getattr(exc, "retry_after", 60)),
            "X-RateLimit-Limit": rate_limit,
        },
    )
