"""
Tests for API key management functionality.
"""

from unittest.mock import AsyncMock, patch

import pytest
from fastapi.testclient import TestClient

from app.core.config import settings

# Admin secret for testing
ADMIN_SECRET = settings.API_KEY_ADMIN_SECRET


@pytest.fixture
def client():
    """Create test client."""
    from app.main import app

    return TestClient(app)


@pytest.fixture
def admin_headers():
    """Headers with admin authentication."""
    return {"X-Admin-Secret": ADMIN_SECRET}


@pytest.fixture
def mock_redis():
    """Mock Redis client."""
    mock = AsyncMock()
    return mock


def test_create_api_key_success(client, admin_headers, mock_redis):
    """Test creating a new API key."""
    with patch("app.api.routes.api_keys.get_redis_client") as mock_get_redis:
        mock_get_redis.return_value = mock_redis

        response = client.post(
            "/api/v1/api-keys",
            json={"name": "Test Key", "description": "Test API key"},
            headers=admin_headers,
        )

        assert response.status_code == 201
        data = response.json()

        # Check response structure
        assert "key" in data
        assert "name" in data
        assert "created_at" in data
        assert data["name"] == "Test Key"

        # Key should be a secure token
        assert len(data["key"]) > 20
        assert isinstance(data["key"], str)


def test_create_api_key_minimal(client, admin_headers, mock_redis):
    """Test creating API key with minimal data (no description)."""
    with patch("app.api.routes.api_keys.get_redis_client") as mock_get_redis:
        mock_get_redis.return_value = mock_redis

        response = client.post(
            "/api/v1/api-keys",
            json={"name": "Minimal Key"},
            headers=admin_headers,
        )

        assert response.status_code == 201
        data = response.json()
        assert data["name"] == "Minimal Key"


def test_create_api_key_validation_error(client, admin_headers):
    """Test validation errors when creating API key."""
    # Missing required name field
    response = client.post("/api/v1/api-keys", json={}, headers=admin_headers)

    assert response.status_code == 422  # Validation error


def test_list_api_keys(client, admin_headers, mock_redis):
    """Test listing API keys."""
    with patch("app.api.routes.api_keys.get_redis_client") as mock_get_redis:
        mock_redis.smembers.return_value = {"hash1", "hash2"}
        mock_redis.hgetall.side_effect = [
            {
                "name": "Key 1",
                "description": "First key",
                "created_at": "2026-01-23T00:00:00Z",
                "last_used": "2026-01-23T01:00:00Z",
                "request_count": "42",
            },
            {
                "name": "Key 2",
                "description": "",
                "created_at": "2026-01-23T00:00:00Z",
                "last_used": "",
                "request_count": "0",
            },
        ]
        mock_get_redis.return_value = mock_redis

        response = client.get("/api/v1/api-keys", headers=admin_headers)

        assert response.status_code == 200
        data = response.json()

        # Should return list of keys
        assert isinstance(data, list)
        assert len(data) == 2

        # Check first key
        assert data[0]["name"] == "Key 1"
        assert data[0]["request_count"] == 42
        assert "key" not in data[0]  # Full key should NOT be in list response

        # Check second key with no last_used
        assert data[1]["name"] == "Key 2"
        assert data[1]["last_used"] is None


def test_revoke_api_key(client, admin_headers, mock_redis):
    """Test revoking an API key."""
    with patch("app.api.routes.api_keys.get_redis_client") as mock_get_redis:
        mock_redis.smembers.return_value = {"abcd123456789012"}
        mock_get_redis.return_value = mock_redis

        response = client.delete("/api/v1/api-keys/abcd12345678", headers=admin_headers)

        assert response.status_code == 204


def test_revoke_nonexistent_key(client, admin_headers, mock_redis):
    """Test revoking a key that doesn't exist."""
    with patch("app.api.routes.api_keys.get_redis_client") as mock_get_redis:
        mock_redis.smembers.return_value = set()  # No keys
        mock_get_redis.return_value = mock_redis

        response = client.delete("/api/v1/api-keys/nonexistent", headers=admin_headers)

        assert response.status_code == 404
        data = response.json()
        assert "not found" in data["detail"].lower()


def test_validate_api_key_success(client):
    """Test successful API key validation."""
    with patch("app.core.security.validate_api_key") as mock_validate:
        mock_validate.return_value = {
            "name": "Valid Key",
            "created_at": "2026-01-23T00:00:00Z",
        }

        # Use API key in a request
        response = client.get("/api/v1/checks", headers={"X-API-Key": "valid_key"})

        # Should succeed (or fail for other reasons, but not auth)
        assert response.status_code != 401


def test_validate_api_key_failure(client):
    """Test API key validation failure."""
    with patch("app.core.security.validate_api_key") as mock_validate:
        mock_validate.return_value = None  # Invalid key

        response = client.get("/api/v1/checks", headers={"X-API-Key": "invalid_key"})

        assert response.status_code == 401
        data = response.json()
        assert "Invalid" in data["detail"] and "API key" in data["detail"]


def test_hash_api_key_for_lookup():
    """Test API key lookup hash uses PBKDF2 (deterministic KDF)."""
    from app.core.security import hash_api_key_for_lookup

    key1 = "test_key_123"
    key2 = "test_key_123"
    key3 = "different_key"

    # Same input should produce same hash (deterministic)
    hash1 = hash_api_key_for_lookup(key1)
    hash2 = hash_api_key_for_lookup(key2)
    assert hash1 == hash2

    # Different input should produce different hash
    hash3 = hash_api_key_for_lookup(key3)
    assert hash1 != hash3

    # Hash should be 64 hex characters (PBKDF2 dklen=32)
    assert len(hash1) == 64
    assert all(c in "0123456789abcdef" for c in hash1)

    # Must differ from plain SHA256 (KDF output is not the same)
    import hashlib

    plain_sha256 = hashlib.sha256(key1.encode()).hexdigest()
    assert hash1 != plain_sha256


def test_anonymous_access_allowed(client):
    """Test that endpoints work without API key (anonymous access)."""
    # Should work without X-API-Key header
    response = client.get("/api/v1/checks")

    # Should succeed or fail for non-auth reasons
    assert response.status_code != 401


# --- F9/F10: no PBKDF2 or per-key bucket for unvalidated API keys --------------------


def _plain_request(api_key=None, peer="203.0.113.9"):
    from starlette.requests import Request

    headers = [(b"x-api-key", api_key.encode())] if api_key else []
    return Request({"type": "http", "headers": headers, "client": (peer, 1)})


def test_rate_limit_key_ignores_unvalidated_api_key():
    from app.core.rate_limit import get_rate_limit_key

    with patch("app.core.security.hash_api_key_for_lookup") as h:
        assert get_rate_limit_key(_plain_request("random-1")) == "203.0.113.9"
        assert get_rate_limit_key(_plain_request("random-2")) == "203.0.113.9"
    h.assert_not_called()


def test_rate_limit_key_uses_validated_key():
    from app.core import security
    from app.core.rate_limit import get_rate_limit_key

    security._VALIDATED_KEY_HASHES["good-key"] = "abc123"
    try:
        assert get_rate_limit_key(_plain_request("good-key")) == "apikey:abc123"
    finally:
        security._VALIDATED_KEY_HASHES.pop("good-key", None)


@pytest.mark.asyncio
async def test_get_data_scope_does_not_hash_raw_header():
    from app.core.session import get_data_scope

    req = _plain_request("whatever")
    req._cookies = {"chemaudit_sid": "sess"}
    with patch("app.core.security.hash_api_key_for_lookup") as h:
        assert await get_data_scope(req) == ("sess", None)
    h.assert_not_called()


@pytest.mark.asyncio
async def test_get_data_scope_uses_validated_hash():
    from app.core.session import get_data_scope

    req = _plain_request("whatever")
    req.state.api_key_hash = "h1"
    assert await get_data_scope(req) == (None, "h1")


@pytest.mark.asyncio
async def test_oversized_key_rejected_without_hashing():
    from fastapi import HTTPException

    from app.core.security import get_api_key

    with patch("app.core.security.hash_api_key_for_lookup") as h:
        with pytest.raises(HTTPException) as exc:
            await get_api_key(_plain_request(), "k" * 300)
    assert exc.value.status_code == 401
    h.assert_not_called()


@pytest.mark.asyncio
async def test_validate_hashes_off_the_event_loop():
    import threading

    from app.core import security

    seen = {}

    def fake_hash(key):
        seen["thread"] = threading.current_thread()
        return "h"

    redis = AsyncMock()
    redis.hgetall.return_value = {}
    with (
        patch.object(security, "hash_api_key_for_lookup", fake_hash),
        patch.object(security, "get_redis_client", AsyncMock(return_value=redis)),
    ):
        assert await security.validate_api_key("never-seen-key") is None
    assert seen["thread"] is not threading.main_thread()


@pytest.mark.asyncio
async def test_valid_key_sets_scope_and_registers_hash():
    from app.core import security

    req = _plain_request()
    with (
        patch.object(security, "validate_api_key", AsyncMock(return_value={"name": "k"})),
        patch.object(security, "hash_api_key_for_lookup", return_value="h9"),
        patch.object(security, "_update_usage_stats", AsyncMock()),
    ):
        assert await security.get_api_key(req, "valid-key") == "valid-key"
    assert req.state.api_key_hash == "h9"
    assert security._VALIDATED_KEY_HASHES.pop("valid-key") == "h9"


@pytest.mark.asyncio
async def test_unvalidated_keys_throttled_per_ip_before_hashing():
    """Final review #1: decorated routes check limits after dependencies, so
    get_api_key itself must cap bad-key attempts per IP before the slow hash."""
    from fastapi import HTTPException
    from limits.storage import MemoryStorage
    from limits.strategies import FixedWindowRateLimiter

    from app.core import security
    from app.core.rate_limit import limiter

    validate = AsyncMock(return_value=None)
    with (
        patch.object(limiter, "enabled", True),
        patch.object(limiter, "_limiter", FixedWindowRateLimiter(MemoryStorage())),
        patch.object(security, "validate_api_key", validate),
        patch("app.core.rate_limit.record_rate_limit_violation") as violation,
    ):
        codes = []
        for i in range(11):
            with pytest.raises(HTTPException) as exc:
                await security.get_api_key(_plain_request(peer="198.51.100.4"), f"bad-{i}")
            codes.append(exc.value.status_code)
        # Another IP is not affected.
        with pytest.raises(HTTPException) as other:
            await security.get_api_key(_plain_request(peer="198.51.100.5"), "bad-x")
    assert codes == [401] * 10 + [429]
    assert validate.await_count == 11  # 10 from the first IP, 1 from the second
    assert other.value.status_code == 401
    violation.assert_called_once_with("198.51.100.4")


def test_routes_do_not_hash_api_keys_on_the_event_loop():
    """Final review #2: routes must take the validated hash from get_data_scope, not
    rerun the slow PBKDF2 synchronously (only key creation hashes a new key)."""
    from pathlib import Path

    routes = Path(__file__).parents[1] / "app" / "api" / "routes"
    offenders = [
        p.name
        for p in routes.glob("*.py")
        if p.name != "api_keys.py" and "hash_api_key_for_lookup(" in p.read_text()
    ]
    assert offenders == []
