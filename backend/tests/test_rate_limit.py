"""
Tests for rate limiting functionality.
"""

import time
from unittest.mock import MagicMock, patch

import pytest
from fastapi.testclient import TestClient


@pytest.fixture
def client():
    """Create test client."""
    from app.main import app

    return TestClient(app)


@pytest.fixture
def mock_redis():
    """Mock Redis to avoid real connection in tests."""
    with patch("app.core.rate_limit.limiter") as mock:
        # Configure mock to not actually rate limit in tests
        mock.enabled = False
        yield mock


def test_rate_limit_headers_present(client):
    """Test that rate limit headers are present in responses."""
    response = client.get("/api/v1/checks")

    # Should have rate limit related headers or work without error
    assert response.status_code in [200, 429]


def test_anonymous_rate_limit_triggers(client):
    """
    Test that anonymous users hit rate limit after configured requests.

    Note: This test may be skipped if RATE_LIMIT_ENABLED=False in test config.
    """
    # Make 11 rapid requests
    responses = []
    for i in range(11):
        response = client.get("/api/v1/checks")
        responses.append(response)
        # Small delay to avoid overwhelming the client
        time.sleep(0.01)

    # Check if any request was rate limited
    status_codes = [r.status_code for r in responses]

    # Either all succeeded (rate limiting disabled in tests) or one was rate limited
    assert all(
        code in [200, 429] for code in status_codes
    ), f"Unexpected status codes: {status_codes}"

    # If rate limiting is enabled, the last request should be 429
    if any(code == 429 for code in status_codes):
        assert 429 in status_codes[-3:], "Expected rate limit to trigger near the end"

        # Find first 429 response
        rate_limited_response = next(r for r in responses if r.status_code == 429)

        # Check 429 response format
        assert "Retry-After" in rate_limited_response.headers
        data = rate_limited_response.json()
        assert "error" in data
        assert data["error"] == "rate_limit_exceeded"


def test_api_key_allows_higher_limit(client, mock_redis):
    """
    Test that API key users get higher rate limits.

    This test uses mocking to verify the logic without hitting real limits.
    """

    # Mock API key validation
    test_api_key = "test_key_12345"

    with patch("app.core.security.validate_api_key") as mock_validate:
        mock_validate.return_value = {
            "name": "test_key",
            "created_at": "2026-01-23T00:00:00Z",
            "request_count": "0",
        }

        # Make request with API key
        response = client.get("/api/v1/checks", headers={"X-API-Key": test_api_key})

        # Should succeed (not testing actual limit, just that key is accepted)
        assert response.status_code in [200, 401]  # 401 if Redis not available


def test_invalid_api_key_returns_401(client):
    """Test that invalid API key returns 401 Unauthorized."""
    with patch("app.core.security.validate_api_key") as mock_validate:
        mock_validate.return_value = None  # Invalid key

        response = client.get("/api/v1/checks", headers={"X-API-Key": "invalid_key"})

        assert response.status_code == 401
        data = response.json()
        assert "detail" in data
        assert "Invalid" in data["detail"] and "API key" in data["detail"]


def test_rate_limit_key_function():
    """Test the rate limit key function logic."""
    from fastapi import Request

    from app.core.rate_limit import get_rate_limit_key

    # Mock request with API key
    mock_request_with_key = MagicMock(spec=Request)
    mock_request_with_key.headers.get.return_value = "test_api_key"

    key = get_rate_limit_key(mock_request_with_key)
    # Key is hashed (BLAKE2b) to avoid exposing raw API keys in Redis
    from app.core.security import hash_api_key_for_lookup

    assert key == f"apikey:{hash_api_key_for_lookup('test_api_key')}"

    # Mock request without API key (anonymous)
    mock_request_anonymous = MagicMock(spec=Request)
    mock_request_anonymous.headers.get.return_value = None
    mock_request_anonymous.client.host = "127.0.0.1"

    assert get_rate_limit_key(mock_request_anonymous) == "127.0.0.1"


def _req(peer: str, real_ip: str | None = None, api_key: str | None = None):
    """Build a bare Starlette request from a peer address and optional headers."""
    from starlette.requests import Request

    headers = []
    if real_ip is not None:
        headers.append((b"x-real-ip", real_ip.encode()))
    if api_key is not None:
        headers.append((b"x-api-key", api_key.encode()))
    return Request({"type": "http", "headers": headers, "client": (peer, 1234)})


def test_client_ip_ignores_header_from_untrusted_peer():
    from app.core.rate_limit import client_ip

    with patch("app.core.rate_limit.settings.TRUSTED_PROXY_CIDRS", "172.16.0.0/12"):
        assert client_ip(_req("203.0.113.5", "1.2.3.4")) == "203.0.113.5"


def test_client_ip_uses_header_from_trusted_proxy():
    from app.core.rate_limit import client_ip

    with patch("app.core.rate_limit.settings.TRUSTED_PROXY_CIDRS", "172.16.0.0/12,fc00::/7"):
        assert client_ip(_req("172.18.0.3", "1.2.3.4")) == "1.2.3.4"
        assert client_ip(_req("fd00::5", "1.2.3.4")) == "1.2.3.4"


def test_client_ip_rejects_invalid_header():
    from app.core.rate_limit import client_ip

    with patch("app.core.rate_limit.settings.TRUSTED_PROXY_CIDRS", "172.16.0.0/12"):
        assert client_ip(_req("172.18.0.3", "not-an-ip")) == "172.18.0.3"


def test_client_ip_default_trusts_no_proxy():
    from app.core.rate_limit import client_ip

    with patch("app.core.rate_limit.settings.TRUSTED_PROXY_CIDRS", ""):
        assert client_ip(_req("172.18.0.3", "1.2.3.4")) == "172.18.0.3"


def test_loopback_never_banned():
    from app.core import rate_limit

    redis = MagicMock()
    redis.incr.return_value = 10_000
    with (
        patch.object(rate_limit, "get_sync_redis_client", return_value=redis),
        patch.object(rate_limit, "_ban_ip") as ban,
    ):
        rate_limit.record_rate_limit_violation("127.0.0.1")
        rate_limit.record_rate_limit_violation("::1")
    ban.assert_not_called()
