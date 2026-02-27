"""
Tests for Phase 2: Batch job ownership enforcement.

Verifies that session-based ownership correctly controls access to batch jobs.
"""

from unittest.mock import MagicMock, patch

import pytest
from fastapi import HTTPException

from app.core.ownership import store_job_owner, verify_job_access


class TestStoreJobOwner:
    """Test storing job ownership in Redis."""

    @patch("app.core.ownership._get_redis")
    def test_stores_owner(self, mock_get_redis):
        """Should store session_id for job_id in Redis."""
        mock_redis = MagicMock()
        mock_get_redis.return_value = mock_redis

        store_job_owner("job-123", "session-abc")

        mock_redis.set.assert_called_once()
        args = mock_redis.set.call_args
        assert "job-123" in args[0][0]
        assert args[0][1] == "session-abc"

    @patch("app.core.ownership._get_redis")
    def test_skips_if_no_session(self, mock_get_redis):
        """Should not store anything when session_id is None."""
        store_job_owner("job-123", None)
        mock_get_redis.assert_not_called()

    @patch("app.core.ownership._get_redis")
    def test_graceful_redis_failure(self, mock_get_redis):
        """Should not raise when Redis is unavailable."""
        mock_get_redis.side_effect = Exception("Redis down")
        # Should not raise
        store_job_owner("job-123", "session-abc")


class TestVerifyJobAccess:
    """Test batch job access verification."""

    def _make_request(self, session_id=None, api_key=None):
        """Create a mock Request object."""
        request = MagicMock()
        cookies = {}
        if session_id:
            cookies["chemaudit_sid"] = session_id
        request.cookies = MagicMock()
        request.cookies.get.return_value = session_id
        headers = {}
        if api_key:
            headers["X-API-Key"] = api_key
        request.headers = MagicMock()
        request.headers.get.side_effect = lambda key: headers.get(key)
        return request

    def test_allows_api_key_users(self):
        """API-key authenticated requests should bypass ownership check."""
        request = self._make_request(api_key="test-key")
        # Should not raise regardless of ownership
        verify_job_access(request, "any-job-id")

    def test_allows_no_session_cookie(self):
        """Requests without session cookie should be allowed (graceful degradation)."""
        request = self._make_request()
        verify_job_access(request, "any-job-id")

    @patch("app.core.ownership._get_redis")
    def test_allows_matching_owner(self, mock_get_redis):
        """Should allow when session matches stored owner."""
        mock_redis = MagicMock()
        mock_redis.get.return_value = "session-abc"
        mock_get_redis.return_value = mock_redis

        request = self._make_request(session_id="session-abc")
        verify_job_access(request, "job-123")

    @patch("app.core.ownership._get_redis")
    def test_denies_wrong_session(self, mock_get_redis):
        """Should return 404 when session doesn't match owner."""
        mock_redis = MagicMock()
        mock_redis.get.return_value = "session-abc"
        mock_get_redis.return_value = mock_redis

        request = self._make_request(session_id="session-xyz")
        with pytest.raises(HTTPException) as exc_info:
            verify_job_access(request, "job-123")
        # Uses 404 (not 403) to prevent job ID enumeration
        assert exc_info.value.status_code == 404

    @patch("app.core.ownership._get_redis")
    def test_allows_no_owner_recorded(self, mock_get_redis):
        """Should allow access when no owner is recorded (legacy job)."""
        mock_redis = MagicMock()
        mock_redis.get.return_value = None
        mock_get_redis.return_value = mock_redis

        request = self._make_request(session_id="session-abc")
        verify_job_access(request, "job-123")

    @patch("app.core.ownership._get_redis")
    def test_graceful_redis_failure(self, mock_get_redis):
        """Should allow access when Redis is unavailable."""
        mock_get_redis.side_effect = Exception("Redis down")

        request = self._make_request(session_id="session-abc")
        # Should not raise
        verify_job_access(request, "job-123")
