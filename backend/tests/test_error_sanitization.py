"""
Tests for Phase 4: Error sanitization.

Verifies that internal exception details are not leaked to API consumers.
"""

from app.core.error_sanitizer import safe_error_detail


class TestSafeErrorDetail:
    """Test error sanitization utility."""

    def test_returns_user_message(self):
        """Should return the user_message, not the exception details."""
        error = ValueError("secret internal path /app/data/file.txt")
        result = safe_error_detail(error, "Something went wrong")
        assert result == "Something went wrong"
        assert "/app/data" not in result

    def test_default_message(self):
        """Should use default message when no user_message provided."""
        error = RuntimeError("internal details")
        result = safe_error_detail(error)
        assert result == "An internal error occurred"
        assert "internal details" not in result

    def test_never_returns_exception_str(self):
        """Should never include exception string in output."""
        errors = [
            ValueError("SQL syntax error near 'DROP TABLE'"),
            FileNotFoundError("/app/secrets/key.pem"),
            ConnectionError("redis://admin:password@host:6379"),
            RuntimeError("Traceback (most recent call last): ..."),
        ]
        for error in errors:
            result = safe_error_detail(error, "Error occurred")
            assert str(error) not in result
