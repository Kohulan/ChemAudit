"""
Tests for Webhook and Email Notification Services

Tests HMAC signature computation, webhook payload shape, email template
rendering, and notification opt-in behavior.
"""

import hashlib
import hmac
import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from app.services.notifications.webhook import (
    build_webhook_payload,
    compute_hmac_signature,
)


class TestWebhookHMAC:
    """Tests for HMAC-SHA256 webhook signatures."""

    def test_hmac_signature_correct(self):
        """Verify HMAC is computed correctly for known input."""
        secret = "test-secret-key"
        payload = {"event": "batch.complete", "job_id": "abc123"}
        body = json.dumps(payload, sort_keys=True).encode("utf-8")

        result = compute_hmac_signature(secret, body)

        # Compute expected manually
        expected_sig = hmac.new(
            secret.encode("utf-8"), body, hashlib.sha256
        ).hexdigest()
        assert result == f"sha256={expected_sig}"

    def test_hmac_signature_changes_with_secret(self):
        """Different secrets produce different signatures."""
        payload = json.dumps({"test": "data"}).encode("utf-8")
        sig1 = compute_hmac_signature("secret1", payload)
        sig2 = compute_hmac_signature("secret2", payload)
        assert sig1 != sig2

    def test_hmac_signature_changes_with_payload(self):
        """Different payloads produce different signatures."""
        secret = "shared-secret"
        sig1 = compute_hmac_signature(secret, b'{"a": 1}')
        sig2 = compute_hmac_signature(secret, b'{"b": 2}')
        assert sig1 != sig2


class TestWebhookPayload:
    """Tests for webhook payload shape."""

    def test_payload_has_required_keys(self):
        """Verify payload contains all required fields."""
        payload = build_webhook_payload(
            job_id="test-job-123",
            status="complete",
            molecule_count=100,
            pass_count=90,
            fail_count=10,
            avg_score=85.3,
        )

        required_keys = {
            "event", "job_id", "status", "molecule_count",
            "pass_count", "fail_count", "avg_score", "report_url", "timestamp",
        }
        assert set(payload.keys()) == required_keys

    def test_payload_values(self):
        """Verify payload values are correct."""
        payload = build_webhook_payload(
            job_id="abc",
            status="complete",
            molecule_count=50,
            pass_count=45,
            fail_count=5,
            avg_score=92.789,
        )
        assert payload["event"] == "batch.complete"
        assert payload["job_id"] == "abc"
        assert payload["molecule_count"] == 50
        assert payload["avg_score"] == 92.79  # Rounded to 2 decimals
        assert "/batch/abc" in payload["report_url"]

    def test_webhook_retry_on_failure(self):
        """Verify task raises on HTTP error (enabling Celery retry)."""
        with patch("app.services.notifications.webhook.httpx.Client") as mock_client:
            mock_response = MagicMock()
            mock_response.raise_for_status.side_effect = Exception("500 Server Error")
            mock_client.return_value.__enter__ = MagicMock(return_value=MagicMock())
            mock_client.return_value.__enter__.return_value.post.return_value = mock_response

            # The task function should propagate the exception for Celery to retry
            from app.services.notifications.webhook import send_webhook

            with pytest.raises(Exception, match="500 Server Error"):
                send_webhook(
                    "https://example.com/hook",
                    "secret",
                    {"event": "batch.complete"},
                )


class TestEmailTemplate:
    """Tests for email template rendering."""

    def test_template_renders_with_stats(self):
        """Email template renders with job stats included."""
        from jinja2 import Environment, FileSystemLoader

        template_dir = (
            Path(__file__).parent.parent.parent
            / "app"
            / "templates"
            / "emails"
        )
        env = Environment(loader=FileSystemLoader(str(template_dir)))
        template = env.get_template("batch_complete.html")
        html = template.render(
            job_id="abcdef12-3456-7890",
            base_url="http://localhost:3002",
            molecule_count=100,
            pass_count=90,
            fail_count=10,
            avg_score=85.3,
        )

        assert "abcdef12" in html
        assert "100" in html
        assert "90" in html
        assert "10" in html
        assert "85.3" in html
        assert "View Results" in html

    def test_email_send_mock(self):
        """Email task calls SMTP with correct args."""
        with (
            patch("app.services.notifications.email.smtplib.SMTP") as mock_smtp,
            patch("app.services.notifications.email.settings") as mock_settings,
        ):
            mock_settings.SMTP_HOST = "smtp.example.com"
            mock_settings.SMTP_PORT = 587
            mock_settings.SMTP_TLS = True
            mock_settings.SMTP_USER = "user@example.com"
            mock_settings.SMTP_PASSWORD = "password"
            mock_settings.SMTP_FROM = "noreply@chemaudit.local"
            mock_settings.BASE_URL = "http://localhost:3002"

            mock_server = MagicMock()
            mock_smtp.return_value.__enter__ = MagicMock(return_value=mock_server)
            mock_smtp.return_value.__exit__ = MagicMock(return_value=False)

            from app.services.notifications.email import send_batch_complete_email

            send_batch_complete_email(
                "user@example.com",
                "test-job-id",
                {"molecule_count": 50, "pass_count": 45, "fail_count": 5, "avg_score": 90.0},
            )

            mock_smtp.assert_called_once_with("smtp.example.com", 587)
            mock_server.starttls.assert_called_once()
            mock_server.login.assert_called_once_with("user@example.com", "password")
            mock_server.sendmail.assert_called_once()
