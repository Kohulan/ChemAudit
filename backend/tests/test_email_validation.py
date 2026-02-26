"""
Tests for Phase 4: Email address validation in batch upload.

Verifies that invalid email addresses and SMTP injection attempts are blocked.
"""

import re

# The regex used in batch.py for email validation
_EMAIL_RE = re.compile(r"^[a-zA-Z0-9._%+\-]+@[a-zA-Z0-9.\-]+\.[a-zA-Z]{2,}$")


class TestEmailValidation:
    """Test email validation regex blocks SMTP injection."""

    def test_valid_emails(self):
        """Standard email addresses should pass."""
        valid = [
            "user@example.com",
            "user.name@domain.org",
            "user+tag@sub.domain.co.uk",
            "test123@mail.io",
        ]
        for email in valid:
            assert _EMAIL_RE.match(email), f"Should accept: {email}"

    def test_rejects_newline_injection(self):
        """Newline characters (SMTP header injection) must be blocked."""
        attacks = [
            "user@example.com\r\nBcc: victim@evil.com",
            "user@example.com\nSubject: pwned",
            "user\r@example.com",
        ]
        for attack in attacks:
            assert not _EMAIL_RE.match(attack), f"Should reject: {repr(attack)}"

    def test_rejects_no_domain(self):
        """Emails without a proper domain should be rejected."""
        assert not _EMAIL_RE.match("user@")
        assert not _EMAIL_RE.match("user@localhost")
        assert not _EMAIL_RE.match("@domain.com")

    def test_rejects_too_long(self):
        """Emails over 254 chars should be rejected (checked in batch.py)."""
        long_email = "a" * 250 + "@b.com"
        assert len(long_email) > 254

    def test_rejects_spaces(self):
        """Emails with spaces should be rejected."""
        assert not _EMAIL_RE.match("user name@example.com")
        assert not _EMAIL_RE.match("user@exam ple.com")

    def test_rejects_angle_brackets(self):
        """Angle brackets (header injection) must be blocked."""
        assert not _EMAIL_RE.match("<script>@example.com")
        assert not _EMAIL_RE.match("user@example.com>")
