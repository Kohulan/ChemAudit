"""Admin auth: backward-compatible static secret + opt-in signed, time-bound tokens."""

import hashlib
import hmac
import time

import app.core.security as sec
from app.core.config import settings


def test_static_secret_accepted_by_default():
    assert sec.verify_admin_secret(settings.API_KEY_ADMIN_SECRET) is True


def test_wrong_secret_rejected():
    assert sec.verify_admin_secret("definitely-not-the-secret") is False
    assert sec.verify_admin_secret("") is False


def test_signed_token_round_trip_accepted():
    token = sec.generate_admin_token()
    assert "." in token
    assert sec.verify_admin_secret(token) is True


def test_expired_signed_token_rejected():
    old_ts = str(int(time.time()) - 99999)
    sig = hmac.new(
        settings.API_KEY_ADMIN_SECRET.encode(), old_ts.encode(), hashlib.sha256
    ).hexdigest()
    assert sec.verify_admin_secret(f"{old_ts}.{sig}") is False


def test_tampered_signed_token_rejected():
    ts = str(int(time.time()))
    assert sec.verify_admin_secret(f"{ts}.deadbeef") is False


def test_strict_mode_rejects_static_secret_but_accepts_signed(monkeypatch):
    monkeypatch.setattr(settings, "ADMIN_AUTH_REQUIRE_SIGNED", True)
    # Static secret no longer accepted in strict mode...
    assert sec.verify_admin_secret(settings.API_KEY_ADMIN_SECRET) is False
    # ...but a fresh signed token is.
    assert sec.verify_admin_secret(sec.generate_admin_token()) is True
