"""A failed OPSIN init must not be retried (JPype JVM state may be inconsistent)."""

import pytest

import app.services.iupac.converter as conv


@pytest.fixture(autouse=True)
def _reset_opsin_state():
    """Snapshot and restore module-level OPSIN state around each test."""
    saved_nts = conv._nts
    saved_failed = getattr(conv, "_init_failed", False)
    conv._nts = None
    conv._init_failed = False
    yield
    conv._nts = saved_nts
    conv._init_failed = saved_failed


def test_failed_init_sets_flag_and_skips_reinit(monkeypatch):
    monkeypatch.setattr(conv.settings, "OPSIN_JAR_PATH", "/nonexistent/path/opsin.jar")

    # First attempt surfaces the failure so startup can log it.
    with pytest.raises(FileNotFoundError):
        conv.init_opsin()

    assert conv._init_failed is True
    assert conv.is_opsin_available() is False

    # Second attempt must short-circuit — no re-raise, no JVM re-attempt.
    conv.init_opsin()
    assert conv.is_opsin_available() is False
