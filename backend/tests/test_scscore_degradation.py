"""SCScore is optional — its absence must degrade gracefully, never crash."""

import app.services.profiler.sa_comparison as sac


def test_load_scscore_returns_none_when_unavailable(monkeypatch):
    monkeypatch.setattr(sac, "_SCSCORE_LOAD_ATTEMPTED", False)
    monkeypatch.setattr(sac, "_SCSCORE_SCORER", None)
    # scscore is not a PyPI dependency; in its absence the loader returns None.
    assert sac._load_scscore() is None


def test_compute_scscore_reports_unavailable(monkeypatch):
    monkeypatch.setattr(sac, "_SCSCORE_LOAD_ATTEMPTED", False)
    monkeypatch.setattr(sac, "_SCSCORE_SCORER", None)
    result = sac._compute_scscore("CCO")
    assert result["available"] is False
    assert "error" in result


def test_compute_scscore_uses_loaded_scorer(monkeypatch):
    """When a scorer is available, _compute_scscore returns a real score."""

    class FakeScorer:
        def get_score_from_smi(self, smi):
            return smi, 2.5

        def apply(self, smi):
            return smi, 2.5

    monkeypatch.setattr(sac, "_SCSCORE_LOAD_ATTEMPTED", True)
    monkeypatch.setattr(sac, "_SCSCORE_SCORER", FakeScorer())
    result = sac._compute_scscore("CCO")
    assert result["available"] is True
    assert result["score"] == 2.5
    assert result["classification"] == "moderate"  # 2 <= 2.5 < 3
