"""SCScore is computed from vendored weights via a self-contained reimplementation.

These tests confirm (a) SCScore is available and numerically correct when the
vendored weights are present, and (b) it degrades gracefully — never crashes —
when the weight file is absent.
"""

import pytest

import app.services.profiler.sa_comparison as sac


@pytest.fixture(autouse=True)
def _reset_scscore_cache(monkeypatch):
    """Reset the module-level weight cache so each test loads fresh."""
    monkeypatch.setattr(sac, "_SCSCORE_LOAD_ATTEMPTED", False)
    monkeypatch.setattr(sac, "_SCSCORE_WEIGHTS", None)


def test_weights_load_as_twelve_arrays():
    weights = sac._load_scscore_weights()
    assert weights is not None, "vendored SCScore weights should be present"
    # 6-layer MLP -> 6 (weight, bias) pairs = 12 arrays.
    assert len(weights) == 12
    assert weights[0].shape == (1024, 300)
    assert weights[-1].shape == (1,)


@pytest.mark.parametrize(
    "smiles,expected_score,expected_class",
    [
        ("CCO", 1.0, "easy"),  # ethanol — trivial
        ("CC(=O)Oc1ccccc1C(=O)O", 1.59, "easy"),  # aspirin
        ("CN1C=NC2=C1C(=O)N(C)C(=O)N2C", 2.79, "moderate"),  # caffeine
        (
            "CCC1OC(=O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(OC2OC(C)CC(N(C)C)C2O)"
            "C(C)(O)CC(C)C(=O)C(C)C(O)C1(C)O",
            4.1,
            "difficult",
        ),  # erythromycin — complex
    ],
)
def test_compute_scscore_matches_reference(smiles, expected_score, expected_class):
    """Scores reproduce the upstream standalone model (validated to <1e-6 raw)."""
    result = sac._compute_scscore(smiles)
    assert result["available"] is True
    assert result["scale"] == "1-5"
    assert result["score"] == expected_score
    assert result["classification"] == expected_class


def test_compute_scscore_invalid_smiles_degrades():
    result = sac._compute_scscore("this-is-not-a-molecule")
    assert result["available"] is False
    assert "error" in result


def test_degrades_gracefully_when_weights_missing(monkeypatch):
    """If the vendored weight file is absent, SCScore reports unavailable, never crashes."""
    monkeypatch.setattr(sac, "_SCSCORE_WEIGHTS_PATH", "/nonexistent/scscore_weights.npz")
    assert sac._load_scscore_weights() is None
    result = sac._compute_scscore("CCO")
    assert result["available"] is False
    assert "error" in result
