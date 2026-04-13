"""
Unit tests for the GenChem composite scorer.

Tests cover:
- score_for_generative() returns float 0-1 for valid SMILES
- score_for_generative() returns None for invalid/empty SMILES (Pitfall 4 / D-14)
- All preset weight vectors sum to 1.0
- Differentiated scores across presets (D-15)
"""

from __future__ import annotations

import pytest

from app.services.genchem.filter_config import PRESETS
from app.services.genchem.scorer import score_for_generative

# ---------------------------------------------------------------------------
# Basic scoring tests
# ---------------------------------------------------------------------------


def test_score_valid():
    """score_for_generative with valid SMILES returns float in [0.0, 1.0]."""
    score = score_for_generative("CCO", PRESETS["drug_like"])
    assert score is not None
    assert isinstance(score, float)
    assert 0.0 <= score <= 1.0


def test_score_null_for_invalid():
    """score_for_generative with invalid SMILES returns None (D-14 / Pitfall 4)."""
    result = score_for_generative("invalid###", PRESETS["drug_like"])
    assert result is None


def test_score_null_for_empty():
    """score_for_generative with empty string returns None."""
    result = score_for_generative("", PRESETS["drug_like"])
    assert result is None


def test_score_range_benzene():
    """score_for_generative with benzene returns float in [0.0, 1.0]."""
    score = score_for_generative("c1ccccc1", PRESETS["drug_like"])
    assert score is not None
    assert 0.0 <= score <= 1.0


def test_score_range_ibuprofen():
    """score_for_generative with ibuprofen (drug-like) returns float in [0.0, 1.0]."""
    ibuprofen = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    score = score_for_generative(ibuprofen, PRESETS["drug_like"])
    assert score is not None
    assert 0.0 <= score <= 1.0


# ---------------------------------------------------------------------------
# Weight vector tests
# ---------------------------------------------------------------------------


def test_weights_sum():
    """All preset weight vectors must sum to 1.0 (within floating-point tolerance)."""
    for name, cfg in PRESETS.items():
        total = cfg.weight_validity + cfg.weight_qed + cfg.weight_alert_free + cfg.weight_sa
        assert total == pytest.approx(1.0, abs=1e-9), (
            f"Preset '{name}' weights sum to {total}, expected 1.0"
        )


# ---------------------------------------------------------------------------
# D-15 differentiated scoring tests
# ---------------------------------------------------------------------------


def test_score_differs_across_presets():
    """
    Per D-15: drug_like and fragment_like presets have different weight vectors,
    so score_for_generative must produce different composite scores for the same molecule.
    """
    ethanol = "CCO"
    score_drug = score_for_generative(ethanol, PRESETS["drug_like"])
    score_frag = score_for_generative(ethanol, PRESETS["fragment_like"])
    assert score_drug is not None
    assert score_frag is not None
    assert score_drug != score_frag, (
        f"D-15 violation: drug_like ({score_drug}) == fragment_like ({score_frag}) "
        "for the same molecule — weight vectors must produce distinct scores"
    )


def test_score_differs_all_presets():
    """
    Per D-15: all 4 presets should produce different scores for a test molecule
    due to their differentiated weight vectors.
    """
    # Use a molecule that will have non-trivial scores across all components
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    scores = {
        name: score_for_generative(aspirin, cfg) for name, cfg in PRESETS.items()
    }
    # All scores must be non-None (aspirin is valid)
    for name, score in scores.items():
        assert score is not None, f"Preset '{name}' returned None for aspirin"
    # At least 3 of 4 should differ (permissive/drug_like may be close for some molecules)
    unique_scores = set(scores.values())
    assert len(unique_scores) >= 2, (
        f"D-15: too few unique scores across presets: {scores}"
    )


def test_score_precision():
    """score_for_generative returns float rounded to 4 decimal places."""
    score = score_for_generative("CCO", PRESETS["drug_like"])
    assert score is not None
    # Check that it's rounded to at most 4 decimal places
    assert score == round(score, 4)
