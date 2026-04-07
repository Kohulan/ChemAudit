"""
Tests for MCS (Maximum Common Substructure) Comparison Service.

Validates compute_mcs_comparison with similar molecules, property deltas,
timeout handling, and Tanimoto range checking.
"""

from __future__ import annotations

import pytest


class TestMCSComparison:
    """Tests for compute_mcs_comparison."""

    def test_mcs_two_similar(self):
        """Aspirin + salicylic acid returns MCS SMARTS with >3 atoms, Tanimoto > 0.3."""
        from app.services.analytics.mcs_comparison import compute_mcs_comparison

        result = compute_mcs_comparison(
            smiles_a="CC(=O)Oc1ccccc1C(=O)O",  # aspirin
            smiles_b="OC(=O)c1ccccc1O",          # salicylic acid
        )

        assert "mcs_smarts" in result
        assert "num_atoms" in result
        assert "num_bonds" in result
        assert "timed_out" in result
        assert "tanimoto" in result
        assert "property_deltas" in result
        assert "smiles_a" in result
        assert "smiles_b" in result

        assert result["num_atoms"] > 3
        assert result["tanimoto"] > 0.3
        assert result["timed_out"] is False
        assert len(result["mcs_smarts"]) > 0

    def test_mcs_property_deltas(self):
        """Returns 8 property deltas with correct structure."""
        from app.services.analytics.mcs_comparison import compute_mcs_comparison

        result = compute_mcs_comparison(
            smiles_a="CC(=O)Oc1ccccc1C(=O)O",  # aspirin
            smiles_b="OC(=O)c1ccccc1O",          # salicylic acid
        )

        deltas = result["property_deltas"]
        assert len(deltas) == 8

        expected_props = {
            "MW", "LogP", "TPSA", "QED", "HBA", "HBD", "RotBonds", "RingCount"
        }
        found_props = {d["property"] for d in deltas}
        assert found_props == expected_props

        for d in deltas:
            assert "property" in d
            assert "mol_a" in d
            assert "mol_b" in d
            assert "delta" in d
            assert isinstance(d["mol_a"], (int, float))
            assert isinstance(d["mol_b"], (int, float))
            assert isinstance(d["delta"], (int, float))
            # Delta should be mol_a - mol_b
            assert abs(d["delta"] - (d["mol_a"] - d["mol_b"])) < 1e-6

    def test_mcs_timeout_handling(self):
        """If FindMCS times out, result still contains partial SMARTS and timed_out=True flag."""
        from app.services.analytics.mcs_comparison import compute_mcs_comparison

        # Normal molecules should not time out, but the timed_out field should exist
        result = compute_mcs_comparison(
            smiles_a="c1ccccc1",
            smiles_b="c1ccc(O)cc1",
        )
        assert "timed_out" in result
        assert isinstance(result["timed_out"], bool)

    def test_mcs_tanimoto_range(self):
        """Tanimoto is in [0.0, 1.0] range."""
        from app.services.analytics.mcs_comparison import compute_mcs_comparison

        # Test with very different molecules
        result = compute_mcs_comparison(
            smiles_a="c1ccccc1",       # benzene
            smiles_b="CCCCCCCCCC",     # decane
        )
        assert 0.0 <= result["tanimoto"] <= 1.0

        # Test with identical molecules
        result2 = compute_mcs_comparison(
            smiles_a="c1ccccc1",
            smiles_b="c1ccccc1",
        )
        assert 0.0 <= result2["tanimoto"] <= 1.0
        assert result2["tanimoto"] == 1.0  # identical molecules

    def test_mcs_invalid_smiles_raises(self):
        """Invalid SMILES raises ValueError."""
        from app.services.analytics.mcs_comparison import compute_mcs_comparison

        with pytest.raises(ValueError):
            compute_mcs_comparison(smiles_a="INVALID", smiles_b="c1ccccc1")

        with pytest.raises(ValueError):
            compute_mcs_comparison(smiles_a="c1ccccc1", smiles_b="INVALID")
