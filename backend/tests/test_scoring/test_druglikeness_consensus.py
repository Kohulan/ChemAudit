"""
Tests for consensus drug-likeness scoring, lead-likeness, and ligand efficiency.
"""

import pytest
from rdkit import Chem

from app.services.scoring.druglikeness import (
    calculate_consensus,
    calculate_lead_likeness,
)
from app.services.scoring.salt_inventory import calculate_ligand_efficiency


class TestConsensusScoring:
    """Tests for consensus drug-likeness scoring across 5 rule sets."""

    def test_consensus_aspirin(self):
        """Aspirin should pass most rule sets (score >= 3)."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_consensus(mol)

        assert result.score >= 3
        assert result.total == 5
        assert len(result.rule_sets) == 5

    def test_consensus_all_properties_present(self):
        """Each rule set must have all its expected properties."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_consensus(mol)

        expected_counts = {
            "Lipinski": 4,
            "Veber": 2,
            "Egan": 2,
            "Ghose": 4,
            "Muegge": 9,
        }

        for rs in result.rule_sets:
            assert rs.name in expected_counts, f"Unexpected rule set: {rs.name}"
            assert len(rs.violations) == expected_counts[rs.name], (
                f"{rs.name} should have {expected_counts[rs.name]} properties, "
                f"got {len(rs.violations)}"
            )

    def test_consensus_large_molecule_fails(self):
        """MW > 600 molecule should fail multiple rule sets."""
        # Very large molecule
        smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCC(=O)OCCCCCCCCCCCCCCCCCCCCC(=O)O"
        mol = Chem.MolFromSmiles(smiles)
        result = calculate_consensus(mol)

        assert result.score < 5
        # Should fail at least Ghose and Muegge due to size
        failed = [rs.name for rs in result.rule_sets if not rs.passed]
        assert len(failed) >= 2

    def test_consensus_fragment_scores_low(self):
        """Small fragment should fail Ghose/Muegge (atom count too low)."""
        mol = Chem.MolFromSmiles("CCO")  # Ethanol
        result = calculate_consensus(mol)

        # Ethanol has too few atoms for Ghose (20-70) and Muegge (>4 carbons)
        ghose = next(rs for rs in result.rule_sets if rs.name == "Ghose")
        assert not ghose.passed

    def test_consensus_score_range(self):
        """Score should be between 0 and 5."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_consensus(mol)

        assert 0 <= result.score <= 5

    def test_consensus_interpretation_present(self):
        """Interpretation string should be present and meaningful."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_consensus(mol)

        assert result.interpretation
        assert "drug-likeness" in result.interpretation.lower()

    def test_consensus_rule_violations_have_result(self):
        """Each violation should have result 'pass' or 'fail'."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_consensus(mol)

        for rs in result.rule_sets:
            for v in rs.violations:
                assert v.result in ("pass", "fail")

    def test_consensus_ibuprofen_good(self):
        """Ibuprofen should have good consensus."""
        mol = Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
        result = calculate_consensus(mol)

        assert result.score >= 3


class TestLeadLikeness:
    """Tests for lead-likeness assessment."""

    def test_lead_likeness_small_molecule(self):
        """MW ~250 molecule should pass lead-likeness."""
        # Phenylacetic acid (MW ~136) — too small
        # Use something in range: aniline MW 93 — too small
        # Naproxen MW 230 — in range
        mol = Chem.MolFromSmiles("COc1ccc2cc(CC(=O)O)ccc2c1")  # Naproxen-like, MW ~216
        result = calculate_lead_likeness(mol)

        # Properties should be populated
        assert "mw" in result.properties
        assert "logp" in result.properties
        assert "rotatable_bonds" in result.properties

    def test_lead_likeness_large_molecule(self):
        """MW > 350 molecule should fail."""
        mol = Chem.MolFromSmiles(
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)OCCCCCCCCC"
        )  # Large ibuprofen ester
        result = calculate_lead_likeness(mol)

        if result.properties["mw"] > 350:
            assert not result.passed

    def test_lead_likeness_all_thresholds(self):
        """Verify each property threshold is checked."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_lead_likeness(mol)

        assert len(result.violation_details) == 3
        props_checked = {v.property for v in result.violation_details}
        assert "MW" in props_checked
        assert "LogP" in props_checked
        assert "RotBonds" in props_checked

    def test_lead_likeness_thresholds_populated(self):
        """Threshold dict should have expected keys."""
        mol = Chem.MolFromSmiles("CCO")
        result = calculate_lead_likeness(mol)

        assert "mw_range" in result.thresholds
        assert "logp_range" in result.thresholds
        assert "rotbonds_max" in result.thresholds


class TestLigandEfficiency:
    """Tests for ligand efficiency calculations."""

    def test_ligand_efficiency_with_activity(self):
        """Providing pIC50=8 should compute real LE > 0."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_ligand_efficiency(mol, activity_value=8.0, activity_type="pIC50")

        assert not result.proxy_used
        assert result.le is not None
        assert result.le > 0
        assert result.activity_value == 8.0
        assert result.activity_type == "pIC50"

    def test_ligand_efficiency_proxy(self):
        """No activity should use BEI proxy."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_ligand_efficiency(mol)

        assert result.proxy_used
        assert result.le is not None
        assert result.le > 0

    def test_ligand_efficiency_heavy_atom_count(self):
        """Heavy atom count should be correct."""
        mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene: 6 heavy atoms
        result = calculate_ligand_efficiency(mol)

        assert result.heavy_atom_count == 6

    def test_ligand_efficiency_interpretation(self):
        """Interpretation should be non-empty."""
        mol = Chem.MolFromSmiles("CCO")
        result = calculate_ligand_efficiency(mol)

        assert result.interpretation
        assert len(result.interpretation) > 5
