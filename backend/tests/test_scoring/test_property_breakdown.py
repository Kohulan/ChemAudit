"""
Tests for property breakdown scoring module.

Tests TPSA per-atom, LogP per-atom, Bertz detail, and Fsp3 detail.
"""

import pytest
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

from app.services.scoring.property_breakdown import (
    PropertyBreakdownScorer,
    calculate_bertz_detail,
    calculate_fsp3_detail,
    calculate_logp_breakdown,
    calculate_tpsa_breakdown,
)


class TestTPSABreakdown:
    """Tests for TPSA per-atom breakdown."""

    def test_tpsa_aspirin_total(self):
        """Sum of per-atom TPSA contributions should match Descriptors.TPSA."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_tpsa_breakdown(mol)
        ref = Descriptors.TPSA(mol)

        assert abs(result.total - ref) < 0.01

    def test_tpsa_atom_count(self):
        """Number of contributions should equal number of heavy atoms."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_tpsa_breakdown(mol)

        assert len(result.atom_contributions) == mol.GetNumAtoms()

    def test_tpsa_nonpolar_zero(self):
        """Carbon atoms should have 0.0 TPSA contribution."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_tpsa_breakdown(mol)

        for ac in result.atom_contributions:
            if ac.symbol == "C":
                assert ac.contribution == 0.0

    def test_tpsa_functional_groups_present(self):
        """At least one functional group should be in the summary."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_tpsa_breakdown(mol)

        assert len(result.functional_group_summary) > 0

    def test_tpsa_ethane_zero(self):
        """Fully nonpolar molecule should have TPSA near 0."""
        mol = Chem.MolFromSmiles("CC")
        result = calculate_tpsa_breakdown(mol)

        assert result.total < 0.1

    def test_tpsa_caffeine(self):
        """Caffeine TPSA breakdown should match reference."""
        mol = Chem.MolFromSmiles("Cn1c(=O)c2c(ncn2C)n(C)c1=O")
        result = calculate_tpsa_breakdown(mol)
        ref = Descriptors.TPSA(mol)

        assert abs(result.total - ref) < 0.01

    def test_tpsa_functional_group_indices(self):
        """Functional group atom indices should be valid."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_tpsa_breakdown(mol)

        for fg in result.functional_group_summary:
            for idx in fg.atom_indices:
                assert 0 <= idx < mol.GetNumAtoms()


class TestLogPBreakdown:
    """Tests for LogP per-atom breakdown."""

    def test_logp_aspirin_total(self):
        """Sum of per-atom LogP contributions should match MolLogP."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_logp_breakdown(mol)
        ref = Crippen.MolLogP(mol)

        assert abs(result.total - ref) < 0.02

    def test_logp_atom_count(self):
        """Contributions should be for heavy atoms only (Hs folded in)."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_logp_breakdown(mol)

        assert len(result.atom_contributions) == mol.GetNumAtoms()

    def test_logp_h_folding(self):
        """H contributions should be folded into parent heavy atoms."""
        mol = Chem.MolFromSmiles("O")  # Water: O with 2 Hs
        result = calculate_logp_breakdown(mol)
        ref = Crippen.MolLogP(mol)

        # Only 1 heavy atom (oxygen), but H contributions folded in
        assert len(result.atom_contributions) == 1
        assert abs(result.total - ref) < 0.02

    def test_logp_functional_groups_present(self):
        """At least one functional group should be in the summary."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_logp_breakdown(mol)

        assert len(result.functional_group_summary) > 0

    def test_logp_benzene(self):
        """Benzene LogP breakdown total should match reference."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = calculate_logp_breakdown(mol)
        ref = Crippen.MolLogP(mol)

        assert abs(result.total - ref) < 0.02

    def test_logp_ethanol(self):
        """Ethanol LogP breakdown total should match reference."""
        mol = Chem.MolFromSmiles("CCO")
        result = calculate_logp_breakdown(mol)
        ref = Crippen.MolLogP(mol)

        assert abs(result.total - ref) < 0.02


class TestBertzDetail:
    """Tests for Bertz complexity detail."""

    def test_bertz_benzene(self):
        """Benzene should have known complexity and ring_complexity > 0."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = calculate_bertz_detail(mol)

        assert result.bertz_ct > 0
        assert result.ring_complexity > 0
        assert result.num_rings == 1
        assert result.num_aromatic_rings == 1

    def test_bertz_methane(self):
        """Methane should have low complexity."""
        mol = Chem.MolFromSmiles("C")
        result = calculate_bertz_detail(mol)

        assert result.bertz_ct >= 0
        assert result.num_bonds == 0
        assert result.num_atoms == 1

    def test_bertz_interpretation(self):
        """Interpretation should match complexity range."""
        # Simple molecule
        mol = Chem.MolFromSmiles("CC")
        result = calculate_bertz_detail(mol)
        assert "simple" in result.interpretation.lower()

        # Complex molecule
        mol = Chem.MolFromSmiles("CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O")
        result = calculate_bertz_detail(mol)
        assert result.bertz_ct > 200

    def test_bertz_ring_complexity_range(self):
        """Ring complexity should be between 0 and 1."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_bertz_detail(mol)

        assert 0 <= result.ring_complexity <= 1

    def test_bertz_atom_count(self):
        """Atom count should match mol.GetNumAtoms()."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = calculate_bertz_detail(mol)

        assert result.num_atoms == mol.GetNumAtoms()


class TestFsp3Detail:
    """Tests for Fsp3 per-carbon detail."""

    def test_fsp3_ethanol(self):
        """Ethanol should have all sp3 carbons -> fsp3 = 1.0."""
        mol = Chem.MolFromSmiles("CCO")
        result = calculate_fsp3_detail(mol)

        assert result.fsp3 == 1.0
        assert result.total_carbons == 2
        assert result.sp3_count == 2
        assert result.sp2_count == 0

    def test_fsp3_benzene(self):
        """Benzene should have all sp2 carbons -> fsp3 = 0.0."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = calculate_fsp3_detail(mol)

        assert result.fsp3 == 0.0
        assert result.total_carbons == 6
        assert result.sp3_count == 0
        assert result.sp2_count == 6

    def test_fsp3_aspirin(self):
        """Aspirin should have mixed -> 0 < fsp3 < 1."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_fsp3_detail(mol)

        assert 0 < result.fsp3 < 1

    def test_fsp3_per_carbon_count(self):
        """per_carbon list length should equal total_carbons."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_fsp3_detail(mol)

        assert len(result.per_carbon) == result.total_carbons

    def test_fsp3_hybridization_values(self):
        """All hybridization values should be in valid set."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_fsp3_detail(mol)

        valid_hyb = {"sp", "sp2", "sp3", "other"}
        for ch in result.per_carbon:
            assert ch.hybridization in valid_hyb

    def test_fsp3_interpretation(self):
        """Interpretation should reflect 3D character."""
        # Flat molecule
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = calculate_fsp3_detail(mol)
        assert "flat" in result.interpretation.lower()

        # 3D molecule
        mol = Chem.MolFromSmiles("CCCCCC")
        result = calculate_fsp3_detail(mol)
        assert "3d" in result.interpretation.lower()

    def test_fsp3_carbon_symbol(self):
        """All per_carbon entries should have symbol 'C'."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_fsp3_detail(mol)

        for ch in result.per_carbon:
            assert ch.symbol == "C"


class TestPropertyBreakdownScorer:
    """Tests for the PropertyBreakdownScorer class."""

    def test_scorer_instantiation(self):
        """Scorer should instantiate."""
        scorer = PropertyBreakdownScorer()
        assert scorer is not None

    def test_scorer_reuse(self):
        """Scorer should be reusable for multiple molecules."""
        scorer = PropertyBreakdownScorer()

        mol1 = Chem.MolFromSmiles("CCO")
        mol2 = Chem.MolFromSmiles("c1ccccc1")

        tpsa1 = scorer.calculate_tpsa_breakdown(mol1)
        tpsa2 = scorer.calculate_tpsa_breakdown(mol2)

        assert tpsa1.total != tpsa2.total
