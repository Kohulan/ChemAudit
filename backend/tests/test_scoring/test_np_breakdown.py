"""
Tests for NP-Likeness Fragment Breakdown

Tests the NPLikenessScorer.calculate_breakdown method and the
calculate_np_breakdown convenience function.
"""

import pytest
from rdkit import Chem

from app.services.scoring.np_likeness import (
    NPBreakdownResult,
    NPFragment,
    NPLikenessScorer,
    calculate_np_breakdown,
)


@pytest.fixture
def scorer():
    """Create a scorer instance."""
    return NPLikenessScorer()


class TestNPBreakdownResult:
    """Tests for NPBreakdownResult dataclass integrity."""

    def test_breakdown_returns_correct_type(self, scorer):
        """calculate_breakdown should return NPBreakdownResult."""
        mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene
        result = scorer.calculate_breakdown(mol)
        assert isinstance(result, NPBreakdownResult)

    def test_breakdown_has_score(self, scorer):
        """Result should include a numeric NP score."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")  # Aspirin
        result = scorer.calculate_breakdown(mol)
        assert isinstance(result.score, float)

    def test_breakdown_has_confidence(self, scorer):
        """Result should include a confidence value between 0 and 1."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = scorer.calculate_breakdown(mol)
        assert 0.0 <= result.confidence <= 1.0

    def test_breakdown_has_interpretation(self, scorer):
        """Result should include interpretation text."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = scorer.calculate_breakdown(mol)
        assert isinstance(result.interpretation, str)
        assert len(result.interpretation) > 0


class TestNPBreakdownFragments:
    """Tests for fragment decomposition in NP breakdown."""

    def test_fragment_count_fields(self, scorer):
        """Fragment counts should be non-negative integers."""
        mol = Chem.MolFromSmiles("CN1CCC23C4=C5C=CC(O)=C4OC2C(O)C=CC3C1C5")  # Morphine
        result = scorer.calculate_breakdown(mol)
        assert isinstance(result.total_fragments, int)
        assert result.total_fragments >= 0
        assert isinstance(result.np_fragment_count, int)
        assert result.np_fragment_count >= 0
        assert isinstance(result.synthetic_fragment_count, int)
        assert result.synthetic_fragment_count >= 0

    def test_fragments_are_np_fragment_type(self, scorer):
        """Each fragment should be an NPFragment instance."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")  # Aspirin
        result = scorer.calculate_breakdown(mol)
        for f in result.fragments:
            assert isinstance(f, NPFragment)

    def test_fragment_classification_values(self, scorer):
        """Fragment classification should be one of three valid values."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = scorer.calculate_breakdown(mol)
        valid_classes = {"np_characteristic", "synthetic_characteristic", "neutral"}
        for f in result.fragments:
            assert f.classification in valid_classes

    def test_fragment_has_bit_id(self, scorer):
        """Each fragment should have a bit_id."""
        mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")  # Naphthalene
        result = scorer.calculate_breakdown(mol)
        for f in result.fragments:
            assert isinstance(f.bit_id, int)

    def test_fragment_sorted_by_abs_contribution(self, scorer):
        """Fragments should be sorted by |contribution| descending."""
        mol = Chem.MolFromSmiles("CN1CCC23C4=C5C=CC(O)=C4OC2C(O)C=CC3C1C5")
        result = scorer.calculate_breakdown(mol)
        if len(result.fragments) >= 2:
            for i in range(len(result.fragments) - 1):
                assert abs(result.fragments[i].contribution) >= abs(
                    result.fragments[i + 1].contribution
                )


class TestNPBreakdownEdgeCases:
    """Tests for edge cases in NP breakdown."""

    def test_invalid_molecule_returns_zero(self, scorer):
        """None mol should return a zero-score result."""
        result = scorer.calculate_breakdown(None)
        assert result.score == 0.0
        assert result.confidence == 0.0

    def test_small_molecule(self, scorer):
        """Very small molecules should still produce a result."""
        mol = Chem.MolFromSmiles("CC")  # Ethane
        result = scorer.calculate_breakdown(mol)
        assert isinstance(result, NPBreakdownResult)
        assert isinstance(result.score, float)

    def test_complex_natural_product(self, scorer):
        """Complex NPs should produce fragment decomposition."""
        # Testosterone
        mol = Chem.MolFromSmiles("CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O")
        result = scorer.calculate_breakdown(mol)
        assert isinstance(result, NPBreakdownResult)
        # Should have at least some fragments
        assert result.total_fragments >= 0

    def test_charged_molecule(self, scorer):
        """Charged molecules should still produce a result."""
        mol = Chem.MolFromSmiles("CC(=O)[O-]")  # Acetate
        result = scorer.calculate_breakdown(mol)
        assert isinstance(result, NPBreakdownResult)


class TestConvenienceFunction:
    """Tests for calculate_np_breakdown module-level function."""

    def test_returns_breakdown_result(self):
        """Convenience function should return NPBreakdownResult."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = calculate_np_breakdown(mol)
        assert isinstance(result, NPBreakdownResult)

    def test_matches_scorer_output(self, scorer):
        """Convenience function should match scorer instance output."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        fn_result = calculate_np_breakdown(mol)
        scorer_result = scorer.calculate_breakdown(mol)
        assert fn_result.score == scorer_result.score
        assert fn_result.total_fragments == scorer_result.total_fragments

    def test_caffeine(self):
        """Test breakdown of caffeine."""
        mol = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        result = calculate_np_breakdown(mol)
        assert isinstance(result.score, float)
        assert isinstance(result.interpretation, str)
