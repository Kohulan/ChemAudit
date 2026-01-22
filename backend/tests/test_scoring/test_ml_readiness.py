"""
Tests for ML-Readiness Scoring

Tests the MLReadinessScorer for various molecule types.
"""
import pytest
from rdkit import Chem

from app.services.scoring.ml_readiness import (
    MLReadinessScorer,
    calculate_ml_readiness,
    MLReadinessResult,
)


@pytest.fixture
def scorer():
    """Create a scorer instance."""
    return MLReadinessScorer()


class TestMLReadinessScorer:
    """Tests for MLReadinessScorer class."""

    def test_aspirin_scores_high(self, scorer):
        """Aspirin should score well (80+) as it's a simple, well-behaved molecule."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")  # Aspirin
        result = scorer.score(mol)

        assert result.score >= 80, f"Aspirin scored {result.score}, expected 80+"
        assert result.interpretation.startswith("Excellent")
        assert result.breakdown.descriptors_score > 0
        assert result.breakdown.fingerprints_score > 0
        assert result.breakdown.size_score > 0

    def test_caffeine_scores_well(self, scorer):
        """Caffeine should score well as a drug-like molecule."""
        mol = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")  # Caffeine
        result = scorer.score(mol)

        assert result.score >= 70, f"Caffeine scored {result.score}, expected 70+"
        assert "morgan" in result.breakdown.fingerprints_successful
        assert "maccs" in result.breakdown.fingerprints_successful

    def test_ethanol_small_molecule(self, scorer):
        """Ethanol is very small, should have lower size score."""
        mol = Chem.MolFromSmiles("CCO")  # Ethanol
        result = scorer.score(mol)

        # Very small molecule, size should be out of optimal range
        assert result.breakdown.num_atoms == 3
        assert result.breakdown.molecular_weight < 100
        # Size category should be acceptable or out_of_range
        assert result.breakdown.size_category in ["acceptable", "out_of_range"]

    def test_large_molecule_size_penalty(self, scorer):
        """Very large molecules should have size penalties."""
        # Create a large molecule (long alkyl chain)
        large_smiles = "C" * 200  # Very long chain
        mol = Chem.MolFromSmiles(large_smiles)
        result = scorer.score(mol)

        assert result.breakdown.size_category == "out_of_range"
        assert result.breakdown.size_score < 20  # Not optimal

    def test_descriptor_breakdown_populated(self, scorer):
        """Check that descriptor breakdown is properly populated."""
        mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene
        result = scorer.score(mol)

        assert result.breakdown.descriptors_total > 0
        assert result.breakdown.descriptors_successful >= 0
        assert result.breakdown.descriptors_score >= 0
        assert result.breakdown.descriptors_max == 40.0

    def test_fingerprint_breakdown_populated(self, scorer):
        """Check that fingerprint breakdown is properly populated."""
        mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene
        result = scorer.score(mol)

        assert result.breakdown.fingerprints_max == 40.0
        # At least some fingerprints should succeed
        assert len(result.breakdown.fingerprints_successful) > 0
        assert result.breakdown.fingerprints_score > 0

    def test_optimal_size_molecule(self, scorer):
        """Molecules in optimal MW/atom range should get full size points."""
        # Ibuprofen - ~206 Da, good size
        mol = Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
        result = scorer.score(mol)

        assert result.breakdown.size_category == "optimal"
        assert result.breakdown.size_score == 20.0

    def test_acceptable_size_molecule(self, scorer):
        """Molecules in acceptable but not optimal range."""
        # Methane - very small
        mol = Chem.MolFromSmiles("C")
        result = scorer.score(mol)

        assert result.breakdown.size_category in ["acceptable", "out_of_range"]

    def test_interpretation_varies_with_score(self, scorer):
        """Test that interpretation text changes based on score."""
        # High-scoring molecule
        high_mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")  # Aspirin
        high_result = scorer.score(high_mol)

        # Low-scoring (very small)
        small_mol = Chem.MolFromSmiles("C")
        small_result = scorer.score(small_mol)

        # Interpretations should differ
        assert high_result.interpretation != small_result.interpretation


class TestCalculateMLReadinessFunction:
    """Tests for the convenience function."""

    def test_calculate_ml_readiness_returns_result(self):
        """Test that calculate_ml_readiness returns correct type."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_ml_readiness(mol)

        assert isinstance(result, MLReadinessResult)
        assert 0 <= result.score <= 100

    def test_calculate_ml_readiness_aspirin(self):
        """Test ML-readiness for aspirin via convenience function."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_ml_readiness(mol)

        assert result.score >= 80
        assert result.breakdown is not None
        assert result.interpretation is not None


class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_molecule_with_uncommon_elements(self):
        """Test molecules with less common elements."""
        # Molecule with silicon
        mol = Chem.MolFromSmiles("[Si](C)(C)(C)C")
        result = calculate_ml_readiness(mol)

        # Should still produce a valid score
        assert 0 <= result.score <= 100
        assert result.breakdown is not None

    def test_charged_molecule(self):
        """Test charged molecules."""
        # Sodium acetate
        mol = Chem.MolFromSmiles("CC(=O)[O-].[Na+]")
        result = calculate_ml_readiness(mol)

        assert 0 <= result.score <= 100

    def test_molecule_with_stereocenters(self):
        """Test molecules with defined stereochemistry."""
        # L-alanine
        mol = Chem.MolFromSmiles("C[C@H](N)C(=O)O")
        result = calculate_ml_readiness(mol)

        assert 0 <= result.score <= 100
        assert "morgan" in result.breakdown.fingerprints_successful

    def test_aromatic_heterocycle(self):
        """Test aromatic heterocycles."""
        # Pyridine
        mol = Chem.MolFromSmiles("c1ccncc1")
        result = calculate_ml_readiness(mol)

        assert result.score >= 60
        assert result.breakdown.fingerprints_score > 0
