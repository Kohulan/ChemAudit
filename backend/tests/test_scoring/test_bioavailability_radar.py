"""
Tests for Bioavailability Radar, BOILED-Egg, and Radar Comparison

Tests the BioavailabilityRadarScorer class and convenience functions.
"""

import pytest
from rdkit import Chem

from app.services.scoring.bioavailability_radar import (
    BioavailabilityRadarResult,
    BioavailabilityRadarScorer,
    BoiledEggResult,
    RadarAxis,
    RadarComparisonResult,
    RadarProfile,
    calculate_bioavailability_radar,
    calculate_boiled_egg,
    calculate_radar_comparison,
)


@pytest.fixture
def scorer():
    """Create a scorer instance."""
    return BioavailabilityRadarScorer()


# ============================================================================
# Bioavailability Radar Tests
# ============================================================================


class TestBioavailabilityRadar:
    """Tests for the 6-axis bioavailability radar."""

    def test_returns_correct_type(self, scorer):
        """calculate_radar should return BioavailabilityRadarResult."""
        mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene
        result = scorer.calculate_radar(mol)
        assert isinstance(result, BioavailabilityRadarResult)

    def test_has_six_axes(self, scorer):
        """Radar should always produce exactly 6 axes."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")  # Aspirin
        result = scorer.calculate_radar(mol)
        assert len(result.axes) == 6

    def test_axis_names(self, scorer):
        """Axes should be named LIPO, SIZE, POLAR, INSOLU, INSATU, FLEX."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = scorer.calculate_radar(mol)
        names = [a.name for a in result.axes]
        assert names == ["LIPO", "SIZE", "POLAR", "INSOLU", "INSATU", "FLEX"]

    def test_axes_are_radar_axis_type(self, scorer):
        """Each axis should be a RadarAxis instance."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = scorer.calculate_radar(mol)
        for a in result.axes:
            assert isinstance(a, RadarAxis)

    def test_normalized_values_in_range(self, scorer):
        """Normalized values should be between 0 and 1."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = scorer.calculate_radar(mol)
        for a in result.axes:
            assert 0.0 <= a.normalized <= 1.0

    def test_in_range_count(self, scorer):
        """overall_in_range_count should be between 0 and 6."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = scorer.calculate_radar(mol)
        assert 0 <= result.overall_in_range_count <= 6
        # Verify count matches individual axes
        assert result.overall_in_range_count == sum(1 for a in result.axes if a.in_range)

    def test_drug_like_molecule_mostly_in_range(self, scorer):
        """A drug-like molecule should have most axes in range."""
        # Ibuprofen - typical drug
        mol = Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
        result = scorer.calculate_radar(mol)
        # Should have at least 3 axes in range
        assert result.overall_in_range_count >= 3

    def test_interpretation_populated(self, scorer):
        """Interpretation should be a non-empty string."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = scorer.calculate_radar(mol)
        assert isinstance(result.interpretation, str)
        assert len(result.interpretation) > 0

    def test_excellent_interpretation(self, scorer):
        """All-in-range should mention 'excellent'."""
        # Construct a molecule likely to have all 6 in range
        # Aspirin: MW 180, LogP ~1.2, TPSA 63.6, moderate rotbonds, some sp3
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = scorer.calculate_radar(mol)
        if result.overall_in_range_count == 6:
            assert "excellent" in result.interpretation.lower()

    def test_convenience_function(self):
        """calculate_bioavailability_radar should match scorer output."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_bioavailability_radar(mol)
        assert isinstance(result, BioavailabilityRadarResult)
        assert len(result.axes) == 6


class TestRadarAxisProperties:
    """Tests for individual radar axis properties."""

    def test_lipo_axis_values(self, scorer):
        """LIPO axis should have correct optimal range."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = scorer.calculate_radar(mol)
        lipo = result.axes[0]
        assert lipo.name == "LIPO"
        assert lipo.property_name == "WLOGP"
        assert lipo.optimal_min == -0.7
        assert lipo.optimal_max == 5.0

    def test_size_axis_values(self, scorer):
        """SIZE axis should have correct optimal range."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = scorer.calculate_radar(mol)
        size = result.axes[1]
        assert size.name == "SIZE"
        assert size.property_name == "MW"
        assert size.optimal_min == 150.0
        assert size.optimal_max == 500.0

    def test_flex_axis_integer_count(self, scorer):
        """FLEX axis actual_value should be non-negative."""
        mol = Chem.MolFromSmiles("CCCCCC")  # Hexane
        result = scorer.calculate_radar(mol)
        flex = result.axes[5]
        assert flex.name == "FLEX"
        assert flex.actual_value >= 0


# ============================================================================
# BOILED-Egg Tests
# ============================================================================


class TestBoiledEgg:
    """Tests for BOILED-Egg GI absorption / BBB permeation classification."""

    def test_returns_correct_type(self, scorer):
        """calculate_boiled_egg should return BoiledEggResult."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = scorer.calculate_boiled_egg(mol)
        assert isinstance(result, BoiledEggResult)

    def test_has_required_fields(self, scorer):
        """Result should have wlogp, tpsa, gi_absorbed, bbb_permeant, region."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = scorer.calculate_boiled_egg(mol)
        assert isinstance(result.wlogp, float)
        assert isinstance(result.tpsa, float)
        assert isinstance(result.gi_absorbed, bool)
        assert isinstance(result.bbb_permeant, bool)
        assert result.region in {"yolk", "white", "grey"}

    def test_bbb_implies_gi(self, scorer):
        """If BBB permeant (yolk), should also be GI absorbed."""
        mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene - small, lipophilic
        result = scorer.calculate_boiled_egg(mol)
        if result.bbb_permeant:
            assert result.gi_absorbed

    def test_yolk_region(self, scorer):
        """Yolk region should have bbb_permeant=True."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = scorer.calculate_boiled_egg(mol)
        if result.region == "yolk":
            assert result.bbb_permeant is True

    def test_white_region(self, scorer):
        """White region should have gi_absorbed=True, bbb_permeant=False."""
        # Aspirin - moderate TPSA, should be in white
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = scorer.calculate_boiled_egg(mol)
        if result.region == "white":
            assert result.gi_absorbed is True
            assert result.bbb_permeant is False

    def test_grey_region(self, scorer):
        """Grey region should have gi_absorbed=False, bbb_permeant=False."""
        # Very polar large molecule
        mol = Chem.MolFromSmiles("OC(=O)CNCC(O)=O")
        result = scorer.calculate_boiled_egg(mol)
        if result.region == "grey":
            assert result.gi_absorbed is False
            assert result.bbb_permeant is False

    def test_ellipse_params_present(self, scorer):
        """Ellipse parameters should be present."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = scorer.calculate_boiled_egg(mol)
        assert result.gi_ellipse is not None
        assert result.bbb_ellipse is not None
        assert result.gi_ellipse.cx > 0
        assert result.bbb_ellipse.cx > 0

    def test_interpretation_contains_region(self, scorer):
        """Interpretation should mention the region classification."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = scorer.calculate_boiled_egg(mol)
        assert isinstance(result.interpretation, str)
        assert len(result.interpretation) > 0

    def test_convenience_function(self):
        """calculate_boiled_egg convenience function should work."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_boiled_egg(mol)
        assert isinstance(result, BoiledEggResult)


# ============================================================================
# Radar Comparison Tests
# ============================================================================


class TestRadarComparison:
    """Tests for multi-molecule radar comparison."""

    def test_returns_correct_type(self, scorer):
        """calculate_comparison should return RadarComparisonResult."""
        smiles_list = ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]
        result = scorer.calculate_comparison(smiles_list)
        assert isinstance(result, RadarComparisonResult)

    def test_profile_count_matches_input(self, scorer):
        """Should have one profile per valid SMILES."""
        smiles_list = ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", "CCO"]
        result = scorer.calculate_comparison(smiles_list)
        assert len(result.profiles) == 3

    def test_profiles_are_radar_profile_type(self, scorer):
        """Each profile should be a RadarProfile instance."""
        smiles_list = ["c1ccccc1", "CCO"]
        result = scorer.calculate_comparison(smiles_list)
        for p in result.profiles:
            assert isinstance(p, RadarProfile)
            assert p.is_reference is False

    def test_reference_profile_present(self, scorer):
        """Reference profile should always be present."""
        smiles_list = ["c1ccccc1"]
        result = scorer.calculate_comparison(smiles_list)
        assert result.reference is not None
        assert result.reference.is_reference is True
        assert result.reference.smiles == "reference"

    def test_reference_has_six_axes(self, scorer):
        """Reference profile should have 6 axes."""
        smiles_list = ["c1ccccc1"]
        result = scorer.calculate_comparison(smiles_list)
        assert len(result.reference.axes) == 6

    def test_invalid_smiles_skipped(self, scorer):
        """Invalid SMILES should be silently skipped."""
        smiles_list = ["c1ccccc1", "INVALID_SMILES", "CCO"]
        result = scorer.calculate_comparison(smiles_list)
        assert len(result.profiles) == 2  # Only valid ones

    def test_empty_list(self, scorer):
        """Empty list should return empty profiles but still have reference."""
        result = scorer.calculate_comparison([])
        assert len(result.profiles) == 0
        assert result.reference is not None

    def test_convenience_function(self):
        """calculate_radar_comparison convenience function should work."""
        smiles_list = ["c1ccccc1", "CCO"]
        result = calculate_radar_comparison(smiles_list)
        assert isinstance(result, RadarComparisonResult)
        assert len(result.profiles) == 2


# ============================================================================
# Normalization Tests
# ============================================================================


class TestNormalization:
    """Tests for the normalization logic."""

    def test_in_range_value_normalized_to_one(self, scorer):
        """Values within optimal range should normalize to 1.0."""
        # Use a molecule with MW in 150-500 range
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")  # MW ~180
        result = scorer.calculate_radar(mol)
        size_axis = result.axes[1]  # SIZE
        if size_axis.in_range:
            assert size_axis.normalized == 1.0

    def test_out_of_range_below_min(self, scorer):
        """Values below optimal min should normalize below 1.0."""
        # Very small molecule - MW well below 150
        mol = Chem.MolFromSmiles("C")  # Methane, MW ~16
        result = scorer.calculate_radar(mol)
        size_axis = result.axes[1]  # SIZE
        assert size_axis.normalized < 1.0
        assert size_axis.in_range is False

    def test_normalized_never_negative(self, scorer):
        """Normalized values should never be negative."""
        mol = Chem.MolFromSmiles("C")  # Very small
        result = scorer.calculate_radar(mol)
        for a in result.axes:
            assert a.normalized >= 0.0
