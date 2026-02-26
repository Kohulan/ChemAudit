"""
Tests for ML-Readiness Scoring (4-dimension redesign)

Tests the MLReadinessScorer across structural quality, property profile,
complexity & feasibility, and representation quality dimensions.
Includes unit tests for scoring math, boundary conditions, and caveats.
"""

import pytest
from rdkit import Chem

from app.services.scoring.ml_readiness import (
    MLReadinessResult,
    MLReadinessScorer,
    calculate_ml_readiness,
)


@pytest.fixture
def scorer():
    """Create a scorer instance."""
    return MLReadinessScorer()


# =============================================================================
# Top-level smoke tests
# =============================================================================


class TestMLReadinessScorer:
    """High-level scorer tests."""

    def test_aspirin_scores_high(self, scorer):
        """Aspirin should score well as a simple, well-behaved molecule."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = scorer.score(mol)

        assert result.score >= 70, f"Aspirin scored {result.score}, expected 70+"
        assert result.label in ("Excellent", "Good")
        assert result.breakdown.structural_quality.score > 0
        assert result.breakdown.property_profile.score > 0
        assert result.breakdown.representation_quality.score > 0

    def test_caffeine_scores_well(self, scorer):
        """Caffeine should score well as a drug-like molecule."""
        mol = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        result = scorer.score(mol)

        assert result.score >= 50, f"Caffeine scored {result.score}, expected 50+"
        for dim_name in ["structural_quality", "property_profile",
                         "complexity_feasibility", "representation_quality"]:
            dim = getattr(result.breakdown, dim_name)
            assert dim.score >= 0

    def test_ethanol_lower_than_aspirin(self, scorer):
        """Ethanol (tiny, non-drug-like) should score lower than aspirin."""
        aspirin = scorer.score(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"))
        ethanol = scorer.score(Chem.MolFromSmiles("CCO"))

        assert ethanol.score < aspirin.score

    def test_label_and_color_always_present(self, scorer):
        """Label and color should always be populated."""
        for smi in ["c1ccccc1", "C", "CC(=O)Oc1ccccc1C(=O)O", "[Na+].[Cl-]"]:
            result = scorer.score(Chem.MolFromSmiles(smi))
            assert result.label in ("Excellent", "Good", "Moderate", "Limited", "Poor")
            assert result.color in ("green", "teal", "amber", "orange", "red")


# =============================================================================
# Dimension 1: Structural Quality (20 pts)
# =============================================================================


class TestStructuralQuality:
    """Tests for binary structural checks."""

    def test_clean_organic_full_score(self, scorer):
        """Clean organic molecule → 20/20."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        sq = scorer.score(mol).breakdown.structural_quality

        assert sq.score == 20.0
        assert sq.max_score == 20.0

    def test_five_items_present(self, scorer):
        """Should always have 5 check items."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        sq = scorer.score(mol).breakdown.structural_quality

        assert len(sq.items) == 5
        names = {i.name for i in sq.items}
        assert "Single component" in names
        assert "Standard organic elements" in names
        assert "No radicals" in names
        assert "Reasonable charge" in names
        assert "No dummy atoms" in names

    def test_all_items_passed_for_clean_mol(self, scorer):
        """All 5 items should pass for benzene."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        sq = scorer.score(mol).breakdown.structural_quality

        for item in sq.items:
            assert item.passed is True, f"{item.name} should pass for benzene"

    def test_mixture_loses_5pts(self, scorer):
        """Mixture → loses single-component 5pts."""
        mol = Chem.MolFromSmiles("CC.CC")
        sq = scorer.score(mol).breakdown.structural_quality

        single_comp = next(i for i in sq.items if i.name == "Single component")
        assert single_comp.passed is False
        assert single_comp.score == 0.0
        assert sq.score == 15.0  # 20 - 5

    def test_inorganic_loses_5pts(self, scorer):
        """Ferrocene-like organometallic → loses organic-elements 5pts."""
        mol = Chem.MolFromSmiles("[Fe]")
        if mol is not None:
            sq = scorer.score(mol).breakdown.structural_quality
            organic = next(i for i in sq.items if i.name == "Standard organic elements")
            assert organic.passed is False
            assert organic.score == 0.0

    def test_radical_loses_3pts(self, scorer):
        """Radical species → loses 3pts."""
        mol = Chem.MolFromSmiles("[CH2]C")
        if mol is not None:
            sq = scorer.score(mol).breakdown.structural_quality
            radical = next(i for i in sq.items if i.name == "No radicals")
            assert radical.passed is False
            assert radical.score == 0.0

    def test_high_charge_loses_3pts(self, scorer):
        """Net charge > 2 → loses charge 3pts."""
        # Triply charged cation
        mol = Chem.MolFromSmiles("[NH4+].[NH4+].[NH4+]")
        if mol is not None:
            sq = scorer.score(mol).breakdown.structural_quality
            charge = next(i for i in sq.items if i.name == "Reasonable charge")
            # |net| = 3 > 2, so should fail
            assert charge.passed is False

    def test_reasonable_charge_passes(self, scorer):
        """Net charge <=2 should pass."""
        mol = Chem.MolFromSmiles("[NH4+]")  # |net| = 1
        sq = scorer.score(mol).breakdown.structural_quality
        charge = next(i for i in sq.items if i.name == "Reasonable charge")
        assert charge.passed is True

    def test_dummy_atom_loses_4pts(self, scorer):
        """Molecule with dummy atom (*) → loses 4pts."""
        mol = Chem.MolFromSmiles("*CC*")
        if mol is not None:
            sq = scorer.score(mol).breakdown.structural_quality
            dummy = next(i for i in sq.items if i.name == "No dummy atoms")
            assert dummy.passed is False
            assert dummy.score == 0.0


# =============================================================================
# Dimension 2: Property Profile (35 pts)
# =============================================================================


class TestPropertyProfile:
    """Tests for desirability-scored property dimension."""

    def test_seven_items(self, scorer):
        """Should have exactly 7 property items."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        pp = scorer.score(mol).breakdown.property_profile
        assert len(pp.items) == 7

    def test_item_names(self, scorer):
        """Check expected item names."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        pp = scorer.score(mol).breakdown.property_profile
        names = [i.name for i in pp.items]
        assert "MW fitness" in names
        assert "LogP fitness" in names
        assert "TPSA fitness" in names
        assert "HBD fitness" in names
        assert "HBA fitness" in names
        assert "RotBonds fitness" in names
        assert "AromRings fitness" in names

    def test_ibuprofen_high_property_score(self, scorer):
        """Ibuprofen (MW 206, LogP ~3.5) → high property profile."""
        mol = Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
        pp = scorer.score(mol).breakdown.property_profile
        assert pp.score > 20  # > 57% of 35

    def test_methane_lower_than_drug(self, scorer):
        """Methane should score lower than a drug-like molecule on properties."""
        methane = scorer.score(Chem.MolFromSmiles("C")).breakdown.property_profile
        ibuprofen = scorer.score(
            Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
        ).breakdown.property_profile
        assert methane.score < ibuprofen.score

    def test_details_contain_all_props(self, scorer):
        """Details dict should have all 7 property values."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        details = scorer.score(mol).breakdown.property_profile.details
        for key in ["mw", "logp", "tpsa", "hbd", "hba", "rotatable_bonds", "aromatic_rings"]:
            assert key in details

    def test_mw_in_range_gets_full_points(self, scorer):
        """MW 200-500 should get ~6/6 for MW item."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")  # MW ~180
        pp = scorer.score(mol).breakdown.property_profile
        mw_item = next(i for i in pp.items if i.name == "MW fitness")
        # Aspirin MW ~180 is slightly below 200, so partial credit
        assert mw_item.max_score == 6.0

    def test_very_large_mw_low_score(self, scorer):
        """Very high MW → low MW item score."""
        # Long chain → high MW
        mol = Chem.MolFromSmiles("C" * 100)
        pp = scorer.score(mol).breakdown.property_profile
        mw_item = next(i for i in pp.items if i.name == "MW fitness")
        assert mw_item.score < mw_item.max_score


# =============================================================================
# Dimension 3: Complexity & Feasibility (25 pts)
# =============================================================================


class TestComplexityFeasibility:
    """Tests for QED, SA Score, Fsp3, stereocenters."""

    def test_four_items(self, scorer):
        """Should have QED, SA Score, Fsp3, Stereocenters."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        cf = scorer.score(mol).breakdown.complexity_feasibility
        assert len(cf.items) == 4
        names = [i.name for i in cf.items]
        assert "QED contribution" in names
        assert "Synthesizability" in names
        assert "Fsp3 fitness" in names
        assert "Stereo. manageability" in names

    def test_details_populated(self, scorer):
        """Details should contain QED, SA score, Fsp3, stereo counts."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        details = scorer.score(mol).breakdown.complexity_feasibility.details
        assert 0 <= details["qed"] <= 1
        assert 1 <= details["sa_score"] <= 10
        assert 0 <= details["fsp3"] <= 1
        assert details["total_stereocenters"] >= 0
        assert details["undefined_stereocenters"] >= 0

    def test_defined_stereo_gets_full_points(self, scorer):
        """Molecule with 2 defined stereocenters → full stereo points."""
        mol = Chem.MolFromSmiles("[C@@H](O)(F)Cl")  # 1 defined center
        cf = scorer.score(mol).breakdown.complexity_feasibility
        stereo = next(i for i in cf.items if i.name == "Stereo. manageability")
        # 1 center, 0 undefined → 5 pts
        assert stereo.score == 5.0

    def test_no_stereocenters_full_points(self, scorer):
        """Zero stereocenters → full 5 pts."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        cf = scorer.score(mol).breakdown.complexity_feasibility
        stereo = next(i for i in cf.items if i.name == "Stereo. manageability")
        assert stereo.score == 5.0


class TestSAScoreToPoints:
    """Unit tests for _sa_score_to_points static method."""

    @pytest.mark.parametrize("sa,expected", [
        (1.0, 8.0),   # Very easy
        (2.0, 8.0),   # Easy
        (3.0, 8.0),   # Boundary of easy
        (4.0, 6.0),   # Mid-range: 8 - (4-3)*2 = 6
        (5.0, 4.0),   # Boundary: 8 - (5-3)*2 = 4
        (6.0, 2.0),   # 4 - (6-5)*2 = 2
        (7.0, 0.0),   # Boundary: 4 - (7-5)*2 = 0
        (8.0, 0.0),   # Very difficult
        (10.0, 0.0),  # Maximum difficulty
    ])
    def test_sa_score_to_points(self, sa, expected):
        """Test SA score conversion at key boundaries."""
        result = MLReadinessScorer._sa_score_to_points(sa)
        assert abs(result - expected) < 0.01, f"SA={sa}: expected {expected}, got {result}"


class TestStereocenterPoints:
    """Unit tests for _stereocenter_points static method."""

    @pytest.mark.parametrize("total,undefined,expected", [
        (0, 0, 5.0),    # No centers → full
        (1, 0, 5.0),    # 1 defined → full
        (4, 0, 5.0),    # 4 defined → full
        (5, 0, 4.25),   # 5 centers: 5 - (5-4)*0.75 = 4.25
        (8, 0, 2.0),    # 8 centers: 5 - (8-4)*0.75 = 2
        (9, 0, 0.0),    # 9 centers → 0
        (20, 0, 0.0),   # Many centers → 0
        (4, 3, 2.5),    # 4 total, 3 undefined (75% > 50%) → 5*0.5 = 2.5
        (4, 1, 5.0),    # 4 total, 1 undefined (25% < 50%) → no penalty
        (4, 2, 5.0),    # 4 total, 2 undefined (50% = 50%) → no penalty (>50% required)
        (6, 4, 1.75),   # 6 total, 4 undefined → base 3.5 * 0.5 = 1.75
    ])
    def test_stereocenter_points(self, total, undefined, expected):
        """Test stereocenter scoring at boundaries."""
        result = MLReadinessScorer._stereocenter_points(total, undefined)
        assert abs(result - expected) < 0.05, (
            f"total={total}, undef={undefined}: expected {expected}, got {result}"
        )


# =============================================================================
# Dimension 4: Representation Quality (20 pts)
# =============================================================================


class TestRepresentationQuality:
    """Tests for descriptor/FP/conformer dimension."""

    def test_four_items(self, scorer):
        """Should have 4 items."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        rq = scorer.score(mol).breakdown.representation_quality
        assert len(rq.items) == 4
        names = [i.name for i in rq.items]
        assert "Descriptor completeness" in names
        assert "Fingerprint generation" in names
        assert "Fingerprint informativeness" in names
        assert "Conformer generation" in names

    def test_aspirin_all_fps_succeed(self, scorer):
        """Aspirin should generate all 7 FP types."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        details = scorer.score(mol).breakdown.representation_quality.details
        assert len(details["fingerprints_successful"]) == 7
        assert len(details["fingerprints_failed"]) == 0

    def test_descriptor_completeness_high_for_drug(self, scorer):
        """Drug-like molecule should have near-perfect descriptor completeness."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        rq = scorer.score(mol).breakdown.representation_quality
        desc = next(i for i in rq.items if i.name == "Descriptor completeness")
        assert desc.score >= 4.5  # Near-perfect out of 5

    def test_conformer_generation_succeeds_for_3d_mol(self, scorer):
        """Standard molecule should succeed at conformer generation."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        rq = scorer.score(mol).breakdown.representation_quality
        conf = next(i for i in rq.items if i.name == "Conformer generation")
        assert conf.score > 0  # Should succeed (5 or 3 for fallback)

    def test_fp_informativeness_for_drug(self, scorer):
        """Drug molecule should have good bit density → 5/5."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        rq = scorer.score(mol).breakdown.representation_quality
        info = next(i for i in rq.items if i.name == "Fingerprint informativeness")
        assert info.score >= 3.0  # Reasonably informative


class TestFingerprintInformativeness:
    """Unit tests for _score_fingerprint_informativeness."""

    def test_empty_fps_returns_zero(self, scorer):
        """No fingerprints → 0 pts."""
        assert scorer._score_fingerprint_informativeness({}) == 0.0

    def test_ideal_density_returns_5(self, scorer):
        """Fingerprints with 1-30% density → 5 pts."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        _, _, _, fps = scorer._score_fingerprints(mol)
        pts = scorer._score_fingerprint_informativeness(fps)
        assert pts == 5.0


# =============================================================================
# Caveats
# =============================================================================


class TestCaveats:
    """Tests for caveat generation from structural checks."""

    def test_mixture_caveat(self):
        """Multi-component → caveat."""
        mol = Chem.MolFromSmiles("CC.CC")
        result = calculate_ml_readiness(mol)
        assert any("Multi-component" in c for c in result.caveats)

    def test_inorganic_caveat(self):
        """Inorganic compound → caveat."""
        mol = Chem.MolFromSmiles("[Na+].[Cl-]")
        result = calculate_ml_readiness(mol)
        caveats_text = " ".join(result.caveats)
        assert "Inorganic" in caveats_text or "non-standard" in caveats_text

    def test_radical_caveat(self):
        """Radical → caveat."""
        mol = Chem.MolFromSmiles("[CH2]C")
        if mol is not None:
            result = calculate_ml_readiness(mol)
            assert any("Radical" in c for c in result.caveats)

    def test_isotope_caveat(self):
        """Deuterium-labeled → caveat."""
        mol = Chem.MolFromSmiles("[2H]C([2H])([2H])O")
        if mol is not None:
            result = calculate_ml_readiness(mol)
            assert any("Isotope" in c for c in result.caveats)

    def test_trivial_molecule_caveat(self):
        """Single heavy atom → trivial caveat."""
        mol = Chem.MolFromSmiles("C")
        result = calculate_ml_readiness(mol)
        assert any("Trivial" in c for c in result.caveats)

    def test_clean_molecule_no_caveats(self):
        """Clean drug-like molecule → zero caveats."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_ml_readiness(mol)
        assert len(result.caveats) == 0


# =============================================================================
# Supplementary data
# =============================================================================


class TestSupplementary:
    """Tests for supplementary (non-scored) data."""

    def test_lipinski_present(self):
        """Supplementary should have Lipinski info."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_ml_readiness(mol)
        assert "lipinski_violations" in result.supplementary
        assert "lipinski_passed" in result.supplementary

    def test_veber_present(self):
        """Supplementary should have Veber info."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_ml_readiness(mol)
        assert "veber_passed" in result.supplementary

    def test_lipinski_violations_correct_type(self):
        """Lipinski violations should be an integer."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_ml_readiness(mol)
        assert isinstance(result.supplementary["lipinski_violations"], int)

    def test_aspirin_passes_lipinski(self):
        """Aspirin should pass Lipinski."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_ml_readiness(mol)
        assert result.supplementary["lipinski_passed"] is True


# =============================================================================
# Label/color thresholds
# =============================================================================


class TestLabelColorThresholds:
    """Test that label/color mapping follows design spec."""

    def test_excellent_threshold(self):
        """Score 85+ → Excellent/green."""
        from app.services.scoring.ml_readiness import _label_and_color
        assert _label_and_color(85) == ("Excellent", "green")
        assert _label_and_color(100) == ("Excellent", "green")

    def test_good_threshold(self):
        """Score 70-84 → Good/teal."""
        from app.services.scoring.ml_readiness import _label_and_color
        assert _label_and_color(70) == ("Good", "teal")
        assert _label_and_color(84) == ("Good", "teal")

    def test_moderate_threshold(self):
        """Score 50-69 → Moderate/amber."""
        from app.services.scoring.ml_readiness import _label_and_color
        assert _label_and_color(50) == ("Moderate", "amber")
        assert _label_and_color(69) == ("Moderate", "amber")

    def test_limited_threshold(self):
        """Score 30-49 → Limited/orange."""
        from app.services.scoring.ml_readiness import _label_and_color
        assert _label_and_color(30) == ("Limited", "orange")
        assert _label_and_color(49) == ("Limited", "orange")

    def test_poor_threshold(self):
        """Score 0-29 → Poor/red."""
        from app.services.scoring.ml_readiness import _label_and_color
        assert _label_and_color(0) == ("Poor", "red")
        assert _label_and_color(29) == ("Poor", "red")


# =============================================================================
# Score distribution and math
# =============================================================================


class TestScoreDistribution:
    """Tests for dimensional budget (20+35+25+20=100)."""

    def test_max_scores_sum_to_100(self, scorer):
        """Dimension max scores should sum to 100."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        b = scorer.score(mol).breakdown
        total = (
            b.structural_quality.max_score
            + b.property_profile.max_score
            + b.complexity_feasibility.max_score
            + b.representation_quality.max_score
        )
        assert total == 100.0

    def test_score_equals_dimension_sum(self, scorer):
        """Total score should equal sum of dimension scores (within rounding)."""
        for smi in ["CC(=O)Oc1ccccc1C(=O)O", "C", "c1ccncc1", "CCCCCCCCCCCCCCCC"]:
            result = scorer.score(Chem.MolFromSmiles(smi))
            b = result.breakdown
            dim_sum = (
                b.structural_quality.score
                + b.property_profile.score
                + b.complexity_feasibility.score
                + b.representation_quality.score
            )
            assert abs(dim_sum - result.score) < 2, (
                f"{smi}: dims sum to {dim_sum} but score is {result.score}"
            )

    def test_each_dimension_bounded(self, scorer):
        """Each dimension score should be 0 <= score <= max_score."""
        for smi in ["CC(=O)Oc1ccccc1C(=O)O", "C", "[Na+].[Cl-]", "[CH2]C"]:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            b = scorer.score(mol).breakdown
            for dim in [b.structural_quality, b.property_profile,
                        b.complexity_feasibility, b.representation_quality]:
                assert 0 <= dim.score <= dim.max_score, (
                    f"{smi}: {dim.name} score {dim.score} out of bounds [0, {dim.max_score}]"
                )

    def test_item_scores_bounded(self, scorer):
        """Each item score should be 0 <= score <= max_score."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        b = scorer.score(mol).breakdown
        for dim in [b.structural_quality, b.property_profile,
                    b.complexity_feasibility, b.representation_quality]:
            for item in dim.items:
                assert 0 <= item.score <= item.max_score, (
                    f"{dim.name}/{item.name}: {item.score} out of [{0}, {item.max_score}]"
                )


# =============================================================================
# Convenience function
# =============================================================================


class TestCalculateMLReadinessFunction:
    """Tests for the module-level convenience function."""

    def test_returns_correct_type(self):
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_ml_readiness(mol)
        assert isinstance(result, MLReadinessResult)

    def test_score_in_range(self):
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = calculate_ml_readiness(mol)
        assert 0 <= result.score <= 100


# =============================================================================
# Edge cases
# =============================================================================


class TestEdgeCases:
    """Tests for unusual molecules and error resilience."""

    def test_silicon_compound(self):
        """Organosilicon → valid score, likely caveats."""
        mol = Chem.MolFromSmiles("[Si](C)(C)(C)C")
        result = calculate_ml_readiness(mol)
        assert 0 <= result.score <= 100

    def test_salt_form(self):
        """Salt (multi-component, charged) → valid but lower score."""
        mol = Chem.MolFromSmiles("CC(=O)[O-].[Na+]")
        result = calculate_ml_readiness(mol)
        assert 0 <= result.score <= 100
        assert result.breakdown.structural_quality.score < 20  # Loses mixture pts

    def test_defined_stereochemistry(self):
        """L-alanine with defined stereo → valid score."""
        mol = Chem.MolFromSmiles("C[C@H](N)C(=O)O")
        result = calculate_ml_readiness(mol)
        assert 0 <= result.score <= 100

    def test_long_alkyl_chain(self):
        """Very long chain → low property profile (high MW, no rings)."""
        mol = Chem.MolFromSmiles("C" * 50)
        result = calculate_ml_readiness(mol)
        assert result.breakdown.property_profile.score < 15  # < 43% of 35

    def test_aromatic_heterocycle(self):
        """Pyridine → decent score."""
        mol = Chem.MolFromSmiles("c1ccncc1")
        result = calculate_ml_readiness(mol)
        assert result.score >= 40


# =============================================================================
# Diverse molecules — discrimination ability
# =============================================================================


class TestScoreDiscrimination:
    """Test that the scorer differentiates structurally diverse molecules."""

    def test_drug_vs_inorganic(self, scorer):
        """Drug-like molecule should score much higher than inorganic."""
        drug = scorer.score(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"))  # Aspirin
        inorg = scorer.score(Chem.MolFromSmiles("[Na+].[Cl-]"))  # NaCl

        assert drug.score > inorg.score + 20  # Meaningful gap

    def test_drug_vs_methane(self, scorer):
        """Drug-like molecule should score higher than methane."""
        drug = scorer.score(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"))
        tiny = scorer.score(Chem.MolFromSmiles("C"))

        assert drug.score > tiny.score + 10

    def test_drug_vs_huge_polymer_like(self, scorer):
        """Drug-like molecule should score higher than huge alkane chain."""
        drug = scorer.score(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"))
        huge = scorer.score(Chem.MolFromSmiles("C" * 100))

        assert drug.score > huge.score
