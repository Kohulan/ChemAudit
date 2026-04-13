"""
Tests for the Audit Column Registry (audit_columns.py)

Covers:
- All 6 sections and their column counts
- Every extractor helper (happy path + missing/None/False)
- Pass/Fail conversion for boolean fields
- Custom one-off extractors
- Public output helpers: get_identity_row, get_flat_headers,
  extract_flat_row, extract_nested, extract_by_section
- Total column count (~91)
"""

from collections import OrderedDict

import pytest

from app.services.export.audit_columns import (
    AUDIT_SECTIONS,
    AuditColumn,
    AuditSection,
    extract_by_section,
    extract_flat_row,
    extract_nested,
    get_flat_headers,
    get_identity_row,
)

# ---------------------------------------------------------------------------
# Shared rich fixture
# ---------------------------------------------------------------------------

FULL_RESULT: dict = {
    "smiles": "CCO",
    "name": "Ethanol",
    "standardized_smiles": "CCO",
    "validation": {
        "canonical_smiles": "CCO",
        "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "overall_score": 95,
        "molecule_info": {
            "inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
            "molecular_formula": "C2H6O",
            "molecular_weight": 46.07,
        },
        "all_checks": [
            {"check_name": "parsability", "passed": True},
            {"check_name": "sanitization", "passed": True},
            {"check_name": "valence", "passed": True},
            {"check_name": "aromaticity", "passed": True},
            {"check_name": "connectivity", "passed": True},
            {"check_name": "undefined_stereocenters", "passed": True},
            {"check_name": "undefined_doublebond_stereo", "passed": True},
            {"check_name": "conflicting_stereo", "passed": True},
            {"check_name": "smiles_roundtrip", "passed": True},
            {"check_name": "inchi_generation", "passed": True},
            {"check_name": "inchi_roundtrip", "passed": True},
            # deep checks
            {"check_name": "stereoisomer_enumeration", "passed": True},
            {"check_name": "tautomer_detection", "passed": True},
            {"check_name": "aromatic_system_validation", "passed": True},
            {"check_name": "coordinate_dimension", "passed": True},
            {"check_name": "mixture_detection", "passed": True},
            {"check_name": "solvent_contamination", "passed": True},
            {"check_name": "inorganic_filter", "passed": True},
            {"check_name": "radical_detection", "passed": True},
            {"check_name": "isotope_label_detection", "passed": True},
            {"check_name": "trivial_molecule", "passed": False},
            {"check_name": "hypervalent_atoms", "passed": True},
            {"check_name": "polymer_detection", "passed": True},
            {"check_name": "ring_strain", "passed": True},
            {"check_name": "macrocycle_detection", "passed": True},
            {"check_name": "charged_species", "passed": True},
            {"check_name": "explicit_hydrogen_audit", "passed": True},
        ],
    },
    "scoring": {
        "ml_readiness_score": 88,
        "ml_readiness": {"label": "Good"},
        "druglikeness": {
            "qed_score": 0.72,
            "lipinski": {
                "passed": True,
                "violations": 0,
                "mw": 46.07,
                "logp": -0.31,
                "hbd": 1,
                "hba": 1,
            },
            "veber": {"passed": True, "rotatable_bonds": 0, "tpsa": 20.23},
            "ro3": {"passed": True},
            "ghose": {"passed": False},
            "egan": {"passed": True},
            "muegge": {"passed": True},
        },
        "np_likeness_score": -0.5,
        "consensus": {"score": 3.2},
        "lead_likeness": {"passed": True},
        "admet": {
            "synthetic_accessibility": {"score": 1.5, "classification": "Easy"},
            "solubility": {"log_s": -1.2, "classification": "Soluble"},
            "complexity": {"fsp3": 0.33},
            "cns_mpo": {"score": 4.5},
            "bioavailability": {"oral_absorption_likely": True},
            "pfizer_rule": {"passed": True},
            "gsk_rule": {"passed": True},
            "golden_triangle": {"in_golden_triangle": True},
        },
        "aggregator": {"likelihood": "Low", "risk_score": 0.05},
        "scaffold": {"scaffold_smiles": ""},
        "boiled_egg": {"gi_absorbed": True, "bbb_permeant": False, "region": "HIA"},
        "safety_filters": {
            "pains": {"passed": True, "alert_count": 0},
            "brenk": {"passed": True, "alert_count": 0},
            "nih": {"passed": True},
            "zinc": {"passed": True},
            "chembl": {"passed": True},
            "all_passed": True,
            "total_alerts": 0,
        },
    },
    "safety_assessment": {
        "herg": {"herg_risk": "low", "risk_score": 1},
        "bro5": {"passed": True},
        "reos": {"passed": True, "n_violations": 0},
        "cyp_softspots": {"n_sites": 2},
        "complexity": {"n_outliers": 0},
    },
    "profiling": {
        "pfi": {"pfi": 3.5, "risk": "low"},
        "stars": {"stars": 5},
        "abbott": {"abbott_score": 0.85, "probability_pct": 85.0},
        "consensus_logp": {"consensus_logp": -0.31},
        "skin_permeation": {"log_kp": -2.5, "classification": "Medium"},
        "sa_comparison": {
            "sa_score": {"score": 1.5},
            "scscore": {"score": 1.2},
            "syba": {"score": 90.0},
        },
        "cns_mpo": {"score": 3.8},
    },
    "standardization": {
        "result": {
            "standardized_smiles": "CCO",
            "success": True,
            "steps_applied": [
                {"name": "cleanup", "applied": True},
                {"name": "normalize", "applied": True},
                {"name": "reionize", "applied": False},
            ],
            "excluded_fragments": ["[Na+]", "[Cl-]"],
            "stereo_comparison": {"lost": 0, "gained": 1},
            "mass_change_percent": 0.0,
        }
    },
}

EMPTY_RESULT: dict = {}


# ---------------------------------------------------------------------------
# Section registry tests
# ---------------------------------------------------------------------------


class TestAuditSectionsRegistry:
    def test_six_sections_exist(self):
        assert len(AUDIT_SECTIONS) == 6

    def test_section_names(self):
        names = [s.name for s in AUDIT_SECTIONS]
        assert names == [
            "Validation",
            "Deep Validation",
            "Scoring",
            "Safety",
            "Compound Profile",
            "Standardization",
        ]

    def test_section_prefixes(self):
        prefixes = [s.prefix for s in AUDIT_SECTIONS]
        assert prefixes == [
            "[Validation]",
            "[Deep Validation]",
            "[Scoring]",
            "[Safety]",
            "[Compound Profile]",
            "[Standardization]",
        ]

    def test_section_types(self):
        for section in AUDIT_SECTIONS:
            assert isinstance(section, AuditSection)
            for col in section.columns:
                assert isinstance(col, AuditColumn)

    def test_total_column_count_at_least_91(self):
        total = sum(len(s.columns) for s in AUDIT_SECTIONS)
        assert total >= 91, f"Expected >=91 columns, got {total}"

    def test_validation_section_has_17_columns(self):
        val_section = next(s for s in AUDIT_SECTIONS if s.name == "Validation")
        assert len(val_section.columns) == 17

    def test_deep_validation_section_has_16_columns(self):
        deep = next(s for s in AUDIT_SECTIONS if s.name == "Deep Validation")
        assert len(deep.columns) == 16

    def test_scoring_section_has_35_columns(self):
        scoring = next(s for s in AUDIT_SECTIONS if s.name == "Scoring")
        assert len(scoring.columns) == 35

    def test_safety_section_has_16_columns(self):
        safety = next(s for s in AUDIT_SECTIONS if s.name == "Safety")
        assert len(safety.columns) == 16

    def test_compound_profile_section_has_12_columns(self):
        cp = next(s for s in AUDIT_SECTIONS if s.name == "Compound Profile")
        assert len(cp.columns) == 12

    def test_standardization_section_has_6_columns(self):
        std = next(s for s in AUDIT_SECTIONS if s.name == "Standardization")
        assert len(std.columns) == 6

    def test_all_keys_are_unique(self):
        keys = [col.key for s in AUDIT_SECTIONS for col in s.columns]
        assert len(keys) == len(set(keys)), "Duplicate column keys found"

    def test_all_headers_are_strings(self):
        for section in AUDIT_SECTIONS:
            for col in section.columns:
                assert isinstance(col.header, str)
                assert col.header.strip() != ""

    def test_all_extractors_are_callable(self):
        for section in AUDIT_SECTIONS:
            for col in section.columns:
                assert callable(col.extractor), f"Extractor for {col.key!r} is not callable"


# ---------------------------------------------------------------------------
# Validation section extractors
# ---------------------------------------------------------------------------


class TestValidationExtractors:
    @pytest.fixture(autouse=True)
    def _section(self):
        self.section = next(s for s in AUDIT_SECTIONS if s.name == "Validation")
        self.cols = {col.key: col for col in self.section.columns}

    def test_canonical_smiles(self):
        assert self.cols["canonical_smiles"].extractor(FULL_RESULT) == "CCO"

    def test_inchikey(self):
        assert self.cols["inchikey"].extractor(FULL_RESULT) == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"

    def test_inchi(self):
        val = self.cols["inchi"].extractor(FULL_RESULT)
        assert val.startswith("InChI=")

    def test_molecular_formula(self):
        assert self.cols["molecular_formula"].extractor(FULL_RESULT) == "C2H6O"

    def test_molecular_weight(self):
        assert self.cols["molecular_weight"].extractor(FULL_RESULT) == 46.07

    def test_overall_score(self):
        assert self.cols["overall_score"].extractor(FULL_RESULT) == 95

    def test_parsability_pass(self):
        assert self.cols["parsability_passed"].extractor(FULL_RESULT) == "Pass"

    def test_trivial_molecule_fail(self):
        # trivial_molecule is set to passed=False in FULL_RESULT
        # But trivial_molecule is in deep validation; basic check test below
        # Test that a failed check returns "Fail"
        result = {
            "validation": {
                "all_checks": [{"check_name": "parsability", "passed": False}]
            }
        }
        assert self.cols["parsability_passed"].extractor(result) == "Fail"

    def test_missing_check_returns_na(self):
        assert self.cols["parsability_passed"].extractor(EMPTY_RESULT) == "N/A"

    def test_none_passed_value_returns_na(self):
        result = {
            "validation": {
                "all_checks": [{"check_name": "parsability", "passed": None}]
            }
        }
        assert self.cols["parsability_passed"].extractor(result) == "N/A"

    def test_empty_result_defaults(self):
        assert self.cols["canonical_smiles"].extractor(EMPTY_RESULT) == ""
        assert self.cols["inchikey"].extractor(EMPTY_RESULT) == ""
        assert self.cols["inchi"].extractor(EMPTY_RESULT) == ""
        assert self.cols["molecular_formula"].extractor(EMPTY_RESULT) == ""
        assert self.cols["molecular_weight"].extractor(EMPTY_RESULT) == "N/A"
        assert self.cols["overall_score"].extractor(EMPTY_RESULT) == "N/A"

    def test_all_basic_checks_present(self):
        expected_keys = [
            "parsability_passed",
            "sanitization_passed",
            "valence_passed",
            "aromaticity_passed",
            "connectivity_passed",
            "undefined_stereocenters_passed",
            "undefined_doublebond_stereo_passed",
            "conflicting_stereo_passed",
            "smiles_roundtrip_passed",
            "inchi_generation_passed",
            "inchi_roundtrip_passed",
        ]
        for key in expected_keys:
            assert key in self.cols, f"Missing column {key!r}"
            assert self.cols[key].extractor(FULL_RESULT) == "Pass"


# ---------------------------------------------------------------------------
# Deep Validation section extractors
# ---------------------------------------------------------------------------


class TestDeepValidationExtractors:
    @pytest.fixture(autouse=True)
    def _section(self):
        self.section = next(s for s in AUDIT_SECTIONS if s.name == "Deep Validation")
        self.cols = {col.key: col for col in self.section.columns}

    def test_trivial_molecule_fail(self):
        assert self.cols["trivial_molecule_passed"].extractor(FULL_RESULT) == "Fail"

    def test_stereoisomer_enumeration_pass(self):
        assert self.cols["stereoisomer_enumeration_passed"].extractor(FULL_RESULT) == "Pass"

    def test_missing_deep_check_returns_na(self):
        assert self.cols["polymer_detection_passed"].extractor(EMPTY_RESULT) == "N/A"

    def test_all_16_deep_check_keys_present(self):
        expected = [
            "stereoisomer_enumeration_passed",
            "tautomer_detection_passed",
            "aromatic_system_validation_passed",
            "coordinate_dimension_passed",
            "mixture_detection_passed",
            "solvent_contamination_passed",
            "inorganic_filter_passed",
            "radical_detection_passed",
            "isotope_label_detection_passed",
            "trivial_molecule_passed",
            "hypervalent_atoms_passed",
            "polymer_detection_passed",
            "ring_strain_passed",
            "macrocycle_detection_passed",
            "charged_species_passed",
            "explicit_hydrogen_audit_passed",
        ]
        for key in expected:
            assert key in self.cols, f"Missing deep validation column {key!r}"


# ---------------------------------------------------------------------------
# Scoring section extractors
# ---------------------------------------------------------------------------


class TestScoringExtractors:
    @pytest.fixture(autouse=True)
    def _section(self):
        self.section = next(s for s in AUDIT_SECTIONS if s.name == "Scoring")
        self.cols = {col.key: col for col in self.section.columns}

    def test_ml_readiness_score(self):
        assert self.cols["ml_readiness_score"].extractor(FULL_RESULT) == 88

    def test_ml_readiness_label(self):
        assert self.cols["ml_readiness_label"].extractor(FULL_RESULT) == "Good"

    def test_qed_score(self):
        assert self.cols["qed_score"].extractor(FULL_RESULT) == 0.72

    def test_lipinski_passed(self):
        assert self.cols["lipinski_passed"].extractor(FULL_RESULT) == "Pass"

    def test_lipinski_violations(self):
        assert self.cols["lipinski_violations"].extractor(FULL_RESULT) == 0

    def test_lipinski_mw(self):
        assert self.cols["lipinski_mw"].extractor(FULL_RESULT) == 46.07

    def test_lipinski_logp(self):
        assert self.cols["lipinski_logp"].extractor(FULL_RESULT) == -0.31

    def test_lipinski_hbd(self):
        assert self.cols["lipinski_hbd"].extractor(FULL_RESULT) == 1

    def test_lipinski_hba(self):
        assert self.cols["lipinski_hba"].extractor(FULL_RESULT) == 1

    def test_veber_passed(self):
        assert self.cols["veber_passed"].extractor(FULL_RESULT) == "Pass"

    def test_veber_rotatable_bonds(self):
        assert self.cols["veber_rotatable_bonds"].extractor(FULL_RESULT) == 0

    def test_veber_tpsa(self):
        assert self.cols["veber_tpsa"].extractor(FULL_RESULT) == 20.23

    def test_ro3_passed(self):
        assert self.cols["ro3_passed"].extractor(FULL_RESULT) == "Pass"

    def test_ghose_passed_false(self):
        assert self.cols["ghose_passed"].extractor(FULL_RESULT) == "Fail"

    def test_egan_passed(self):
        assert self.cols["egan_passed"].extractor(FULL_RESULT) == "Pass"

    def test_muegge_passed(self):
        assert self.cols["muegge_passed"].extractor(FULL_RESULT) == "Pass"

    def test_np_likeness_score(self):
        assert self.cols["np_likeness_score"].extractor(FULL_RESULT) == -0.5

    def test_consensus_score(self):
        assert self.cols["consensus_score"].extractor(FULL_RESULT) == 3.2

    def test_lead_likeness_passed(self):
        assert self.cols["lead_likeness_passed"].extractor(FULL_RESULT) == "Pass"

    def test_sa_score_uses_correct_path(self):
        # path is scoring.admet.synthetic_accessibility.score NOT scoring.admet.sa_score
        assert self.cols["sa_score"].extractor(FULL_RESULT) == 1.5

    def test_sa_classification(self):
        assert self.cols["sa_classification"].extractor(FULL_RESULT) == "Easy"

    def test_esol_log_s(self):
        assert self.cols["esol_log_s"].extractor(FULL_RESULT) == -1.2

    def test_esol_classification(self):
        assert self.cols["esol_classification"].extractor(FULL_RESULT) == "Soluble"

    def test_fsp3(self):
        assert self.cols["fsp3"].extractor(FULL_RESULT) == 0.33

    def test_cns_mpo(self):
        assert self.cols["cns_mpo"].extractor(FULL_RESULT) == 4.5

    def test_oral_absorption(self):
        assert self.cols["oral_absorption"].extractor(FULL_RESULT) == "Pass"

    def test_pfizer_rule(self):
        assert self.cols["pfizer_rule"].extractor(FULL_RESULT) == "Pass"

    def test_gsk_rule(self):
        assert self.cols["gsk_rule"].extractor(FULL_RESULT) == "Pass"

    def test_golden_triangle(self):
        assert self.cols["golden_triangle"].extractor(FULL_RESULT) == "Pass"

    def test_aggregator_likelihood(self):
        assert self.cols["aggregator_likelihood"].extractor(FULL_RESULT) == "Low"

    def test_aggregator_risk(self):
        assert self.cols["aggregator_risk"].extractor(FULL_RESULT) == 0.05

    def test_scaffold_smiles_empty_string(self):
        assert self.cols["scaffold_smiles"].extractor(FULL_RESULT) == ""

    def test_boiled_egg_gi(self):
        assert self.cols["boiled_egg_gi"].extractor(FULL_RESULT) == "Pass"

    def test_boiled_egg_bbb_false(self):
        assert self.cols["boiled_egg_bbb"].extractor(FULL_RESULT) == "Fail"

    def test_boiled_egg_region(self):
        assert self.cols["boiled_egg_region"].extractor(FULL_RESULT) == "HIA"

    def test_missing_scoring_returns_na(self):
        assert self.cols["ml_readiness_score"].extractor(EMPTY_RESULT) == "N/A"
        assert self.cols["qed_score"].extractor(EMPTY_RESULT) == "N/A"
        assert self.cols["lipinski_passed"].extractor(EMPTY_RESULT) == "N/A"

    def test_false_bool_returns_fail(self):
        result = {
            "scoring": {
                "druglikeness": {"lipinski": {"passed": False}}
            }
        }
        assert self.cols["lipinski_passed"].extractor(result) == "Fail"

    def test_none_bool_returns_na(self):
        result = {
            "scoring": {
                "druglikeness": {"lipinski": {"passed": None}}
            }
        }
        assert self.cols["lipinski_passed"].extractor(result) == "N/A"


# ---------------------------------------------------------------------------
# Safety section extractors
# ---------------------------------------------------------------------------


class TestSafetyExtractors:
    @pytest.fixture(autouse=True)
    def _section(self):
        self.section = next(s for s in AUDIT_SECTIONS if s.name == "Safety")
        self.cols = {col.key: col for col in self.section.columns}

    def test_pains_passed(self):
        assert self.cols["pains_passed"].extractor(FULL_RESULT) == "Pass"

    def test_pains_count(self):
        assert self.cols["pains_count"].extractor(FULL_RESULT) == 0

    def test_brenk_passed(self):
        assert self.cols["brenk_passed"].extractor(FULL_RESULT) == "Pass"

    def test_brenk_count(self):
        assert self.cols["brenk_count"].extractor(FULL_RESULT) == 0

    def test_nih_passed(self):
        assert self.cols["nih_passed"].extractor(FULL_RESULT) == "Pass"

    def test_zinc_passed(self):
        assert self.cols["zinc_passed"].extractor(FULL_RESULT) == "Pass"

    def test_chembl_passed(self):
        assert self.cols["chembl_passed"].extractor(FULL_RESULT) == "Pass"

    def test_safety_all_passed_direct_field(self):
        # all_passed is directly under safety_filters, not under a sub-filter
        assert self.cols["safety_all_passed"].extractor(FULL_RESULT) == "Pass"

    def test_total_alerts(self):
        assert self.cols["total_alerts"].extractor(FULL_RESULT) == 0

    def test_herg_risk(self):
        assert self.cols["herg_risk"].extractor(FULL_RESULT) == "low"

    def test_herg_risk_score(self):
        assert self.cols["herg_risk_score"].extractor(FULL_RESULT) == 1

    def test_bro5_passed(self):
        assert self.cols["bro5_passed"].extractor(FULL_RESULT) == "Pass"

    def test_reos_passed(self):
        assert self.cols["reos_passed"].extractor(FULL_RESULT) == "Pass"

    def test_reos_violations(self):
        assert self.cols["reos_violations"].extractor(FULL_RESULT) == 0

    def test_cyp_softspot_count(self):
        assert self.cols["cyp_softspot_count"].extractor(FULL_RESULT) == 2

    def test_complexity_outliers(self):
        assert self.cols["complexity_outliers"].extractor(FULL_RESULT) == 0

    def test_pains_failed_converts_bool(self):
        result = {
            "scoring": {
                "safety_filters": {
                    "pains": {"passed": False, "alert_count": 3}
                }
            }
        }
        assert self.cols["pains_passed"].extractor(result) == "Fail"
        assert self.cols["pains_count"].extractor(result) == 3

    def test_bro5_false_returns_fail(self):
        result = {"safety_assessment": {"bro5": {"passed": False}}}
        assert self.cols["bro5_passed"].extractor(result) == "Fail"

    def test_reos_false_returns_fail(self):
        result = {"safety_assessment": {"reos": {"passed": False}}}
        assert self.cols["reos_passed"].extractor(result) == "Fail"

    def test_safety_all_passed_false(self):
        result = {"scoring": {"safety_filters": {"all_passed": False}}}
        assert self.cols["safety_all_passed"].extractor(result) == "Fail"

    def test_missing_safety_returns_na(self):
        assert self.cols["pains_passed"].extractor(EMPTY_RESULT) == "N/A"
        assert self.cols["herg_risk"].extractor(EMPTY_RESULT) == "N/A"
        assert self.cols["bro5_passed"].extractor(EMPTY_RESULT) == "N/A"
        assert self.cols["reos_passed"].extractor(EMPTY_RESULT) == "N/A"


# ---------------------------------------------------------------------------
# Compound Profile section extractors
# ---------------------------------------------------------------------------


class TestCompoundProfileExtractors:
    @pytest.fixture(autouse=True)
    def _section(self):
        self.section = next(s for s in AUDIT_SECTIONS if s.name == "Compound Profile")
        self.cols = {col.key: col for col in self.section.columns}

    def test_pfi_score(self):
        assert self.cols["pfi_score"].extractor(FULL_RESULT) == 3.5

    def test_pfi_risk(self):
        assert self.cols["pfi_risk"].extractor(FULL_RESULT) == "low"

    def test_stars_count(self):
        assert self.cols["stars_count"].extractor(FULL_RESULT) == 5

    def test_abbott_score(self):
        assert self.cols["abbott_score"].extractor(FULL_RESULT) == 0.85

    def test_abbott_probability(self):
        assert self.cols["abbott_probability"].extractor(FULL_RESULT) == 85.0

    def test_consensus_logp(self):
        assert self.cols["consensus_logp"].extractor(FULL_RESULT) == -0.31

    def test_skin_perm_logkp(self):
        assert self.cols["skin_perm_logkp"].extractor(FULL_RESULT) == -2.5

    def test_skin_perm_class(self):
        assert self.cols["skin_perm_class"].extractor(FULL_RESULT) == "Medium"

    def test_sa_comparison_sa(self):
        assert self.cols["sa_comparison_sa"].extractor(FULL_RESULT) == 1.5

    def test_sa_comparison_scscore(self):
        assert self.cols["sa_comparison_scscore"].extractor(FULL_RESULT) == 1.2

    def test_sa_comparison_syba(self):
        assert self.cols["sa_comparison_syba"].extractor(FULL_RESULT) == 90.0

    def test_profile_cns_mpo(self):
        assert self.cols["profile_cns_mpo"].extractor(FULL_RESULT) == 3.8

    def test_missing_profiling_returns_na(self):
        assert self.cols["pfi_score"].extractor(EMPTY_RESULT) == "N/A"
        assert self.cols["sa_comparison_sa"].extractor(EMPTY_RESULT) == "N/A"

    def test_partial_profiling_returns_na(self):
        result = {"profiling": {"pfi": {"pfi": 5.0}}}  # missing risk
        assert self.cols["pfi_score"].extractor(result) == 5.0
        assert self.cols["pfi_risk"].extractor(result) == "N/A"


# ---------------------------------------------------------------------------
# Standardization section extractors
# ---------------------------------------------------------------------------


class TestStandardizationExtractors:
    @pytest.fixture(autouse=True)
    def _section(self):
        self.section = next(s for s in AUDIT_SECTIONS if s.name == "Standardization")
        self.cols = {col.key: col for col in self.section.columns}

    def test_standardized_smiles_prefers_top_level(self):
        # FULL_RESULT has top-level standardized_smiles="CCO"
        assert self.cols["standardized_smiles"].extractor(FULL_RESULT) == "CCO"

    def test_standardized_smiles_falls_back_to_nested(self):
        result = {
            "standardization": {
                "result": {"standardized_smiles": "c1ccccc1"}
            }
        }
        assert self.cols["standardized_smiles"].extractor(result) == "c1ccccc1"

    def test_standardized_smiles_empty_when_missing(self):
        assert self.cols["standardized_smiles"].extractor(EMPTY_RESULT) == ""

    def test_std_success_pass(self):
        assert self.cols["std_success"].extractor(FULL_RESULT) == "Pass"

    def test_std_success_fail(self):
        result = {"standardization": {"result": {"success": False}}}
        assert self.cols["std_success"].extractor(result) == "Fail"

    def test_std_success_missing_returns_na(self):
        assert self.cols["std_success"].extractor(EMPTY_RESULT) == "N/A"

    def test_std_steps_count_counts_applied_true(self):
        # FULL_RESULT has 3 steps: 2 applied=True, 1 applied=False
        assert self.cols["std_steps_count"].extractor(FULL_RESULT) == 2

    def test_std_steps_count_missing_returns_na(self):
        assert self.cols["std_steps_count"].extractor(EMPTY_RESULT) == "N/A"

    def test_std_steps_count_all_false(self):
        result = {
            "standardization": {
                "result": {
                    "steps_applied": [
                        {"name": "a", "applied": False},
                        {"name": "b", "applied": False},
                    ]
                }
            }
        }
        assert self.cols["std_steps_count"].extractor(result) == 0

    def test_std_excluded_fragments_count(self):
        # FULL_RESULT has 2 excluded fragments
        assert self.cols["std_excluded_fragments"].extractor(FULL_RESULT) == 2

    def test_std_excluded_fragments_missing_returns_na(self):
        assert self.cols["std_excluded_fragments"].extractor(EMPTY_RESULT) == "N/A"

    def test_std_stereo_changes_format(self):
        assert self.cols["std_stereo_changes"].extractor(FULL_RESULT) == "lost: 0, gained: 1"

    def test_std_stereo_changes_missing_returns_na(self):
        assert self.cols["std_stereo_changes"].extractor(EMPTY_RESULT) == "N/A"

    def test_std_mass_change(self):
        assert self.cols["std_mass_change"].extractor(FULL_RESULT) == 0.0

    def test_std_mass_change_missing_returns_na(self):
        assert self.cols["std_mass_change"].extractor(EMPTY_RESULT) == "N/A"


# ---------------------------------------------------------------------------
# Public output helper tests
# ---------------------------------------------------------------------------


class TestGetIdentityRow:
    def test_basic_fields(self):
        row = get_identity_row(0, FULL_RESULT)
        assert isinstance(row, OrderedDict)
        assert row["index"] == 1
        assert row["name"] == "Ethanol"
        assert row["input_smiles"] == "CCO"

    def test_index_is_one_based(self):
        assert get_identity_row(0, FULL_RESULT)["index"] == 1
        assert get_identity_row(4, FULL_RESULT)["index"] == 5

    def test_missing_name_defaults_empty(self):
        row = get_identity_row(0, {"smiles": "C"})
        assert row["name"] == ""

    def test_missing_smiles_defaults_empty(self):
        row = get_identity_row(0, {"name": "X"})
        assert row["input_smiles"] == ""

    def test_empty_result(self):
        row = get_identity_row(0, EMPTY_RESULT)
        assert row["index"] == 1
        assert row["name"] == ""
        assert row["input_smiles"] == ""

    def test_key_order(self):
        row = get_identity_row(0, FULL_RESULT)
        assert list(row.keys()) == ["index", "name", "input_smiles"]


class TestGetFlatHeaders:
    def test_returns_list(self):
        headers = get_flat_headers()
        assert isinstance(headers, list)

    def test_count_matches_total_columns(self):
        total_cols = sum(len(s.columns) for s in AUDIT_SECTIONS)
        assert len(get_flat_headers()) == total_cols

    def test_headers_have_prefix(self):
        headers = get_flat_headers()
        assert headers[0].startswith("[Validation]")

    def test_all_sections_represented(self):
        headers = get_flat_headers()
        for section in AUDIT_SECTIONS:
            assert any(h.startswith(section.prefix) for h in headers)

    def test_first_header(self):
        headers = get_flat_headers()
        assert headers[0] == "[Validation] Canonical SMILES"

    def test_no_duplicate_headers(self):
        headers = get_flat_headers()
        assert len(headers) == len(set(headers)), "Duplicate headers found"


class TestExtractFlatRow:
    def test_returns_ordered_dict(self):
        row = extract_flat_row(FULL_RESULT)
        assert isinstance(row, OrderedDict)

    def test_key_count_matches_headers(self):
        row = extract_flat_row(FULL_RESULT)
        headers = get_flat_headers()
        assert list(row.keys()) == headers

    def test_first_key_canonical_smiles(self):
        row = extract_flat_row(FULL_RESULT)
        first_key = list(row.keys())[0]
        assert first_key == "[Validation] Canonical SMILES"
        assert row[first_key] == "CCO"

    def test_empty_result_no_exceptions(self):
        row = extract_flat_row(EMPTY_RESULT)
        assert isinstance(row, OrderedDict)
        assert len(row) == len(get_flat_headers())

    def test_scoring_sa_score_path(self):
        row = extract_flat_row(FULL_RESULT)
        assert row["[Scoring] SA Score (1-10)"] == 1.5

    def test_safety_all_passed_direct_path(self):
        row = extract_flat_row(FULL_RESULT)
        assert row["[Safety] All Safety Filters (Pass/Fail)"] == "Pass"

    def test_standardization_steps_count(self):
        row = extract_flat_row(FULL_RESULT)
        assert row["[Standardization] Steps Applied"] == 2

    def test_standardization_stereo_changes(self):
        row = extract_flat_row(FULL_RESULT)
        assert row["[Standardization] Stereo Changes"] == "lost: 0, gained: 1"


class TestExtractNested:
    def test_returns_dict(self):
        nested = extract_nested(FULL_RESULT)
        assert isinstance(nested, dict)

    def test_has_six_section_keys(self):
        nested = extract_nested(FULL_RESULT)
        assert set(nested.keys()) == {
            "validation",
            "deep_validation",
            "scoring",
            "safety",
            "compound_profile",
            "standardization",
        }

    def test_validation_canonical_smiles(self):
        nested = extract_nested(FULL_RESULT)
        assert nested["validation"]["canonical_smiles"] == "CCO"

    def test_scoring_sa_score(self):
        nested = extract_nested(FULL_RESULT)
        assert nested["scoring"]["sa_score"] == 1.5

    def test_safety_bro5(self):
        nested = extract_nested(FULL_RESULT)
        assert nested["safety"]["bro5_passed"] == "Pass"

    def test_compound_profile_pfi(self):
        nested = extract_nested(FULL_RESULT)
        assert nested["compound_profile"]["pfi_score"] == 3.5

    def test_standardization_steps(self):
        nested = extract_nested(FULL_RESULT)
        assert nested["standardization"]["std_steps_count"] == 2

    def test_empty_result_no_exceptions(self):
        nested = extract_nested(EMPTY_RESULT)
        assert isinstance(nested, dict)
        assert len(nested) == 6

    def test_inner_dicts_keyed_by_col_key(self):
        nested = extract_nested(FULL_RESULT)
        val_keys = list(nested["validation"].keys())
        # Should use col.key, not col.header
        assert "canonical_smiles" in val_keys
        assert "Canonical SMILES" not in val_keys


class TestExtractBySection:
    def test_returns_dict(self):
        by_sec = extract_by_section(FULL_RESULT)
        assert isinstance(by_sec, dict)

    def test_has_six_section_names(self):
        by_sec = extract_by_section(FULL_RESULT)
        assert set(by_sec.keys()) == {
            "Validation",
            "Deep Validation",
            "Scoring",
            "Safety",
            "Compound Profile",
            "Standardization",
        }

    def test_each_section_is_ordered_dict(self):
        by_sec = extract_by_section(FULL_RESULT)
        for section_name, data in by_sec.items():
            assert isinstance(data, OrderedDict), f"{section_name!r} is not OrderedDict"

    def test_validation_keyed_by_header(self):
        by_sec = extract_by_section(FULL_RESULT)
        assert "Canonical SMILES" in by_sec["Validation"]
        assert by_sec["Validation"]["Canonical SMILES"] == "CCO"

    def test_scoring_sa_score_header(self):
        by_sec = extract_by_section(FULL_RESULT)
        assert "SA Score (1-10)" in by_sec["Scoring"]
        assert by_sec["Scoring"]["SA Score (1-10)"] == 1.5

    def test_empty_result_no_exceptions(self):
        by_sec = extract_by_section(EMPTY_RESULT)
        assert len(by_sec) == 6


# ---------------------------------------------------------------------------
# Edge case / defensive tests
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_nested_helper_handles_non_dict_midpath(self):
        # If a path segment is a non-dict value, return default
        result = {"scoring": "not_a_dict"}
        cols = {col.key: col for s in AUDIT_SECTIONS for col in s.columns}
        assert cols["ml_readiness_score"].extractor(result) == "N/A"

    def test_pass_fail_helper_handles_non_dict_midpath(self):
        result = {"scoring": {"druglikeness": "oops"}}
        cols = {col.key: col for s in AUDIT_SECTIONS for col in s.columns}
        assert cols["lipinski_passed"].extractor(result) == "N/A"

    def test_safety_filter_handles_missing_sub_filter(self):
        result = {"scoring": {"safety_filters": {}}}
        cols = {col.key: col for s in AUDIT_SECTIONS for col in s.columns}
        assert cols["pains_passed"].extractor(result) == "N/A"

    def test_safety_assess_handles_missing_category(self):
        result = {"safety_assessment": {}}
        cols = {col.key: col for s in AUDIT_SECTIONS for col in s.columns}
        assert cols["herg_risk"].extractor(result) == "N/A"

    def test_profile_handles_missing_key(self):
        result = {"profiling": {"pfi": {}}}
        cols = {col.key: col for s in AUDIT_SECTIONS for col in s.columns}
        assert cols["pfi_score"].extractor(result) == "N/A"

    def test_std_steps_count_non_list_returns_na(self):
        result = {"standardization": {"result": {"steps_applied": "not_a_list"}}}
        cols = {col.key: col for s in AUDIT_SECTIONS for col in s.columns}
        assert cols["std_steps_count"].extractor(result) == "N/A"

    def test_std_excluded_fragments_non_list_returns_na(self):
        result = {"standardization": {"result": {"excluded_fragments": None}}}
        cols = {col.key: col for s in AUDIT_SECTIONS for col in s.columns}
        assert cols["std_excluded_fragments"].extractor(result) == "N/A"

    def test_std_stereo_changes_non_dict_returns_na(self):
        result = {"standardization": {"result": {"stereo_comparison": "n/a"}}}
        cols = {col.key: col for s in AUDIT_SECTIONS for col in s.columns}
        assert cols["std_stereo_changes"].extractor(result) == "N/A"

    def test_check_passed_with_empty_all_checks(self):
        result = {"validation": {"all_checks": []}}
        cols = {col.key: col for s in AUDIT_SECTIONS for col in s.columns}
        assert cols["parsability_passed"].extractor(result) == "N/A"

    def test_safety_all_passed_not_nested_under_sub_filter(self):
        # Confirm that all_passed is read from top of safety_filters, not a sub-dict
        result = {
            "scoring": {
                "safety_filters": {
                    "all_passed": True,
                    "pains": {"passed": False},  # individual filter can differ
                }
            }
        }
        cols = {col.key: col for s in AUDIT_SECTIONS for col in s.columns}
        assert cols["safety_all_passed"].extractor(result) == "Pass"
        assert cols["pains_passed"].extractor(result) == "Fail"

    def test_extract_flat_row_and_nested_same_count(self):
        flat = extract_flat_row(FULL_RESULT)
        nested = extract_nested(FULL_RESULT)
        flat_count = len(flat)
        nested_count = sum(len(v) for v in nested.values())
        assert flat_count == nested_count
