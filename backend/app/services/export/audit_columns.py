"""
Declarative Audit Column Registry

Defines ~91 audit columns across 6 sections (Validation, Deep Validation,
Scoring, Safety, Compound Profile, Standardization). All exporters (CSV,
Excel, JSON, SDF, PDF) import from this module as the single source of truth.
"""

from collections import OrderedDict
from dataclasses import dataclass
from typing import Any, Callable, Dict, List


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class AuditColumn:
    """A single audit column with key, display header, and extractor callable."""

    key: str
    header: str
    extractor: Callable[[Dict[str, Any]], Any]


@dataclass
class AuditSection:
    """A group of related audit columns with a display name and CSV/Excel prefix."""

    name: str
    prefix: str
    columns: List[AuditColumn]


# ---------------------------------------------------------------------------
# Private extractor factory helpers
# ---------------------------------------------------------------------------


def _check_passed(check_name: str) -> Callable[[Dict[str, Any]], str]:
    """Return an extractor that looks up *check_name* in validation.all_checks."""

    def extractor(result: Dict[str, Any]) -> str:
        validation = result.get("validation") or {}
        all_checks = validation.get("all_checks", [])
        for check in all_checks:
            if check.get("check_name") == check_name:
                passed = check.get("passed")
                if passed is True:
                    return "Pass"
                elif passed is False:
                    return "Fail"
                return "N/A"
        return "N/A"

    return extractor


def _nested(*path: str, default: Any = "N/A") -> Callable[[Dict[str, Any]], Any]:
    """Return an extractor that safely traverses nested dict keys."""

    def extractor(result: Dict[str, Any]) -> Any:
        node: Any = result
        for key in path:
            if not isinstance(node, dict):
                return default
            node = node.get(key)
            if node is None:
                return default
        return node

    return extractor


def _pass_fail(*path: str) -> Callable[[Dict[str, Any]], str]:
    """Like _nested but converts True->'Pass', False->'Fail', None/missing->'N/A'."""

    def extractor(result: Dict[str, Any]) -> str:
        node: Any = result
        for key in path:
            if not isinstance(node, dict):
                return "N/A"
            node = node.get(key)
            if node is None:
                return "N/A"
        if node is True:
            return "Pass"
        elif node is False:
            return "Fail"
        return "N/A"

    return extractor


def _safety_filter(filter_name: str, field: str) -> Callable[[Dict[str, Any]], Any]:
    """Extract result["scoring"]["safety_filters"][filter_name][field]; bool->Pass/Fail."""

    def extractor(result: Dict[str, Any]) -> Any:
        scoring = result.get("scoring") or {}
        safety_filters = scoring.get("safety_filters") or {}
        sub = safety_filters.get(filter_name) or {}
        if not isinstance(sub, dict):
            return "N/A"
        value = sub.get(field)
        if value is None:
            return "N/A"
        if isinstance(value, bool):
            return "Pass" if value else "Fail"
        return value

    return extractor


def _safety_assess(category: str, field: str) -> Callable[[Dict[str, Any]], Any]:
    """Extract result["safety_assessment"][category][field]; bool->Pass/Fail."""

    def extractor(result: Dict[str, Any]) -> Any:
        safety = result.get("safety_assessment") or {}
        sub = safety.get(category) or {}
        if not isinstance(sub, dict):
            return "N/A"
        value = sub.get(field)
        if value is None:
            return "N/A"
        if isinstance(value, bool):
            return "Pass" if value else "Fail"
        return value

    return extractor


def _profile(*path: str) -> Callable[[Dict[str, Any]], Any]:
    """Extract from result["profiling"][path...]."""

    def extractor(result: Dict[str, Any]) -> Any:
        node: Any = result.get("profiling") or {}
        for key in path:
            if not isinstance(node, dict):
                return "N/A"
            node = node.get(key)
            if node is None:
                return "N/A"
        return node

    return extractor


def _std_field(field: str) -> Callable[[Dict[str, Any]], Any]:
    """Extract result["standardization"]["result"][field]."""

    def extractor(result: Dict[str, Any]) -> Any:
        std = result.get("standardization") or {}
        res = std.get("result") or {}
        value = res.get(field)
        return value if value is not None else "N/A"

    return extractor


def _pass_fail_std(field: str) -> Callable[[Dict[str, Any]], str]:
    """Like _std_field but converts bool to Pass/Fail."""

    def extractor(result: Dict[str, Any]) -> str:
        std = result.get("standardization") or {}
        res = std.get("result") or {}
        value = res.get(field)
        if value is True:
            return "Pass"
        elif value is False:
            return "Fail"
        return "N/A"

    return extractor


# ---------------------------------------------------------------------------
# Custom extractors (one-off logic)
# ---------------------------------------------------------------------------


def _extract_standardized_smiles(result: Dict[str, Any]) -> str:
    """Prefer top-level standardized_smiles (set by batch), fall back to nested."""
    top = result.get("standardized_smiles")
    if top:
        return str(top)
    std = result.get("standardization") or {}
    res = std.get("result") or {}
    return res.get("standardized_smiles", "")


def _extract_std_steps_count(result: Dict[str, Any]) -> Any:
    """Count standardization steps where applied=True."""
    std = result.get("standardization") or {}
    res = std.get("result") or {}
    steps = res.get("steps_applied")
    if not isinstance(steps, list):
        return "N/A"
    return sum(1 for s in steps if isinstance(s, dict) and s.get("applied") is True)


def _extract_std_excluded_fragments(result: Dict[str, Any]) -> Any:
    """Count excluded fragments list length."""
    std = result.get("standardization") or {}
    res = std.get("result") or {}
    frags = res.get("excluded_fragments")
    if not isinstance(frags, list):
        return "N/A"
    return len(frags)


def _extract_std_stereo_changes(result: Dict[str, Any]) -> str:
    """Summarize stereo_comparison as 'lost: X, gained: Y' or 'N/A'."""
    std = result.get("standardization") or {}
    res = std.get("result") or {}
    stereo = res.get("stereo_comparison")
    if not isinstance(stereo, dict):
        return "N/A"
    lost = stereo.get("lost", 0)
    gained = stereo.get("gained", 0)
    return f"lost: {lost}, gained: {gained}"


def _extract_bro5_passed(result: Dict[str, Any]) -> str:
    safety = result.get("safety_assessment") or {}
    bro5 = safety.get("bro5") or {}
    val = bro5.get("passed")
    if val is True:
        return "Pass"
    elif val is False:
        return "Fail"
    return "N/A"


def _extract_reos_passed(result: Dict[str, Any]) -> str:
    safety = result.get("safety_assessment") or {}
    reos = safety.get("reos") or {}
    val = reos.get("passed")
    if val is True:
        return "Pass"
    elif val is False:
        return "Fail"
    return "N/A"


# ---------------------------------------------------------------------------
# Section 1: Validation (~17 columns)
# ---------------------------------------------------------------------------

_VALIDATION_COLUMNS: List[AuditColumn] = [
    AuditColumn(
        key="canonical_smiles",
        header="Canonical SMILES",
        extractor=_nested("validation", "canonical_smiles", default=""),
    ),
    AuditColumn(
        key="inchikey",
        header="InChIKey",
        extractor=_nested("validation", "inchikey", default=""),
    ),
    AuditColumn(
        key="inchi",
        header="InChI",
        extractor=_nested("validation", "molecule_info", "inchi", default=""),
    ),
    AuditColumn(
        key="molecular_formula",
        header="Molecular Formula",
        extractor=_nested("validation", "molecule_info", "molecular_formula", default=""),
    ),
    AuditColumn(
        key="molecular_weight",
        header="Molecular Weight",
        extractor=_nested("validation", "molecule_info", "molecular_weight"),
    ),
    AuditColumn(
        key="overall_score",
        header="Overall Score (0-100)",
        extractor=_nested("validation", "overall_score"),
    ),
    # 11 basic checks
    AuditColumn(
        key="parsability_passed",
        header="Parsability (Pass/Fail)",
        extractor=_check_passed("parsability"),
    ),
    AuditColumn(
        key="sanitization_passed",
        header="Sanitization (Pass/Fail)",
        extractor=_check_passed("sanitization"),
    ),
    AuditColumn(
        key="valence_passed",
        header="Valence (Pass/Fail)",
        extractor=_check_passed("valence"),
    ),
    AuditColumn(
        key="aromaticity_passed",
        header="Aromaticity (Pass/Fail)",
        extractor=_check_passed("aromaticity"),
    ),
    AuditColumn(
        key="connectivity_passed",
        header="Connectivity (Pass/Fail)",
        extractor=_check_passed("connectivity"),
    ),
    AuditColumn(
        key="undefined_stereocenters_passed",
        header="Undefined Stereocenters (Pass/Fail)",
        extractor=_check_passed("undefined_stereocenters"),
    ),
    AuditColumn(
        key="undefined_doublebond_stereo_passed",
        header="Undefined Doublebond Stereo (Pass/Fail)",
        extractor=_check_passed("undefined_doublebond_stereo"),
    ),
    AuditColumn(
        key="conflicting_stereo_passed",
        header="Conflicting Stereo (Pass/Fail)",
        extractor=_check_passed("conflicting_stereo"),
    ),
    AuditColumn(
        key="smiles_roundtrip_passed",
        header="SMILES Roundtrip (Pass/Fail)",
        extractor=_check_passed("smiles_roundtrip"),
    ),
    AuditColumn(
        key="inchi_generation_passed",
        header="InChI Generation (Pass/Fail)",
        extractor=_check_passed("inchi_generation"),
    ),
    AuditColumn(
        key="inchi_roundtrip_passed",
        header="InChI Roundtrip (Pass/Fail)",
        extractor=_check_passed("inchi_roundtrip"),
    ),
]

# ---------------------------------------------------------------------------
# Section 2: Deep Validation (~16 columns)
# ---------------------------------------------------------------------------

_DEEP_VALIDATION_COLUMNS: List[AuditColumn] = [
    AuditColumn(
        key="stereoisomer_enumeration_passed",
        header="Stereoisomer Enumeration (Pass/Fail)",
        extractor=_check_passed("stereoisomer_enumeration"),
    ),
    AuditColumn(
        key="tautomer_detection_passed",
        header="Tautomer Detection (Pass/Fail)",
        extractor=_check_passed("tautomer_detection"),
    ),
    AuditColumn(
        key="aromatic_system_validation_passed",
        header="Aromatic System Validation (Pass/Fail)",
        extractor=_check_passed("aromatic_system_validation"),
    ),
    AuditColumn(
        key="coordinate_dimension_passed",
        header="Coordinate Dimension (Pass/Fail)",
        extractor=_check_passed("coordinate_dimension"),
    ),
    AuditColumn(
        key="mixture_detection_passed",
        header="Mixture Detection (Pass/Fail)",
        extractor=_check_passed("mixture_detection"),
    ),
    AuditColumn(
        key="solvent_contamination_passed",
        header="Solvent Contamination (Pass/Fail)",
        extractor=_check_passed("solvent_contamination"),
    ),
    AuditColumn(
        key="inorganic_filter_passed",
        header="Inorganic Filter (Pass/Fail)",
        extractor=_check_passed("inorganic_filter"),
    ),
    AuditColumn(
        key="radical_detection_passed",
        header="Radical Detection (Pass/Fail)",
        extractor=_check_passed("radical_detection"),
    ),
    AuditColumn(
        key="isotope_label_detection_passed",
        header="Isotope Label Detection (Pass/Fail)",
        extractor=_check_passed("isotope_label_detection"),
    ),
    AuditColumn(
        key="trivial_molecule_passed",
        header="Trivial Molecule (Pass/Fail)",
        extractor=_check_passed("trivial_molecule"),
    ),
    AuditColumn(
        key="hypervalent_atoms_passed",
        header="Hypervalent Atoms (Pass/Fail)",
        extractor=_check_passed("hypervalent_atoms"),
    ),
    AuditColumn(
        key="polymer_detection_passed",
        header="Polymer Detection (Pass/Fail)",
        extractor=_check_passed("polymer_detection"),
    ),
    AuditColumn(
        key="ring_strain_passed",
        header="Ring Strain (Pass/Fail)",
        extractor=_check_passed("ring_strain"),
    ),
    AuditColumn(
        key="macrocycle_detection_passed",
        header="Macrocycle Detection (Pass/Fail)",
        extractor=_check_passed("macrocycle_detection"),
    ),
    AuditColumn(
        key="charged_species_passed",
        header="Charged Species (Pass/Fail)",
        extractor=_check_passed("charged_species"),
    ),
    AuditColumn(
        key="explicit_hydrogen_audit_passed",
        header="Explicit Hydrogen Audit (Pass/Fail)",
        extractor=_check_passed("explicit_hydrogen_audit"),
    ),
]

# ---------------------------------------------------------------------------
# Section 3: Scoring (~35 columns)
# ---------------------------------------------------------------------------

_SCORING_COLUMNS: List[AuditColumn] = [
    AuditColumn(
        key="ml_readiness_score",
        header="ML-Readiness Score (0-100)",
        extractor=_nested("scoring", "ml_readiness_score"),
    ),
    AuditColumn(
        key="ml_readiness_label",
        header="ML-Readiness Label",
        extractor=_nested("scoring", "ml_readiness", "label"),
    ),
    AuditColumn(
        key="qed_score",
        header="QED (0-1)",
        extractor=_nested("scoring", "druglikeness", "qed_score"),
    ),
    AuditColumn(
        key="lipinski_passed",
        header="Lipinski (Pass/Fail)",
        extractor=_pass_fail("scoring", "druglikeness", "lipinski", "passed"),
    ),
    AuditColumn(
        key="lipinski_violations",
        header="Lipinski Violations (0-4)",
        extractor=_nested("scoring", "druglikeness", "lipinski", "violations"),
    ),
    AuditColumn(
        key="lipinski_mw",
        header="Lipinski MW",
        extractor=_nested("scoring", "druglikeness", "lipinski", "mw"),
    ),
    AuditColumn(
        key="lipinski_logp",
        header="Lipinski LogP",
        extractor=_nested("scoring", "druglikeness", "lipinski", "logp"),
    ),
    AuditColumn(
        key="lipinski_hbd",
        header="Lipinski HBD",
        extractor=_nested("scoring", "druglikeness", "lipinski", "hbd"),
    ),
    AuditColumn(
        key="lipinski_hba",
        header="Lipinski HBA",
        extractor=_nested("scoring", "druglikeness", "lipinski", "hba"),
    ),
    AuditColumn(
        key="veber_passed",
        header="Veber (Pass/Fail)",
        extractor=_pass_fail("scoring", "druglikeness", "veber", "passed"),
    ),
    AuditColumn(
        key="veber_rotatable_bonds",
        header="Veber Rotatable Bonds",
        extractor=_nested("scoring", "druglikeness", "veber", "rotatable_bonds"),
    ),
    AuditColumn(
        key="veber_tpsa",
        header="Veber TPSA",
        extractor=_nested("scoring", "druglikeness", "veber", "tpsa"),
    ),
    AuditColumn(
        key="ro3_passed",
        header="Ro3 (Pass/Fail)",
        extractor=_pass_fail("scoring", "druglikeness", "ro3", "passed"),
    ),
    AuditColumn(
        key="ghose_passed",
        header="Ghose (Pass/Fail)",
        extractor=_pass_fail("scoring", "druglikeness", "ghose", "passed"),
    ),
    AuditColumn(
        key="egan_passed",
        header="Egan (Pass/Fail)",
        extractor=_pass_fail("scoring", "druglikeness", "egan", "passed"),
    ),
    AuditColumn(
        key="muegge_passed",
        header="Muegge (Pass/Fail)",
        extractor=_pass_fail("scoring", "druglikeness", "muegge", "passed"),
    ),
    AuditColumn(
        key="np_likeness_score",
        header="NP-Likeness (-5 to +5)",
        extractor=_nested("scoring", "np_likeness_score"),
    ),
    AuditColumn(
        key="consensus_score",
        header="Consensus Drug-Likeness (0-5)",
        extractor=_nested("scoring", "consensus", "score"),
    ),
    AuditColumn(
        key="lead_likeness_passed",
        header="Lead-Likeness (Pass/Fail)",
        extractor=_pass_fail("scoring", "lead_likeness", "passed"),
    ),
    AuditColumn(
        key="sa_score",
        header="SA Score (1-10)",
        extractor=_nested("scoring", "admet", "synthetic_accessibility", "score"),
    ),
    AuditColumn(
        key="sa_classification",
        header="SA Classification",
        extractor=_nested("scoring", "admet", "synthetic_accessibility", "classification"),
    ),
    AuditColumn(
        key="esol_log_s",
        header="ESOL Log S",
        extractor=_nested("scoring", "admet", "solubility", "log_s"),
    ),
    AuditColumn(
        key="esol_classification",
        header="ESOL Classification",
        extractor=_nested("scoring", "admet", "solubility", "classification"),
    ),
    AuditColumn(
        key="fsp3",
        header="Fsp3 (0-1)",
        extractor=_nested("scoring", "admet", "complexity", "fsp3"),
    ),
    AuditColumn(
        key="cns_mpo",
        header="CNS MPO (0-6)",
        extractor=_nested("scoring", "admet", "cns_mpo", "score"),
    ),
    AuditColumn(
        key="oral_absorption",
        header="Oral Absorption Likely (Pass/Fail)",
        extractor=_pass_fail("scoring", "admet", "bioavailability", "oral_absorption_likely"),
    ),
    AuditColumn(
        key="pfizer_rule",
        header="Pfizer Rule (Pass/Fail)",
        extractor=_pass_fail("scoring", "admet", "pfizer_rule", "passed"),
    ),
    AuditColumn(
        key="gsk_rule",
        header="GSK Rule (Pass/Fail)",
        extractor=_pass_fail("scoring", "admet", "gsk_rule", "passed"),
    ),
    AuditColumn(
        key="golden_triangle",
        header="Golden Triangle (Pass/Fail)",
        extractor=_pass_fail("scoring", "admet", "golden_triangle", "in_golden_triangle"),
    ),
    AuditColumn(
        key="aggregator_likelihood",
        header="Aggregator Likelihood",
        extractor=_nested("scoring", "aggregator", "likelihood"),
    ),
    AuditColumn(
        key="aggregator_risk",
        header="Aggregator Risk Score (0-1)",
        extractor=_nested("scoring", "aggregator", "risk_score"),
    ),
    AuditColumn(
        key="scaffold_smiles",
        header="Murcko Scaffold SMILES",
        extractor=_nested("scoring", "scaffold", "scaffold_smiles", default=""),
    ),
    AuditColumn(
        key="boiled_egg_gi",
        header="GI Absorption (Pass/Fail)",
        extractor=_pass_fail("scoring", "boiled_egg", "gi_absorbed"),
    ),
    AuditColumn(
        key="boiled_egg_bbb",
        header="BBB Permeant (Pass/Fail)",
        extractor=_pass_fail("scoring", "boiled_egg", "bbb_permeant"),
    ),
    AuditColumn(
        key="boiled_egg_region",
        header="BOILED-Egg Region",
        extractor=_nested("scoring", "boiled_egg", "region"),
    ),
]

# ---------------------------------------------------------------------------
# Section 4: Safety (~16 columns)
# ---------------------------------------------------------------------------

_SAFETY_COLUMNS: List[AuditColumn] = [
    AuditColumn(
        key="pains_passed",
        header="PAINS (Pass/Fail)",
        extractor=_safety_filter("pains", "passed"),
    ),
    AuditColumn(
        key="pains_count",
        header="PAINS Alert Count",
        extractor=_safety_filter("pains", "alert_count"),
    ),
    AuditColumn(
        key="brenk_passed",
        header="BRENK (Pass/Fail)",
        extractor=_safety_filter("brenk", "passed"),
    ),
    AuditColumn(
        key="brenk_count",
        header="BRENK Alert Count",
        extractor=_safety_filter("brenk", "alert_count"),
    ),
    AuditColumn(
        key="nih_passed",
        header="NIH (Pass/Fail)",
        extractor=_safety_filter("nih", "passed"),
    ),
    AuditColumn(
        key="zinc_passed",
        header="ZINC (Pass/Fail)",
        extractor=_safety_filter("zinc", "passed"),
    ),
    AuditColumn(
        key="chembl_passed",
        header="ChEMBL (Pass/Fail)",
        extractor=_pass_fail("scoring", "safety_filters", "chembl", "passed"),
    ),
    AuditColumn(
        key="safety_all_passed",
        header="All Safety Filters (Pass/Fail)",
        extractor=_pass_fail("scoring", "safety_filters", "all_passed"),
    ),
    AuditColumn(
        key="total_alerts",
        header="Total Structural Alerts",
        extractor=_nested("scoring", "safety_filters", "total_alerts"),
    ),
    AuditColumn(
        key="herg_risk",
        header="hERG Risk (low/moderate/high)",
        extractor=_safety_assess("herg", "herg_risk"),
    ),
    AuditColumn(
        key="herg_risk_score",
        header="hERG Risk Score (0-4)",
        extractor=_safety_assess("herg", "risk_score"),
    ),
    AuditColumn(
        key="bro5_passed",
        header="bRo5 (Pass/Fail)",
        extractor=_extract_bro5_passed,
    ),
    AuditColumn(
        key="reos_passed",
        header="REOS (Pass/Fail)",
        extractor=_extract_reos_passed,
    ),
    AuditColumn(
        key="reos_violations",
        header="REOS Violations",
        extractor=_safety_assess("reos", "n_violations"),
    ),
    AuditColumn(
        key="cyp_softspot_count",
        header="CYP Soft Spots",
        extractor=_safety_assess("cyp_softspots", "n_sites"),
    ),
    AuditColumn(
        key="complexity_outliers",
        header="Complexity Outliers",
        extractor=_safety_assess("complexity", "n_outliers"),
    ),
]

# ---------------------------------------------------------------------------
# Section 5: Compound Profile (~12 columns)
# ---------------------------------------------------------------------------

_COMPOUND_PROFILE_COLUMNS: List[AuditColumn] = [
    AuditColumn(
        key="pfi_score",
        header="PFI (0-15+)",
        extractor=_profile("pfi", "pfi"),
    ),
    AuditColumn(
        key="pfi_risk",
        header="PFI Risk (low/moderate/high)",
        extractor=_profile("pfi", "risk"),
    ),
    AuditColumn(
        key="stars_count",
        header="Stars (0-8)",
        extractor=_profile("stars", "stars"),
    ),
    AuditColumn(
        key="abbott_score",
        header="Abbott Bioavailability Score",
        extractor=_profile("abbott", "abbott_score"),
    ),
    AuditColumn(
        key="abbott_probability",
        header="Abbott Probability (%)",
        extractor=_profile("abbott", "probability_pct"),
    ),
    AuditColumn(
        key="consensus_logp",
        header="Consensus LogP",
        extractor=_profile("consensus_logp", "consensus_logp"),
    ),
    AuditColumn(
        key="skin_perm_logkp",
        header="Skin Permeation Log Kp",
        extractor=_profile("skin_permeation", "log_kp"),
    ),
    AuditColumn(
        key="skin_perm_class",
        header="Skin Permeation Class",
        extractor=_profile("skin_permeation", "classification"),
    ),
    AuditColumn(
        key="sa_comparison_sa",
        header="SA Score Comparison (1-10)",
        extractor=_profile("sa_comparison", "sa_score", "score"),
    ),
    AuditColumn(
        key="sa_comparison_scscore",
        header="SCScore Comparison",
        extractor=_profile("sa_comparison", "scscore", "score"),
    ),
    AuditColumn(
        key="sa_comparison_syba",
        header="SYBA Comparison",
        extractor=_profile("sa_comparison", "syba", "score"),
    ),
    AuditColumn(
        key="profile_cns_mpo",
        header="Profile CNS MPO (0-4)",
        extractor=_profile("cns_mpo", "score"),
    ),
]

# ---------------------------------------------------------------------------
# Section 6: Standardization (~6 columns)
# ---------------------------------------------------------------------------

_STANDARDIZATION_COLUMNS: List[AuditColumn] = [
    AuditColumn(
        key="standardized_smiles",
        header="Standardized SMILES",
        extractor=_extract_standardized_smiles,
    ),
    AuditColumn(
        key="std_success",
        header="Success (Pass/Fail)",
        extractor=_pass_fail_std("success"),
    ),
    AuditColumn(
        key="std_steps_count",
        header="Steps Applied",
        extractor=_extract_std_steps_count,
    ),
    AuditColumn(
        key="std_excluded_fragments",
        header="Excluded Fragments",
        extractor=_extract_std_excluded_fragments,
    ),
    AuditColumn(
        key="std_stereo_changes",
        header="Stereo Changes",
        extractor=_extract_std_stereo_changes,
    ),
    AuditColumn(
        key="std_mass_change",
        header="Mass Change (%)",
        extractor=_std_field("mass_change_percent"),
    ),
]

# ---------------------------------------------------------------------------
# Master registry: AUDIT_SECTIONS
# ---------------------------------------------------------------------------

AUDIT_SECTIONS: List[AuditSection] = [
    AuditSection(
        name="Validation",
        prefix="[Validation]",
        columns=_VALIDATION_COLUMNS,
    ),
    AuditSection(
        name="Deep Validation",
        prefix="[Deep Validation]",
        columns=_DEEP_VALIDATION_COLUMNS,
    ),
    AuditSection(
        name="Scoring",
        prefix="[Scoring]",
        columns=_SCORING_COLUMNS,
    ),
    AuditSection(
        name="Safety",
        prefix="[Safety]",
        columns=_SAFETY_COLUMNS,
    ),
    AuditSection(
        name="Compound Profile",
        prefix="[Compound Profile]",
        columns=_COMPOUND_PROFILE_COLUMNS,
    ),
    AuditSection(
        name="Standardization",
        prefix="[Standardization]",
        columns=_STANDARDIZATION_COLUMNS,
    ),
]

# ---------------------------------------------------------------------------
# Public output helper functions
# ---------------------------------------------------------------------------


def get_identity_row(idx: int, result: Dict[str, Any]) -> OrderedDict:
    """Return identity columns: index (1-based), name, input_smiles.

    Args:
        idx: 0-based index of the result in the batch.
        result: Single result dictionary from batch processing.

    Returns:
        OrderedDict with keys: index, name, input_smiles.
    """
    return OrderedDict(
        [
            ("index", idx + 1),
            ("name", result.get("name", "")),
            ("input_smiles", result.get("smiles", "")),
        ]
    )


def get_flat_headers() -> List[str]:
    """Return a list of all prefixed column headers in section order.

    Returns:
        List of strings like '[Validation] Canonical SMILES'.
    """
    headers: List[str] = []
    for section in AUDIT_SECTIONS:
        for col in section.columns:
            headers.append(f"{section.prefix} {col.header}")
    return headers


def extract_flat_row(result: Dict[str, Any]) -> OrderedDict:
    """Extract all audit columns as a flat OrderedDict with prefixed headers.

    Args:
        result: Single result dictionary from batch processing.

    Returns:
        OrderedDict keyed by prefixed headers in section order.
    """
    row: OrderedDict = OrderedDict()
    for section in AUDIT_SECTIONS:
        for col in section.columns:
            header = f"{section.prefix} {col.header}"
            row[header] = col.extractor(result)
    return row


def extract_nested(result: Dict[str, Any]) -> Dict[str, Any]:
    """Extract all audit columns as a nested dict keyed by section snake_case name.

    Args:
        result: Single result dictionary from batch processing.

    Returns:
        Dict[section_snake_case, Dict[col_key, value]].
    """
    nested: Dict[str, Any] = {}
    for section in AUDIT_SECTIONS:
        section_key = section.name.lower().replace(" ", "_")
        nested[section_key] = {}
        for col in section.columns:
            nested[section_key][col.key] = col.extractor(result)
    return nested


def extract_by_section(result: Dict[str, Any]) -> Dict[str, OrderedDict]:
    """Extract audit columns grouped by section name, with headers as keys.

    Args:
        result: Single result dictionary from batch processing.

    Returns:
        Dict[section_name, OrderedDict[header, value]].
    """
    by_section: Dict[str, OrderedDict] = {}
    for section in AUDIT_SECTIONS:
        data: OrderedDict = OrderedDict()
        for col in section.columns:
            data[col.header] = col.extractor(result)
        by_section[section.name] = data
    return by_section
