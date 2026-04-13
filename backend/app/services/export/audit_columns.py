"""
Declarative Audit Column Registry

Defines ~78 audit columns across 6 sections (Validation, Deep Validation,
Scoring, Safety, Compound Profile, Standardization). All exporters (CSV,
Excel, JSON, SDF, PDF) import from this module as the single source of truth.

Extractor paths match the batch result dict built by
``app.services.batch.tasks._process_single_molecule()``.
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


def _safety_filter_flat(field: str) -> Callable[[Dict[str, Any]], Any]:
    """Extract result['scoring']['safety_filters'][field]; bool->Pass/Fail.

    The batch processor stores safety filter results as flat keys
    (e.g. ``pains_passed``, ``brenk_passed``) not sub-dicts.
    """

    def extractor(result: Dict[str, Any]) -> Any:
        scoring = result.get("scoring") or {}
        sf = scoring.get("safety_filters") or {}
        if not isinstance(sf, dict):
            return "N/A"
        value = sf.get(field)
        if value is None:
            return "N/A"
        if isinstance(value, bool):
            return "Pass" if value else "Fail"
        return value

    return extractor


def _alert_count(catalog_upper: str) -> Callable[[Dict[str, Any]], Any]:
    """Count alerts from result['alerts']['alerts'] matching a catalog name."""

    def extractor(result: Dict[str, Any]) -> Any:
        alerts = result.get("alerts") or {}
        alert_list = alerts.get("alerts") or []
        if not isinstance(alert_list, list):
            return 0
        return sum(
            1
            for a in alert_list
            if isinstance(a, dict) and str(a.get("catalog", "")).upper() == catalog_upper
        )

    return extractor


def _safety_assess(category: str, field: str) -> Callable[[Dict[str, Any]], Any]:
    """Extract result['safety_assessment'][category][field]; bool->Pass/Fail."""

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
    """Extract from result['profiling'][path...]."""

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


# ---------------------------------------------------------------------------
# Custom extractors (one-off logic)
# ---------------------------------------------------------------------------


def _extract_standardized_smiles(result: Dict[str, Any]) -> str:
    """Get standardized SMILES from standardization dict (no 'result' sub-key)."""
    std = result.get("standardization") or {}
    return std.get("standardized_smiles") or ""


def _extract_std_steps_count(result: Dict[str, Any]) -> Any:
    """Count standardization steps where applied=True (no 'result' sub-key)."""
    std = result.get("standardization") or {}
    steps = std.get("steps_applied")
    if not isinstance(steps, list):
        return "N/A"
    return sum(1 for s in steps if isinstance(s, dict) and s.get("applied") is True)


def _extract_std_excluded_fragments(result: Dict[str, Any]) -> Any:
    """Count excluded fragments (no 'result' sub-key)."""
    std = result.get("standardization") or {}
    frags = std.get("excluded_fragments")
    if not isinstance(frags, list):
        return "N/A"
    return len(frags)


_extract_bro5_passed = _pass_fail("safety_assessment", "bro5", "passed")
_extract_reos_passed = _pass_fail("safety_assessment", "reos", "passed")


def _extract_cyp_count(result: Dict[str, Any]) -> Any:
    """Count CYP soft-spot matches (cyp_softspots is a list, not a dict)."""
    safety = result.get("safety_assessment") or {}
    cyp = safety.get("cyp_softspots")
    if not isinstance(cyp, list):
        return "N/A"
    return len(cyp)


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
# Section 3: Scoring (~17 columns)
#
# Only includes fields the batch processor actually computes.
# Batch uses flat keys (lipinski_passed, mw, sa_score, etc.) not nested
# sub-dicts.  Extended druglikeness (ro3/ghose/egan/muegge), NP-likeness,
# consensus, lead-likeness, aggregator, scaffold, boiled-egg, and advanced
# ADMET (cns_mpo, pfizer_rule, gsk_rule, golden_triangle, bioavailability)
# are not computed in batch mode and are therefore excluded.
# ---------------------------------------------------------------------------

_SCORING_COLUMNS: List[AuditColumn] = [
    AuditColumn(
        key="ml_readiness_score",
        header="ML-Readiness Score (0-100)",
        extractor=_nested("scoring", "ml_readiness", "score"),
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
        extractor=_pass_fail("scoring", "druglikeness", "lipinski_passed"),
    ),
    AuditColumn(
        key="lipinski_violations",
        header="Lipinski Violations (0-4)",
        extractor=_nested("scoring", "druglikeness", "lipinski_violations"),
    ),
    AuditColumn(
        key="lipinski_mw",
        header="Lipinski MW",
        extractor=_nested("scoring", "druglikeness", "mw"),
    ),
    AuditColumn(
        key="lipinski_logp",
        header="Lipinski LogP",
        extractor=_nested("scoring", "druglikeness", "logp"),
    ),
    AuditColumn(
        key="lipinski_hbd",
        header="Lipinski HBD",
        extractor=_nested("scoring", "druglikeness", "hbd"),
    ),
    AuditColumn(
        key="lipinski_hba",
        header="Lipinski HBA",
        extractor=_nested("scoring", "druglikeness", "hba"),
    ),
    AuditColumn(
        key="veber_passed",
        header="Veber (Pass/Fail)",
        extractor=_pass_fail("scoring", "druglikeness", "veber_passed"),
    ),
    AuditColumn(
        key="veber_rotatable_bonds",
        header="Veber Rotatable Bonds",
        extractor=_nested("scoring", "druglikeness", "rotatable_bonds"),
    ),
    AuditColumn(
        key="veber_tpsa",
        header="Veber TPSA",
        extractor=_nested("scoring", "druglikeness", "tpsa"),
    ),
    AuditColumn(
        key="aromatic_rings",
        header="Aromatic Rings",
        extractor=_nested("scoring", "druglikeness", "aromatic_rings"),
    ),
    AuditColumn(
        key="sa_score",
        header="SA Score (1-10)",
        extractor=_nested("scoring", "admet", "sa_score"),
    ),
    AuditColumn(
        key="sa_classification",
        header="SA Classification",
        extractor=_nested("scoring", "admet", "sa_classification"),
    ),
    AuditColumn(
        key="solubility_class",
        header="Solubility Classification",
        extractor=_nested("scoring", "admet", "solubility_class"),
    ),
    AuditColumn(
        key="fsp3",
        header="Fsp3 (0-1)",
        extractor=_nested("scoring", "admet", "fsp3"),
    ),
]

# ---------------------------------------------------------------------------
# Section 4: Safety (~16 columns)
#
# Safety filter results use flat keys (pains_passed, brenk_passed, etc.)
# not sub-dicts.  Alert counts are derived from result["alerts"]["alerts"].
# Safety assessment (hERG, bRo5, REOS, CYP, complexity) comes from the
# optional "safety_assessment" top-level key.
# ---------------------------------------------------------------------------

_SAFETY_COLUMNS: List[AuditColumn] = [
    AuditColumn(
        key="pains_passed",
        header="PAINS (Pass/Fail)",
        extractor=_safety_filter_flat("pains_passed"),
    ),
    AuditColumn(
        key="pains_count",
        header="PAINS Alert Count",
        extractor=_alert_count("PAINS"),
    ),
    AuditColumn(
        key="brenk_passed",
        header="BRENK (Pass/Fail)",
        extractor=_safety_filter_flat("brenk_passed"),
    ),
    AuditColumn(
        key="brenk_count",
        header="BRENK Alert Count",
        extractor=_alert_count("BRENK"),
    ),
    AuditColumn(
        key="nih_passed",
        header="NIH (Pass/Fail)",
        extractor=_safety_filter_flat("nih_passed"),
    ),
    AuditColumn(
        key="zinc_passed",
        header="ZINC (Pass/Fail)",
        extractor=_safety_filter_flat("zinc_passed"),
    ),
    AuditColumn(
        key="chembl_passed",
        header="ChEMBL (Pass/Fail)",
        extractor=_safety_filter_flat("chembl_passed"),
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
        extractor=_extract_cyp_count,
    ),
    AuditColumn(
        key="complexity_outliers",
        header="Complexity Outliers",
        extractor=_safety_assess("complexity", "n_outliers"),
    ),
]

# ---------------------------------------------------------------------------
# Section 5: Compound Profile (~8 columns)
#
# ``compute_full_profile()`` provides: pfi, stars, abbott, consensus_logp,
# skin_permeation, and a druglikeness sub-dict.
# ``sa_comparison`` and ``cns_mpo`` are NOT computed and therefore excluded.
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
]

# ---------------------------------------------------------------------------
# Section 6: Standardization (~4 columns)
#
# The batch processor stores standardization fields directly on the
# ``standardization`` dict — there is NO ``result`` sub-key.
# ``stereo_comparison`` and ``mass_change_percent`` are not stored in batch
# and therefore excluded.
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
        extractor=_pass_fail("standardization", "success"),
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
