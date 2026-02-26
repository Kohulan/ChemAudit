"""
Scoring Schemas

Pydantic schemas for scoring requests and responses.
Includes drug-likeness, safety filters, and ADMET predictions.
"""

from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field, field_validator

# =============================================================================
# Drug-likeness Schemas
# =============================================================================


class LipinskiSchema(BaseModel):
    """Lipinski's Rule of Five results."""

    passed: bool = Field(description="Whether molecule passes Ro5 (<=1 violation)")
    violations: int = Field(ge=0, le=4, description="Number of rule violations")
    mw: float = Field(description="Molecular weight")
    logp: float = Field(description="Calculated LogP")
    hbd: int = Field(description="Hydrogen bond donors")
    hba: int = Field(description="Hydrogen bond acceptors")
    details: Dict[str, bool] = Field(
        default_factory=dict, description="Individual rule pass/fail"
    )


class QEDSchema(BaseModel):
    """QED score results."""

    score: float = Field(ge=0, le=1, description="QED score (0-1)")
    properties: Dict[str, Any] = Field(
        default_factory=dict, description="Component property values"
    )
    interpretation: str = Field(description="Human-readable interpretation")


class VeberSchema(BaseModel):
    """Veber rules results."""

    passed: bool = Field(description="Whether molecule passes Veber rules")
    rotatable_bonds: int = Field(description="Number of rotatable bonds")
    tpsa: float = Field(description="Topological polar surface area")


class RuleOfThreeSchema(BaseModel):
    """Rule of Three (fragment-likeness) results."""

    passed: bool = Field(description="Whether molecule passes Ro3")
    violations: int = Field(description="Number of rule violations")
    mw: float = Field(description="Molecular weight")
    logp: float = Field(description="Calculated LogP")
    hbd: int = Field(description="Hydrogen bond donors")
    hba: int = Field(description="Hydrogen bond acceptors")
    rotatable_bonds: int = Field(description="Number of rotatable bonds")
    tpsa: float = Field(description="Topological polar surface area")


class GhoseSchema(BaseModel):
    """Ghose filter results."""

    passed: bool = Field(description="Whether molecule passes Ghose filter")
    violations: int = Field(description="Number of violations")
    mw: float = Field(description="Molecular weight")
    logp: float = Field(description="Calculated LogP")
    atom_count: int = Field(description="Heavy atom count")
    molar_refractivity: float = Field(description="Molar refractivity")


class EganSchema(BaseModel):
    """Egan filter results."""

    passed: bool = Field(description="Whether molecule passes Egan filter")
    logp: float = Field(description="Calculated LogP")
    tpsa: float = Field(description="Topological polar surface area")


class MueggeSchema(BaseModel):
    """Muegge filter results."""

    passed: bool = Field(description="Whether molecule passes Muegge filter")
    violations: int = Field(description="Number of violations")
    details: Dict[str, bool] = Field(
        default_factory=dict, description="Individual criterion pass/fail"
    )


class DrugLikenessResultSchema(BaseModel):
    """Complete drug-likeness scoring results."""

    lipinski: LipinskiSchema = Field(description="Lipinski's Rule of Five")
    qed: QEDSchema = Field(description="QED score")
    veber: VeberSchema = Field(description="Veber rules")
    ro3: RuleOfThreeSchema = Field(description="Rule of Three (fragment-likeness)")
    ghose: Optional[GhoseSchema] = Field(None, description="Ghose filter")
    egan: Optional[EganSchema] = Field(None, description="Egan filter")
    muegge: Optional[MueggeSchema] = Field(None, description="Muegge filter")
    interpretation: str = Field(description="Overall drug-likeness interpretation")


# =============================================================================
# Safety Filters Schemas
# =============================================================================


class AlertDetailSchema(BaseModel):
    """Enriched detail for a single matched alert pattern."""

    name: str = Field(description="Pattern name from the filter catalog")
    category: str = Field(
        description="Functional group category: Reactive Group, Toxicophore, "
        "Metabolic Liability, Assay Interference, Physicochemical, or Unwanted Functionality"
    )


class FilterAlertSchema(BaseModel):
    """Result for a single filter category."""

    passed: bool = Field(description="Whether molecule passed this filter")
    alerts: List[str] = Field(default_factory=list, description="List of alert names")
    alert_details: List[AlertDetailSchema] = Field(
        default_factory=list,
        description="Enriched alert details with pattern name and category",
    )
    alert_count: int = Field(default=0, description="Number of alerts triggered")


class ChEMBLAlertsSchema(BaseModel):
    """ChEMBL structural alerts results."""

    passed: bool = Field(description="Whether all ChEMBL filters passed")
    total_alerts: int = Field(default=0, description="Total ChEMBL alerts")
    bms: Optional[FilterAlertSchema] = Field(None, description="BMS filter results")
    dundee: Optional[FilterAlertSchema] = Field(
        None, description="Dundee filter results"
    )
    glaxo: Optional[FilterAlertSchema] = Field(None, description="Glaxo filter results")
    inpharmatica: Optional[FilterAlertSchema] = Field(
        None, description="Inpharmatica filter results"
    )
    lint: Optional[FilterAlertSchema] = Field(None, description="LINT filter results")
    mlsmr: Optional[FilterAlertSchema] = Field(None, description="MLSMR filter results")
    schembl: Optional[FilterAlertSchema] = Field(
        None, description="SureChEMBL filter results"
    )


class SafetyFilterResultSchema(BaseModel):
    """Complete safety filter results."""

    pains: FilterAlertSchema = Field(description="PAINS filter results")
    brenk: FilterAlertSchema = Field(description="Brenk filter results")
    nih: Optional[FilterAlertSchema] = Field(None, description="NIH filter results")
    zinc: Optional[FilterAlertSchema] = Field(None, description="ZINC filter results")
    chembl: Optional[ChEMBLAlertsSchema] = Field(
        None, description="ChEMBL structural alerts"
    )
    all_passed: bool = Field(description="Whether all filters passed")
    total_alerts: int = Field(default=0, description="Total number of alerts")
    interpretation: str = Field(description="Safety assessment interpretation")


# =============================================================================
# ADMET Schemas
# =============================================================================


class SyntheticAccessibilitySchema(BaseModel):
    """Synthetic accessibility score result."""

    score: float = Field(ge=1, le=10, description="SA score (1=easy, 10=difficult)")
    classification: str = Field(description="easy, moderate, or difficult")
    interpretation: str = Field(description="Human-readable interpretation")


class SolubilitySchema(BaseModel):
    """ESOL solubility prediction result."""

    log_s: float = Field(description="Log solubility (mol/L)")
    solubility_mg_ml: float = Field(description="Solubility in mg/mL")
    classification: str = Field(
        description="highly_soluble, soluble, moderate, poor, or insoluble"
    )
    interpretation: str = Field(description="Human-readable interpretation")


class ComplexitySchema(BaseModel):
    """Molecular complexity metrics."""

    fsp3: float = Field(ge=0, le=1, description="Fraction of sp3 carbons")
    num_stereocenters: int = Field(description="Number of stereocenters")
    num_rings: int = Field(description="Total number of rings")
    num_aromatic_rings: int = Field(description="Number of aromatic rings")
    bertz_ct: float = Field(description="Bertz complexity index")
    classification: str = Field(description="flat, moderate, or 3d")
    interpretation: str = Field(description="Human-readable interpretation")


class CNSMPOSchema(BaseModel):
    """CNS MPO score result."""

    score: float = Field(ge=0, le=6, description="CNS MPO score (0-6)")
    components: Dict[str, float] = Field(
        default_factory=dict, description="Individual component scores"
    )
    cns_penetrant: bool = Field(description="Predicted CNS penetration")
    interpretation: str = Field(description="Human-readable interpretation")


class BioavailabilitySchema(BaseModel):
    """Bioavailability indicators."""

    tpsa: float = Field(description="Topological polar surface area")
    rotatable_bonds: int = Field(description="Number of rotatable bonds")
    hbd: int = Field(description="Hydrogen bond donors")
    hba: int = Field(description="Hydrogen bond acceptors")
    mw: float = Field(description="Molecular weight")
    logp: float = Field(description="Calculated LogP")
    oral_absorption_likely: bool = Field(description="Predicted oral absorption")
    cns_penetration_likely: bool = Field(description="Predicted CNS penetration")
    interpretation: str = Field(description="Human-readable interpretation")


class PfizerRuleSchema(BaseModel):
    """Pfizer 3/75 Rule result."""

    passed: bool = Field(description="Whether molecule passes Pfizer 3/75 rule")
    logp: float = Field(description="Calculated LogP")
    tpsa: float = Field(description="Topological polar surface area")
    interpretation: str = Field(description="Rule interpretation")


class GSKRuleSchema(BaseModel):
    """GSK 4/400 Rule result."""

    passed: bool = Field(description="Whether molecule passes GSK 4/400 rule")
    mw: float = Field(description="Molecular weight")
    logp: float = Field(description="Calculated LogP")
    interpretation: str = Field(description="Rule interpretation")


class GoldenTriangleSchema(BaseModel):
    """Golden Triangle (Abbott) analysis."""

    in_golden_triangle: bool = Field(
        description="Whether molecule is in Golden Triangle"
    )
    mw: float = Field(description="Molecular weight")
    logd: float = Field(description="LogD (using LogP as proxy)")
    interpretation: str = Field(description="Analysis interpretation")


class ADMETResultSchema(BaseModel):
    """Complete ADMET prediction results."""

    synthetic_accessibility: SyntheticAccessibilitySchema = Field(
        description="Synthetic accessibility score"
    )
    solubility: SolubilitySchema = Field(description="Solubility prediction")
    complexity: ComplexitySchema = Field(description="Molecular complexity metrics")
    cns_mpo: Optional[CNSMPOSchema] = Field(None, description="CNS MPO score")
    bioavailability: BioavailabilitySchema = Field(
        description="Bioavailability indicators"
    )
    pfizer_rule: Optional[PfizerRuleSchema] = Field(
        None, description="Pfizer 3/75 rule assessment"
    )
    gsk_rule: Optional[GSKRuleSchema] = Field(
        None, description="GSK 4/400 rule assessment"
    )
    golden_triangle: Optional[GoldenTriangleSchema] = Field(
        None, description="Golden Triangle analysis"
    )
    molar_refractivity: Optional[float] = Field(None, description="Molar refractivity")
    interpretation: str = Field(description="Overall ADMET interpretation")


# =============================================================================
# Aggregator Likelihood Schema
# =============================================================================


class AggregatorLikelihoodSchema(BaseModel):
    """Aggregator likelihood prediction result."""

    likelihood: str = Field(description="Aggregation likelihood: low, moderate, high")
    risk_score: float = Field(ge=0, le=1, description="Risk score (0-1)")
    logp: float = Field(description="Calculated LogP")
    tpsa: float = Field(description="Topological polar surface area")
    mw: float = Field(description="Molecular weight")
    aromatic_rings: int = Field(description="Number of aromatic rings")
    risk_factors: List[str] = Field(
        default_factory=list, description="Identified risk factors"
    )
    interpretation: str = Field(description="Interpretation of aggregation risk")
    confidence: float = Field(
        default=0.0, ge=0, le=1, description="Confidence score (0-1)"
    )
    evidence: List[Dict[str, Any]] = Field(
        default_factory=list, description="Evidence detail per indicator"
    )


# =============================================================================
# Consensus Drug-Likeness Schemas
# =============================================================================


class RuleViolationSchema(BaseModel):
    """Per-property result within a rule set."""

    property: str = Field(description="Property name")
    value: float = Field(description="Actual property value")
    threshold: str = Field(description="Threshold expression (e.g., '<=500', '160-480')")
    result: str = Field(description="'pass' or 'fail'")


class RuleSetDetailSchema(BaseModel):
    """Detail for a single rule set in consensus scoring."""

    name: str = Field(description="Rule set name (e.g., 'Lipinski', 'Veber')")
    passed: bool = Field(description="Whether molecule passes this rule set")
    violations: List[RuleViolationSchema] = Field(
        default_factory=list, description="Per-property results"
    )


class ConsensusScoreSchema(BaseModel):
    """Consensus drug-likeness score across 5 rule sets."""

    score: int = Field(ge=0, le=5, description="Number of rule sets passed (0-5)")
    total: int = Field(description="Total rule sets (always 5)")
    rule_sets: List[RuleSetDetailSchema] = Field(
        default_factory=list, description="Per-rule-set detail"
    )
    interpretation: str = Field(description="Human-readable interpretation")


# =============================================================================
# Lead-Likeness Schema
# =============================================================================


class LeadLikenessSchema(BaseModel):
    """Lead-likeness assessment result."""

    passed: bool = Field(description="Whether molecule is lead-like")
    violations: int = Field(description="Number of violations")
    properties: Dict[str, float] = Field(
        default_factory=dict, description="Actual property values"
    )
    thresholds: Dict[str, str] = Field(
        default_factory=dict, description="Threshold expressions"
    )
    violation_details: List[RuleViolationSchema] = Field(
        default_factory=list, description="Per-property violation detail"
    )


# =============================================================================
# Salt Inventory Schemas
# =============================================================================


class SaltFragmentSchema(BaseModel):
    """A single fragment from a salt-form molecule."""

    smiles: str = Field(description="Fragment SMILES")
    name: str = Field(description="Fragment name")
    category: str = Field(
        description="Category: counterion, salt, solvent, drug, unknown"
    )
    mw: float = Field(description="Fragment molecular weight")
    heavy_atom_count: int = Field(description="Number of heavy atoms")


class SaltInventorySchema(BaseModel):
    """Salt/counterion inventory result."""

    has_salts: bool = Field(description="Whether salt forms were detected")
    parent_smiles: str = Field(description="Parent compound SMILES")
    fragments: List[SaltFragmentSchema] = Field(
        default_factory=list, description="Classified fragments"
    )
    total_fragments: int = Field(description="Total fragment count")
    interpretation: str = Field(description="Inventory interpretation")


# =============================================================================
# Ligand Efficiency Schema
# =============================================================================


class LigandEfficiencySchema(BaseModel):
    """Ligand efficiency result."""

    le: Optional[float] = Field(None, description="Ligand efficiency value")
    heavy_atom_count: int = Field(description="Number of heavy atoms")
    activity_value: Optional[float] = Field(
        None, description="Activity value used for LE"
    )
    activity_type: Optional[str] = Field(
        None, description="Activity type (e.g., pIC50, pKi)"
    )
    proxy_used: bool = Field(description="Whether BEI proxy was used")
    interpretation: str = Field(description="Efficiency interpretation")


# =============================================================================
# Property Breakdown Schemas
# =============================================================================


class AtomContributionSchema(BaseModel):
    """Per-atom property contribution."""

    atom_index: int = Field(description="Atom index")
    symbol: str = Field(description="Atom symbol")
    contribution: float = Field(description="Property contribution value")


class FunctionalGroupContributionSchema(BaseModel):
    """Per-functional-group property contribution."""

    group_name: str = Field(description="Functional group name")
    contribution: float = Field(description="Total group contribution")
    atom_indices: List[int] = Field(
        default_factory=list, description="Atom indices in this group"
    )


class TPSABreakdownSchema(BaseModel):
    """TPSA per-atom breakdown."""

    total: float = Field(description="Total TPSA")
    atom_contributions: List[AtomContributionSchema] = Field(
        default_factory=list, description="Per-atom TPSA contributions"
    )
    functional_group_summary: List[FunctionalGroupContributionSchema] = Field(
        default_factory=list, description="Per-group TPSA summary"
    )


class LogPBreakdownSchema(BaseModel):
    """LogP per-atom breakdown."""

    total: float = Field(description="Total LogP")
    atom_contributions: List[AtomContributionSchema] = Field(
        default_factory=list, description="Per-atom LogP contributions"
    )
    functional_group_summary: List[FunctionalGroupContributionSchema] = Field(
        default_factory=list, description="Per-group LogP summary"
    )


class BertzDetailSchema(BaseModel):
    """Bertz complexity detail."""

    bertz_ct: float = Field(description="Bertz complexity index")
    num_bonds: int = Field(description="Number of bonds")
    num_atoms: int = Field(description="Number of heavy atoms")
    num_rings: int = Field(description="Number of rings")
    num_aromatic_rings: int = Field(description="Number of aromatic rings")
    ring_complexity: float = Field(description="Fraction of bonds in rings")
    interpretation: str = Field(description="Complexity interpretation")


class CarbonHybridizationSchema(BaseModel):
    """Per-carbon hybridization data."""

    atom_index: int = Field(description="Carbon atom index")
    symbol: str = Field(default="C", description="Atom symbol (always C)")
    hybridization: str = Field(description="Hybridization: sp, sp2, sp3, other")


class Fsp3DetailSchema(BaseModel):
    """Fsp3 per-carbon detail."""

    fsp3: float = Field(ge=0, le=1, description="Fraction of sp3 carbons")
    total_carbons: int = Field(description="Total carbon count")
    sp3_count: int = Field(description="Number of sp3 carbons")
    sp2_count: int = Field(description="Number of sp2 carbons")
    sp_count: int = Field(description="Number of sp carbons")
    per_carbon: List[CarbonHybridizationSchema] = Field(
        default_factory=list, description="Per-carbon hybridization"
    )
    interpretation: str = Field(description="3D character interpretation")


# =============================================================================
# NP-Likeness Breakdown Schemas
# =============================================================================


class NPFragmentSchema(BaseModel):
    """A single NP fragment contribution."""

    smiles: str = Field(description="Fragment SMILES")
    contribution: float = Field(description="Contribution to NP score")
    bit_id: int = Field(description="Morgan FP bit ID")
    radius: int = Field(description="Morgan FP radius")
    center_atom_idx: int = Field(description="Center atom index")
    classification: str = Field(
        description="np_characteristic, synthetic_characteristic, or neutral"
    )


class NPBreakdownSchema(BaseModel):
    """NP-likeness fragment breakdown."""

    score: float = Field(description="NP-likeness score")
    confidence: float = Field(ge=0, le=1, description="Confidence (0-1)")
    fragments: List[NPFragmentSchema] = Field(
        default_factory=list, description="Per-fragment contributions"
    )
    total_fragments: int = Field(description="Total fragment count")
    np_fragment_count: int = Field(
        description="Number of NP-characteristic fragments"
    )
    synthetic_fragment_count: int = Field(
        description="Number of synthetic-characteristic fragments"
    )
    interpretation: str = Field(description="NP-likeness interpretation")


# =============================================================================
# Bioavailability Radar & BOILED-Egg Schemas
# =============================================================================


class RadarAxisSchema(BaseModel):
    """A single axis of the bioavailability radar."""

    name: str = Field(description="Axis name (LIPO, SIZE, POLAR, etc.)")
    actual_value: float = Field(description="Actual property value")
    normalized: float = Field(ge=0, le=1, description="Normalized value (0-1)")
    optimal_min: float = Field(description="Optimal range minimum")
    optimal_max: float = Field(description="Optimal range maximum")
    in_range: bool = Field(description="Whether value is in optimal range")
    property_name: str = Field(description="Property name (WLOGP, MW, etc.)")
    unit: str = Field(description="Property unit")


class BioavailabilityRadarSchema(BaseModel):
    """Bioavailability radar result."""

    axes: List[RadarAxisSchema] = Field(
        default_factory=list, description="6 radar axes"
    )
    overall_in_range_count: int = Field(
        description="Number of axes in optimal range"
    )
    interpretation: str = Field(description="Bioavailability interpretation")


class EllipseParamsSchema(BaseModel):
    """Ellipse parameters for BOILED-Egg model."""

    cx: float = Field(description="Center x (TPSA axis)")
    cy: float = Field(description="Center y (WLOGP axis)")
    a: float = Field(description="Semi-axis x (TPSA)")
    b: float = Field(description="Semi-axis y (WLOGP)")


class BoiledEggSchema(BaseModel):
    """BOILED-Egg classification result."""

    wlogp: float = Field(description="Wildman-Crippen LogP")
    tpsa: float = Field(description="Topological polar surface area")
    gi_absorbed: bool = Field(description="Predicted GI absorption")
    bbb_permeant: bool = Field(description="Predicted BBB permeation")
    region: str = Field(description="Classification: yolk, white, or grey")
    gi_ellipse: Optional[EllipseParamsSchema] = Field(
        None, description="GI (white) ellipse parameters"
    )
    bbb_ellipse: Optional[EllipseParamsSchema] = Field(
        None, description="BBB (yolk) ellipse parameters"
    )
    interpretation: str = Field(description="Classification interpretation")


class RadarProfileSchema(BaseModel):
    """A single molecule's radar profile."""

    smiles: str = Field(description="Molecule SMILES")
    axes: List[RadarAxisSchema] = Field(
        default_factory=list, description="Radar axes"
    )
    is_reference: bool = Field(description="Whether this is the reference profile")


class RadarComparisonSchema(BaseModel):
    """Multi-molecule radar comparison."""

    profiles: List[RadarProfileSchema] = Field(
        default_factory=list, description="Molecule profiles"
    )
    reference: Optional[RadarProfileSchema] = Field(
        None, description="Reference profile"
    )


class ComparisonRequest(BaseModel):
    """Request for molecule comparison."""

    smiles_list: List[str] = Field(
        ...,
        min_length=1,
        max_length=10,
        description="List of SMILES to compare (max 10)",
    )


# =============================================================================
# ML-Readiness Schemas (4-dimension redesign)
# =============================================================================


class MLDimensionItemSchema(BaseModel):
    """A single scored item within a dimension."""

    name: str = Field(description="Item name")
    score: float = Field(description="Points earned")
    max_score: float = Field(description="Maximum possible points")
    passed: Optional[bool] = Field(None, description="Pass/fail (for binary checks)")
    subtitle: Optional[str] = Field(None, description="Contextual subtitle (e.g. actual value)")
    tooltip: Optional[str] = Field(None, description="Explanatory tooltip text")


class MLDimensionSchema(BaseModel):
    """A scored dimension with constituent items."""

    name: str = Field(description="Dimension name")
    score: float = Field(description="Dimension score")
    max_score: float = Field(description="Maximum dimension score")
    items: List[MLDimensionItemSchema] = Field(
        default_factory=list, description="Scored items within this dimension"
    )
    details: Dict[str, Any] = Field(
        default_factory=dict, description="Additional dimension-specific data"
    )


class MLReadinessBreakdownSchema(BaseModel):
    """Breakdown of ML-readiness score into 4 dimensions."""

    structural_quality: MLDimensionSchema = Field(
        description="Structural Quality (20 pts)"
    )
    property_profile: MLDimensionSchema = Field(
        description="Property Profile (35 pts)"
    )
    complexity_feasibility: MLDimensionSchema = Field(
        description="Complexity & Feasibility (25 pts)"
    )
    representation_quality: MLDimensionSchema = Field(
        description="Representation Quality (20 pts)"
    )


class MLReadinessResultSchema(BaseModel):
    """Result of ML-readiness scoring."""

    score: int = Field(ge=0, le=100, description="Overall ML-readiness score (0-100)")
    label: str = Field(description="Score label (Excellent/Good/Moderate/Limited/Poor)")
    color: str = Field(description="Color key (green/teal/amber/orange/red)")
    breakdown: MLReadinessBreakdownSchema = Field(
        description="Score breakdown by 4 dimensions"
    )
    caveats: List[str] = Field(
        default_factory=list, description="Relevant structural warnings"
    )
    supplementary: Dict[str, Any] = Field(
        default_factory=dict,
        description="Informational data not used in scoring (e.g. Lipinski, Veber)",
    )
    interpretation: str = Field(
        description="Human-readable interpretation of the score"
    )


class NPLikenessResultSchema(BaseModel):
    """Result of NP-likeness scoring."""

    score: float = Field(description="NP-likeness score (typically -5 to +5)")
    interpretation: str = Field(description="Human-readable interpretation")
    caveats: List[str] = Field(
        default_factory=list, description="Warnings or limitations"
    )
    details: Dict[str, Any] = Field(
        default_factory=dict, description="Additional scoring details"
    )


class ScaffoldResultSchema(BaseModel):
    """Result of scaffold extraction."""

    scaffold_smiles: str = Field(description="SMILES of the Murcko scaffold")
    generic_scaffold_smiles: str = Field(
        description="SMILES of the generic (framework) scaffold"
    )
    has_scaffold: bool = Field(
        description="Whether a scaffold was found (False for acyclic molecules)"
    )
    message: str = Field(description="Status message about the extraction")
    details: Dict[str, Any] = Field(
        default_factory=dict, description="Additional extraction details"
    )


class ScoringRequest(BaseModel):
    """Request for molecule scoring."""

    molecule: str = Field(
        ...,
        min_length=1,
        max_length=10000,
        description="Molecule string (SMILES, InChI, or MOL block)",
    )
    format: str = Field(
        default="auto", pattern="^(auto|smiles|inchi|mol)$", description="Input format"
    )
    include: List[str] = Field(
        default=[
            "ml_readiness",
            "np_likeness",
            "scaffold",
            "druglikeness",
            "safety_filters",
            "admet",
            "aggregator",
        ],
        description="Scoring types to include",
    )

    @field_validator("molecule")
    @classmethod
    def sanitize_molecule_input(cls, v: str) -> str:
        """Sanitize molecule input."""
        dangerous = ["<", ">", "&", ";", "|", "$", "`"]
        if any(c in v for c in dangerous):
            raise ValueError("Invalid characters in molecule string")
        return v.strip()

    @field_validator("include")
    @classmethod
    def validate_include(cls, v: List[str]) -> List[str]:
        """Validate include list."""
        valid_types = {
            "ml_readiness",
            "np_likeness",
            "scaffold",
            "druglikeness",
            "safety_filters",
            "admet",
            "aggregator",
            "consensus",
            "lead_likeness",
            "salt_inventory",
            "ligand_efficiency",
            "tpsa_breakdown",
            "logp_breakdown",
            "bertz_detail",
            "fsp3_detail",
            "np_breakdown",
            "bioavailability_radar",
            "boiled_egg",
        }
        for item in v:
            if item not in valid_types:
                raise ValueError(
                    f"Invalid scoring type: {item}. Valid types: {valid_types}"
                )
        return v


class MoleculeInfoSchema(BaseModel):
    """Basic molecule information."""

    input_string: str = Field(description="Original input string")
    canonical_smiles: Optional[str] = Field(None, description="Canonical SMILES")
    molecular_formula: Optional[str] = Field(None, description="Molecular formula")
    molecular_weight: Optional[float] = Field(None, description="Molecular weight")


class ScoringResponse(BaseModel):
    """Response containing scoring results."""

    status: str = Field(default="completed", description="Request status")
    molecule_info: MoleculeInfoSchema = Field(description="Basic molecule information")
    ml_readiness: Optional[MLReadinessResultSchema] = Field(
        None, description="ML-readiness scoring results"
    )
    np_likeness: Optional[NPLikenessResultSchema] = Field(
        None, description="NP-likeness scoring results"
    )
    scaffold: Optional[ScaffoldResultSchema] = Field(
        None, description="Scaffold extraction results"
    )
    druglikeness: Optional[DrugLikenessResultSchema] = Field(
        None, description="Drug-likeness scoring results"
    )
    safety_filters: Optional[SafetyFilterResultSchema] = Field(
        None, description="Safety filter results"
    )
    admet: Optional[ADMETResultSchema] = Field(
        None, description="ADMET prediction results"
    )
    aggregator: Optional[AggregatorLikelihoodSchema] = Field(
        None, description="Aggregator likelihood prediction"
    )
    consensus: Optional[ConsensusScoreSchema] = Field(
        None, description="Consensus drug-likeness score (0-5)"
    )
    lead_likeness: Optional[LeadLikenessSchema] = Field(
        None, description="Lead-likeness assessment"
    )
    salt_inventory: Optional[SaltInventorySchema] = Field(
        None, description="Salt/counterion inventory"
    )
    ligand_efficiency: Optional[LigandEfficiencySchema] = Field(
        None, description="Ligand efficiency"
    )
    tpsa_breakdown: Optional[TPSABreakdownSchema] = Field(
        None, description="TPSA per-atom breakdown"
    )
    logp_breakdown: Optional[LogPBreakdownSchema] = Field(
        None, description="LogP per-atom breakdown"
    )
    bertz_detail: Optional[BertzDetailSchema] = Field(
        None, description="Bertz complexity detail"
    )
    fsp3_detail: Optional[Fsp3DetailSchema] = Field(
        None, description="Fsp3 per-carbon detail"
    )
    np_breakdown: Optional[NPBreakdownSchema] = Field(
        None, description="NP-likeness fragment breakdown"
    )
    bioavailability_radar: Optional[BioavailabilityRadarSchema] = Field(
        None, description="Bioavailability radar (6 axes)"
    )
    boiled_egg: Optional[BoiledEggSchema] = Field(
        None, description="BOILED-Egg GI/BBB classification"
    )
    execution_time_ms: int = Field(description="Execution time in milliseconds")
