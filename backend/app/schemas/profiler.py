"""
Profiler Schemas

Pydantic v2 request/response models for the compound profiler endpoints.
Covers: PFI, #stars, Abbott bioavailability, consensus LogP, skin permeation,
3D shape descriptors, CNS MPO, custom MPO, and extended ligand efficiency.
"""

from typing import Any, Optional

from pydantic import BaseModel, Field

# ---------------------------------------------------------------------------
# Request Models
# ---------------------------------------------------------------------------


class ProfileRequest(BaseModel):
    """Request for full compound profile computation."""

    smiles: str = Field(..., min_length=1, description="SMILES string to profile")


class LERequest(BaseModel):
    """Request for extended ligand efficiency computation."""

    smiles: str = Field(..., min_length=1, description="SMILES string")
    activity_value: float = Field(..., description="Activity measurement value")
    activity_type: str = Field(
        ...,
        pattern="^(IC50_nM|IC50_uM|Ki_nM|pIC50|pKd)$",
        description="Activity type: IC50_nM | IC50_uM | Ki_nM | pIC50 | pKd",
    )


class MPOProperty(BaseModel):
    """A single property entry in a custom MPO profile."""

    property: str = Field(..., description="Property name (MW, LogP, TPSA, HBD, etc.)")
    low: float = Field(..., description="Lower boundary of the favorable range")
    high: float = Field(..., description="Upper boundary of the favorable range")
    weight: float = Field(default=1.0, ge=0, description="Relative importance weight")
    shape: str = Field(
        default="sigmoid",
        pattern="^(sigmoid|ramp|step)$",
        description="Desirability curve shape: sigmoid | ramp | step",
    )


class MPORequest(BaseModel):
    """Request for custom MPO score computation."""

    smiles: str = Field(..., min_length=1, description="SMILES string")
    profile: list[MPOProperty] = Field(
        ...,
        min_length=1,
        description="List of MPO property desirability definitions",
    )


class Shape3DRequest(BaseModel):
    """Request for 3D shape descriptor computation."""

    smiles: str = Field(..., min_length=1, description="SMILES string")


class SAComparisonRequest(BaseModel):
    """Request for SA Score / SCScore / SYBA comparison."""

    smiles: str = Field(..., min_length=1, description="SMILES string")


# ---------------------------------------------------------------------------
# Response Models
# ---------------------------------------------------------------------------


class PFIResult(BaseModel):
    """Property Forecast Index result."""

    pfi: float = Field(description="PFI score (cLogP + aromatic rings)")
    clogp: float = Field(description="Wildman-Crippen cLogP")
    aromatic_rings: int = Field(description="Number of aromatic rings")
    risk: str = Field(description="Risk classification: low | moderate | high")


class StarsPropertyDetail(BaseModel):
    """Per-property detail for #stars."""

    property: str = Field(description="Property name")
    value: float = Field(description="Computed property value")
    range_low: float = Field(description="95th-percentile lower bound")
    range_high: float = Field(description="95th-percentile upper bound")
    violated: bool = Field(description="True if value is outside the range")


class StarsResult(BaseModel):
    """#Stars outlier count result."""

    stars: int = Field(ge=0, description="Number of 95th-percentile violations")
    details: list[StarsPropertyDetail] = Field(
        description="Per-property detail for all 8 properties"
    )


class AbbottResult(BaseModel):
    """Abbott Bioavailability Score result."""

    abbott_score: float = Field(description="Bioavailability probability (0.11/0.17/0.56/0.85)")
    probability_pct: int = Field(description="Bioavailability probability as integer percent")
    tpsa: float = Field(description="Topological polar surface area used in calculation")
    lipinski_violations: int = Field(description="Number of Lipinski Ro5 violations")


class ConsensusLogPResult(BaseModel):
    """Consensus LogP result."""

    consensus_logp: float = Field(description="Average of Wildman-Crippen and XLOGP3 approx")
    wildman_crippen: float = Field(description="Wildman-Crippen cLogP (RDKit)")
    xlogp3_approx: float = Field(
        description="XLOGP3 approximation (0.92 * Wildman-Crippen + 0.12)"
    )
    xlogp3_is_approximation: bool = Field(
        description="Always True — XLOGP3 value is approximated, not computed by XLogP3 itself"
    )


class SkinPermeationResult(BaseModel):
    """Skin permeation coefficient result."""

    log_kp: float = Field(description="Skin permeation coefficient (log Kp, cm/h)")
    classification: str = Field(description="Permeation class: low | moderate | high")


class Shape3DResult(BaseModel):
    """3D shape descriptor result."""

    pmi1: float = Field(description="Principal Moment of Inertia 1 (smallest)")
    pmi2: float = Field(description="Principal Moment of Inertia 2")
    pmi3: float = Field(description="Principal Moment of Inertia 3 (largest)")
    npr1: float = Field(ge=0, le=1, description="Normalized PMI ratio 1 (rod axis)")
    npr2: float = Field(ge=0, le=1, description="Normalized PMI ratio 2 (disc axis)")
    pbf: float = Field(description="Plane of Best Fit (PBF)")
    shape_class: str = Field(description="Shape classification: sphere | disc | rod")


class CNSMPOResult(BaseModel):
    """CNS MPO score result (Wager 2010, 4-component RDKit version)."""

    score: float = Field(ge=0, le=4, description="CNS MPO score (0-4 for 4-component version)")
    max_score: float = Field(description="Maximum possible score (4.0 for 4-component)")
    components: dict[str, float] = Field(
        description="Per-component desirability values: clogp, tpsa, mw, hbd"
    )


class MPOComponentDetail(BaseModel):
    """Per-property detail for custom MPO."""

    property: str = Field(description="Property name")
    value: float = Field(description="Computed property value")
    desirability: float = Field(ge=0, le=1, description="Desirability (0-1)")
    weight: float = Field(description="Property weight")
    contribution: float = Field(description="Weighted contribution")


class CustomMPOResult(BaseModel):
    """Custom MPO score result."""

    score: float = Field(description="Weighted sum of desirabilities")
    max_score: float = Field(description="Sum of all weights (maximum possible score)")
    normalized: float = Field(ge=0, le=1, description="score / max_score (0-1)")
    components: list[MPOComponentDetail] = Field(description="Per-property breakdown")


class LEResult(BaseModel):
    """Extended ligand efficiency result."""

    pIC50: float = Field(description="Converted pIC50 value")  # noqa: N815
    LE: float = Field(description="Ligand Efficiency (1.4 * pIC50 / HA)")
    LLE: float = Field(description="Lipophilic Ligand Efficiency (pIC50 - LogP)")
    LELP: float = Field(description="Ligand Efficiency Lipophilicity Price (LogP / LE)")
    BEI: float = Field(description="Binding Efficiency Index (pIC50 * 1000 / MW)")
    SEI: float = Field(description="Surface Efficiency Index (pIC50 * 100 / TPSA)")


class ProfileResponse(BaseModel):
    """Composite compound profile response (full endpoint)."""

    pfi: PFIResult = Field(description="PFI score and risk classification")
    stars: StarsResult = Field(description="#Stars outlier count with property details")
    abbott: AbbottResult = Field(description="Abbott Bioavailability Score")
    consensus_logp: ConsensusLogPResult = Field(description="Consensus LogP")
    skin_permeation: SkinPermeationResult = Field(description="Skin permeation (log Kp)")
    druglikeness: Optional[dict[str, Any]] = Field(
        None, description="Drug-likeness rules summary (Lipinski, QED, Veber)"
    )
