"""
Standardization Schemas

Pydantic schemas for standardization requests and responses.
"""

from typing import List, Literal, Optional

from pydantic import BaseModel, Field

from app.schemas.validation import MoleculeInfo


# ---------------------------------------------------------------------------
# Provenance schemas (Phase 02, Plan 01 — Milestone 2.1)
# ---------------------------------------------------------------------------


class ChargeChange(BaseModel):
    """Atom-level formal charge change during standardization (STD-02, STD-03)."""

    atom_idx: int = Field(..., description="Atom index in the molecule")
    element: str = Field(..., description="Element symbol (e.g. 'N', 'O', 'S')")
    before_charge: int = Field(..., description="Formal charge before standardization")
    after_charge: int = Field(..., description="Formal charge after standardization")
    rule_name: str = Field(..., description="Name of the normalization rule applied")
    smarts: str = Field(default="", description="SMARTS pattern of the rule (if known)")


class BondChange(BaseModel):
    """Bond-level type change during standardization (STD-03)."""

    bond_idx: int = Field(..., description="Bond index in the molecule")
    atom1_idx: int = Field(..., description="Index of the first atom")
    atom2_idx: int = Field(..., description="Index of the second atom")
    before_type: str = Field(..., description="Bond type before standardization")
    after_type: str = Field(..., description="Bond type after standardization")
    rule_name: str = Field(..., description="Name of the normalization rule applied")


class RadicalChange(BaseModel):
    """Atom-level radical electron change during standardization (STD-03)."""

    atom_idx: int = Field(..., description="Atom index in the molecule")
    element: str = Field(..., description="Element symbol")
    before_radicals: int = Field(..., description="Radical electrons before standardization")
    after_radicals: int = Field(..., description="Radical electrons after standardization")


class FragmentRemoval(BaseModel):
    """A fragment removed during parent extraction (STD-04)."""

    smiles: str = Field(..., description="Canonical SMILES of the removed fragment")
    name: Optional[str] = Field(None, description="Human-readable name (if known)")
    role: Literal["salt", "solvent", "counterion", "unknown"] = Field(
        ..., description="Role classification of the fragment"
    )
    mw: float = Field(..., description="Exact molecular weight of the fragment")


class TautomerProvenance(BaseModel):
    """Tautomer canonicalization provenance record (STD-01)."""

    input_smiles: str = Field(..., description="SMILES before tautomer canonicalization")
    canonical_smiles: str = Field(..., description="Canonical tautomer SMILES")
    num_tautomers_enumerated: int = Field(
        ..., description="Number of tautomers enumerated"
    )
    modified_atoms: List[int] = Field(
        default_factory=list, description="Atom indices modified during tautomerization"
    )
    modified_bonds: List[int] = Field(
        default_factory=list, description="Bond indices modified during tautomerization"
    )
    stereo_stripped: bool = Field(
        ..., description="Whether stereochemistry was lost during tautomerization"
    )
    complexity_flag: bool = Field(
        default=False,
        description="True when > 100 tautomers enumerated (high-complexity molecule)",
    )


class ProvStageRecord(BaseModel):
    """A single pipeline stage provenance record."""

    stage_name: str = Field(..., description="Pipeline stage name")
    input_smiles: str = Field(..., description="Input SMILES for this stage")
    output_smiles: str = Field(..., description="Output SMILES from this stage")
    applied: bool = Field(..., description="Whether the stage made any changes")
    charge_changes: List[ChargeChange] = Field(
        default_factory=list, description="Atom-level charge changes (STD-02, STD-03)"
    )
    bond_changes: List[BondChange] = Field(
        default_factory=list, description="Bond-level type changes (STD-03)"
    )
    radical_changes: List[RadicalChange] = Field(
        default_factory=list, description="Radical electron changes (STD-03)"
    )
    fragment_removals: List[FragmentRemoval] = Field(
        default_factory=list, description="Fragments removed during parent extraction (STD-04)"
    )
    dval_cross_refs: List[str] = Field(
        default_factory=list,
        description="Cross-references to deep validation check IDs",
    )


class StereoProvenance(BaseModel):
    """Stereochemistry summary across the full standardization pipeline."""

    stereo_stripped: bool = Field(
        ..., description="Whether any stereochemistry was lost"
    )
    centers_lost: int = Field(..., description="Number of chiral centers lost")
    bonds_lost: int = Field(..., description="Number of E/Z stereo bonds lost")
    per_center: List[dict] = Field(
        default_factory=list,
        description="Per-center detail (reserved for STD-06 in Plan 02)",
    )


class StandardizationProvenance(BaseModel):
    """Full provenance record for a standardization run."""

    stages: List[ProvStageRecord] = Field(
        ..., description="Per-stage provenance records"
    )
    tautomer: Optional[TautomerProvenance] = Field(
        None, description="Tautomer provenance (only when include_tautomer=True)"
    )
    stereo_summary: Optional[StereoProvenance] = Field(
        None, description="Stereochemistry summary across all stages"
    )


# ---------------------------------------------------------------------------
# Existing schemas (unchanged)
# ---------------------------------------------------------------------------


class StandardizationOptions(BaseModel):
    """Options for the standardization pipeline."""

    include_tautomer: bool = Field(
        default=False,
        description="Include tautomer canonicalization (WARNING: may lose E/Z stereochemistry)",
    )
    preserve_stereo: bool = Field(
        default=True,
        description="Attempt to preserve stereochemistry during standardization",
    )
    include_provenance: bool = Field(
        default=False,
        description="Include detailed provenance records for each pipeline stage",
    )


class StandardizationStep(BaseModel):
    """A single step in the standardization pipeline."""

    step_name: str = Field(..., description="Name of the pipeline step")
    applied: bool = Field(..., description="Whether the step was successfully applied")
    description: str = Field(
        ..., description="Human-readable description of what this step does"
    )
    changes: str = Field(
        default="", description="Description of changes made by this step"
    )


class StereoComparison(BaseModel):
    """Comparison of stereochemistry before and after standardization."""

    before_count: int = Field(
        default=0, description="Number of defined stereocenters before"
    )
    after_count: int = Field(
        default=0, description="Number of defined stereocenters after"
    )
    lost: int = Field(default=0, description="Number of stereocenters lost")
    gained: int = Field(default=0, description="Number of stereocenters gained")
    double_bond_stereo_lost: int = Field(
        default=0, description="Number of E/Z stereo bonds lost"
    )
    warning: Optional[str] = Field(
        None, description="Warning message if stereochemistry was lost"
    )


class StructureComparisonSchema(BaseModel):
    """Comparison between original and standardized structures."""

    original_atom_count: int = Field(default=0)
    standardized_atom_count: int = Field(default=0)
    original_formula: Optional[str] = None
    standardized_formula: Optional[str] = None
    original_mw: Optional[float] = None
    standardized_mw: Optional[float] = None
    mass_change_percent: float = Field(default=0.0)
    is_identical: bool = Field(default=False)
    diff_summary: List[str] = Field(default_factory=list)


class CheckerIssue(BaseModel):
    """An issue found by the ChEMBL structure checker."""

    penalty_score: int = Field(
        ..., description="Penalty score for this issue (higher = worse)"
    )
    message: str = Field(..., description="Description of the issue")


class StandardizationResult(BaseModel):
    """Result of the standardization pipeline."""

    original_smiles: str = Field(..., description="Original input SMILES")
    standardized_smiles: Optional[str] = Field(None, description="Standardized SMILES")

    success: bool = Field(
        default=False, description="Whether standardization was successful"
    )
    error_message: Optional[str] = Field(None, description="Error message if failed")

    steps_applied: List[StandardizationStep] = Field(
        default_factory=list, description="List of pipeline steps and their results"
    )

    checker_issues: List[CheckerIssue] = Field(
        default_factory=list, description="Issues found by ChEMBL structure checker"
    )

    excluded_fragments: List[str] = Field(
        default_factory=list,
        description="SMILES of fragments removed (salts, solvents)",
    )

    stereo_comparison: Optional[StereoComparison] = Field(
        None, description="Stereochemistry comparison before/after"
    )

    structure_comparison: Optional[StructureComparisonSchema] = Field(
        None, description="Structure comparison before/after"
    )

    mass_change_percent: float = Field(
        default=0.0, description="Percentage change in molecular weight"
    )

    provenance: Optional[StandardizationProvenance] = Field(
        None,
        description="Detailed provenance records when include_provenance=True",
    )


class StandardizeRequest(BaseModel):
    """Request to standardize a molecule."""

    molecule: str = Field(
        ...,
        min_length=1,
        max_length=10000,
        description="Molecule string (SMILES, InChI, or MOL block)",
    )
    format: str = Field(
        default="auto", pattern="^(auto|smiles|inchi|mol)$", description="Input format"
    )
    options: StandardizationOptions = Field(
        default_factory=StandardizationOptions, description="Standardization options"
    )


class StandardizeResponse(BaseModel):
    """Response from standardization endpoint."""

    molecule_info: MoleculeInfo = Field(
        ..., description="Information about the original molecule"
    )
    result: StandardizationResult = Field(..., description="Standardization result")
    execution_time_ms: int = Field(..., description="Execution time in milliseconds")
