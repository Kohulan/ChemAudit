"""
Scoring Schemas

Pydantic schemas for scoring requests and responses.
"""
from pydantic import BaseModel, Field, field_validator
from typing import List, Optional, Dict, Any


class MLReadinessBreakdownSchema(BaseModel):
    """Breakdown of ML-readiness score components."""

    descriptors_score: float = Field(description="Score for descriptor calculability (0-40)")
    descriptors_max: float = Field(default=40.0, description="Maximum descriptor score")
    descriptors_successful: int = Field(description="Number of successfully calculated descriptors")
    descriptors_total: int = Field(description="Total number of descriptors attempted")

    fingerprints_score: float = Field(description="Score for fingerprint generation (0-40)")
    fingerprints_max: float = Field(default=40.0, description="Maximum fingerprint score")
    fingerprints_successful: List[str] = Field(description="Successfully generated fingerprint types")
    fingerprints_failed: List[str] = Field(default_factory=list, description="Failed fingerprint types")

    size_score: float = Field(description="Score for size constraints (0-20)")
    size_max: float = Field(default=20.0, description="Maximum size score")
    molecular_weight: Optional[float] = Field(None, description="Molecular weight")
    num_atoms: Optional[int] = Field(None, description="Number of atoms")
    size_category: str = Field(description="Size category: optimal, acceptable, or out_of_range")


class MLReadinessResultSchema(BaseModel):
    """Result of ML-readiness scoring."""

    score: int = Field(ge=0, le=100, description="Overall ML-readiness score (0-100)")
    breakdown: MLReadinessBreakdownSchema = Field(description="Score breakdown by component")
    interpretation: str = Field(description="Human-readable interpretation of the score")
    failed_descriptors: List[str] = Field(default_factory=list, description="List of failed descriptor names")


class NPLikenessResultSchema(BaseModel):
    """Result of NP-likeness scoring."""

    score: float = Field(description="NP-likeness score (typically -5 to +5)")
    interpretation: str = Field(description="Human-readable interpretation")
    caveats: List[str] = Field(default_factory=list, description="Warnings or limitations")
    details: Dict[str, Any] = Field(default_factory=dict, description="Additional scoring details")


class ScaffoldResultSchema(BaseModel):
    """Result of scaffold extraction."""

    scaffold_smiles: str = Field(description="SMILES of the Murcko scaffold")
    generic_scaffold_smiles: str = Field(description="SMILES of the generic (framework) scaffold")
    has_scaffold: bool = Field(description="Whether a scaffold was found (False for acyclic molecules)")
    message: str = Field(description="Status message about the extraction")
    details: Dict[str, Any] = Field(default_factory=dict, description="Additional extraction details")


class ScoringRequest(BaseModel):
    """Request for molecule scoring."""

    molecule: str = Field(
        ...,
        min_length=1,
        max_length=10000,
        description="Molecule string (SMILES, InChI, or MOL block)"
    )
    format: str = Field(
        default="auto",
        pattern="^(auto|smiles|inchi|mol)$",
        description="Input format"
    )
    include: List[str] = Field(
        default=["ml_readiness", "np_likeness", "scaffold"],
        description="Scoring types to include"
    )

    @field_validator('molecule')
    @classmethod
    def sanitize_molecule_input(cls, v: str) -> str:
        """Sanitize molecule input."""
        dangerous = ['<', '>', '&', ';', '|', '$', '`']
        if any(c in v for c in dangerous):
            raise ValueError("Invalid characters in molecule string")
        return v.strip()

    @field_validator('include')
    @classmethod
    def validate_include(cls, v: List[str]) -> List[str]:
        """Validate include list."""
        valid_types = {"ml_readiness", "np_likeness", "scaffold"}
        for item in v:
            if item not in valid_types:
                raise ValueError(f"Invalid scoring type: {item}. Valid types: {valid_types}")
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
    execution_time_ms: int = Field(description="Execution time in milliseconds")
