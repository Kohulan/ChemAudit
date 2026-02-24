"""
Validation Schemas

Pydantic schemas for validation requests and responses.
"""

from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field, field_validator

from app.schemas.common import Severity


class CheckResultSchema(BaseModel):
    """Schema for a single validation check result."""

    check_name: str
    passed: bool
    severity: Severity
    message: str
    affected_atoms: List[int] = Field(default_factory=list)
    details: Dict[str, Any] = Field(default_factory=dict)


class MoleculeInfo(BaseModel):
    """Schema for molecule information."""

    input_smiles: str
    canonical_smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    num_atoms: Optional[int] = None


class SimilarityRequest(BaseModel):
    """Schema for pairwise ECFP4 Tanimoto similarity request."""

    smiles_a: str = Field(..., min_length=1, max_length=10000, description="First SMILES string")
    smiles_b: str = Field(..., min_length=1, max_length=10000, description="Second SMILES string")


class SimilarityResponse(BaseModel):
    """Schema for pairwise ECFP4 Tanimoto similarity response."""

    tanimoto_similarity: float = Field(ge=0.0, le=1.0, description="Tanimoto similarity (0-1)")
    fingerprint_type: str = "ECFP4"
    radius: int = 2
    n_bits: int = 2048
    common_bits: int = Field(description="Number of bits set in both fingerprints")
    bits_a: int = Field(description="Number of bits set in fingerprint A")
    bits_b: int = Field(description="Number of bits set in fingerprint B")


class InputInterpretation(BaseModel):
    """Schema for input type detection and IUPAC conversion transparency."""

    detected_as: str = Field(description="Input detected as: 'smiles' or 'iupac'")
    original_input: str = Field(description="Original input string before conversion")
    converted_smiles: Optional[str] = Field(
        None, description="Converted SMILES if input was IUPAC"
    )
    conversion_source: Optional[str] = Field(
        None, description="Conversion engine used: 'opsin' or None"
    )


class ValidationRequest(BaseModel):
    """Schema for validation request."""

    molecule: str = Field(
        ...,
        min_length=1,
        max_length=10000,
        description="Molecule string (SMILES, InChI, MOL block, or IUPAC name)",
    )
    format: str = Field(
        default="auto", pattern="^(auto|smiles|inchi|mol|iupac)$", description="Input format"
    )
    input_type: Optional[str] = Field(
        default=None,
        description="Override input type detection: 'smiles', 'iupac', or 'auto' (default: auto)",
    )
    checks: List[str] = Field(default=["all"], description="List of checks to run")
    preserve_aromatic: bool = Field(
        default=False,
        description="Preserve aromatic notation in canonical SMILES output",
    )

    @field_validator("molecule")
    @classmethod
    def sanitize_molecule_input(cls, v: str) -> str:
        """Sanitize molecule input to prevent injection attacks."""
        dangerous = ["<", ">", "&", ";", "|", "$", "`"]
        if any(c in v for c in dangerous):
            raise ValueError("Invalid characters in molecule string")
        return v.strip()


class ValidationResponse(BaseModel):
    """Schema for validation response."""

    status: str = "completed"
    molecule_info: MoleculeInfo
    overall_score: int = Field(
        ge=0, le=100, description="Overall validation score (0-100)"
    )
    issues: List[CheckResultSchema] = Field(
        default_factory=list, description="Failed checks"
    )
    all_checks: List[CheckResultSchema] = Field(
        default_factory=list, description="All check results"
    )
    execution_time_ms: int = Field(description="Execution time in milliseconds")
    cached: bool = Field(
        default=False, description="Whether result was served from cache"
    )
    input_interpretation: Optional[InputInterpretation] = Field(
        None, description="Input type detection result (present when IUPAC conversion occurred)"
    )
