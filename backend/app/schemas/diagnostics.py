"""
Pydantic v2 schemas for Structure Quality Diagnostics (Phase 9).

Provides request and response models for all five diagnostic tools:
- DIAG-01: SMILES error diagnostics
- DIAG-02: InChI layer-by-layer diff comparison
- DIAG-03: Format round-trip lossiness checker
- DIAG-04: Cross-pipeline standardization comparison
- DIAG-05: File pre-validation (SDF and CSV)
"""

from typing import Any, Optional

from pydantic import BaseModel, Field, field_validator

# =============================================================================
# Request Models
# =============================================================================


class SMILESDiagnosticsRequest(BaseModel):
    """Request for SMILES error diagnostics (DIAG-01)."""

    smiles: str = Field(
        ...,
        min_length=1,
        max_length=10000,
        description="Raw SMILES string to diagnose",
    )

    @field_validator("smiles")
    @classmethod
    def sanitize_smiles_input(cls, v: str) -> str:
        """Sanitize SMILES input to prevent injection attacks."""
        dangerous = ["<", ">", "&", ";", "|", "$", "`"]
        if any(c in v for c in dangerous):
            raise ValueError("Invalid characters in SMILES string")
        return v.strip()


class InChIDiffRequest(BaseModel):
    """Request for InChI layer-by-layer comparison (DIAG-02)."""

    inchi_a: str = Field(
        ...,
        min_length=1,
        max_length=10000,
        description="First InChI string",
    )
    inchi_b: str = Field(
        ...,
        min_length=1,
        max_length=10000,
        description="Second InChI string",
    )

    @field_validator("inchi_a", "inchi_b")
    @classmethod
    def sanitize_inchi_input(cls, v: str) -> str:
        """Sanitize InChI input to prevent injection attacks."""
        dangerous = ["<", ">", "&", ";", "|", "$", "`"]
        if any(c in v for c in dangerous):
            raise ValueError("Invalid characters in InChI string")
        return v.strip()


class RoundTripRequest(BaseModel):
    """Request for format round-trip lossiness check (DIAG-03)."""

    smiles: str = Field(
        ...,
        min_length=1,
        max_length=10000,
        description="SMILES string to check for round-trip lossiness",
    )
    route: str = Field(
        default="smiles_inchi_smiles",
        pattern="^(smiles_inchi_smiles|smiles_mol_smiles)$",
        description="Conversion route: smiles_inchi_smiles or smiles_mol_smiles",
    )

    @field_validator("smiles")
    @classmethod
    def sanitize_smiles_input(cls, v: str) -> str:
        """Sanitize SMILES input to prevent injection attacks."""
        dangerous = ["<", ">", "&", ";", "|", "$", "`"]
        if any(c in v for c in dangerous):
            raise ValueError("Invalid characters in SMILES string")
        return v.strip()


class CrossPipelineRequest(BaseModel):
    """Request for cross-pipeline standardization comparison (DIAG-04)."""

    molecule: str = Field(
        ...,
        min_length=1,
        max_length=10000,
        description="Molecule input (SMILES, InChI, or MOL format)",
    )
    format: str = Field(
        default="auto",
        pattern="^(auto|smiles|inchi|mol)$",
        description="Input format: auto, smiles, inchi, or mol",
    )

    @field_validator("molecule")
    @classmethod
    def sanitize_molecule_input(cls, v: str) -> str:
        """Sanitize molecule input to prevent injection attacks."""
        dangerous = ["<", ">", "&", ";", "|", "$", "`"]
        if any(c in v for c in dangerous):
            raise ValueError("Invalid characters in molecule string")
        return v.strip()


# =============================================================================
# Response Sub-Models
# =============================================================================


class FixSuggestion(BaseModel):
    """A single fix suggestion for a SMILES error."""

    description: str = Field(description="Human-readable description of the suggested fix")
    corrected_smiles: Optional[str] = Field(
        default=None,
        description="Corrected SMILES string if auto-correction is possible",
    )
    confidence: float = Field(description="Confidence in this suggestion (0.0 to 1.0)")


class SMILESError(BaseModel):
    """A single error found while diagnosing a SMILES string."""

    raw_message: str = Field(description="Raw error message from RDKit")
    position: Optional[int] = Field(
        default=None,
        description="Zero-indexed character position of the error, if known",
    )
    error_type: str = Field(
        description=(
            "Structured error type: unmatched_bracket, valence_error, "
            "ring_closure_mismatch, unknown_atom_symbol, invalid_charge, or parse_error"
        )
    )
    message: str = Field(description="Human-readable error description")
    suggestions: list[FixSuggestion] = Field(
        default_factory=list,
        description="Fix suggestions sorted by confidence descending",
    )


class SMILESWarning(BaseModel):
    """A single warning from RDKit chemistry problem detection on a valid SMILES."""

    atom_index: Optional[int] = Field(
        default=None,
        description="Zero-indexed atom index associated with the warning, if applicable",
    )
    type: str = Field(description="RDKit warning type string")
    message: str = Field(description="Warning message from RDKit")


class LayerRow(BaseModel):
    """A single layer comparison row in an InChI diff result."""

    layer: str = Field(description="Layer name (e.g., formula, connections, hydrogens)")
    value_a: Optional[str] = Field(
        default=None,
        description="Layer value from the first InChI string",
    )
    value_b: Optional[str] = Field(
        default=None,
        description="Layer value from the second InChI string",
    )
    match: bool = Field(description="True if both values are identical")


class RoundTripLoss(BaseModel):
    """A single type of information lost during a format round-trip."""

    type: str = Field(description="Loss type: stereo, charge, or isotope")
    description: str = Field(description="Human-readable description of what was lost")
    before: int = Field(description="Count/value before round-trip conversion")
    after: int = Field(description="Count/value after round-trip conversion")


class PipelineResult(BaseModel):
    """Standardization result from a single pipeline."""

    name: str = Field(description="Pipeline name (e.g., RDKit MolStandardize, ChEMBL Pipeline)")
    smiles: Optional[str] = Field(default=None, description="Canonical SMILES output")
    inchikey: Optional[str] = Field(default=None, description="InChIKey of standardized molecule")
    mw: Optional[float] = Field(
        default=None, description="Exact molecular weight (4 decimal places)"
    )
    formula: Optional[str] = Field(default=None, description="Molecular formula")
    charge: Optional[int] = Field(default=None, description="Net formal charge")
    stereo_count: Optional[int] = Field(default=None, description="Number of stereo centers")
    error: Optional[str] = Field(
        default=None,
        description="Error message if this pipeline failed",
    )
    highlight_atoms: list[int] = Field(
        default_factory=list,
        description=(
            "Zero-indexed atom indices (relative to the canonical SMILES in `smiles`) "
            "that are not part of the cross-pipeline MCS. Empty when pipelines agree."
        ),
    )
    highlight_bonds: list[int] = Field(
        default_factory=list,
        description=(
            "Zero-indexed bond indices that touch at least one non-MCS atom. "
            "Empty when pipelines agree."
        ),
    )


class PropertyComparison(BaseModel):
    """Per-property agreement comparison across all pipelines."""

    property: str = Field(description="Property name (e.g., smiles, inchikey, mw)")
    values: list[Any] = Field(description="Property values from each pipeline in order")
    agrees: bool = Field(description="True if all pipelines agree on this property")
    structural: bool = Field(
        default=False,
        description="True if this is a structural property (SMILES or InChIKey)",
    )


class FileIssue(BaseModel):
    """A single structural issue found during file pre-validation."""

    block: Optional[int] = Field(
        default=None,
        description="SDF block number (1-indexed) where the issue was found",
    )
    line: Optional[int] = Field(
        default=None,
        description="Absolute file line number (1-indexed) where the issue was found",
    )
    issue_type: str = Field(
        description=(
            "Issue type: missing_m_end, malformed_count_line, encoding_fallback, "
            "encoding_error, suspicious_content, missing_smiles_column, empty_rows, empty_file, "
            "duplicate_columns"
        )
    )
    severity: str = Field(description="Severity level: error, warning, or info")
    description: str = Field(description="Human-readable description of the issue")


# =============================================================================
# Response Models
# =============================================================================


class SMILESDiagnosticsResponse(BaseModel):
    """Response for SMILES error diagnostics (DIAG-01)."""

    valid: bool = Field(description="True if the SMILES string is valid")
    canonical_smiles: Optional[str] = Field(
        default=None,
        description="Canonical SMILES form if valid, else None",
    )
    warnings: list[SMILESWarning] = Field(
        default_factory=list,
        description="Chemistry warnings on valid SMILES",
    )
    errors: list[SMILESError] = Field(
        default_factory=list,
        description="Errors on invalid SMILES with position and fix suggestions",
    )


class InChIDiffResponse(BaseModel):
    """Response for InChI layer-by-layer comparison (DIAG-02)."""

    identical: bool = Field(description="True if all layers in both InChI strings are identical")
    layer_rows: list[LayerRow] = Field(
        description="Per-layer comparison rows sorted by LAYER_DISPLAY_ORDER",
    )
    layers_a: dict[str, str] = Field(description="Parsed layers from the first InChI string")
    layers_b: dict[str, str] = Field(description="Parsed layers from the second InChI string")


class RoundTripResponse(BaseModel):
    """Response for format round-trip lossiness check (DIAG-03)."""

    route: str = Field(description="Conversion route used")
    original_smiles: str = Field(description="Canonical SMILES of the original molecule")
    intermediate: Optional[str] = Field(
        default=None,
        description="Intermediate representation (InChI string or MOL block)",
    )
    roundtrip_smiles: Optional[str] = Field(
        default=None,
        description="Canonical SMILES after round-trip conversion, or None on error",
    )
    lossy: bool = Field(description="True if any information was lost during round-trip")
    losses: list[RoundTripLoss] = Field(
        default_factory=list,
        description="List of detected information losses",
    )
    error: Optional[str] = Field(
        default=None,
        description="Error message if the round-trip could not be completed",
    )


class CrossPipelineResponse(BaseModel):
    """Response for cross-pipeline standardization comparison (DIAG-04)."""

    pipelines: list[PipelineResult] = Field(
        description="Results from each of the three standardization pipelines",
    )
    disagreements: int = Field(
        description="Number of properties where not all pipelines agree",
    )
    structural_disagreements: int = Field(
        description="Number of structural properties (SMILES, InChIKey) that disagree",
    )
    all_agree: bool = Field(
        description="True if all three pipelines agree on all properties",
    )
    property_comparison: list[PropertyComparison] = Field(
        description="Per-property comparison table",
    )


class FilePreValidationResponse(BaseModel):
    """Response for file pre-validation (DIAG-05)."""

    file_type: str = Field(description="Detected file type: sdf or csv")
    total_blocks: Optional[int] = Field(
        default=None,
        description="Total number of SDF molecule blocks (SDF only)",
    )
    total_rows: Optional[int] = Field(
        default=None,
        description="Total number of data rows excluding header (CSV only)",
    )
    encoding: Optional[str] = Field(
        default=None,
        description="Encoding used to decode the file (CSV only)",
    )
    issue_count: int = Field(description="Total number of issues found")
    issues: list[FileIssue] = Field(description="List of all issues found in the file")
    valid: bool = Field(
        description="True if no error-severity issues were found (warnings do not affect validity)",
    )
