"""
Pydantic v2 schemas for QSAR-Ready Pipeline (Phase 10).

Provides request and response models for the QSAR-Ready API surface:
- Single molecule processing
- Batch file upload and status
- Batch results with pagination and summary statistics
"""

from typing import Optional

from pydantic import BaseModel, Field


# =============================================================================
# Config Schema
# =============================================================================


class QSARReadyConfigSchema(BaseModel):
    """Pipeline configuration for QSAR-Ready curation.

    Mirrors QSARReadyConfig dataclass from pipeline.py.
    All boolean toggle fields default to match the default QSARReadyConfig.
    """

    model_config = {"from_attributes": True}

    enable_metals: bool = Field(
        default=True,
        description="Disconnect metal-organic bonds (step 2: MetalDisconnector)",
    )
    enable_desalt: bool = Field(
        default=True,
        description="Keep largest fragment / remove salts (step 3: LargestFragmentChooser)",
    )
    enable_normalize: bool = Field(
        default=True,
        description="Normalize functional groups (step 4: Normalizer)",
    )
    enable_neutralize: bool = Field(
        default=True,
        description="Neutralize charges (step 5: Uncharger)",
    )
    enable_tautomer: bool = Field(
        default=True,
        description="Canonicalize tautomers (step 6: TautomerEnumerator)",
    )
    enable_stereo_strip: bool = Field(
        default=False,
        description="Strip stereochemistry (step 7). True for QSAR-2D, False for QSAR-3D.",
    )
    enable_isotope_strip: bool = Field(
        default=True,
        description="Strip isotope labels (step 8).",
    )
    min_heavy_atoms: int = Field(
        default=3,
        ge=1,
        le=500,
        description="Minimum heavy atom count for composition filter (step 9)",
    )
    max_heavy_atoms: int = Field(
        default=100,
        ge=1,
        le=500,
        description="Maximum heavy atom count for composition filter (step 9)",
    )
    max_mw: float = Field(
        default=1500.0,
        gt=0,
        le=5000.0,
        description="Maximum molecular weight (Da) for composition filter (step 9)",
    )
    remove_inorganics: bool = Field(
        default=True,
        description="Reject molecules with no carbon atoms in composition filter (step 9)",
    )


# =============================================================================
# Step and Result Schemas
# =============================================================================


class QSARStepResultSchema(BaseModel):
    """Result of a single QSAR-Ready pipeline step.

    Mirrors QSARStepResult dataclass from pipeline.py.
    """

    model_config = {"from_attributes": True}

    step_name: str = Field(
        description=(
            "Step identifier: parse, metals, desalt, normalize, neutralize, "
            "tautomer, stereo, isotope, filter, or canonical"
        )
    )
    step_index: int = Field(description="1-based step index (1 through 10)")
    enabled: bool = Field(description="Whether this step was enabled in the config")
    status: str = Field(
        description="Step outcome: 'applied', 'no_change', 'skipped', or 'error'"
    )
    before_smiles: Optional[str] = Field(
        default=None,
        description="SMILES before this step was applied (None for skipped steps)",
    )
    after_smiles: Optional[str] = Field(
        default=None,
        description="SMILES after this step was applied (None for skipped/error steps)",
    )
    detail: Optional[str] = Field(
        default=None,
        description="Human-readable detail: change description, rejection reason, or error message",
    )


class QSARReadyResultSchema(BaseModel):
    """Full result of processing a single molecule through the QSAR-Ready pipeline.

    Mirrors QSARReadyResult dataclass from pipeline.py.
    """

    model_config = {"from_attributes": True}

    original_smiles: str = Field(description="The input SMILES string")
    original_inchikey: Optional[str] = Field(
        default=None,
        description="InChIKey computed from parsed mol BEFORE any pipeline steps (D-14)",
    )
    curated_smiles: Optional[str] = Field(
        default=None,
        description="Canonical SMILES after all pipeline steps (None if rejected or error)",
    )
    standardized_inchikey: Optional[str] = Field(
        default=None,
        description="InChIKey computed from final mol AFTER all pipeline steps (D-14)",
    )
    inchikey_changed: bool = Field(
        default=False,
        description="True if original_inchikey != standardized_inchikey (D-14 locked decision)",
    )
    status: str = Field(
        description="Pipeline outcome: 'ok', 'rejected', 'duplicate', or 'error'"
    )
    rejection_reason: Optional[str] = Field(
        default=None,
        description="Human-readable reason for rejection or duplication",
    )
    steps: list[QSARStepResultSchema] = Field(
        default_factory=list,
        description="Per-step provenance list (10 entries when status != 'error')",
    )


# =============================================================================
# Request Models
# =============================================================================


class QSARSingleRequest(BaseModel):
    """Request body for single-molecule QSAR-Ready processing."""

    model_config = {"from_attributes": True}

    smiles: str = Field(
        ...,
        min_length=1,
        max_length=5000,
        description="SMILES string to process through the QSAR-Ready pipeline",
    )
    config: QSARReadyConfigSchema = Field(
        default_factory=QSARReadyConfigSchema,
        description="Pipeline configuration (defaults to standard QSAR-Ready config)",
    )


# =============================================================================
# Batch Response Models
# =============================================================================


class QSARBatchUploadResponse(BaseModel):
    """Response from submitting a QSAR-Ready batch job."""

    model_config = {"from_attributes": True}

    job_id: str = Field(description="Unique job identifier for tracking progress")
    total_molecules: int = Field(
        description="Number of molecules detected in the uploaded file"
    )
    status: str = Field(default="pending", description="Initial job status")
    message: str = Field(
        default="Job submitted successfully",
        description="Human-readable status message",
    )


class QSARBatchStatusResponse(BaseModel):
    """Current status of a QSAR-Ready batch job."""

    model_config = {"from_attributes": True}

    job_id: str = Field(description="Job identifier")
    status: str = Field(
        description="Job status: 'pending', 'processing', 'complete', or 'failed'"
    )
    progress: int = Field(default=0, description="Progress percentage (0-100)")
    processed: int = Field(default=0, description="Number of molecules processed so far")
    total: int = Field(default=0, description="Total number of molecules in the job")
    eta_seconds: Optional[int] = Field(
        default=None, description="Estimated seconds remaining until job completion"
    )


class QSARBatchSummary(BaseModel):
    """Aggregate statistics for a completed QSAR-Ready batch job."""

    model_config = {"from_attributes": True}

    total: int = Field(description="Total number of input molecules")
    ok: int = Field(description="Molecules successfully curated (status='ok')")
    rejected: int = Field(description="Molecules rejected by pipeline steps (status='rejected')")
    duplicate: int = Field(
        description="Molecules that are InChIKey duplicates of earlier molecules (status='duplicate')"
    )
    error: int = Field(description="Molecules that encountered pipeline errors (status='error')")
    steps_applied_counts: dict[str, int] = Field(
        default_factory=dict,
        description="Per-step count of molecules where the step was applied (status='applied')",
    )


class QSARBatchResultsResponse(BaseModel):
    """Paginated results response for a completed QSAR-Ready batch job."""

    model_config = {"from_attributes": True}

    job_id: str = Field(description="Job identifier")
    status: str = Field(
        description="Job status: 'pending', 'processing', 'complete', or 'failed'"
    )
    config: QSARReadyConfigSchema = Field(
        description="Pipeline configuration used for this batch job"
    )
    summary: QSARBatchSummary = Field(
        description="Aggregate statistics for the completed batch job"
    )
    results: list[QSARReadyResultSchema] = Field(
        default_factory=list,
        description="Paginated list of per-molecule pipeline results",
    )
    page: int = Field(default=1, description="Current page number (1-indexed)")
    per_page: int = Field(default=50, description="Number of results per page")
    total_pages: int = Field(default=1, description="Total number of result pages")
    total_results: int = Field(
        default=0, description="Total number of result records across all pages"
    )
