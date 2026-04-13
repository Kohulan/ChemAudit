"""
Pydantic v2 schemas for Dataset Intelligence API (Phase 12).

Provides request and response models for:
- POST /dataset/upload       -- upload CSV/SDF for async health audit
- GET  /dataset/{job_id}/status   -- poll job progress
- GET  /dataset/{job_id}/results  -- fetch audit results
- POST /dataset/{job_id}/diff     -- compare two datasets
- GET  /dataset/{job_id}/download/report -- download curation report JSON
- GET  /dataset/{job_id}/download/csv    -- download curated CSV
"""

from typing import Optional

from pydantic import BaseModel, Field

# =============================================================================
# Upload response
# =============================================================================


class DatasetUploadResponse(BaseModel):
    """Response from POST /dataset/upload."""

    job_id: str = Field(description="Unique job identifier for async tracking")
    filename: str = Field(description="Original uploaded filename")
    file_type: str = Field(description="Detected file type: 'csv' or 'sdf'")
    status: str = Field(default="pending", description="Initial job status")
    message: str = Field(
        default="Dataset audit job submitted",
        description="Human-readable status message",
    )


# =============================================================================
# Status response
# =============================================================================


class DatasetAuditStatusResponse(BaseModel):
    """Response from GET /dataset/{job_id}/status."""

    job_id: str = Field(description="Job identifier")
    status: str = Field(
        description="Job status: 'pending', 'processing', 'complete', or 'error'"
    )
    progress: float = Field(default=0.0, description="Progress percentage (0-100)")
    current_stage: Optional[str] = Field(
        default=None, description="Name of the current processing stage"
    )
    eta_seconds: Optional[float] = Field(
        default=None, description="Estimated time remaining in seconds"
    )


# =============================================================================
# Health audit result models
# =============================================================================


class SubScoreDetail(BaseModel):
    """Detail for a single health sub-score component."""

    name: str = Field(description="Sub-score name (e.g. 'parsability', 'stereo')")
    score: float = Field(description="Sub-score value (0-1)")
    weight: float = Field(description="Weight applied to this sub-score")
    count: int = Field(description="Number of affected molecules")
    total: int = Field(description="Total molecules evaluated for this sub-score")


class HealthAuditResult(BaseModel):
    """Health audit results from compute_health_score."""

    overall_score: float = Field(description="Composite health score (0-100)")
    sub_scores: list[SubScoreDetail] = Field(
        description="Breakdown of 5 sub-score components"
    )
    weights: dict[str, float] = Field(description="Weight vector used for scoring")
    molecule_count: int = Field(description="Total molecules in dataset")
    issues: list[dict] = Field(
        default_factory=list, description="Per-molecule issue records"
    )
    property_distributions: dict = Field(
        default_factory=dict,
        description="MW/LogP/TPSA histograms with drug-space reference",
    )
    std_pipeline_comparison: dict = Field(
        default_factory=dict,
        description="Standardization pipeline agreement statistics",
    )
    std_sample_size: int = Field(
        default=0, description="Number of molecules sampled for std consistency"
    )
    dedup_groups: list[dict] = Field(
        default_factory=list, description="Groups of duplicate InChIKeys"
    )


# =============================================================================
# Contradictory label result
# =============================================================================


class ContradictoryLabelResult(BaseModel):
    """A single contradictory label group."""

    inchikey: str = Field(description="InChIKey shared by conflicting entries")
    entries: list[dict] = Field(
        description="Individual entries with row_index, smiles, activity"
    )
    fold_difference: float = Field(
        description="Fold-difference between max and min activity"
    )
    entry_count: int = Field(description="Number of entries in this group")
    smiles: str = Field(description="Representative SMILES for this group")


# =============================================================================
# Full audit response
# =============================================================================


class DatasetAuditResponse(BaseModel):
    """Response from GET /dataset/{job_id}/results."""

    job_id: str = Field(description="Job identifier")
    status: str = Field(description="Job status")
    health_audit: Optional[HealthAuditResult] = Field(
        default=None, description="Health audit results (None if not yet complete)"
    )
    contradictions: list[ContradictoryLabelResult] = Field(
        default_factory=list,
        description="Contradictory label groups detected",
    )
    numeric_columns: list[dict] = Field(
        default_factory=list,
        description="Detected numeric columns with priority ranking",
    )
    curation_report: Optional[dict] = Field(
        default=None, description="Full curation report (JSON-serializable dict)"
    )
    curated_csv_available: bool = Field(
        default=False,
        description="Whether a curated CSV download is available",
    )


# =============================================================================
# Diff request/response
# =============================================================================


class DatasetDiffRequest(BaseModel):
    """Request body for POST /dataset/{job_id}/diff."""

    job_id: str = Field(description="Primary dataset job ID")


class DiffMolecule(BaseModel):
    """A molecule entry in a dataset diff result."""

    inchikey: str = Field(description="InChIKey identifier")
    smiles: str = Field(description="SMILES representation")
    row_index: int = Field(description="Row index in source dataset")
    properties: dict = Field(
        default_factory=dict, description="Molecule properties"
    )
    changes: list[dict] = Field(
        default_factory=list,
        description="Property changes (only for modified molecules)",
    )


class DatasetDiffResponse(BaseModel):
    """Response from POST /dataset/{job_id}/diff."""

    added: list[DiffMolecule] = Field(
        description="Molecules in comparison but not in primary"
    )
    removed: list[DiffMolecule] = Field(
        description="Molecules in primary but not in comparison"
    )
    modified: list[DiffMolecule] = Field(
        description="Molecules in both with property changes"
    )
    added_count: int = Field(description="Number of added molecules")
    removed_count: int = Field(description="Number of removed molecules")
    modified_count: int = Field(description="Number of modified molecules")
    unchanged_count: int = Field(description="Number of unchanged molecules")
    unique_columns_primary: int = Field(
        default=0, description="Columns unique to the primary dataset"
    )
    unique_columns_comparison: int = Field(
        default=0, description="Columns unique to the comparison dataset"
    )
