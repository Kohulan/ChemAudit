"""
Pydantic v2 schemas for the QSAR-Ready Pipeline API.

These schemas mirror the backend dataclasses (QSARReadyConfig, QSARStepResult,
QSARReadyResult) for use in FastAPI request/response models.
"""

from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field


# =============================================================================
# Config Schema
# =============================================================================


class QSARReadyConfigSchema(BaseModel):
    """Mirrors QSARReadyConfig dataclass — configures which pipeline steps to run."""

    enable_metals: bool = True
    enable_desalt: bool = True
    enable_normalize: bool = True
    enable_neutralize: bool = True
    enable_tautomer: bool = True
    enable_stereo_strip: bool = False
    enable_isotope_strip: bool = True
    min_heavy_atoms: int = Field(default=3, ge=0)
    max_heavy_atoms: int = Field(default=100, ge=1)
    max_mw: float = Field(default=1500.0, gt=0.0)
    remove_inorganics: bool = True


# =============================================================================
# Result Schemas
# =============================================================================


class QSARStepResultSchema(BaseModel):
    """Schema for a single pipeline step result."""

    step_name: str
    step_index: int
    enabled: bool
    status: str  # "applied" | "no_change" | "skipped" | "error"
    before_smiles: Optional[str] = None
    after_smiles: Optional[str] = None
    detail: Optional[str] = None


class QSARReadyResultSchema(BaseModel):
    """Schema for the result of processing a single molecule through the pipeline."""

    original_smiles: str
    original_inchikey: Optional[str] = None
    curated_smiles: Optional[str] = None
    standardized_inchikey: Optional[str] = None
    inchikey_changed: bool = False
    status: str  # "ok" | "rejected" | "duplicate" | "error"
    rejection_reason: Optional[str] = None
    steps: List[QSARStepResultSchema] = Field(default_factory=list)


# =============================================================================
# Request Schemas
# =============================================================================


class QSARSingleRequest(BaseModel):
    """Request body for POST /api/v1/qsar-ready/single."""

    smiles: str
    config: QSARReadyConfigSchema = Field(default_factory=QSARReadyConfigSchema)


# =============================================================================
# Batch Summary Schema
# =============================================================================


class QSARBatchSummary(BaseModel):
    """Summary statistics for a batch QSAR job."""

    total: int
    ok: int
    rejected: int
    duplicate: int
    error: int
    steps_applied_counts: Dict[str, int] = Field(default_factory=dict)


# =============================================================================
# Batch Response Schemas
# =============================================================================


class QSARBatchUploadResponse(BaseModel):
    """Response for POST /api/v1/qsar-ready/batch/upload."""

    job_id: str
    total_molecules: int
    status: str  # "pending"
    message: str


class QSARBatchStatusResponse(BaseModel):
    """Response for GET /api/v1/qsar-ready/batch/{job_id}/status."""

    job_id: str
    status: str  # "pending" | "processing" | "complete" | "failed" | "cancelled"
    progress: int  # 0-100 percentage
    processed: int
    total: int
    eta_seconds: Optional[int] = None


class QSARBatchResultsResponse(BaseModel):
    """Response for GET /api/v1/qsar-ready/batch/{job_id}/results."""

    job_id: str
    status: str
    config: QSARReadyConfigSchema
    summary: QSARBatchSummary
    results: List[QSARReadyResultSchema]
    page: int
    per_page: int
    total_pages: int
    total_results: int
