"""
Batch Processing Schemas

Pydantic models for batch upload, job status, and results.
"""
from typing import Optional, List, Dict, Any, Literal
from pydantic import BaseModel, Field


class BatchUploadResponse(BaseModel):
    """Response from batch file upload."""

    job_id: str = Field(..., description="Unique job identifier for tracking")
    status: str = Field(default="pending", description="Initial job status")
    total_molecules: int = Field(..., description="Number of molecules detected in file")
    message: str = Field(default="Job submitted successfully")


class BatchJobStatus(BaseModel):
    """Current status of a batch job."""

    job_id: str
    status: Literal["pending", "processing", "complete", "failed", "cancelled"]
    progress: int = Field(..., ge=0, le=100, description="Progress percentage")
    processed: int = Field(..., description="Number of molecules processed")
    total: int = Field(..., description="Total number of molecules")
    eta_seconds: Optional[int] = Field(None, description="Estimated time remaining in seconds")
    error_message: Optional[str] = Field(None, description="Error message if failed")


class BatchResultItem(BaseModel):
    """Single molecule result in batch."""

    smiles: str
    name: Optional[str] = None
    index: int
    status: Literal["success", "error"]
    error: Optional[str] = None
    validation: Optional[Dict[str, Any]] = Field(
        None, description="Validation results with overall_score and issues"
    )
    alerts: Optional[Dict[str, Any]] = Field(
        None, description="Structural alert screening results"
    )
    scoring: Optional[Dict[str, Any]] = Field(
        None, description="ML-readiness scoring results"
    )


class BatchStatistics(BaseModel):
    """Aggregate statistics for a batch job."""

    total: int = Field(..., description="Total molecules processed")
    successful: int = Field(..., description="Molecules processed successfully")
    errors: int = Field(..., description="Molecules with errors")
    avg_validation_score: Optional[float] = Field(
        None, description="Average validation score"
    )
    avg_ml_readiness_score: Optional[float] = Field(
        None, description="Average ML-readiness score"
    )
    score_distribution: Dict[str, int] = Field(
        default_factory=dict,
        description="Score distribution buckets (excellent, good, moderate, poor)",
    )
    alert_summary: Dict[str, int] = Field(
        default_factory=dict, description="Alert counts by catalog type"
    )
    processing_time_seconds: Optional[float] = Field(
        None, description="Total processing time in seconds"
    )


class BatchResultsResponse(BaseModel):
    """Paginated batch results response."""

    job_id: str
    status: str
    statistics: Optional[BatchStatistics] = None
    results: List[BatchResultItem] = Field(default_factory=list)
    page: int = Field(..., ge=1)
    page_size: int = Field(..., ge=1, le=100)
    total_results: int
    total_pages: int


class CSVColumnsResponse(BaseModel):
    """Response with detected CSV columns for SMILES selection."""

    columns: List[str] = Field(..., description="Available column names")
    suggested_smiles: Optional[str] = Field(
        None, description="Suggested SMILES column based on name"
    )
    row_count_estimate: int = Field(..., description="Estimated number of rows")


class BatchProgressMessage(BaseModel):
    """WebSocket progress message format."""

    job_id: str
    status: str
    progress: int = Field(..., ge=0, le=100)
    processed: int
    total: int
    eta_seconds: Optional[int] = None
    error_message: Optional[str] = None
