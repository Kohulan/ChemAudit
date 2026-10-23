"""
ChemAudit Client Models

Pydantic models for API requests and responses.

Models mirror the backend schemas in ``backend/app/schemas``. Unknown fields are
kept (``extra="allow"``) so new server fields never break older clients.
"""

from typing import List, Optional, Dict, Any, Literal
from pydantic import BaseModel, ConfigDict, Field
from enum import Enum


class _APIModel(BaseModel):
    """Base model that keeps fields added by newer server versions."""

    model_config = ConfigDict(extra="allow")


class Severity(str, Enum):
    """Severity level for validation checks."""

    CRITICAL = "critical"
    ERROR = "error"
    WARNING = "warning"
    INFO = "info"
    PASS = "pass"


class CheckResult(_APIModel):
    """Result from a single validation check."""

    check_name: str
    passed: bool
    severity: Severity
    message: str
    affected_atoms: List[int] = Field(default_factory=list)
    details: Dict[str, Any] = Field(default_factory=dict)


class MoleculeInfo(_APIModel):
    """Molecule information."""

    input_smiles: str
    canonical_smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    num_atoms: Optional[int] = None


class ValidationResult(_APIModel):
    """Complete validation result for a molecule."""

    status: str = "completed"
    molecule_info: MoleculeInfo
    overall_score: int = Field(ge=0, le=100)
    issues: List[CheckResult] = Field(default_factory=list)
    all_checks: List[CheckResult] = Field(default_factory=list)
    execution_time_ms: int
    cached: bool = False


class AlertResult(_APIModel):
    """A single structural alert match."""

    pattern_name: str
    description: str
    severity: str
    matched_atoms: List[int] = Field(default_factory=list)
    catalog_source: str
    smarts: Optional[str] = None
    category: Optional[str] = None


class AlertScreeningResult(_APIModel):
    """Complete alert screening result."""

    status: str = "completed"
    molecule_info: MoleculeInfo
    alerts: List[AlertResult] = Field(default_factory=list)
    total_alerts: int = 0
    screened_catalogs: List[str] = Field(default_factory=list)
    has_critical: bool = False
    has_warning: bool = False
    execution_time_ms: int
    educational_note: Optional[str] = None


class ScoreResult(_APIModel):
    """Scoring results. Each score block is None unless it was requested."""

    status: str = "completed"
    molecule_info: MoleculeInfo
    ml_readiness: Optional[Dict[str, Any]] = None
    np_likeness: Optional[Dict[str, Any]] = None
    scaffold: Optional[Dict[str, Any]] = None
    druglikeness: Optional[Dict[str, Any]] = None
    safety_filters: Optional[Dict[str, Any]] = None
    admet: Optional[Dict[str, Any]] = None
    aggregator: Optional[Dict[str, Any]] = None
    consensus: Optional[Dict[str, Any]] = None
    lead_likeness: Optional[Dict[str, Any]] = None
    execution_time_ms: int


class StandardizationDetails(_APIModel):
    """Outcome of the ChEMBL standardization pipeline."""

    original_smiles: str
    standardized_smiles: Optional[str] = None
    success: bool = True
    error_message: Optional[str] = None
    steps_applied: List[Dict[str, Any]] = Field(default_factory=list)
    checker_issues: List[Dict[str, Any]] = Field(default_factory=list)
    excluded_fragments: List[Any] = Field(default_factory=list)
    stereo_comparison: Optional[Dict[str, Any]] = None
    mass_change_percent: float = 0.0


class StandardizationResult(_APIModel):
    """Molecule standardization result."""

    molecule_info: MoleculeInfo
    result: StandardizationDetails
    execution_time_ms: int


class BatchJobStatus(str, Enum):
    """Batch job status."""

    PENDING = "pending"
    PROCESSING = "processing"
    COMPLETE = "complete"
    FAILED = "failed"
    CANCELLED = "cancelled"


class BatchJob(_APIModel):
    """Batch job status information."""

    job_id: str
    status: BatchJobStatus
    progress: int = Field(ge=0, le=100)
    processed: int
    total: int
    eta_seconds: Optional[int] = None
    error_message: Optional[str] = None


class BatchResultItem(_APIModel):
    """Single molecule result in batch."""

    smiles: str
    name: Optional[str] = None
    index: int
    status: Literal["success", "error"]
    error: Optional[str] = None
    validation: Optional[Dict[str, Any]] = None
    alerts: Optional[Dict[str, Any]] = None
    scoring: Optional[Dict[str, Any]] = None
    standardization: Optional[Dict[str, Any]] = None
    profiling: Optional[Dict[str, Any]] = None
    safety_assessment: Optional[Dict[str, Any]] = None


class BatchStatistics(_APIModel):
    """Aggregate statistics for batch job."""

    total: int
    successful: int
    errors: int
    avg_validation_score: Optional[float] = None
    avg_ml_readiness_score: Optional[float] = None
    score_distribution: Dict[str, int] = Field(default_factory=dict)
    alert_summary: Dict[str, int] = Field(default_factory=dict)
    processing_time_seconds: Optional[float] = None


class BatchResult(_APIModel):
    """Batch processing results with pagination."""

    job_id: str
    status: BatchJobStatus
    statistics: Optional[BatchStatistics] = None
    results: List[BatchResultItem] = Field(default_factory=list)
    page: int = Field(ge=1)
    page_size: int = Field(ge=1, le=100)
    total_results: int
    total_pages: int


class BatchUploadResponse(_APIModel):
    """Response from batch upload."""

    job_id: str
    status: str = "pending"
    total_molecules: int
    message: str = "Job submitted successfully"
