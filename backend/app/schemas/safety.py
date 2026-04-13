"""
Safety Assessment Schemas

Pydantic schemas for CYP soft-spot, hERG, bRo5, REOS, and complexity
percentile assessment requests and responses.
"""

from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field, field_validator

from app.schemas.validation import MoleculeInfo


class SafetyAssessRequest(BaseModel):
    """Request for full safety assessment."""

    molecule: str = Field(..., min_length=1, max_length=10000)
    format: str = Field(default="auto", pattern="^(auto|smiles|inchi|mol)$")

    @field_validator("molecule")
    @classmethod
    def sanitize_molecule_input(cls, v: str) -> str:
        """Sanitize molecule input to prevent injection attacks."""
        dangerous = ["<", ">", "&", ";", "|", "$", "`"]
        if any(c in v for c in dangerous):
            raise ValueError("Invalid characters in molecule string")
        return v.strip()


class CypSiteSchema(BaseModel):
    """A single CYP metabolism soft-spot site match."""

    site_name: str
    reaction_type: str
    matched_atoms: List[int]


class CypResultSchema(BaseModel):
    """CYP soft-spot prediction result."""

    sites: List[CypSiteSchema] = Field(default_factory=list)
    n_sites: int = Field(default=0)


class HergResultSchema(BaseModel):
    """hERG liability risk assessment result."""

    herg_risk: str  # "low", "moderate", "high"
    risk_score: int
    max_score: int = 4
    flags: List[str] = Field(default_factory=list)
    descriptors: Dict[str, Any] = Field(default_factory=dict)


class Bro5ViolationSchema(BaseModel):
    """A single bRo5 property violation."""

    property: str
    value: float
    threshold: float
    direction: str


class Bro5ResultSchema(BaseModel):
    """Beyond-Rule-of-5 assessment result."""

    applicable: bool
    passed: bool
    message: Optional[str] = None
    violations: List[Bro5ViolationSchema] = Field(default_factory=list)
    values: Dict[str, float] = Field(default_factory=dict)


class ReosViolationSchema(BaseModel):
    """A single REOS property violation."""

    property: str
    value: float
    range: List[float]
    exceeded: bool


class ReosResultSchema(BaseModel):
    """REOS filter assessment result."""

    passed: bool
    violations: List[ReosViolationSchema] = Field(default_factory=list)
    n_violations: int = 0
    descriptors: Dict[str, float] = Field(default_factory=dict)


class ComplexityPropertySchema(BaseModel):
    """Assessment for a single complexity property."""

    value: float
    p5: float
    p95: float
    outlier: bool
    direction: Optional[str] = None  # "below", "above", or None


class ComplexityResultSchema(BaseModel):
    """Complexity percentile filter assessment result."""

    properties: Dict[str, ComplexityPropertySchema] = Field(default_factory=dict)
    n_outliers: int = 0
    outlier_properties: List[str] = Field(default_factory=list)
    within_range: bool = True


class SafetyAssessResponse(BaseModel):
    """Full safety assessment response."""

    status: str = Field(default="completed")
    molecule_info: MoleculeInfo
    cyp_softspots: CypResultSchema
    herg: HergResultSchema
    bro5: Bro5ResultSchema
    reos: ReosResultSchema
    complexity: ComplexityResultSchema
    execution_time_ms: int = Field(default=0)


class SafetySummaryResponse(BaseModel):
    """Lightweight safety summary for SV/Profiler badges."""

    status: str = Field(default="completed")
    total_alerts: int = Field(default=0)
    has_critical: bool = Field(default=False)
    cyp_status: str = Field(default="default")  # always "default" (informational)
    herg_status: str = Field(default="success")  # "success", "warning", "error"
    bro5_status: str = Field(default="success")  # "success", "error", "default" (N/A)
    reos_status: str = Field(default="success")  # "success", "warning", "error"
    complexity_outliers: int = Field(default=0)
