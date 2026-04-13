"""
Pydantic v2 schemas for Structure Filter API (Phase 11).

Provides request and response models for:
- POST /structure-filter/filter        — multi-stage funnel pipeline
- POST /structure-filter/score         — composite 0-1 scorer per SMILES
- POST /structure-filter/reinvent-score — REINVENT 4-compatible scoring API
- Batch upload, status, results, download
"""

from typing import Optional

from pydantic import BaseModel, Field

# =============================================================================
# Config Schema
# =============================================================================


class FilterConfigSchema(BaseModel):
    """Pydantic mirror of FilterConfig dataclass for API serialization.

    Specifies property thresholds, alert catalog selection, novelty settings,
    and composite scorer weight vectors for the structure filter pipeline.
    """

    model_config = {"from_attributes": True}

    # Property thresholds (per D-12)
    min_mw: float = Field(default=200.0, description="Minimum molecular weight (Da)")
    max_mw: float = Field(default=500.0, description="Maximum molecular weight (Da)")
    min_logp: float = Field(default=-1.0, description="Minimum LogP")
    max_logp: float = Field(default=5.0, description="Maximum LogP")
    max_tpsa: float = Field(
        default=140.0, description="Maximum topological polar surface area (Å²)"
    )
    max_rot_bonds: int = Field(default=10, description="Maximum rotatable bonds")
    max_rings: Optional[int] = Field(
        default=None, description="Maximum ring count (None = no limit)"
    )
    max_sa_score: float = Field(default=5.0, description="Maximum SA Score (1=easy, 10=hard)")
    # Alert catalog selection
    use_pains: bool = Field(default=True, description="Screen against PAINS catalog")
    use_brenk: bool = Field(default=True, description="Screen against Brenk catalog")
    use_kazius: bool = Field(default=True, description="Screen against Kazius mutagenicity rules")
    use_nibr: bool = Field(default=False, description="Screen against NIBR Novartis deck filters")
    # Novelty settings
    enable_novelty: bool = Field(
        default=False, description="Enable ChEMBL Tanimoto novelty filter"
    )
    novelty_threshold: float = Field(
        default=0.85,
        description="Maximum Tanimoto similarity to ChEMBL compounds (lower = more novel)",
    )
    # Composite scorer weights (per D-15)
    weight_validity: float = Field(default=0.3, description="Weight for validity component")
    weight_qed: float = Field(default=0.3, description="Weight for QED component")
    weight_alert_free: float = Field(
        default=0.2, description="Weight for alert-free binary component"
    )
    weight_sa: float = Field(default=0.2, description="Weight for SA Score component")


# =============================================================================
# Filter endpoint schemas
# =============================================================================


class FilterRequest(BaseModel):
    """Request body for /structure-filter/filter endpoint.

    Accepts a SMILES list with either a named preset or an explicit config.
    If preset is provided, it overrides config.
    """

    smiles_list: list[str] = Field(description="List of SMILES strings to filter")
    config: Optional[FilterConfigSchema] = Field(
        default=None, description="Explicit filter configuration (overridden by preset)"
    )
    preset: Optional[str] = Field(
        default=None,
        description="Named preset: 'drug_like', 'lead_like', 'fragment_like', 'permissive'",
    )


class StageResultSchema(BaseModel):
    """Result summary for a single pipeline stage."""

    model_config = {"from_attributes": True}

    stage_name: str = Field(description="Stage identifier")
    stage_index: int = Field(description="1-based stage index")
    input_count: int = Field(description="Molecules entering this stage")
    passed_count: int = Field(description="Molecules passing this stage")
    rejected_count: int = Field(description="Molecules rejected at this stage")
    enabled: bool = Field(default=True, description="Whether this stage was enabled")


class MoleculeResultSchema(BaseModel):
    """Per-molecule result from the filter pipeline."""

    model_config = {"from_attributes": True}

    smiles: str = Field(description="Input SMILES string")
    status: str = Field(description="Outcome: 'passed', 'rejected', 'duplicate', or 'error'")
    failed_at: Optional[str] = Field(
        default=None, description="Stage name where molecule was rejected (None if passed)"
    )
    rejection_reason: Optional[str] = Field(
        default=None, description="Human-readable rejection reason"
    )


class FilterResponse(BaseModel):
    """Response from /structure-filter/filter endpoint."""

    model_config = {"from_attributes": True}

    input_count: int = Field(description="Total input SMILES count")
    output_count: int = Field(description="Molecules passing all stages")
    stages: list[StageResultSchema] = Field(description="Per-stage funnel results (6 or 7 stages)")
    molecules: list[MoleculeResultSchema] = Field(description="Per-molecule results")


# =============================================================================
# Score endpoint schemas
# =============================================================================


class ScoreRequest(BaseModel):
    """Request body for /structure-filter/score endpoint."""

    smiles_list: list[str] = Field(description="List of SMILES strings to score")
    preset: Optional[str] = Field(
        default="drug_like",
        description="Named preset for weight vector: 'drug_like', 'lead_like', etc.",
    )


class ScoreResponse(BaseModel):
    """Response from /structure-filter/score endpoint.

    Scores are in [0, 1] range. None (JSON null) for unparseable SMILES (D-14).
    """

    scores: list[Optional[float]] = Field(
        description="Composite scores (0-1) in input order. null for invalid SMILES."
    )


# =============================================================================
# REINVENT endpoint schemas
# =============================================================================


class REINVENTInput(BaseModel):
    """Single input item for the REINVENT 4 scoring endpoint."""

    input_string: str = Field(description="SMILES string to score")
    query_id: str = Field(description="Identifier for this scoring request")


class REINVENTSuccessItem(BaseModel):
    """Single item in the successes_list of a REINVENT scoring response.

    Only valid SMILES appear here — invalid SMILES are omitted (Pitfall 4).
    """

    query_id: str = Field(description="Identifier matching the input query_id")
    output_value: float = Field(description="Composite score in [0, 1] range")


class REINVENTOutput(BaseModel):
    """Inner output wrapper for REINVENT 4 response format."""

    successes_list: list[REINVENTSuccessItem] = Field(
        description="Scored items. Items with invalid SMILES are omitted, not scored 0.0."
    )


class REINVENTResponse(BaseModel):
    """Response from /structure-filter/reinvent-score endpoint.

    Matches the REINVENT 4 component contract exactly:
    {output: {successes_list: [{query_id, output_value}]}}
    """

    output: REINVENTOutput = Field(description="REINVENT 4 output wrapper")


# =============================================================================
# Batch endpoint schemas
# =============================================================================


class StructureFilterBatchUploadResponse(BaseModel):
    """Response from submitting a structure filter async batch job."""

    job_id: str = Field(description="Unique job identifier for tracking progress")
    total_molecules: int = Field(description="Number of molecules submitted")
    status: str = Field(default="pending", description="Initial job status")


class StructureFilterBatchStatusResponse(BaseModel):
    """Current status of a structure filter async batch job."""

    job_id: str = Field(description="Job identifier")
    status: str = Field(
        description="Job status: 'pending', 'processing', 'complete', or 'failed'"
    )
    progress: Optional[float] = Field(
        default=None, description="Progress percentage (0-100)"
    )
    current_stage: Optional[str] = Field(
        default=None, description="Name of the currently processing stage"
    )


class StructureFilterBatchResultsResponse(BaseModel):
    """Results response for a completed structure filter async batch job."""

    job_id: str = Field(description="Job identifier")
    status: str = Field(description="Job status")
    result: Optional[FilterResponse] = Field(
        default=None, description="Full filter result (None if not yet complete)"
    )
