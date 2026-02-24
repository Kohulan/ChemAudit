"""
Batch Processing API Routes

Endpoints for batch file upload, job status, results retrieval, and analytics.
"""

import uuid
from typing import Any, List, Optional

from fastapi import (
    APIRouter,
    Body,
    Depends,
    File,
    Form,
    HTTPException,
    Query,
    Request,
    UploadFile,
)
from pydantic import BaseModel

from app.core.config import settings
from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.schemas.analytics import (
    AnalysisStatus,
    AnalyticsTriggerResponse,
    BatchAnalyticsResponse,
)
from app.schemas.batch import (
    BatchJobStatus,
    BatchResultItem,
    BatchResultsResponse,
    BatchStatistics,
    BatchUploadResponse,
    CSVColumnsResponse,
)
from app.services.analytics.storage import analytics_storage
from app.services.batch.analytics_tasks import run_expensive_analytics
from app.services.batch.file_parser import (
    detect_csv_columns,
    parse_csv,
    parse_sdf,
    validate_file_content_type,
)
from app.services.batch.progress_tracker import progress_tracker
from app.services.batch.result_aggregator import result_storage
from app.services.batch.subset_actions import (
    export_subset,
    rescore_subset,
    revalidate_subset,
)
from app.services.batch.tasks import cancel_batch_job, process_batch_job

router = APIRouter()


@router.post("/batch/upload", response_model=BatchUploadResponse)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def upload_batch(
    request: Request,
    file: UploadFile = File(..., description="SDF or CSV file to process"),
    smiles_column: Optional[str] = Form(
        default="SMILES",
        description="Column name containing SMILES (for CSV files)",
    ),
    name_column: Optional[str] = Form(
        default=None,
        description="Column name containing molecule names/IDs (for CSV files)",
    ),
    include_extended_safety: bool = Form(
        default=False,
        description="Include NIH and ZINC structural alert filters",
    ),
    include_chembl_alerts: bool = Form(
        default=False,
        description="Include ChEMBL pharma company structural alert filters",
    ),
    include_standardization: bool = Form(
        default=False,
        description="Run ChEMBL standardization pipeline on each molecule",
    ),
    notification_email: Optional[str] = Form(
        default=None,
        description="Email address for batch completion notification (overrides global setting)",
    ),
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Upload a file for batch processing.

    Accepts SDF or CSV files with up to 1,000,000 molecules.
    Returns a job_id for tracking progress and retrieving results.

    - **file**: SDF (.sdf) or CSV (.csv) file
    - **smiles_column**: For CSV files, the column containing SMILES strings
    - **name_column**: For CSV files, the column containing molecule names/IDs

    Returns job_id immediately. Use `/batch/{job_id}` to check status
    and `/ws/batch/{job_id}` for real-time progress via WebSocket.
    """
    # Validate file extension
    filename = file.filename or ""
    filename_lower = filename.lower()

    # Supported text-based formats: .csv, .tsv, .txt (all parsed as delimited text)
    text_formats = (".csv", ".tsv", ".txt")
    if not (filename_lower.endswith(".sdf") or filename_lower.endswith(text_formats)):
        raise HTTPException(
            status_code=400,
            detail="Invalid file type. Supported formats: .sdf, .csv, .tsv, .txt",
        )

    # Read file content
    try:
        content = await file.read()
    except Exception:
        raise HTTPException(
            status_code=400,
            detail="Failed to read file",
        )

    # Security: Check file size
    file_size_mb = len(content) / (1024 * 1024)
    if file_size_mb > settings.MAX_FILE_SIZE_MB:
        raise HTTPException(
            status_code=400,
            detail=f"File too large: {file_size_mb:.1f}MB exceeds limit of {settings.MAX_FILE_SIZE_MB}MB",
        )

    # Security: Validate file content matches extension
    # .csv, .tsv, and .txt are all treated as delimited text files
    expected_type = "sdf" if filename_lower.endswith(".sdf") else "csv"
    is_valid, error_msg = validate_file_content_type(content, expected_type, filename)
    if not is_valid:
        raise HTTPException(
            status_code=400,
            detail=error_msg or "Invalid file content",
        )

    # Parse file
    try:
        if filename_lower.endswith(".sdf"):
            molecules = parse_sdf(content, max_file_size_mb=settings.MAX_FILE_SIZE_MB)
        else:  # .csv, .tsv, .txt - all parsed as delimited text
            molecules = parse_csv(
                content,
                smiles_column=smiles_column or "SMILES",
                name_column=name_column,
                max_file_size_mb=settings.MAX_FILE_SIZE_MB,
            )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception:
        raise HTTPException(
            status_code=400,
            detail="Failed to parse file. Please check the file format.",
        )

    # Validate molecule count
    if len(molecules) == 0:
        raise HTTPException(
            status_code=400,
            detail="No valid molecules found in file",
        )

    if len(molecules) > settings.MAX_BATCH_SIZE:
        raise HTTPException(
            status_code=400,
            detail=f"File contains {len(molecules)} molecules. "
            f"Maximum allowed is {settings.MAX_BATCH_SIZE}.",
        )

    # Generate job ID and start processing
    job_id = str(uuid.uuid4())

    # Convert MoleculeData objects to dicts for Celery serialization
    mol_dicts = [
        {
            "smiles": m.smiles,
            "name": m.name,
            "index": m.index,
            "properties": m.properties,
            "parse_error": m.parse_error,
        }
        for m in molecules
    ]

    # Build safety screening options
    safety_options = {
        "include_extended": include_extended_safety,
        "include_chembl": include_chembl_alerts,
        "include_standardization": include_standardization,
    }

    # Start batch processing
    process_batch_job(job_id, mol_dicts, safety_options=safety_options)

    # Store notification email for this job (if provided or globally configured)
    email_target = notification_email or settings.NOTIFICATION_EMAIL
    if email_target:
        from app.services.batch.progress_tracker import progress_tracker as _pt

        r = _pt._get_redis()
        r.set(f"batch:email:{job_id}", email_target, ex=settings.BATCH_RESULT_TTL)

    return BatchUploadResponse(
        job_id=job_id,
        status="pending",
        total_molecules=len(molecules),
        message=f"Job submitted. Processing {len(molecules)} molecules.",
    )


@router.get("/batch/{job_id}", response_model=BatchResultsResponse)
@limiter.limit("60/minute", key_func=get_rate_limit_key)
async def get_batch_results(
    request: Request,
    job_id: str,
    page: int = Query(default=1, ge=1, description="Page number"),
    page_size: int = Query(default=50, ge=1, le=100, description="Results per page"),
    status_filter: Optional[str] = Query(
        default=None,
        description="Filter by status (success, error)",
    ),
    min_score: Optional[int] = Query(
        default=None, ge=0, le=100, description="Minimum validation score"
    ),
    max_score: Optional[int] = Query(
        default=None, ge=0, le=100, description="Maximum validation score"
    ),
    sort_by: Optional[str] = Query(
        default=None,
        description="Sort field (index, name, smiles, score, qed, safety, status, issues)",
    ),
    sort_dir: Optional[str] = Query(
        default=None,
        description="Sort direction (asc, desc)",
    ),
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Get batch job status and results.

    - **job_id**: Job identifier from upload response
    - **page**: Page number for results (1-indexed)
    - **page_size**: Results per page (max 100)
    - **status_filter**: Filter results by 'success' or 'error'
    - **min_score**: Filter by minimum validation score
    - **max_score**: Filter by maximum validation score

    Returns job status and paginated results with statistics.
    """
    # Get job status
    progress_info = progress_tracker.get_progress(job_id)

    if not progress_info:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    # Get statistics
    stats_data = result_storage.get_statistics(job_id)
    statistics = None
    if stats_data:
        statistics = BatchStatistics(
            total=stats_data.total,
            successful=stats_data.successful,
            errors=stats_data.errors,
            avg_validation_score=stats_data.avg_validation_score,
            avg_ml_readiness_score=stats_data.avg_ml_readiness_score,
            avg_qed_score=stats_data.avg_qed_score,
            avg_sa_score=stats_data.avg_sa_score,
            lipinski_pass_rate=stats_data.lipinski_pass_rate,
            safety_pass_rate=stats_data.safety_pass_rate,
            score_distribution=stats_data.score_distribution,
            alert_summary=stats_data.alert_summary,
            issue_summary=stats_data.issue_summary,
            processing_time_seconds=stats_data.processing_time_seconds,
        )

    # Get paginated results
    result_data = result_storage.get_results(
        job_id=job_id,
        page=page,
        page_size=page_size,
        status_filter=status_filter,
        min_score=min_score,
        max_score=max_score,
        sort_by=sort_by,
        sort_dir=sort_dir,
    )

    # Convert to response model
    results = [
        BatchResultItem(
            smiles=r.get("smiles", ""),
            name=r.get("name"),
            index=r.get("index", 0),
            status=r.get("status", "error"),
            error=r.get("error"),
            validation=r.get("validation"),
            alerts=r.get("alerts"),
            scoring=r.get("scoring"),
            standardization=r.get("standardization"),
        )
        for r in result_data.get("results", [])
    ]

    return BatchResultsResponse(
        job_id=job_id,
        status=progress_info.status,
        statistics=statistics,
        results=results,
        page=result_data.get("page", page),
        page_size=result_data.get("page_size", page_size),
        total_results=result_data.get("total_results", 0),
        total_pages=result_data.get("total_pages", 0),
    )


@router.get("/batch/{job_id}/status", response_model=BatchJobStatus)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def get_batch_status(
    request: Request, job_id: str, api_key: Optional[str] = Depends(get_api_key)
):
    """
    Get current status of a batch job.

    Returns progress information without results data.
    Lighter endpoint for polling status.
    """
    progress_info = progress_tracker.get_progress(job_id)

    if not progress_info:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    return BatchJobStatus(
        job_id=progress_info.job_id,
        status=progress_info.status,
        progress=progress_info.progress,
        processed=progress_info.processed,
        total=progress_info.total,
        eta_seconds=progress_info.eta_seconds,
        error_message=progress_info.error_message,
    )


@router.get("/batch/{job_id}/stats", response_model=BatchStatistics)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def get_batch_stats(
    request: Request, job_id: str, api_key: Optional[str] = Depends(get_api_key)
):
    """
    Get statistics for a completed batch job.

    Returns aggregate statistics without individual results.
    """
    stats_data = result_storage.get_statistics(job_id)

    if not stats_data:
        raise HTTPException(
            status_code=404,
            detail=f"Statistics not found for job {job_id}. Job may still be processing.",
        )

    return BatchStatistics(
        total=stats_data.total,
        successful=stats_data.successful,
        errors=stats_data.errors,
        avg_validation_score=stats_data.avg_validation_score,
        avg_ml_readiness_score=stats_data.avg_ml_readiness_score,
        avg_qed_score=stats_data.avg_qed_score,
        avg_sa_score=stats_data.avg_sa_score,
        lipinski_pass_rate=stats_data.lipinski_pass_rate,
        safety_pass_rate=stats_data.safety_pass_rate,
        score_distribution=stats_data.score_distribution,
        alert_summary=stats_data.alert_summary,
        issue_summary=stats_data.issue_summary,
        processing_time_seconds=stats_data.processing_time_seconds,
    )


@router.delete("/batch/{job_id}")
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def cancel_batch(
    request: Request, job_id: str, api_key: Optional[str] = Depends(get_api_key)
):
    """
    Cancel a batch job.

    Stops processing and marks job as cancelled.
    Already completed results are retained.
    """
    progress_info = progress_tracker.get_progress(job_id)

    if not progress_info:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    if progress_info.status in ("complete", "cancelled", "failed"):
        return {
            "job_id": job_id,
            "status": progress_info.status,
            "message": f"Job already {progress_info.status}",
        }

    # Cancel the job
    cancel_batch_job.delay(job_id)

    return {
        "job_id": job_id,
        "status": "cancelling",
        "message": "Job cancellation requested",
    }


@router.get("/batch/{job_id}/analytics", response_model=BatchAnalyticsResponse)
@limiter.limit("120/minute", key_func=get_rate_limit_key)
async def get_batch_analytics(
    request: Request,
    job_id: str,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Get analytics status and cached results for a completed batch job.

    Returns status for all analytics types and any results that are ready.
    Analytics are computed asynchronously after batch processing completes.

    - **job_id**: Batch job identifier from the upload response.
    """
    status_data = analytics_storage.get_status(job_id)
    if status_data is None:
        raise HTTPException(
            status_code=404,
            detail=f"Analytics not available for job {job_id}. "
            "Job may still be processing or results have expired.",
        )

    # Build validated status dict
    status_dict: dict[str, AnalysisStatus] = {
        analysis_type: AnalysisStatus(**entry)
        for analysis_type, entry in status_data.items()
    }

    # Fetch completed results
    _RESULT_FIELD_MAP = {
        "deduplication": "deduplication",
        "scaffold": "scaffold",
        "chemical_space": "chemical_space",
        "similarity_search": "similarity_matrix",
        "mmp": "mmp",
        "statistics": "statistics",
    }

    result_kwargs: dict[str, Any] = {}
    for analysis_type, field_name in _RESULT_FIELD_MAP.items():
        entry = status_dict.get(analysis_type)
        if entry and entry.status == "complete":
            raw = analytics_storage.get_result(job_id, analysis_type)
            if raw is not None:
                result_kwargs[field_name] = raw

    return BatchAnalyticsResponse(
        job_id=job_id,
        status=status_dict,
        **result_kwargs,
    )


_ALLOWED_EXPENSIVE_ANALYSES = {
    "scaffold",
    "chemical_space",
    "mmp",
    "similarity_search",
    "rgroup",
}


@router.post(
    "/batch/{job_id}/analytics/{analysis_type}",
    response_model=AnalyticsTriggerResponse,
)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def trigger_batch_analytics(
    request: Request,
    job_id: str,
    analysis_type: str,
    params: Optional[dict[str, Any]] = Body(default=None),
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Trigger an expensive analytics computation for a completed batch job.

    Supported analysis types:
    - **scaffold**: Murcko scaffold diversity analysis.
    - **chemical_space**: PCA or t-SNE 2-D chemical space embedding.
    - **mmp**: Matched molecular pair analysis (requires activity_column param).
    - **similarity_search**: Nearest-neighbor similarity search.
    - **rgroup**: R-group decomposition (requires core_smarts param).

    Optional params body (JSON object):
    - `method`: "pca" or "tsne" (for chemical_space)
    - `activity_column`: property key for MMP cliff detection
    - `query_smiles`: query SMILES for similarity_search
    - `query_index`: query molecule index for similarity_search
    - `top_k`: number of neighbors to return (similarity_search)
    - `core_smarts`: SMARTS pattern for R-group core
    """
    if analysis_type not in _ALLOWED_EXPENSIVE_ANALYSES:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown analysis type '{analysis_type}'. "
            f"Allowed: {sorted(_ALLOWED_EXPENSIVE_ANALYSES)}",
        )

    run_expensive_analytics.delay(job_id, analysis_type, params)

    return AnalyticsTriggerResponse(
        job_id=job_id,
        analysis_type=analysis_type,
        status="queued",
    )


@router.post("/batch/detect-columns", response_model=CSVColumnsResponse)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def detect_columns(
    request: Request,
    file: UploadFile = File(..., description="CSV file to analyze"),
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Detect columns in a delimited text file for SMILES and Name/ID selection.

    Returns list of column names, suggested columns, and sample values.
    Accepts .csv, .tsv, and .txt files.
    """
    filename = file.filename or ""
    text_formats = (".csv", ".tsv", ".txt")
    if not filename.lower().endswith(text_formats):
        raise HTTPException(
            status_code=400,
            detail="This endpoint only accepts delimited text files (.csv, .tsv, .txt)",
        )

    try:
        content = await file.read()

        # Security: Check file size before processing
        file_size_mb = len(content) / (1024 * 1024)
        if file_size_mb > settings.MAX_FILE_SIZE_MB:
            raise HTTPException(
                status_code=400,
                detail=f"File too large: {file_size_mb:.1f}MB exceeds limit",
            )

        result = detect_csv_columns(content, max_file_size_mb=settings.MAX_FILE_SIZE_MB)
        return CSVColumnsResponse(
            columns=result["columns"],
            suggested_smiles=result.get("suggested_smiles"),
            suggested_name=result.get("suggested_name"),
            column_samples=result.get("column_samples", {}),
            row_count_estimate=result.get("row_count_estimate", 0),
            file_size_mb=result.get("file_size_mb", 0),
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(
            status_code=400,
            detail=f"Failed to analyze CSV: {str(e)}",
        )


# =============================================================================
# Batch Subset Actions
# =============================================================================


class SubsetRequest(BaseModel):
    """Request body for subset actions."""

    indices: List[int]


class SubsetRescoreRequest(BaseModel):
    """Request body for subset re-score action."""

    indices: List[int]
    profile_id: Optional[int] = None


class SubsetExportRequest(BaseModel):
    """Request body for subset export action."""

    indices: List[int]
    format: str = "csv"


@router.post("/batch/{job_id}/subset/revalidate")
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def subset_revalidate(
    request: Request,
    job_id: str,
    body: SubsetRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Re-validate a subset of molecules from an existing batch job.

    Creates a new batch job with the selected molecules.
    Returns the new job_id for tracking progress.
    """
    try:
        new_job_id = revalidate_subset(job_id, body.indices)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))

    return {
        "job_id": new_job_id,
        "source_job_id": job_id,
        "molecule_count": len(body.indices),
        "action": "revalidate",
    }


@router.post("/batch/{job_id}/subset/rescore")
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def subset_rescore(
    request: Request,
    job_id: str,
    body: SubsetRescoreRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Re-score a subset of molecules, optionally with a custom scoring profile.

    Creates a new batch job for the selected molecules.
    Returns the new job_id for tracking progress.
    """
    try:
        new_job_id = rescore_subset(job_id, body.indices, body.profile_id)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))

    return {
        "job_id": new_job_id,
        "source_job_id": job_id,
        "molecule_count": len(body.indices),
        "action": "rescore",
        "profile_id": body.profile_id,
    }


@router.post("/batch/{job_id}/subset/export")
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def subset_export(
    request: Request,
    job_id: str,
    body: SubsetExportRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Export a subset of molecules in the specified format.

    Returns a file download with the exported data.
    """
    from datetime import datetime

    from fastapi.responses import StreamingResponse

    try:
        export_buffer = export_subset(job_id, body.indices, body.format)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Export failed: {str(e)}")

    # Determine media type and extension
    from app.services.export.base import ExporterFactory, ExportFormat

    export_format = ExportFormat(body.format)
    exporter = ExporterFactory.create(export_format)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"subset_{job_id[:8]}_{timestamp}.{exporter.file_extension}"

    def iterfile():
        export_buffer.seek(0)
        while chunk := export_buffer.read(1024 * 1024):
            yield chunk

    return StreamingResponse(
        iterfile(),
        media_type=exporter.media_type,
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )
