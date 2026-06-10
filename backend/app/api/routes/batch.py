"""
Batch Processing API Routes

Endpoints for batch file upload, job status, results retrieval, and analytics.
"""

import logging
import re
import uuid
from typing import Any, Optional

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

from app.core.config import settings
from app.core.error_sanitizer import safe_error_detail
from app.core.ownership import store_job_owner, verify_job_access
from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.core.session import get_session_id
from app.schemas.analytics import (
    AnalysisStatus,
    AnalyticsTriggerResponse,
    BatchAnalyticsResponse,
    MCSCompareRequest,
    MCSComparisonResult,
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
    detect_suspicious_content,
    parse_csv,
    parse_sdf,
    validate_file_content_type,
)
from app.services.batch.parsing import parse_json_field
from app.services.batch.progress_tracker import progress_tracker
from app.services.batch.result_aggregator import result_storage
from app.services.batch.tasks import cancel_batch_job, process_batch_job

logger = logging.getLogger(__name__)

router = APIRouter()


@router.post("/batch/upload", response_model=BatchUploadResponse)
@limiter.limit("3/minute", key_func=get_rate_limit_key)
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
    profile_id: Optional[int] = Form(
        default=None,
        description="Scoring profile ID for profile-based desirability scoring",
    ),
    include_profiling: bool = Form(
        default=False,
        description="Compute compound profiling (PFI, stars, bioavailability, etc.) for each molecule",
    ),
    include_safety_assessment: bool = Form(
        default=False,
        description="Run safety assessment (CYP, hERG, bRo5, REOS, complexity) for each molecule",
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

    # Security: Audit-log ALL suspicious patterns before strict validation.
    # validate_file_content_type below will block on the first pattern match
    # with a 400 response, but the audit log should capture every finding so
    # post-hoc investigation has the full picture for rejected uploads too.
    suspicious_findings = detect_suspicious_content(content)
    if suspicious_findings:
        client_ip = request.client.host if request.client else "unknown"
        logger.warning(
            "Suspicious patterns detected in batch upload: filename=%s "
            "client_ip=%s size_bytes=%d findings=%s",
            filename,
            client_ip,
            len(content),
            suspicious_findings,
        )

    # Security: Validate file content matches extension and block on malicious
    # patterns. .csv, .tsv, and .txt are all treated as delimited text files.
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
    except ValueError:
        raise HTTPException(
            status_code=400,
            detail="Invalid file content. Please check the file format and column names.",
        )
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

    # Store ownership for session-based access control
    store_job_owner(job_id, get_session_id(request))

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
        "include_profiling": include_profiling,
        "include_safety_assessment": include_safety_assessment,
    }

    # If profile selected, fetch thresholds and attach to safety_options
    if profile_id is not None:
        from app.db import async_session
        from app.services.profiles.service import ProfileService

        async with async_session() as db:
            profile = await ProfileService().get(db, profile_id)
        if profile is None:
            raise HTTPException(status_code=404, detail=f"Profile {profile_id} not found")
        safety_options["profile_id"] = profile.id
        safety_options["profile_name"] = profile.name
        safety_options["profile_thresholds"] = parse_json_field(profile.thresholds)
        safety_options["profile_weights"] = parse_json_field(profile.weights)

    # Validate notification email before dispatching job
    email_target = notification_email or settings.NOTIFICATION_EMAIL
    if email_target:
        email_re = re.compile(r"^[a-zA-Z0-9._%+\-]+@[a-zA-Z0-9.\-]+\.[a-zA-Z]{2,}$")
        if len(email_target) > 254 or not email_re.match(email_target):
            raise HTTPException(
                status_code=400,
                detail="Invalid notification email address",
            )

    # Start batch processing
    process_batch_job(job_id, mol_dicts, safety_options=safety_options)

    # Store notification email for this job (if provided or globally configured)
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
    issue_filter: Optional[str] = Query(
        default=None,
        max_length=200,
        description="Filter to molecules with a specific failed validation check name",
    ),
    alert_filter: Optional[str] = Query(
        default=None,
        max_length=200,
        description="Filter to molecules with alerts from a specific catalog",
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
    - **issue_filter**: Filter to molecules with a specific failed validation check

    Returns job status and paginated results with statistics.
    """
    verify_job_access(request, job_id)

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
        issue_filter=issue_filter,
        alert_filter=alert_filter,
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
            profiling=r.get("profiling"),
            safety_assessment=r.get("safety_assessment"),
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
    verify_job_access(request, job_id)

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
    verify_job_access(request, job_id)

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
    verify_job_access(request, job_id)

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
    verify_job_access(request, job_id)

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
    result_field_map = {
        "deduplication": "deduplication",
        "scaffold": "scaffold",
        "chemical_space": "chemical_space",
        "similarity_search": "similarity_matrix",
        "mmp": "mmp",
        "statistics": "statistics",
        "clustering": "clustering",
        "taxonomy": "taxonomy",
        "registration": "registration",
    }

    result_kwargs: dict[str, Any] = {}
    for analysis_type, field_name in result_field_map.items():
        entry = status_dict.get(analysis_type)
        if entry and entry.status == "complete":
            raw = analytics_storage.get_result(job_id, analysis_type)
            if raw is not None:
                result_kwargs[field_name] = raw

    # Backfill smiles_map for clustering results that predate the field
    clustering_data = result_kwargs.get("clustering")
    if clustering_data and not clustering_data.get("smiles_map"):
        all_results = result_storage.get_all_results(job_id)
        if all_results:
            smap: dict[str, str] = {}
            for r in all_results:
                idx = r.get("index", 0)
                std = r.get("standardization")
                std_smi = std.get("standardized_smiles") if isinstance(std, dict) else None
                smi = std_smi or r.get("smiles", "")
                if smi:
                    smap[str(idx)] = smi
            clustering_data["smiles_map"] = smap

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
    "clustering",
    "taxonomy",
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
    - **clustering**: Butina sphere-exclusion clustering (optional cutoff param).
    - **taxonomy**: SMARTS-based chemotype classification.

    Optional params body (JSON object):
    - `method`: "pca" or "tsne" (for chemical_space)
    - `activity_column`: property key for MMP cliff detection
    - `query_smiles`: query SMILES for similarity_search
    - `query_index`: query molecule index for similarity_search
    - `top_k`: number of neighbors to return (similarity_search)
    - `core_smarts`: SMARTS pattern for R-group core
    - `cutoff`: Tanimoto distance threshold for clustering (default 0.35)
    """
    verify_job_access(request, job_id)

    if analysis_type not in _ALLOWED_EXPENSIVE_ANALYSES:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown analysis type '{analysis_type}'. "
            f"Allowed: {sorted(_ALLOWED_EXPENSIVE_ANALYSES)}",
        )

    # --- Idempotency check: avoid unnecessary re-computation ----------------
    # Without this, every trigger call unconditionally sets "computing" and
    # dispatches a Celery task — even when the result is already cached.
    # This defeats the Celery task's own idempotency guard (which checks for
    # "complete" status that the route already overwrote to "computing").
    status_data = analytics_storage.get_status(job_id)
    if status_data:
        current_entry = status_data.get(analysis_type, {})
        current_status_val = (
            current_entry.get("status") if isinstance(current_entry, dict) else None
        )

        if current_status_val == "computing":
            # Already computing — don't dispatch a duplicate task.
            return AnalyticsTriggerResponse(
                job_id=job_id,
                analysis_type=analysis_type,
                status="already_computing",
            )

        if current_status_val == "complete":
            # Allow re-computation when params differ from stored result.
            needs_recompute = False
            if analysis_type == "chemical_space" and params and params.get("method"):
                stored = analytics_storage.get_result(job_id, "chemical_space")
                if stored and stored.get("method") != params["method"]:
                    needs_recompute = True
            elif analysis_type == "clustering" and params and params.get("cutoff"):
                stored = analytics_storage.get_result(job_id, "clustering")
                stored_cutoff = stored.get("distance_cutoff") if stored else None
                new_cutoff = float(params["cutoff"])
                if stored_cutoff is None or abs(stored_cutoff - new_cutoff) > 1e-9:
                    needs_recompute = True

            if not needs_recompute:
                return AnalyticsTriggerResponse(
                    job_id=job_id,
                    analysis_type=analysis_type,
                    status="already_complete",
                )

    # --- Clustering cap check (D-06: hard cap at 1,000 molecules) ---
    if analysis_type == "clustering":
        results = result_storage.get_all_results(job_id)
        if results and len(results) > 1000:
            raise HTTPException(
                status_code=400,
                detail=f"Clustering is limited to 1,000 molecules. "
                f"This batch has {len(results)} molecules. "
                "Filter or subsample before clustering.",
            )

    # Set status to "computing" immediately so frontend polling sees the
    # transition before the Celery worker picks up the task.
    analytics_storage.update_status(job_id, analysis_type, "computing")

    run_expensive_analytics.delay(job_id, analysis_type, params)

    return AnalyticsTriggerResponse(
        job_id=job_id,
        analysis_type=analysis_type,
        status="queued",
    )


@router.post("/batch/{job_id}/mcs", response_model=MCSComparisonResult)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def compute_batch_mcs(
    request: Request,
    job_id: str,
    body: MCSCompareRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Compute Maximum Common Substructure between two molecules from a batch.

    Runs synchronously (timeout=10s max). Returns MCS SMARTS, matched
    atom/bond counts, Tanimoto similarity, and property deltas.

    - **job_id**: Batch job identifier.
    - **body.index_a**: Index of the first molecule.
    - **body.index_b**: Index of the second molecule.
    """
    verify_job_access(request, job_id)

    results = result_storage.get_all_results(job_id)
    if not results:
        raise HTTPException(
            status_code=404,
            detail=f"Batch results not found for job {job_id}. Results may have expired.",
        )

    # Find molecules by index
    mol_a = next((r for r in results if r.get("index") == body.index_a), None)
    mol_b = next((r for r in results if r.get("index") == body.index_b), None)

    if mol_a is None or mol_b is None:
        missing = []
        if mol_a is None:
            missing.append(str(body.index_a))
        if mol_b is None:
            missing.append(str(body.index_b))
        raise HTTPException(
            status_code=404,
            detail=f"Molecule(s) at index {', '.join(missing)} not found in batch.",
        )

    smiles_a = mol_a.get("smiles", "")
    smiles_b = mol_b.get("smiles", "")

    if not smiles_a or not smiles_b:
        raise HTTPException(
            status_code=400,
            detail="Both molecules must have valid SMILES.",
        )

    try:
        from app.services.analytics.mcs_comparison import compute_mcs_comparison

        result = compute_mcs_comparison(smiles_a, smiles_b)
        return MCSComparisonResult(**result)
    except ValueError:
        raise HTTPException(
            status_code=400, detail="Invalid input for MCS comparison"
        )
    except Exception as exc:
        raise HTTPException(
            status_code=500,
            detail=safe_error_detail(exc, "MCS computation failed"),
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
    except ValueError:
        raise HTTPException(
            status_code=400, detail="Invalid CSV file or column structure"
        )
    except Exception as e:
        raise HTTPException(
            status_code=400,
            detail=safe_error_detail(e, "Failed to analyze CSV file"),
        )


