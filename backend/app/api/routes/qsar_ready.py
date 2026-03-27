"""
QSAR-Ready Pipeline API Routes

Endpoints for single-molecule and batch QSAR-Ready processing.

Routes:
  POST /api/v1/qsar-ready/single                           — run single molecule through pipeline
  POST /api/v1/qsar-ready/batch/upload                     — upload file/paste SMILES for batch
  GET  /api/v1/qsar-ready/batch/{job_id}/status            — poll batch job progress
  GET  /api/v1/qsar-ready/batch/{job_id}/results           — paginated batch results
  GET  /api/v1/qsar-ready/batch/{job_id}/download/{format} — download CSV, SDF, or JSON

WebSocket endpoint /ws/qsar/{job_id} is registered in main.py.
"""

import csv
import dataclasses
import importlib
import io
import json
import logging
import uuid
from typing import Optional

from fastapi import (
    APIRouter,
    Depends,
    File,
    Form,
    HTTPException,
    Query,
    Request,
    UploadFile,
)
from fastapi.responses import StreamingResponse

from app.core.config import settings
from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.schemas.qsar_ready import (
    QSARBatchResultsResponse,
    QSARBatchStatusResponse,
    QSARBatchSummary,
    QSARBatchUploadResponse,
    QSARReadyConfigSchema,
    QSARReadyResultSchema,
    QSARSingleRequest,
    QSARStepResultSchema,
)
from app.services.batch.progress_tracker import progress_tracker
from app.services.qsar_ready.pipeline import QSARReadyConfig, qsar_ready_single

logger = logging.getLogger(__name__)

router = APIRouter()


# =============================================================================
# Redis helpers
# =============================================================================


def _get_redis():
    """Get a synchronous Redis client."""
    import redis as sync_redis

    return sync_redis.from_url(settings.REDIS_URL, decode_responses=True)


# =============================================================================
# Endpoint 1: Single molecule
# =============================================================================


@router.post("/qsar-ready/single", response_model=QSARReadyResultSchema)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def qsar_ready_single_endpoint(
    request: Request,
    body: QSARSingleRequest,
    api_key: Optional[str] = Depends(get_api_key),
) -> QSARReadyResultSchema:
    """
    Run a single SMILES string through the configurable QSAR-Ready 10-step pipeline.

    Accepts a SMILES string and a QSARReadyConfig and returns the full pipeline
    result with per-step provenance, InChIKey change detection, and status.

    Args:
        body: QSARSingleRequest with smiles and config fields.

    Returns:
        QSARReadyResultSchema with status, curated_smiles, steps, and InChIKey fields.

    Raises:
        HTTPException: 500 on unexpected server error.
    """
    try:
        config = QSARReadyConfig(**body.config.model_dump())
        result = qsar_ready_single(body.smiles, config)
        return QSARReadyResultSchema(**dataclasses.asdict(result))
    except Exception as exc:
        logger.exception("Error in qsar_ready_single: %s", exc)
        raise HTTPException(
            status_code=500,
            detail={"error": "Pipeline processing failed", "detail": str(exc)},
        ) from exc


# =============================================================================
# Endpoint 2: Batch upload
# =============================================================================


@router.post("/qsar-ready/batch/upload", response_model=QSARBatchUploadResponse)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def qsar_batch_upload(
    request: Request,
    config: str = Form(..., description="JSON-encoded QSARReadyConfig"),
    file: Optional[UploadFile] = File(default=None, description="CSV or SDF file"),
    smiles_text: Optional[str] = Form(
        default=None,
        description="Pasted SMILES, one per line (D-03). Used when no file is uploaded.",
    ),
    api_key: Optional[str] = Depends(get_api_key),
) -> QSARBatchUploadResponse:
    """
    Upload a CSV/SDF file or paste SMILES for batch QSAR-Ready processing.

    Accepts either a file upload or SMILES text input (D-03).
    Auto-runs Phase 9 file pre-validator before processing (D-13, Pitfall 6).
    Dispatches a Celery task and returns the job_id immediately.

    Args:
        config: JSON-encoded QSARReadyConfigSchema.
        file: Optional SDF or CSV file upload.
        smiles_text: Optional SMILES pasted as text, one per line.

    Returns:
        QSARBatchUploadResponse with job_id, total_molecules, status, message.

    Raises:
        HTTPException: 400 for invalid input, 422 for critical pre-validation issues.
    """
    # Parse and validate config
    try:
        config_data = json.loads(config)
        config_schema = QSARReadyConfigSchema(**config_data)
    except (json.JSONDecodeError, Exception) as exc:
        raise HTTPException(
            status_code=400,
            detail=f"Invalid config JSON: {exc}",
        ) from exc

    smiles_list: list[str] = []

    if file is not None:
        # File upload path — read content first
        try:
            content = await file.read()
        except Exception as exc:
            raise HTTPException(status_code=400, detail="Failed to read file") from exc

        # File size check
        file_size_mb = len(content) / (1024 * 1024)
        if file_size_mb > settings.MAX_FILE_SIZE_MB:
            raise HTTPException(
                status_code=400,
                detail=(
                    f"File too large: {file_size_mb:.1f}MB exceeds limit of "
                    f"{settings.MAX_FILE_SIZE_MB}MB"
                ),
            )

        filename = file.filename or ""

        # Phase 9 file pre-validation (D-13, Pitfall 6)
        # importlib bypasses app.services.diagnostics.__init__ which might trigger
        # heavy imports. If the diagnostics service is not yet deployed (Phase 9
        # not merged), degrade gracefully.
        try:
            prevalidator = importlib.import_module("app.services.diagnostics.file_prevalidator")
            filename_lower = filename.lower()
            if filename_lower.endswith(".sdf") or filename_lower.endswith(".sd"):
                prevalidation_result = prevalidator.prevalidate_sdf(content)
            else:
                prevalidation_result = prevalidator.prevalidate_csv(content)
            # If critical issues found — return 422 with pre-validation results
            if prevalidation_result.get("critical_issues_found", False):
                raise HTTPException(
                    status_code=422,
                    detail={
                        "error": "File has critical pre-validation issues",
                        "prevalidation": prevalidation_result,
                    },
                )
        except HTTPException:
            raise
        except Exception:
            # Phase 9 diagnostics not available — skip pre-validation gracefully
            logger.debug("File pre-validation skipped: diagnostics module not available")

        # Parse SMILES from file using existing file_parser
        try:
            from app.services.batch.file_parser import parse_csv, parse_sdf

            filename_lower = filename.lower()
            if filename_lower.endswith(".sdf") or filename_lower.endswith(".sd"):
                molecules = parse_sdf(content, max_file_size_mb=settings.MAX_FILE_SIZE_MB)
            else:
                molecules = parse_csv(
                    content,
                    smiles_column="SMILES",
                    name_column=None,
                    max_file_size_mb=settings.MAX_FILE_SIZE_MB,
                )
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        except Exception as exc:
            raise HTTPException(
                status_code=400,
                detail="Failed to parse file. Please check the file format.",
            ) from exc

        smiles_list = [
            m.smiles for m in molecules if m.smiles and not m.parse_error
        ]

    elif smiles_text:
        # Text paste path (D-03)
        smiles_list = [
            line.strip()
            for line in smiles_text.splitlines()
            if line.strip()
        ]
    else:
        raise HTTPException(
            status_code=400,
            detail="Either a file or smiles_text must be provided",
        )

    if not smiles_list:
        raise HTTPException(
            status_code=400,
            detail="No valid SMILES found in input",
        )

    if len(smiles_list) > settings.MAX_BATCH_SIZE:
        raise HTTPException(
            status_code=400,
            detail=(
                f"Input contains {len(smiles_list)} molecules. "
                f"Maximum allowed is {settings.MAX_BATCH_SIZE}."
            ),
        )

    # Generate job ID
    job_id = str(uuid.uuid4())

    # Serialize config for Celery (must be JSON-serializable)
    config_dict = config_schema.model_dump()

    # Store job metadata in Redis for status/results retrieval
    try:
        r = _get_redis()
        r.set(
            f"qsar:job:{job_id}",
            json.dumps({"config": config_dict, "total": len(smiles_list)}),
            ex=3600,
        )
        # Store SMILES list for results retrieval after completion
        r.set(
            f"qsar:smiles:{job_id}",
            json.dumps(smiles_list),
            ex=3600,
        )
        # Store session ownership for WebSocket access control
        from app.core.session import get_session_id

        session_id = get_session_id(request)
        if session_id:
            r.set(f"qsar:owner:{job_id}", session_id, ex=settings.BATCH_RESULT_TTL)
    except Exception as exc:
        logger.warning("Failed to store QSAR job metadata in Redis: %s", exc)

    # Dispatch Celery task
    try:
        from app.services.batch.qsar_tasks import process_qsar_batch_job

        process_qsar_batch_job.delay(job_id, smiles_list, config_dict)
    except Exception as exc:
        logger.error("Failed to dispatch QSAR batch task for %s: %s", job_id, exc)
        raise HTTPException(
            status_code=500,
            detail="Failed to start batch processing",
        ) from exc

    return QSARBatchUploadResponse(
        job_id=job_id,
        total_molecules=len(smiles_list),
        status="pending",
        message=f"Batch job queued. Processing {len(smiles_list)} molecules.",
    )


# =============================================================================
# Endpoint 3: Batch job status
# =============================================================================


@router.get("/qsar-ready/batch/{job_id}/status", response_model=QSARBatchStatusResponse)
@limiter.limit("60/minute", key_func=get_rate_limit_key)
async def qsar_batch_status(
    request: Request,
    job_id: str,
    api_key: Optional[str] = Depends(get_api_key),
) -> QSARBatchStatusResponse:
    """
    Get the current processing status of a QSAR batch job.

    Reads progress from the Redis key `batch:job:{job_id}` maintained by
    progress_tracker. Returns progress percentage, processed count, and total.

    Args:
        job_id: Batch job identifier (UUID).

    Returns:
        QSARBatchStatusResponse with status, progress, and counts.

    Raises:
        HTTPException: 404 if job not found, 400 if job_id format is invalid.
    """
    _validate_uuid(job_id)

    progress = progress_tracker.get_progress(job_id)
    if progress is None:
        raise HTTPException(
            status_code=404,
            detail=f"Job {job_id} not found",
        )

    return QSARBatchStatusResponse(
        job_id=job_id,
        status=progress.status,
        progress=progress.progress,
        processed=progress.processed,
        total=progress.total,
        eta_seconds=progress.eta_seconds,
    )


# =============================================================================
# Endpoint 4: Batch results
# =============================================================================


@router.get("/qsar-ready/batch/{job_id}/results", response_model=QSARBatchResultsResponse)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def qsar_batch_results(
    request: Request,
    job_id: str,
    page: int = Query(default=1, ge=1),
    per_page: int = Query(default=50, ge=1, le=500),
    api_key: Optional[str] = Depends(get_api_key),
) -> QSARBatchResultsResponse:
    """
    Retrieve paginated results for a completed QSAR batch job.

    Reads results stored in Redis by the Celery task. Returns paginated
    QSARReadyResult records with summary statistics.

    Args:
        job_id: Batch job identifier (UUID).
        page: Page number (1-based).
        per_page: Results per page (1-500).

    Returns:
        QSARBatchResultsResponse with paginated results and summary.

    Raises:
        HTTPException: 404 if job not found or results not yet available.
    """
    _validate_uuid(job_id)

    # Check job progress to get status and config
    progress = progress_tracker.get_progress(job_id)
    if progress is None:
        raise HTTPException(
            status_code=404,
            detail=f"Job {job_id} not found",
        )

    # Read results from Redis
    try:
        r = _get_redis()
        raw_results = r.get(f"qsar:results:{job_id}")
        raw_job = r.get(f"qsar:job:{job_id}")
    except Exception as exc:
        raise HTTPException(
            status_code=503,
            detail="Redis unavailable",
        ) from exc

    if raw_results is None:
        raise HTTPException(
            status_code=404,
            detail=f"Results for job {job_id} not yet available (status: {progress.status})",
        )

    all_results = json.loads(raw_results)
    config_dict = {}
    if raw_job:
        job_meta = json.loads(raw_job)
        config_dict = job_meta.get("config", {})

    # Compute summary from full results
    summary = _compute_summary(all_results)

    # Paginate
    total_results = len(all_results)
    total_pages = max(1, (total_results + per_page - 1) // per_page)
    start = (page - 1) * per_page
    end = start + per_page
    page_results = all_results[start:end]

    # Build result schemas
    result_schemas = [QSARReadyResultSchema(**_normalize_result(r)) for r in page_results]

    # Build config schema (use defaults if not stored)
    try:
        config_schema = QSARReadyConfigSchema(**config_dict) if config_dict else QSARReadyConfigSchema()
    except Exception:
        config_schema = QSARReadyConfigSchema()

    return QSARBatchResultsResponse(
        job_id=job_id,
        status=progress.status,
        config=config_schema,
        summary=summary,
        results=result_schemas,
        page=page,
        per_page=per_page,
        total_pages=total_pages,
        total_results=total_results,
    )


# =============================================================================
# Endpoint 5: Batch download
# =============================================================================


@router.get("/qsar-ready/batch/{job_id}/download/{format}")
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def qsar_batch_download(
    request: Request,
    job_id: str,
    format: str,
    api_key: Optional[str] = Depends(get_api_key),
) -> StreamingResponse:
    """
    Download batch QSAR results in CSV, SDF, or JSON format (D-12).

    CSV columns: original_smiles, curated_smiles, original_inchikey,
                 curated_inchikey, inchikey_changed, status, steps_applied.
    SDF: curated mol blocks from curated_smiles, title = original_inchikey.
    JSON: full provenance dump with all results, summary, config, and duplicate list.

    Args:
        job_id: Batch job identifier (UUID).
        format: Output format — "csv", "sdf", or "json".

    Returns:
        StreamingResponse with appropriate Content-Type and Content-Disposition headers.

    Raises:
        HTTPException: 400 for invalid format, 404 for missing job/results.
    """
    _validate_uuid(job_id)

    if format not in ("csv", "sdf", "json"):
        raise HTTPException(
            status_code=400,
            detail=f"Invalid format '{format}'. Supported: csv, sdf, json",
        )

    # Read results from Redis
    try:
        r = _get_redis()
        raw_results = r.get(f"qsar:results:{job_id}")
        raw_job = r.get(f"qsar:job:{job_id}")
    except Exception as exc:
        raise HTTPException(status_code=503, detail="Redis unavailable") from exc

    if raw_results is None:
        raise HTTPException(
            status_code=404,
            detail=f"Results for job {job_id} not found",
        )

    all_results = json.loads(raw_results)
    config_dict = {}
    if raw_job:
        job_meta = json.loads(raw_job)
        config_dict = job_meta.get("config", {})

    if format == "csv":
        return _build_csv_response(job_id, all_results)
    elif format == "sdf":
        return _build_sdf_response(job_id, all_results)
    else:
        return _build_json_response(job_id, all_results, config_dict)


# =============================================================================
# Download helpers
# =============================================================================


def _build_csv_response(job_id: str, results: list) -> StreamingResponse:
    """Build a CSV StreamingResponse from results (D-12)."""
    output = io.StringIO()
    writer = csv.DictWriter(
        output,
        fieldnames=[
            "original_smiles",
            "curated_smiles",
            "original_inchikey",
            "curated_inchikey",
            "inchikey_changed",
            "status",
            "rejection_reason",
            "steps_applied",
        ],
    )
    writer.writeheader()
    for r in results:
        steps = r.get("steps", [])
        applied_steps = ",".join(
            s["step_name"] for s in steps if s.get("status") == "applied"
        )
        writer.writerow(
            {
                "original_smiles": r.get("original_smiles", ""),
                "curated_smiles": r.get("curated_smiles", ""),
                "original_inchikey": r.get("original_inchikey", ""),
                "curated_inchikey": r.get("standardized_inchikey", ""),
                "inchikey_changed": str(r.get("inchikey_changed", False)).lower(),
                "status": r.get("status", ""),
                "rejection_reason": r.get("rejection_reason", ""),
                "steps_applied": applied_steps,
            }
        )

    output.seek(0)
    return StreamingResponse(
        iter([output.getvalue()]),
        media_type="text/csv",
        headers={
            "Content-Disposition": f"attachment; filename=qsar_results_{job_id[:8]}.csv"
        },
    )


def _build_sdf_response(job_id: str, results: list) -> StreamingResponse:
    """Build an SDF StreamingResponse from curated_smiles (Pitfall 8 from RESEARCH.md)."""
    from rdkit import Chem

    sdf_blocks = []
    for r in results:
        curated_smiles = r.get("curated_smiles")
        if not curated_smiles:
            continue
        mol = Chem.MolFromSmiles(curated_smiles)
        if mol is None:
            continue
        # Use original InChIKey as molecule title
        title = r.get("original_inchikey") or r.get("original_smiles", "")[:60]
        mol.SetProp("_Name", title)
        mol.SetProp("original_smiles", r.get("original_smiles", ""))
        mol.SetProp("curated_smiles", curated_smiles)
        mol.SetProp("status", r.get("status", ""))
        try:
            block = Chem.MolToMolBlock(mol)
            sdf_blocks.append(block)
            sdf_blocks.append("$$$$\n")
        except Exception:
            continue

    sdf_content = "\n".join(sdf_blocks) if sdf_blocks else ""

    return StreamingResponse(
        iter([sdf_content]),
        media_type="chemical/x-mdl-sdfile",
        headers={
            "Content-Disposition": f"attachment; filename=qsar_results_{job_id[:8]}.sdf"
        },
    )


def _build_json_response(job_id: str, results: list, config_dict: dict) -> StreamingResponse:
    """Build a full-provenance JSON StreamingResponse (D-12)."""
    summary = _compute_summary_dict(results)
    duplicates = [
        {
            "original_smiles": r.get("original_smiles"),
            "standardized_inchikey": r.get("standardized_inchikey"),
            "rejection_reason": r.get("rejection_reason"),
        }
        for r in results
        if r.get("status") == "duplicate"
    ]

    payload = {
        "job_id": job_id,
        "summary": summary,
        "config": config_dict,
        "duplicates": duplicates,
        "results": results,
    }
    json_str = json.dumps(payload, indent=2, default=str)
    return StreamingResponse(
        iter([json_str]),
        media_type="application/json",
        headers={
            "Content-Disposition": f"attachment; filename=qsar_results_{job_id[:8]}.json"
        },
    )


# =============================================================================
# Utility helpers
# =============================================================================


def _validate_uuid(job_id: str) -> None:
    """Validate UUID format; raises HTTPException 400 if invalid."""
    import re

    _UUID_RE = re.compile(
        r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$", re.I
    )
    if not _UUID_RE.match(job_id):
        raise HTTPException(
            status_code=400,
            detail="Invalid job ID format — must be a UUID",
        )


def _normalize_result(r: dict) -> dict:
    """Normalize a result dict to match QSARReadyResultSchema field names."""
    # Ensure steps have correct structure
    steps = r.get("steps", [])
    normalized_steps = []
    for s in steps:
        normalized_steps.append(
            {
                "step_name": s.get("step_name", ""),
                "step_index": s.get("step_index", 0),
                "enabled": s.get("enabled", True),
                "status": s.get("status", "skipped"),
                "before_smiles": s.get("before_smiles"),
                "after_smiles": s.get("after_smiles"),
                "detail": s.get("detail"),
            }
        )
    return {
        "original_smiles": r.get("original_smiles", ""),
        "original_inchikey": r.get("original_inchikey"),
        "curated_smiles": r.get("curated_smiles"),
        "standardized_inchikey": r.get("standardized_inchikey"),
        "inchikey_changed": r.get("inchikey_changed", False),
        "status": r.get("status", "error"),
        "rejection_reason": r.get("rejection_reason"),
        "steps": normalized_steps,
    }


def _compute_summary(results: list) -> QSARBatchSummary:
    """Compute summary statistics from a list of result dicts."""
    ok = sum(1 for r in results if r.get("status") == "ok")
    rejected = sum(1 for r in results if r.get("status") == "rejected")
    duplicate = sum(1 for r in results if r.get("status") == "duplicate")
    error = sum(1 for r in results if r.get("status") == "error")

    steps_applied_counts: dict = {}
    for r in results:
        for step in r.get("steps", []):
            if step.get("status") == "applied":
                name = step.get("step_name", "")
                steps_applied_counts[name] = steps_applied_counts.get(name, 0) + 1

    return QSARBatchSummary(
        total=len(results),
        ok=ok,
        rejected=rejected,
        duplicate=duplicate,
        error=error,
        steps_applied_counts=steps_applied_counts,
    )


def _compute_summary_dict(results: list) -> dict:
    """Compute summary statistics as a plain dict."""
    summary = _compute_summary(results)
    return summary.model_dump()
