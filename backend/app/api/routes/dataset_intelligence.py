"""
Dataset Intelligence API Routes (Phase 12).

Endpoints for dataset health auditing, contradictory label detection,
dataset diff comparison, and curation report/CSV downloads.

Routes:
  POST /api/v1/dataset/upload                  -- upload CSV/SDF for async audit
  GET  /api/v1/dataset/{job_id}/status         -- poll job progress
  GET  /api/v1/dataset/{job_id}/results        -- fetch audit results
  POST /api/v1/dataset/{job_id}/diff           -- compare with a second dataset
  GET  /api/v1/dataset/{job_id}/download/report -- download curation report JSON
  GET  /api/v1/dataset/{job_id}/download/csv   -- download curated CSV

WebSocket /ws/dataset/{job_id} is registered in main.py.
"""

import base64
import csv
import io
import json
import logging
import os
import re
import tempfile
from typing import Optional
from uuid import uuid4

from fastapi import APIRouter, Depends, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import Response

from app.core.config import settings
from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.core.session import get_session_id
from app.schemas.dataset_intelligence import (
    DatasetAuditResponse,
    DatasetAuditStatusResponse,
    DatasetDiffResponse,
    DatasetUploadResponse,
    DiffMolecule,
)

logger = logging.getLogger(__name__)

router = APIRouter(tags=["dataset-intelligence"])

_UUID_RE = re.compile(
    r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$", re.I
)

_ALLOWED_EXTENSIONS = {"csv", "sdf"}


def _validate_uuid(job_id: str) -> None:
    """Validate UUID format; raises HTTPException 400 if invalid."""
    if not _UUID_RE.match(job_id):
        raise HTTPException(
            status_code=400,
            detail="Invalid job ID format -- must be a UUID",
        )


def _get_redis():
    """Get a synchronous Redis client."""
    import redis as sync_redis

    return sync_redis.from_url(settings.REDIS_URL, decode_responses=True)


def _get_file_extension(filename: str) -> str:
    """Extract and validate file extension from filename."""
    if not filename:
        return ""
    parts = filename.rsplit(".", 1)
    return parts[-1].lower() if len(parts) > 1 else ""


# =============================================================================
# Endpoint 1: Upload dataset
# =============================================================================


@router.post("/dataset/upload", response_model=DatasetUploadResponse)
@limiter.limit("3/minute", key_func=get_rate_limit_key)
async def upload_dataset(
    request: Request,
    file: UploadFile = File(..., description="CSV or SDF file with molecules"),
    smiles_column: Optional[str] = Form(
        default=None, description="Name of the SMILES column (auto-detected if None)"
    ),
    activity_column: Optional[str] = Form(
        default=None, description="Name of the activity column (auto-detected if None)"
    ),
    api_key: Optional[str] = Depends(get_api_key),
) -> DatasetUploadResponse:
    """Upload a CSV or SDF file for dataset health auditing.

    The file is saved to a temp directory and a Celery task is dispatched
    for async processing. Returns a job_id for tracking progress via
    WebSocket or polling.

    Args:
        file: CSV or SDF file containing molecules.
        smiles_column: Optional SMILES column name (auto-detected if None).
        activity_column: Optional activity column name (auto-detected if None).

    Returns:
        DatasetUploadResponse with job_id, filename, file_type.

    Raises:
        HTTPException: 400 for invalid file type, 413 for oversized files.
    """
    filename = file.filename or "unknown"
    ext = _get_file_extension(filename)

    # Validate file extension
    if ext not in _ALLOWED_EXTENSIONS:
        raise HTTPException(
            status_code=400,
            detail="Invalid file type. Please upload a CSV or SDF file.",
        )

    # Read file content
    try:
        content = await file.read()
    except Exception as exc:
        raise HTTPException(
            status_code=400, detail="Failed to read uploaded file"
        ) from exc

    # Validate file size
    file_size_mb = len(content) / (1024 * 1024)
    if file_size_mb > settings.MAX_FILE_SIZE_MB:
        raise HTTPException(
            status_code=413,
            detail="File exceeds maximum size limit.",
        )

    # Encode file content for Celery serialization (avoids cross-container
    # temp-file issues — Celery workers run in separate Docker containers
    # and cannot access the backend's /tmp filesystem).
    file_content_b64 = base64.b64encode(content).decode("ascii")

    # Generate job ID and initialize Redis metadata
    job_id = str(uuid4())
    options = {
        "smiles_column": smiles_column,
        "activity_column": activity_column,
    }

    try:
        r = _get_redis()
        r.hset(f"dataset:meta:{job_id}", mapping={
            "status": "pending",
            "progress": "0",
            "current_stage": "",
            "filename": filename,
        })
        r.expire(f"dataset:meta:{job_id}", 3600)
        # Store session ownership for WebSocket access control
        session_id = get_session_id(request)
        if session_id:
            r.set(f"dataset:owner:{job_id}", session_id, ex=3600)
    except Exception as exc:
        logger.warning("Failed to store dataset job metadata: %s", exc)

    # Dispatch Celery task
    try:
        from app.services.dataset_intelligence.batch_processor import (
            process_dataset_audit,
        )

        process_dataset_audit.delay(file_content_b64, filename, ext, job_id, options)
    except Exception as exc:
        logger.error("Failed to dispatch dataset audit task for %s: %s", job_id, exc)
        raise HTTPException(
            status_code=500,
            detail="Failed to start dataset processing",
        ) from exc

    return DatasetUploadResponse(
        job_id=job_id,
        filename=filename,
        file_type=ext,
    )


# =============================================================================
# Endpoint 2: Job status
# =============================================================================


@router.get("/dataset/{job_id}/status", response_model=DatasetAuditStatusResponse)
@limiter.limit("60/minute", key_func=get_rate_limit_key)
async def get_dataset_status(
    request: Request,
    job_id: str,
    api_key: Optional[str] = Depends(get_api_key),
) -> DatasetAuditStatusResponse:
    """Get the current processing status of a dataset audit job.

    Reads job metadata from Redis hash ``dataset:meta:{job_id}``.

    Args:
        job_id: Dataset audit job identifier (UUID).

    Returns:
        DatasetAuditStatusResponse with status, progress, current_stage.

    Raises:
        HTTPException: 404 if job not found.
    """
    _validate_uuid(job_id)

    try:
        r = _get_redis()
        meta = r.hgetall(f"dataset:meta:{job_id}")
    except Exception as exc:
        raise HTTPException(status_code=503, detail="Redis unavailable") from exc

    if not meta:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    return DatasetAuditStatusResponse(
        job_id=job_id,
        status=meta.get("status", "unknown"),
        progress=float(meta.get("progress", 0)),
        current_stage=meta.get("current_stage") or None,
    )


# =============================================================================
# Endpoint 3: Audit results
# =============================================================================


@router.get("/dataset/{job_id}/results", response_model=DatasetAuditResponse)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def get_dataset_results(
    request: Request,
    job_id: str,
    api_key: Optional[str] = Depends(get_api_key),
):
    """Retrieve results for a completed dataset audit job.

    Reads results from Redis key ``dataset:results:{job_id}``.

    Args:
        job_id: Dataset audit job identifier (UUID).

    Returns:
        DatasetAuditResponse with health_audit, contradictions, numeric_columns,
        curation_report. Returns 202 if job is still processing.

    Raises:
        HTTPException: 404 if job not found, 202 if still processing.
    """
    _validate_uuid(job_id)

    try:
        r = _get_redis()
        meta = r.hgetall(f"dataset:meta:{job_id}")
        raw_results = r.get(f"dataset:results:{job_id}")
    except Exception as exc:
        raise HTTPException(status_code=503, detail="Redis unavailable") from exc

    if not meta:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    status = meta.get("status", "unknown")

    # Return 202 if still processing
    if status in ("pending", "processing"):
        return Response(
            content=json.dumps({
                "job_id": job_id,
                "status": status,
                "message": "Job is still processing",
            }),
            status_code=202,
            media_type="application/json",
        )

    if raw_results is None:
        raise HTTPException(
            status_code=404,
            detail=f"Results for job {job_id} not found",
        )

    result_data = json.loads(raw_results)

    return DatasetAuditResponse(
        job_id=job_id,
        status=status,
        health_audit=result_data.get("health_audit"),
        contradictions=result_data.get("contradictions", []),
        numeric_columns=result_data.get("numeric_columns", []),
        curation_report=result_data.get("curation_report"),
        curated_csv_available=result_data.get("curated_csv_available", False),
    )


# =============================================================================
# Endpoint 4: Dataset diff
# =============================================================================


@router.post("/dataset/{job_id}/diff", response_model=DatasetDiffResponse)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def compute_diff(
    request: Request,
    job_id: str,
    file: UploadFile = File(..., description="Comparison CSV or SDF file"),
    api_key: Optional[str] = Depends(get_api_key),
) -> DatasetDiffResponse:
    """Compare a previously uploaded dataset with a new file.

    Reads the primary dataset molecules from Redis results, parses the
    comparison file, and computes the InChIKey-based diff.

    Args:
        job_id: Primary dataset job ID (must be complete).
        file: Comparison CSV or SDF file.

    Returns:
        DatasetDiffResponse with added, removed, modified, unchanged counts.

    Raises:
        HTTPException: 404 if primary job not found, 400 if job not complete.
    """
    _validate_uuid(job_id)

    # Read primary dataset results
    try:
        r = _get_redis()
        meta = r.hgetall(f"dataset:meta:{job_id}")
        raw_results = r.get(f"dataset:results:{job_id}")
    except Exception as exc:
        raise HTTPException(status_code=503, detail="Redis unavailable") from exc

    if not meta:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    if meta.get("status") != "complete":
        raise HTTPException(
            status_code=400,
            detail="Primary dataset job is not yet complete",
        )

    if raw_results is None:
        raise HTTPException(
            status_code=404,
            detail=f"Results for job {job_id} not found",
        )

    result_data = json.loads(raw_results)
    primary_molecules = result_data.get("molecules", [])

    # Parse comparison file
    comp_filename = file.filename or "comparison"
    comp_ext = _get_file_extension(comp_filename)
    if comp_ext not in _ALLOWED_EXTENSIONS:
        raise HTTPException(
            status_code=400,
            detail="Invalid file type. Please upload a CSV or SDF file.",
        )

    try:
        comp_content = await file.read()
    except Exception as exc:
        raise HTTPException(
            status_code=400, detail="Failed to read comparison file"
        ) from exc

    # Save to temp and parse
    suffix = f".{comp_ext}"
    try:
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
        tmp.write(comp_content)
        tmp.flush()
        tmp.close()
        comp_file_path = tmp.name
    except Exception as exc:
        raise HTTPException(
            status_code=500, detail="Failed to save comparison file"
        ) from exc

    try:
        from app.services.dataset_intelligence.batch_processor import (
            _parse_csv_full,
            _parse_sdf_full,
        )

        if comp_ext == "sdf":
            comparison_molecules, _ = _parse_sdf_full(comp_file_path)
        else:
            comparison_molecules, _ = _parse_csv_full(comp_file_path, None)

        # Strip mol objects for the diff function (it only needs dict fields)
        comparison_mols_clean = [
            {
                "index": mol_dict["index"],
                "smiles": mol_dict["smiles"],
                "inchikey": mol_dict.get("inchikey"),
                "properties": mol_dict.get("properties", {}),
            }
            for mol_dict in comparison_molecules
        ]
    finally:
        try:
            os.unlink(comp_file_path)
        except OSError:
            pass

    # Compute diff
    from app.services.dataset_intelligence.dataset_diff import compute_dataset_diff

    diff_result = compute_dataset_diff(primary_molecules, comparison_mols_clean)

    # Convert to response schema
    def _to_diff_molecule(mol: dict) -> DiffMolecule:
        return DiffMolecule(
            inchikey=mol.get("inchikey", ""),
            smiles=mol.get("smiles", ""),
            row_index=mol.get("row_index", 0),
            properties=mol.get("properties", {}),
            changes=mol.get("changes", []),
        )

    return DatasetDiffResponse(
        added=[_to_diff_molecule(m) for m in diff_result["added"]],
        removed=[_to_diff_molecule(m) for m in diff_result["removed"]],
        modified=[_to_diff_molecule(m) for m in diff_result["modified"]],
        added_count=diff_result["added_count"],
        removed_count=diff_result["removed_count"],
        modified_count=diff_result["modified_count"],
        unchanged_count=diff_result["unchanged_count"],
        unique_columns_primary=diff_result.get("unique_columns_primary", 0),
        unique_columns_comparison=diff_result.get("unique_columns_comparison", 0),
    )


# =============================================================================
# Endpoint 5: Download curation report
# =============================================================================


@router.get("/dataset/{job_id}/download/report")
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def download_report(
    request: Request,
    job_id: str,
    api_key: Optional[str] = Depends(get_api_key),
) -> Response:
    """Download the curation report as a JSON file.

    Args:
        job_id: Dataset audit job identifier.

    Returns:
        JSON file download response.

    Raises:
        HTTPException: 404 if job or results not found.
    """
    _validate_uuid(job_id)

    try:
        r = _get_redis()
        raw_results = r.get(f"dataset:results:{job_id}")
    except Exception as exc:
        raise HTTPException(status_code=503, detail="Redis unavailable") from exc

    if raw_results is None:
        raise HTTPException(
            status_code=404,
            detail=f"Results for job {job_id} not found",
        )

    result_data = json.loads(raw_results)
    curation_report = result_data.get("curation_report")

    if curation_report is None:
        raise HTTPException(
            status_code=404,
            detail="Curation report not available for this job",
        )

    json_bytes = json.dumps(curation_report, indent=2).encode("utf-8")
    return Response(
        content=json_bytes,
        media_type="application/json",
        headers={
            "Content-Disposition": (
                f"attachment; filename=curation_report_{job_id[:8]}.json"
            ),
        },
    )


# =============================================================================
# Endpoint 6: Download curated CSV
# =============================================================================


@router.get("/dataset/{job_id}/download/csv")
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def download_curated_csv(
    request: Request,
    job_id: str,
    api_key: Optional[str] = Depends(get_api_key),
) -> Response:
    """Download the curated dataset as a CSV file.

    Curated CSV includes original properties plus appended columns:
    _health_issues, _is_duplicate, _standardized_smiles, _alert_flags.

    Args:
        job_id: Dataset audit job identifier.

    Returns:
        CSV file download response.

    Raises:
        HTTPException: 404 if job or results not found.
    """
    _validate_uuid(job_id)

    try:
        r = _get_redis()
        raw_results = r.get(f"dataset:results:{job_id}")
    except Exception as exc:
        raise HTTPException(status_code=503, detail="Redis unavailable") from exc

    if raw_results is None:
        raise HTTPException(
            status_code=404,
            detail=f"Results for job {job_id} not found",
        )

    result_data = json.loads(raw_results)
    curated_rows = result_data.get("curated_csv_rows", [])

    if not curated_rows:
        raise HTTPException(
            status_code=404,
            detail="Curated CSV data not available for this job",
        )

    # Collect all field names from all rows
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in curated_rows:
        for key in row.keys():
            if key not in seen:
                fieldnames.append(key)
                seen.add(key)

    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    for row in curated_rows:
        writer.writerow(row)

    csv_str = output.getvalue()

    return Response(
        content=csv_str,
        media_type="text/csv",
        headers={
            "Content-Disposition": f"attachment; filename=curated_{job_id[:8]}.csv",
        },
    )
