"""
GenChem Filter API Routes (Phase 11).

Endpoints for multi-stage generative chemistry filtering, REINVENT-compatible
scoring, and async batch processing.

Routes:
  POST /api/v1/genchem/filter          — multi-stage funnel (sync ≤1000, async >1000)
  POST /api/v1/genchem/score           — composite 0-1 scores per SMILES
  POST /api/v1/genchem/reinvent-score                      — REINVENT 4-compatible scoring endpoint
  POST /api/v1/genchem/batch/upload                        — queue async batch job
  GET  /api/v1/genchem/batch/{job_id}/status               — poll batch progress
  GET  /api/v1/genchem/batch/{job_id}/results              — batch results (Redis-stored)
  GET  /api/v1/genchem/batch/{job_id}/download/{format}    — download passed_txt or full_csv

WebSocket /ws/genchem/{job_id} is registered in main.py.
"""

import csv
import dataclasses
import io
import json
import logging
import uuid
from typing import Optional

from fastapi import APIRouter, Depends, File, Form, HTTPException, Query, Request, UploadFile
from fastapi.responses import StreamingResponse

from app.core.config import settings
from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.schemas.genchem import (
    FilterConfigSchema,
    FilterRequest,
    FilterResponse,
    GenChemBatchResultsResponse,
    GenChemBatchStatusResponse,
    GenChemBatchUploadResponse,
    MoleculeResultSchema,
    REINVENTInput,
    REINVENTOutput,
    REINVENTResponse,
    REINVENTSuccessItem,
    ScoreRequest,
    ScoreResponse,
    StageResultSchema,
)
from app.services.genchem.filter_config import PRESETS, FilterConfig
from app.services.genchem.filter_pipeline import filter_batch
from app.services.genchem.scorer import score_for_generative

logger = logging.getLogger(__name__)

router = APIRouter()

# Sync/async split threshold (D-22)
_SYNC_THRESHOLD = 1000


# =============================================================================
# Config resolution helper
# =============================================================================


def _resolve_config(
    preset: Optional[str], config_schema: Optional[FilterConfigSchema]
) -> FilterConfig:
    """Resolve FilterConfig from a preset name or FilterConfigSchema.

    If preset is provided, it takes priority over config_schema.
    If neither is provided, defaults to the 'drug_like' preset.

    Args:
        preset: Named preset ('drug_like', 'lead_like', 'fragment_like', 'permissive').
        config_schema: Optional explicit FilterConfigSchema.

    Returns:
        FilterConfig dataclass instance.

    Raises:
        HTTPException: 400 if an unknown preset name is provided.
    """
    if preset is not None:
        if preset not in PRESETS:
            raise HTTPException(
                status_code=400,
                detail=(
                    f"Unknown preset '{preset}'. "
                    f"Valid options: {sorted(PRESETS.keys())}"
                ),
            )
        return PRESETS[preset]

    if config_schema is not None:
        return FilterConfig(**config_schema.model_dump())

    # Default
    return PRESETS["drug_like"]


def _get_redis():
    """Get a synchronous Redis client."""
    import redis as sync_redis

    return sync_redis.from_url(settings.REDIS_URL, decode_responses=True)


# =============================================================================
# Endpoint 1: Filter
# =============================================================================


@router.post("/genchem/filter")
@limiter.limit("20/minute", key_func=get_rate_limit_key)
async def genchem_filter(
    request: Request,
    body: FilterRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Run the multi-stage GenChem filter pipeline on a list of SMILES.

    For ≤1,000 SMILES: runs synchronously and returns FilterResponse.
    For >1,000 SMILES: queues a Celery task and returns GenChemBatchUploadResponse
    with a job_id for async tracking (D-22).

    Stages: parse → valence → alerts → property → SA score → dedup (→ novelty if enabled).

    Args:
        body: FilterRequest with smiles_list, optional preset or config.

    Returns:
        FilterResponse (sync) or GenChemBatchUploadResponse (async).

    Raises:
        HTTPException: 400 for unknown preset or invalid input.
    """
    config = _resolve_config(body.preset, body.config)

    if len(body.smiles_list) <= _SYNC_THRESHOLD:
        # Synchronous path
        try:
            result = filter_batch(body.smiles_list, config)
        except Exception as exc:
            logger.exception("Error in genchem filter_batch: %s", exc)
            raise HTTPException(
                status_code=500,
                detail={"error": "Filter pipeline failed", "detail": str(exc)},
            ) from exc

        stages = [StageResultSchema(**dataclasses.asdict(s)) for s in result.stages]
        molecules = [MoleculeResultSchema(**dataclasses.asdict(m)) for m in result.molecules]
        return FilterResponse(
            input_count=result.input_count,
            output_count=result.output_count,
            stages=stages,
            molecules=molecules,
        )

    # Async path (>1000 SMILES)
    job_id = str(uuid.uuid4())
    config_dict = dataclasses.asdict(config)

    try:
        r = _get_redis()
        r.set(
            f"genchem:meta:{job_id}",
            json.dumps({"status": "pending", "total": len(body.smiles_list)}),
            ex=3600,
        )
    except Exception as exc:
        logger.warning("Failed to store genchem job metadata: %s", exc)

    try:
        from app.services.genchem.batch_processor import process_genchem_batch

        process_genchem_batch.delay(body.smiles_list, config_dict, job_id)
    except Exception as exc:
        logger.error("Failed to dispatch genchem batch task for %s: %s", job_id, exc)
        raise HTTPException(
            status_code=500,
            detail="Failed to start batch processing",
        ) from exc

    return GenChemBatchUploadResponse(
        job_id=job_id,
        total_molecules=len(body.smiles_list),
        status="pending",
    )


# =============================================================================
# Endpoint 2: Score
# =============================================================================


@router.post("/genchem/score", response_model=ScoreResponse)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def genchem_score(
    request: Request,
    body: ScoreRequest,
    api_key: Optional[str] = Depends(get_api_key),
) -> ScoreResponse:
    """
    Compute composite 0-1 scores for a list of SMILES.

    Returns a score per SMILES in input order. Invalid or unparseable SMILES
    return null (None) — not 0.0 — per D-14 / Pitfall 4.

    Args:
        body: ScoreRequest with smiles_list and optional preset.

    Returns:
        ScoreResponse with scores list (null for invalid SMILES).

    Raises:
        HTTPException: 400 for unknown preset.
    """
    config = _resolve_config(body.preset, None)

    scores: list[Optional[float]] = []
    for smi in body.smiles_list:
        try:
            score = score_for_generative(smi, config)
        except Exception as exc:
            logger.warning("Error scoring SMILES '%s': %s", smi, exc)
            score = None
        scores.append(score)

    return ScoreResponse(scores=scores)


# =============================================================================
# Endpoint 3: REINVENT Score
# =============================================================================


@router.post("/genchem/reinvent-score", response_model=REINVENTResponse)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def genchem_reinvent_score(
    request: Request,
    body: list[REINVENTInput],
    preset: str = Query(default="drug_like"),
    api_key: Optional[str] = Depends(get_api_key),
) -> REINVENTResponse:
    """
    REINVENT 4-compatible scoring endpoint.

    Accepts a raw list of {input_string, query_id} items (not wrapped in a model).
    Returns {output: {successes_list: [{query_id, output_value}]}} per REINVENT 4 contract.

    Invalid SMILES are OMITTED from successes_list — not scored as 0.0 (Pitfall 4 / D-14).

    Args:
        body: List of REINVENTInput items.
        preset: Named preset for weight vector (default 'drug_like').

    Returns:
        REINVENTResponse with output.successes_list.

    Raises:
        HTTPException: 400 for unknown preset.
    """
    config = _resolve_config(preset, None)

    successes: list[REINVENTSuccessItem] = []
    for item in body:
        try:
            score = score_for_generative(item.input_string, config)
        except Exception as exc:
            logger.warning("Error scoring REINVENT input '%s': %s", item.input_string, exc)
            score = None

        # Only append if score is not None (invalid SMILES are omitted)
        if score is not None:
            successes.append(REINVENTSuccessItem(query_id=item.query_id, output_value=score))

    return REINVENTResponse(output=REINVENTOutput(successes_list=successes))


# =============================================================================
# Endpoint 4: Batch upload
# =============================================================================


@router.post("/genchem/batch/upload", response_model=GenChemBatchUploadResponse)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def genchem_batch_upload(
    request: Request,
    file: UploadFile = File(..., description="CSV or SDF file with SMILES"),
    preset: Optional[str] = Form(default=None, description="Named preset name"),
    config: Optional[str] = Form(default=None, description="JSON-encoded FilterConfigSchema"),
    api_key: Optional[str] = Depends(get_api_key),
) -> GenChemBatchUploadResponse:
    """
    Upload a file for async GenChem batch filtering.

    Parses the file, queues a Celery task on the 'default' queue, and returns
    immediately with a job_id for async tracking.

    Args:
        file: CSV or SDF file containing SMILES.
        preset: Optional preset name (takes priority over config).
        config: Optional JSON-encoded FilterConfigSchema.

    Returns:
        GenChemBatchUploadResponse with job_id, total_molecules, status.

    Raises:
        HTTPException: 400 for missing/invalid input.
    """
    # Resolve config from form fields
    config_schema: Optional[FilterConfigSchema] = None
    if config:
        try:
            config_data = json.loads(config)
            config_schema = FilterConfigSchema(**config_data)
        except Exception as exc:
            raise HTTPException(
                status_code=400,
                detail=f"Invalid config JSON: {exc}",
            ) from exc

    filter_config = _resolve_config(preset, config_schema)

    # Parse file
    try:
        content = await file.read()
    except Exception as exc:
        raise HTTPException(status_code=400, detail="Failed to read file") from exc

    filename = file.filename or ""
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

    smiles_list = [m.smiles for m in molecules if m.smiles and not m.parse_error]

    if not smiles_list:
        raise HTTPException(status_code=400, detail="No valid SMILES found in input")

    job_id = str(uuid.uuid4())
    config_dict = dataclasses.asdict(filter_config)

    # Store metadata in Redis
    try:
        r = _get_redis()
        r.set(
            f"genchem:meta:{job_id}",
            json.dumps({"status": "pending", "total": len(smiles_list)}),
            ex=3600,
        )
    except Exception as exc:
        logger.warning("Failed to store genchem batch metadata in Redis: %s", exc)

    # Dispatch Celery task
    try:
        from app.services.genchem.batch_processor import process_genchem_batch

        process_genchem_batch.delay(smiles_list, config_dict, job_id)
    except Exception as exc:
        logger.error("Failed to dispatch genchem batch task for %s: %s", job_id, exc)
        raise HTTPException(
            status_code=500,
            detail="Failed to start batch processing",
        ) from exc

    return GenChemBatchUploadResponse(
        job_id=job_id,
        total_molecules=len(smiles_list),
        status="pending",
    )


# =============================================================================
# Endpoint 5: Batch status
# =============================================================================


@router.get("/genchem/batch/{job_id}/status", response_model=GenChemBatchStatusResponse)
@limiter.limit("60/minute", key_func=get_rate_limit_key)
async def genchem_batch_status(
    request: Request,
    job_id: str,
    api_key: Optional[str] = Depends(get_api_key),
) -> GenChemBatchStatusResponse:
    """
    Get the current processing status of a GenChem async batch job.

    Reads job metadata from Redis key `genchem:meta:{job_id}`.

    Args:
        job_id: Batch job identifier (UUID).

    Returns:
        GenChemBatchStatusResponse with status, progress, current_stage.

    Raises:
        HTTPException: 404 if job not found.
    """
    _validate_uuid(job_id)

    try:
        r = _get_redis()
        raw_meta = r.get(f"genchem:meta:{job_id}")
    except Exception as exc:
        raise HTTPException(status_code=503, detail="Redis unavailable") from exc

    if raw_meta is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    meta = json.loads(raw_meta)
    return GenChemBatchStatusResponse(
        job_id=job_id,
        status=meta.get("status", "unknown"),
        progress=meta.get("progress"),
        current_stage=meta.get("current_stage"),
    )


# =============================================================================
# Endpoint 6: Batch results
# =============================================================================


@router.get("/genchem/batch/{job_id}/results", response_model=GenChemBatchResultsResponse)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def genchem_batch_results(
    request: Request,
    job_id: str,
    api_key: Optional[str] = Depends(get_api_key),
) -> GenChemBatchResultsResponse:
    """
    Retrieve results for a completed GenChem async batch job.

    Reads results from Redis key `genchem:results:{job_id}` (not Celery result backend).

    Args:
        job_id: Batch job identifier (UUID).

    Returns:
        GenChemBatchResultsResponse with status and optional result.

    Raises:
        HTTPException: 404 if job or results not found.
    """
    _validate_uuid(job_id)

    try:
        r = _get_redis()
        raw_meta = r.get(f"genchem:meta:{job_id}")
        raw_results = r.get(f"genchem:results:{job_id}")
    except Exception as exc:
        raise HTTPException(status_code=503, detail="Redis unavailable") from exc

    if raw_meta is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    meta = json.loads(raw_meta)
    status = meta.get("status", "unknown")

    if raw_results is None:
        return GenChemBatchResultsResponse(job_id=job_id, status=status, result=None)

    result_data = json.loads(raw_results)
    stages = [StageResultSchema(**s) for s in result_data.get("stages", [])]
    molecules = [MoleculeResultSchema(**m) for m in result_data.get("molecules", [])]
    filter_response = FilterResponse(
        input_count=result_data.get("input_count", 0),
        output_count=result_data.get("output_count", 0),
        stages=stages,
        molecules=molecules,
    )
    return GenChemBatchResultsResponse(
        job_id=job_id,
        status=status,
        result=filter_response,
    )


# =============================================================================
# Endpoint 7: Batch download
# =============================================================================


@router.get("/genchem/batch/{job_id}/download/{format}")
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def genchem_batch_download(
    request: Request,
    job_id: str,
    format: str,
    api_key: Optional[str] = Depends(get_api_key),
) -> StreamingResponse:
    """
    Download batch GenChem results in 'passed_txt' or 'full_csv' format.

    passed_txt: one SMILES per line, only molecules with status == 'passed'.
    full_csv: all molecules with columns smiles,status,failed_at,rejection_reason.

    Args:
        job_id: Batch job identifier (UUID).
        format: 'passed_txt' or 'full_csv'.

    Returns:
        StreamingResponse with appropriate Content-Type header.

    Raises:
        HTTPException: 400 for invalid format, 404 for missing results.
    """
    _validate_uuid(job_id)

    if format not in ("passed_txt", "full_csv"):
        raise HTTPException(
            status_code=400,
            detail=f"Invalid format '{format}'. Supported: passed_txt, full_csv",
        )

    try:
        r = _get_redis()
        raw_results = r.get(f"genchem:results:{job_id}")
    except Exception as exc:
        raise HTTPException(status_code=503, detail="Redis unavailable") from exc

    if raw_results is None:
        raise HTTPException(
            status_code=404,
            detail=f"Results for job {job_id} not found",
        )

    result_data = json.loads(raw_results)
    molecules = result_data.get("molecules", [])

    if format == "passed_txt":
        return _build_passed_txt_response(job_id, molecules)
    else:
        return _build_full_csv_response(job_id, molecules)


# =============================================================================
# Download helpers
# =============================================================================


def _build_passed_txt_response(job_id: str, molecules: list) -> StreamingResponse:
    """Build a plain-text StreamingResponse with one passed SMILES per line."""
    passed_smiles = [m["smiles"] for m in molecules if m.get("status") == "passed"]
    content = "\n".join(passed_smiles)
    return StreamingResponse(
        iter([content]),
        media_type="text/plain",
        headers={
            "Content-Disposition": f"attachment; filename=genchem_passed_{job_id[:8]}.txt"
        },
    )


def _build_full_csv_response(job_id: str, molecules: list) -> StreamingResponse:
    """Build a CSV StreamingResponse with all molecule results."""
    output = io.StringIO()
    writer = csv.DictWriter(
        output,
        fieldnames=["smiles", "status", "failed_at", "rejection_reason"],
    )
    writer.writeheader()
    for m in molecules:
        writer.writerow(
            {
                "smiles": m.get("smiles", ""),
                "status": m.get("status", ""),
                "failed_at": m.get("failed_at") or "",
                "rejection_reason": m.get("rejection_reason") or "",
            }
        )
    output.seek(0)
    return StreamingResponse(
        iter([output.getvalue()]),
        media_type="text/csv",
        headers={
            "Content-Disposition": f"attachment; filename=genchem_results_{job_id[:8]}.csv"
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
