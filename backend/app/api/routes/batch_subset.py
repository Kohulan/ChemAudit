"""Batch subset action routes.

Endpoints for acting on a user-selected subset of a batch result
(revalidate, rescore, export, inline score). Split out of batch.py to keep
each route module focused.
"""

import logging
from typing import Any, Dict, List, Optional

from fastapi import (
    APIRouter,
    Depends,
    HTTPException,
    Request,
)
from pydantic import BaseModel

from app.core.error_sanitizer import safe_error_detail
from app.core.ownership import verify_job_access
from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.services.batch.parsing import parse_json_field
from app.services.batch.result_aggregator import result_storage
from app.services.batch.subset_actions import (
    export_subset,
    rescore_subset,
    revalidate_subset,
)

logger = logging.getLogger(__name__)

router = APIRouter()


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
    verify_job_access(request, job_id)

    try:
        new_job_id = revalidate_subset(job_id, body.indices)
    except ValueError:
        raise HTTPException(status_code=404, detail="Job not found or invalid indices")

    return {
        "new_job_id": new_job_id,
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
    verify_job_access(request, job_id)

    # Build safety_options with full profile data if profile_id is given
    safety_options: Dict[str, Any] = {}
    if body.profile_id is not None:
        try:
            from app.db import async_session
            from app.services.profiles.service import ProfileService

            async with async_session() as session:
                profile = await ProfileService().get(session, body.profile_id)
                if profile:
                    safety_options["profile_id"] = profile.id
                    safety_options["profile_name"] = profile.name
                    safety_options["profile_thresholds"] = parse_json_field(
                        profile.thresholds
                    )
                    safety_options["profile_weights"] = parse_json_field(
                        profile.weights
                    )
        except Exception as exc:
            logger.exception("Failed to fetch scoring profile %s", body.profile_id)
            raise HTTPException(
                status_code=500,
                detail="Failed to fetch the requested scoring profile",
            ) from exc

    try:
        new_job_id = rescore_subset(
            job_id, body.indices, safety_options=safety_options
        )
    except ValueError:
        raise HTTPException(status_code=404, detail="Job not found or invalid indices")

    return {
        "new_job_id": new_job_id,
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
    verify_job_access(request, job_id)

    from datetime import datetime

    from fastapi.responses import StreamingResponse

    try:
        export_buffer = export_subset(job_id, body.indices, body.format)
    except ValueError:
        raise HTTPException(status_code=404, detail="Job not found or invalid indices")
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=safe_error_detail(e, "Export failed"),
        )

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


@router.post("/batch/{job_id}/subset/score-inline")
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def score_subset_inline(
    request: Request,
    job_id: str,
    body: SubsetRescoreRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Score a subset of molecules with a profile inline (no Celery).

    Reads existing batch results, applies profile desirability scoring,
    and returns scored molecules synchronously.
    """
    verify_job_access(request, job_id)

    from app.db import async_session
    from app.services.profiles.service import ProfileService
    from app.services.scoring.profile_scoring import compute_profile_result

    if body.profile_id is None:
        raise HTTPException(status_code=400, detail="profile_id is required")

    # Fetch profile
    try:
        async with async_session() as session:
            profile = await ProfileService().get(session, body.profile_id)
    except Exception as exc:
        logger.exception("Failed to fetch scoring profile %s", body.profile_id)
        raise HTTPException(
            status_code=500,
            detail="Failed to fetch the requested scoring profile",
        ) from exc

    if not profile:
        raise HTTPException(status_code=404, detail="Profile not found")

    thresholds = parse_json_field(profile.thresholds)
    weights = parse_json_field(profile.weights)

    # Fetch existing batch results for the selected indices
    result_data = result_storage.get_results(
        job_id=job_id, page=1, page_size=10000
    )
    all_results = result_data.get("results", [])
    indices_set = set(body.indices)
    subset = [r for r in all_results if r.get("index") in indices_set]

    if not subset:
        raise HTTPException(
            status_code=404,
            detail=f"No results found for job {job_id} with given indices",
        )

    # Score each molecule using existing property data
    scored = []
    for mol in subset:
        scoring = mol.get("scoring") or {}
        dl = scoring.get("druglikeness") or {}
        admet = scoring.get("admet") or {}
        mol_properties = {
            "mw": dl.get("mw"),
            "logp": dl.get("logp"),
            "hbd": dl.get("hbd"),
            "hba": dl.get("hba"),
            "tpsa": dl.get("tpsa"),
            "rotatable_bonds": dl.get("rotatable_bonds"),
            "aromatic_rings": dl.get("aromatic_rings"),
            "fsp3": admet.get("fsp3"),
        }
        try:
            profile_result = compute_profile_result(
                properties=mol_properties,
                profile_id=profile.id,
                profile_name=profile.name,
                thresholds=thresholds,
                weights=weights,
            )
        except Exception:
            logger.exception(
                "Profile scoring failed for molecule index=%s in job %s",
                mol.get("index"),
                job_id,
            )
            profile_result = {
                "profile_id": profile.id,
                "profile_name": profile.name,
                "score": None,
                "properties": {},
                "error": "Scoring failed for this molecule",
            }
        scored.append(
            {
                "index": mol.get("index"),
                "name": mol.get("name"),
                "smiles": mol.get("smiles"),
                "profile": profile_result,
            }
        )

    return {
        "profile_name": profile.name,
        "profile_id": profile.id,
        "molecules": scored,
    }
