"""
QSAR-Ready Batch Processor

Processes a list of SMILES strings through the QSAR-Ready pipeline sequentially,
tracking progress via the shared progress_tracker service and deduplicating
results by standardized InChIKey.

Progress is published to the Redis channel `qsar:{job_id}` — distinct from the
regular batch pipeline's `batch:{job_id}` channel (Pitfall 5 from RESEARCH.md).
"""

import logging
from dataclasses import asdict
from typing import Any, Dict, List, Optional

from app.services.batch.progress_tracker import progress_tracker
from app.services.qsar_ready.pipeline import QSARReadyConfig, QSARReadyResult, qsar_ready_single

logger = logging.getLogger(__name__)


def _compute_summary(results: List[QSARReadyResult]) -> Dict[str, Any]:
    """
    Compute summary statistics from a list of QSAR pipeline results.

    Args:
        results: List of QSARReadyResult instances.

    Returns:
        Dict with total, ok, rejected, duplicate, error counts and steps_applied_counts.
    """
    ok = sum(1 for r in results if r.status == "ok")
    rejected = sum(1 for r in results if r.status == "rejected")
    duplicate = sum(1 for r in results if r.status == "duplicate")
    error = sum(1 for r in results if r.status == "error")

    # Aggregate step-level applied counts
    steps_applied_counts: Dict[str, int] = {}
    for result in results:
        for step in result.steps:
            if step.status == "applied":
                steps_applied_counts[step.step_name] = (
                    steps_applied_counts.get(step.step_name, 0) + 1
                )

    return {
        "total": len(results),
        "ok": ok,
        "rejected": rejected,
        "duplicate": duplicate,
        "error": error,
        "steps_applied_counts": steps_applied_counts,
    }


def run_qsar_batch(
    job_id: str,
    smiles_list: List[str],
    config: QSARReadyConfig,
) -> Dict[str, Any]:
    """
    Process a batch of SMILES strings through the QSAR-Ready pipeline.

    Initializes progress tracking, processes each molecule sequentially,
    applies InChIKey-based deduplication, and marks the job complete.

    Progress updates are published to the `qsar:{job_id}` Redis pub/sub channel,
    which the `/ws/qsar/{job_id}` WebSocket endpoint subscribes to.

    Args:
        job_id: Unique job identifier for progress tracking.
        smiles_list: List of SMILES strings to process.
        config: QSARReadyConfig controlling pipeline behavior.

    Returns:
        Dict with "results" (list of serialized results), "summary" (counts dict),
        and "config" (serialized config dict).
    """
    total = len(smiles_list)

    # Initialize progress tracking (publishes to qsar:{job_id} via progress_tracker)
    progress_tracker.init_job(job_id, total)

    results: List[QSARReadyResult] = []
    seen_inchikeys: Dict[str, int] = {}  # InChIKey -> original index

    for i, smiles in enumerate(smiles_list):
        result = qsar_ready_single(smiles, config)

        # Deduplication by standardized InChIKey (after all 10 steps)
        if result.status == "ok" and result.standardized_inchikey:
            key = result.standardized_inchikey
            if key in seen_inchikeys:
                result.status = "duplicate"
                result.rejection_reason = (
                    f"Duplicate of molecule {seen_inchikeys[key]}"
                )
            else:
                seen_inchikeys[key] = i

        results.append(result)

        # Update progress — throttling is handled by progress_tracker internally
        progress_tracker.update_progress(
            job_id=job_id,
            processed=i + 1,
            total=total,
            status="processing",
        )

    # Mark job complete (publishes final status to qsar:{job_id} channel)
    progress_tracker.mark_complete(job_id)

    summary = _compute_summary(results)

    # Serialize results: convert dataclasses to dicts for JSON/Celery serialization
    serialized_results = []
    for r in results:
        r_dict = asdict(r)
        serialized_results.append(r_dict)

    return {
        "results": serialized_results,
        "summary": summary,
        "config": {
            "enable_metals": config.enable_metals,
            "enable_desalt": config.enable_desalt,
            "enable_normalize": config.enable_normalize,
            "enable_neutralize": config.enable_neutralize,
            "enable_tautomer": config.enable_tautomer,
            "enable_stereo_strip": config.enable_stereo_strip,
            "enable_isotope_strip": config.enable_isotope_strip,
            "min_heavy_atoms": config.min_heavy_atoms,
            "max_heavy_atoms": config.max_heavy_atoms,
            "max_mw": config.max_mw,
            "remove_inorganics": config.remove_inorganics,
        },
    }
