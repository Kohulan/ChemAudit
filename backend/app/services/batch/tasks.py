"""
Celery Tasks for Batch Processing

Implements chunk-based molecule processing with progress tracking.
Uses chord pattern: process_molecule_chunk tasks -> aggregate_batch_results.
Supports priority queues for concurrent job handling.
"""

import asyncio
import logging
import time
from dataclasses import asdict
from typing import Any, Dict, List, Optional

from celery import chord, group
from rdkit import Chem
from rdkit.Chem import Lipinski

from app.celery_app import SMALL_JOB_THRESHOLD, celery_app
from app.core.config import settings
from app.services.alerts import alert_manager
from app.services.batch.progress_tracker import progress_tracker
from app.services.batch.result_aggregator import compute_statistics, result_storage
from app.services.scoring.admet import calculate_admet
from app.services.scoring.druglikeness import calculate_druglikeness
from app.services.scoring.ml_readiness import calculate_ml_readiness
from app.services.scoring.profile_scoring import compute_profile_result
from app.services.scoring.safety_filters import calculate_safety_filters
from app.services.standardization import standardize_molecule
from app.services.validation.engine import validation_engine

CHUNK_SIZE = 100  # Process 100 molecules per chunk for progress updates

logger = logging.getLogger(__name__)


# =============================================================================
# Single Molecule Validation Task (Priority Queue)
# =============================================================================


@celery_app.task(
    bind=True,
    queue="high_priority",
    autoretry_for=(Exception,),
    retry_backoff=True,
    max_retries=2,
)
def validate_single_molecule(
    self,
    smiles: str,
    checks: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Validate a single molecule using the priority queue.

    This task runs on the high_priority queue to ensure single molecule
    validations are not blocked by large batch jobs.

    Args:
        smiles: SMILES string of the molecule
        checks: Optional list of specific checks to run

    Returns:
        Dictionary with validation results
    """
    start_time = time.time()

    result = {
        "smiles": smiles,
        "status": "error",
        "error": None,
        "validation": None,
        "alerts": None,
        "scoring": None,
    }

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result["error"] = "Failed to parse SMILES"
            return result

        # Run validation
        try:
            check_results, validation_score = validation_engine.validate(mol, checks)
            result["validation"] = {
                "overall_score": validation_score,
                "issues": [
                    {
                        "check_name": r.check_name,
                        "passed": r.passed,
                        "severity": (
                            r.severity.value
                            if hasattr(r.severity, "value")
                            else str(r.severity)
                        ),
                        "message": r.message,
                        "affected_atoms": r.affected_atoms,
                    }
                    for r in check_results
                    if not r.passed
                ],
            }
        except Exception as e:
            result["validation"] = {"error": str(e)}

        # Run alerts screening
        try:
            alert_result = alert_manager.screen(mol, catalogs=["PAINS", "BRENK"])
            result["alerts"] = {
                "has_alerts": alert_result.total_alerts > 0,
                "alert_count": alert_result.total_alerts,
                "alerts": [
                    {
                        "catalog": a.catalog_source,
                        "rule_name": a.pattern_name,
                        "severity": (
                            a.severity.value
                            if hasattr(a.severity, "value")
                            else str(a.severity)
                        ),
                        "matched_atoms": a.matched_atoms,
                    }
                    for a in alert_result.alerts
                ],
            }
        except Exception as e:
            result["alerts"] = {"error": str(e)}

        # Calculate ML-readiness score
        try:
            ml_result = calculate_ml_readiness(mol)
            result["scoring"] = {
                "ml_readiness": {
                    "score": ml_result.score,
                    "interpretation": ml_result.interpretation,
                }
            }
        except Exception as e:
            result["scoring"] = {"error": str(e)}

        result["status"] = "success"
        result["execution_time_ms"] = int((time.time() - start_time) * 1000)

    except Exception as e:
        result["error"] = f"Processing error: {str(e)}"

    return result


@celery_app.task(
    bind=True,
    autoretry_for=(Exception,),
    retry_backoff=True,
    max_retries=3,
    acks_late=True,
)
def process_molecule_chunk(
    self,
    job_id: str,
    chunk_idx: int,
    molecules: List[Dict[str, Any]],
    total_molecules: int,
    total_chunks: int,
    safety_options: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    """
    Process a chunk of molecules with progress updates.

    Args:
        job_id: Job identifier for progress tracking
        chunk_idx: Index of this chunk (for progress calculation)
        molecules: List of molecule data dicts (smiles, name, index, properties)
        total_molecules: Total molecules in the entire batch
        total_chunks: Total number of chunks
        safety_options: Optional safety screening options

    Returns:
        List of result dictionaries for each molecule
    """
    results = []

    for mol_data in molecules:
        result = _process_single_molecule(mol_data, safety_options=safety_options)
        results.append(result)

        # Atomically increment counter and update progress
        processed_so_far = progress_tracker.increment_processed(job_id)
        progress_tracker.update_progress(
            job_id=job_id,
            processed=processed_so_far,
            total=total_molecules,
            status="processing",
        )

    return results


def _process_single_molecule(
    mol_data: Dict[str, Any],
    safety_options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Process a single molecule through validation, alerts, and scoring.

    Args:
        mol_data: Dictionary with smiles, name, index, properties, and optional parse_error
        safety_options: Optional dict with include_extended and include_chembl flags

    Returns:
        Result dictionary with all validation results or error info
    """
    result = {
        "smiles": mol_data.get("smiles", ""),
        "name": mol_data.get("name", ""),
        "index": mol_data.get("index", 0),
        "status": "error",
        "error": None,
        "validation": None,
        "alerts": None,
        "scoring": None,
        "standardization": None,
    }

    # Check for pre-existing parse error from file parsing
    if mol_data.get("parse_error"):
        result["error"] = mol_data["parse_error"]
        return result

    smiles = mol_data.get("smiles", "")
    if not smiles:
        result["error"] = "Empty SMILES string"
        return result

    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result["error"] = "Failed to parse SMILES"
            return result

        # Run standardization if enabled
        opts = safety_options or {}
        if opts.get("include_standardization"):
            try:
                std_result = standardize_molecule(mol)
                result["standardization"] = {
                    "standardized_smiles": std_result.standardized_smiles,
                    "success": std_result.success,
                    "error": std_result.error_message,
                    "steps_applied": [
                        {
                            "step_name": s.step_name,
                            "applied": s.applied,
                            "description": s.description,
                            "changes": s.changes,
                        }
                        for s in std_result.steps_applied
                        if s.applied
                    ],
                    "excluded_fragments": std_result.excluded_fragments,
                    "changed": (
                        std_result.standardized_smiles is not None
                        and std_result.standardized_smiles != smiles
                    ),
                }
            except Exception as e:
                result["standardization"] = {"error": str(e)}

        # Run validation
        try:
            check_results, validation_score = validation_engine.validate(mol)
            result["validation"] = {
                "overall_score": validation_score,
                "issues": [
                    {
                        "check_name": r.check_name,
                        "passed": r.passed,
                        "severity": (
                            r.severity.value
                            if hasattr(r.severity, "value")
                            else str(r.severity)
                        ),
                        "message": r.message,
                        "affected_atoms": r.affected_atoms,
                    }
                    for r in check_results
                    if not r.passed
                ],
            }
        except Exception as e:
            result["validation"] = {"error": str(e)}

        # Run structural alerts screening
        try:
            opts = safety_options or {}
            catalogs = ["PAINS", "BRENK"]
            if opts.get("include_extended"):
                catalogs.extend(["NIH", "ZINC"])
            if opts.get("include_chembl"):
                catalogs.extend(
                    [
                        "CHEMBL_BMS",
                        "CHEMBL_DUNDEE",
                        "CHEMBL_GLAXO",
                        "CHEMBL_INPHARMATICA",
                        "CHEMBL_LINT",
                        "CHEMBL_MLSMR",
                        "CHEMBL_SURECHEMBL",
                    ]
                )
            alert_result = alert_manager.screen(mol, catalogs=catalogs)
            result["alerts"] = {
                "has_alerts": alert_result.total_alerts > 0,
                "alert_count": alert_result.total_alerts,
                "alerts": [
                    {
                        "catalog": a.catalog_source,
                        "rule_name": a.pattern_name,
                        "severity": (
                            a.severity.value
                            if hasattr(a.severity, "value")
                            else str(a.severity)
                        ),
                        "matched_atoms": a.matched_atoms,
                    }
                    for a in alert_result.alerts
                ],
            }
        except Exception as e:
            result["alerts"] = {"error": str(e)}

        # Calculate ML-readiness score
        try:
            ml_result = calculate_ml_readiness(mol)
            result["scoring"] = {
                "ml_readiness": {
                    "score": ml_result.score,
                    "interpretation": ml_result.interpretation,
                }
            }
        except Exception as e:
            result["scoring"] = {"error": str(e)}

        # Calculate drug-likeness scores
        try:
            dl_result = calculate_druglikeness(mol, include_extended=False)
            result["scoring"]["druglikeness"] = {
                "qed_score": dl_result.qed.score,
                "lipinski_passed": dl_result.lipinski.passed,
                "lipinski_violations": dl_result.lipinski.violations,
                "mw": dl_result.lipinski.mw,
                "logp": dl_result.lipinski.logp,
                "hbd": dl_result.lipinski.hbd,
                "hba": dl_result.lipinski.hba,
                "tpsa": dl_result.veber.tpsa if dl_result.veber else None,
                "rotatable_bonds": (
                    dl_result.veber.rotatable_bonds if dl_result.veber else None
                ),
                "aromatic_rings": Lipinski.NumAromaticRings(mol),
            }
        except Exception as e:
            if "scoring" not in result or result["scoring"] is None:
                result["scoring"] = {}
            result["scoring"]["druglikeness"] = {"error": str(e)}

        # Calculate safety filters
        try:
            opts = safety_options or {}
            sf_result = calculate_safety_filters(
                mol,
                include_extended=opts.get("include_extended", False),
                include_chembl=opts.get("include_chembl", False),
            )
            sf_data: Dict[str, Any] = {
                "all_passed": sf_result.all_passed,
                "total_alerts": sf_result.total_alerts,
                "pains_passed": sf_result.pains.passed,
                "brenk_passed": sf_result.brenk.passed,
            }
            if sf_result.nih is not None:
                sf_data["nih_passed"] = sf_result.nih.passed
            if sf_result.zinc is not None:
                sf_data["zinc_passed"] = sf_result.zinc.passed
            if sf_result.chembl is not None:
                sf_data["chembl_passed"] = sf_result.chembl.passed
            result["scoring"]["safety_filters"] = sf_data
        except Exception as e:
            if "scoring" not in result or result["scoring"] is None:
                result["scoring"] = {}
            result["scoring"]["safety_filters"] = {"error": str(e)}

        # Calculate ADMET predictions
        try:
            admet_result = calculate_admet(mol, include_cns_mpo=False)
            result["scoring"]["admet"] = {
                "sa_score": admet_result.synthetic_accessibility.score,
                "sa_classification": admet_result.synthetic_accessibility.classification,
                "solubility_class": admet_result.solubility.classification,
                "fsp3": admet_result.complexity.fsp3,
            }
        except Exception as e:
            if "scoring" not in result or result["scoring"] is None:
                result["scoring"] = {}
            result["scoring"]["admet"] = {"error": str(e)}

        # Calculate profile score (if profile was selected at upload)
        opts = safety_options or {}
        if opts.get("profile_id") is not None:
            try:
                dl = result.get("scoring", {}).get("druglikeness", {})
                admet = result.get("scoring", {}).get("admet", {})
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
                result["scoring"]["profile"] = compute_profile_result(
                    properties=mol_properties,
                    profile_id=opts["profile_id"],
                    profile_name=opts.get("profile_name", ""),
                    thresholds=opts.get("profile_thresholds", {}),
                    weights=opts.get("profile_weights", {}),
                )
            except Exception as e:
                if "scoring" not in result or result["scoring"] is None:
                    result["scoring"] = {}
                result["scoring"]["profile"] = {"error": str(e)}

        # If we got this far with validation, mark as success
        result["status"] = "success"

    except Exception as e:
        result["error"] = f"Processing error: {str(e)}"

    return result


async def _log_batch_audit(
    job_id: str,
    molecule_count: int,
    pass_count: int,
    fail_count: int,
) -> None:
    """
    Create a DB session and log a batch completion event to the audit trail.

    This async helper is called via asyncio.run() from Celery tasks, which run
    in their own process with no existing event loop.

    Args:
        job_id: Batch job identifier
        molecule_count: Total molecules processed
        pass_count: Number of molecules passing validation
        fail_count: Number of molecules with errors
    """
    from app.db import async_session
    from app.services.audit.service import log_batch_event

    async with async_session() as db:
        await log_batch_event(
            db=db,
            job_id=job_id,
            molecule_count=molecule_count,
            pass_count=pass_count,
            fail_count=fail_count,
        )


@celery_app.task(bind=True)
def aggregate_batch_results(
    self,
    chunk_results: List[List[Dict[str, Any]]],
    job_id: str,
    start_time: float,
) -> Dict[str, Any]:
    """
    Aggregate all chunk results and compute statistics.

    Args:
        chunk_results: List of result lists from each chunk
        job_id: Job identifier
        start_time: Job start timestamp for processing time calculation

    Returns:
        Summary dictionary with job_id, total, and statistics
    """
    # Flatten chunk results
    all_results = [result for chunk in chunk_results for result in chunk]

    # Compute statistics
    stats = compute_statistics(all_results)
    stats.processing_time_seconds = round(time.time() - start_time, 2)

    # Store results
    result_storage.store_results(job_id, all_results, stats)

    # Mark job as complete
    progress_tracker.mark_complete(job_id)

    # Log batch event to audit trail
    try:
        total = len(all_results)
        fail_count = stats.errors
        pass_count = total - fail_count
        asyncio.run(_log_batch_audit(
            job_id=job_id,
            molecule_count=total,
            pass_count=pass_count,
            fail_count=fail_count,
        ))
    except Exception as e:
        logger.warning("Failed to log batch audit event for %s: %s", job_id, e)

    # Dispatch webhook if configured
    try:
        webhook_url = settings.WEBHOOK_URL
        if webhook_url:
            from app.services.notifications.webhook import send_webhook

            total = len(all_results)
            payload = {
                "event": "batch_complete",
                "job_id": job_id,
                "status": "complete",
                "molecule_count": total,
                "summary_url": f"/batch/{job_id}",
            }
            send_webhook.delay(webhook_url, settings.WEBHOOK_SECRET, payload)
    except Exception as e:
        logger.warning("Failed to dispatch webhook for %s: %s", job_id, e)

    # Dispatch email notification if configured
    try:
        r = progress_tracker._get_redis()
        notification_email = r.get(f"batch:email:{job_id}")
        if notification_email:
            if isinstance(notification_email, bytes):
                notification_email = notification_email.decode("utf-8")
            from app.services.notifications.email import send_batch_complete_email

            total = len(all_results)
            fail_count = stats.errors
            pass_count = total - fail_count
            avg_score = getattr(stats, "avg_validation_score", 0) or 0
            email_stats = {
                "molecule_count": total,
                "pass_count": pass_count,
                "fail_count": fail_count,
                "avg_score": round(avg_score, 1),
            }
            send_batch_complete_email.delay(notification_email, job_id, email_stats)
            # Clean up the Redis key
            r.delete(f"batch:email:{job_id}")
    except Exception as e:
        logger.warning("Failed to dispatch email notification for %s: %s", job_id, e)

    # Trigger cheap analytics computation
    try:
        from app.services.batch.analytics_tasks import run_cheap_analytics

        run_cheap_analytics.delay(job_id)
    except Exception:
        logger.warning("Failed to dispatch cheap analytics for %s", job_id)

    return {
        "job_id": job_id,
        "total": len(all_results),
        "statistics": asdict(stats),
    }


# Priority queue versions of tasks for small jobs
@celery_app.task(
    bind=True,
    autoretry_for=(Exception,),
    retry_backoff=True,
    max_retries=3,
    acks_late=True,
    queue="high_priority",
)
def process_molecule_chunk_priority(
    self,
    job_id: str,
    chunk_idx: int,
    molecules: List[Dict[str, Any]],
    total_molecules: int,
    total_chunks: int,
    safety_options: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    """Priority queue version of process_molecule_chunk for small jobs."""
    results = []

    for mol_data in molecules:
        result = _process_single_molecule(mol_data, safety_options=safety_options)
        results.append(result)

        # Atomically increment counter and update progress
        processed_so_far = progress_tracker.increment_processed(job_id)
        progress_tracker.update_progress(
            job_id=job_id,
            processed=processed_so_far,
            total=total_molecules,
            status="processing",
        )

    return results


@celery_app.task(bind=True, queue="high_priority")
def aggregate_batch_results_priority(
    self,
    chunk_results: List[List[Dict[str, Any]]],
    job_id: str,
    start_time: float,
) -> Dict[str, Any]:
    """Priority queue version of aggregate_batch_results for small jobs."""
    all_results = [result for chunk in chunk_results for result in chunk]

    stats = compute_statistics(all_results)
    stats.processing_time_seconds = round(time.time() - start_time, 2)

    result_storage.store_results(job_id, all_results, stats)
    progress_tracker.mark_complete(job_id)

    # Log batch event to audit trail
    try:
        total = len(all_results)
        fail_count = stats.errors
        pass_count = total - fail_count
        asyncio.run(_log_batch_audit(
            job_id=job_id,
            molecule_count=total,
            pass_count=pass_count,
            fail_count=fail_count,
        ))
    except Exception as e:
        logger.warning("Failed to log batch audit event for %s: %s", job_id, e)

    # Dispatch webhook if configured
    try:
        webhook_url = settings.WEBHOOK_URL
        if webhook_url:
            from app.services.notifications.webhook import send_webhook

            total = len(all_results)
            payload = {
                "event": "batch_complete",
                "job_id": job_id,
                "status": "complete",
                "molecule_count": total,
                "summary_url": f"/batch/{job_id}",
            }
            send_webhook.delay(webhook_url, settings.WEBHOOK_SECRET, payload)
    except Exception as e:
        logger.warning("Failed to dispatch webhook for %s: %s", job_id, e)

    # Dispatch email notification if configured
    try:
        r = progress_tracker._get_redis()
        notification_email = r.get(f"batch:email:{job_id}")
        if notification_email:
            if isinstance(notification_email, bytes):
                notification_email = notification_email.decode("utf-8")
            from app.services.notifications.email import send_batch_complete_email

            total = len(all_results)
            fail_count = stats.errors
            pass_count = total - fail_count
            avg_score = getattr(stats, "avg_validation_score", 0) or 0
            email_stats = {
                "molecule_count": total,
                "pass_count": pass_count,
                "fail_count": fail_count,
                "avg_score": round(avg_score, 1),
            }
            send_batch_complete_email.delay(notification_email, job_id, email_stats)
            # Clean up the Redis key
            r.delete(f"batch:email:{job_id}")
    except Exception as e:
        logger.warning("Failed to dispatch email notification for %s: %s", job_id, e)

    # Trigger cheap analytics computation
    try:
        from app.services.batch.analytics_tasks import run_cheap_analytics

        run_cheap_analytics.delay(job_id)
    except Exception:
        logger.warning("Failed to dispatch cheap analytics for %s", job_id)

    return {
        "job_id": job_id,
        "total": len(all_results),
        "statistics": asdict(stats),
    }


def process_batch_job(
    job_id: str,
    molecules: List[Dict[str, Any]],
    safety_options: Optional[Dict[str, Any]] = None,
) -> str:
    """
    Start batch processing job using Celery chord pattern.

    Routes small jobs to high_priority queue for faster processing.
    Large jobs go to default queue to avoid blocking small jobs.

    Args:
        job_id: Unique job identifier
        molecules: List of molecule data dicts from file parsing
        safety_options: Optional safety screening options

    Returns:
        job_id for tracking
    """
    total_molecules = len(molecules)
    start_time = time.time()

    # Determine if this is a small job (use priority queue)
    is_small_job = total_molecules <= SMALL_JOB_THRESHOLD

    # Initialize progress tracking
    progress_tracker.init_job(job_id, total_molecules)

    # Split into chunks
    chunks = []
    for i in range(0, total_molecules, CHUNK_SIZE):
        chunk = molecules[i : i + CHUNK_SIZE]
        # Convert MoleculeData objects to dicts if needed
        chunk_dicts = []
        for m in chunk:
            if hasattr(m, "__dict__"):
                chunk_dicts.append(
                    {
                        "smiles": m.smiles,
                        "name": m.name,
                        "index": m.index,
                        "properties": m.properties if hasattr(m, "properties") else {},
                        "parse_error": (
                            m.parse_error if hasattr(m, "parse_error") else None
                        ),
                    }
                )
            else:
                chunk_dicts.append(m)
        chunks.append(chunk_dicts)

    total_chunks = len(chunks)

    # Select task based on job size
    if is_small_job:
        # Small jobs use priority queue
        chunk_task_func = process_molecule_chunk_priority
        aggregate_task_func = aggregate_batch_results_priority
    else:
        # Large jobs use default queue
        chunk_task_func = process_molecule_chunk
        aggregate_task_func = aggregate_batch_results

    # Create chord: all chunk tasks -> aggregate task
    chunk_tasks = group(
        chunk_task_func.s(
            job_id=job_id,
            chunk_idx=idx,
            molecules=chunk,
            total_molecules=total_molecules,
            total_chunks=total_chunks,
            safety_options=safety_options,
        )
        for idx, chunk in enumerate(chunks)
    )

    aggregate_task = aggregate_task_func.s(job_id=job_id, start_time=start_time)

    # Execute chord
    chord(chunk_tasks)(aggregate_task)

    return job_id


@celery_app.task
def cancel_batch_job(job_id: str) -> Dict[str, Any]:
    """
    Cancel a batch job (revoke pending tasks).

    Args:
        job_id: Job identifier

    Returns:
        Status dictionary
    """
    progress_tracker.mark_cancelled(job_id)
    return {"job_id": job_id, "status": "cancelled"}
