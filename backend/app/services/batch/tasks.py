"""
Celery Tasks for Batch Processing

Implements chunk-based molecule processing with progress tracking.
Uses chord pattern: process_molecule_chunk tasks -> aggregate_batch_results.
"""
import time
from typing import List, Dict, Any, Optional
from dataclasses import asdict

from celery import chord, group
from rdkit import Chem

from app.celery_app import celery_app
from app.services.validation.engine import validation_engine
from app.services.alerts import alert_manager
from app.services.scoring.ml_readiness import calculate_ml_readiness
from app.services.batch.progress_tracker import progress_tracker
from app.services.batch.result_aggregator import (
    compute_statistics,
    result_storage,
)


CHUNK_SIZE = 100  # Process 100 molecules per chunk for progress updates


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
) -> List[Dict[str, Any]]:
    """
    Process a chunk of molecules with progress updates.

    Args:
        job_id: Job identifier for progress tracking
        chunk_idx: Index of this chunk (for progress calculation)
        molecules: List of molecule data dicts (smiles, name, index, properties)
        total_molecules: Total molecules in the entire batch
        total_chunks: Total number of chunks

    Returns:
        List of result dictionaries for each molecule
    """
    results = []
    chunk_start = chunk_idx * CHUNK_SIZE

    for i, mol_data in enumerate(molecules):
        result = _process_single_molecule(mol_data)
        results.append(result)

        # Calculate overall progress and update
        processed_so_far = chunk_start + i + 1
        progress_tracker.update_progress(
            job_id=job_id,
            processed=processed_so_far,
            total=total_molecules,
            status="processing",
        )

    return results


def _process_single_molecule(mol_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process a single molecule through validation, alerts, and scoring.

    Args:
        mol_data: Dictionary with smiles, name, index, properties, and optional parse_error

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

        # Run validation
        try:
            check_results, validation_score = validation_engine.validate(mol)
            result["validation"] = {
                "overall_score": validation_score,
                "issues": [
                    {
                        "check_name": r.check_name,
                        "passed": r.passed,
                        "severity": r.severity.value if hasattr(r.severity, "value") else str(r.severity),
                        "message": r.message,
                    }
                    for r in check_results
                    if not r.passed
                ],
            }
        except Exception as e:
            result["validation"] = {"error": str(e)}

        # Run structural alerts screening
        try:
            alert_result = alert_manager.screen(mol, catalogs=["PAINS", "BRENK"])
            result["alerts"] = {
                "has_alerts": alert_result.has_alerts,
                "alert_count": len(alert_result.alerts),
                "alerts": [
                    {
                        "catalog": a.catalog,
                        "rule_name": a.rule_name,
                        "severity": a.severity.value if hasattr(a.severity, "value") else str(a.severity),
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

        # If we got this far with validation, mark as success
        result["status"] = "success"

    except Exception as e:
        result["error"] = f"Processing error: {str(e)}"

    return result


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

    return {
        "job_id": job_id,
        "total": len(all_results),
        "statistics": asdict(stats),
    }


def process_batch_job(job_id: str, molecules: List[Dict[str, Any]]) -> str:
    """
    Start batch processing job using Celery chord pattern.

    Args:
        job_id: Unique job identifier
        molecules: List of molecule data dicts from file parsing

    Returns:
        job_id for tracking
    """
    total_molecules = len(molecules)
    start_time = time.time()

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
                chunk_dicts.append({
                    "smiles": m.smiles,
                    "name": m.name,
                    "index": m.index,
                    "properties": m.properties if hasattr(m, "properties") else {},
                    "parse_error": m.parse_error if hasattr(m, "parse_error") else None,
                })
            else:
                chunk_dicts.append(m)
        chunks.append(chunk_dicts)

    total_chunks = len(chunks)

    # Create chord: all chunk tasks -> aggregate task
    chunk_tasks = group(
        process_molecule_chunk.s(
            job_id=job_id,
            chunk_idx=idx,
            molecules=chunk,
            total_molecules=total_molecules,
            total_chunks=total_chunks,
        )
        for idx, chunk in enumerate(chunks)
    )

    aggregate_task = aggregate_batch_results.s(job_id=job_id, start_time=start_time)

    # Execute chord
    workflow = chord(chunk_tasks)(aggregate_task)

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
