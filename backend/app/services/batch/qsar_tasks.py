"""
Celery Tasks for QSAR-Ready Batch Processing

Provides the Celery task for async QSAR-Ready pipeline batch processing.
Registered on the `default` queue per D-11 from CONTEXT.md.

Progress is published to the `qsar:{job_id}` Redis pub/sub channel — distinct
from regular batch processing's `batch:{job_id}` channel (Pitfall 5 from RESEARCH.md).
"""

import logging

from app.celery_app import celery_app

logger = logging.getLogger(__name__)


@celery_app.task(
    bind=True,
    queue="default",
    name="app.services.batch.qsar_tasks.process_qsar_batch_job",
)
def process_qsar_batch_job(
    self,
    job_id: str,
    molecules: list,
    config_dict: dict,
) -> dict:
    """
    Process a QSAR batch job asynchronously.

    Accepts a serialized config dict and a list of SMILES strings, runs
    them through the QSAR-Ready pipeline with deduplication and progress
    tracking, and returns the serialized results.

    Publishes progress to qsar:{job_id} Redis channel for real-time WebSocket
    forwarding via /ws/qsar/{job_id}.

    Uses lazy imports inside the task body to avoid circular imports
    (same pattern as existing batch tasks in app.services.batch.tasks).

    Args:
        job_id: Unique job identifier for progress tracking and Redis storage.
        molecules: List of SMILES strings parsed from the uploaded file.
        config_dict: Serialized QSARReadyConfig as a plain dict.

    Returns:
        Dict with "results" (list of serialized QSARReadyResult dicts),
        "summary" (ok/rejected/duplicate/error counts + steps_applied_counts),
        and "config" (serialized config dict).
    """
    from app.services.qsar_ready.batch_processor import run_qsar_batch
    from app.services.qsar_ready.pipeline import QSARReadyConfig

    config = QSARReadyConfig(**config_dict)
    return run_qsar_batch(job_id, molecules, config)
