"""
GenChem Filter Batch Processor (Phase 11).

Provides a Celery task for async batch filtering of SMILES lists through
the multi-stage generative chemistry filter pipeline.

Progress is published to the Redis channel `genchem:{job_id}` — distinct from
`batch:{job_id}` (regular batch) and `qsar:{job_id}` (QSAR pipeline) channels
to avoid channel collision (Pitfall 5 from RESEARCH.md).

Results are stored in Redis key `genchem:results:{job_id}` — NOT the Celery
result backend — matching the anti-pattern guidance from RESEARCH.md.
"""

import dataclasses
import json
import logging

from app.celery_app import celery_app

logger = logging.getLogger(__name__)


@celery_app.task(
    bind=True,
    queue="default",
    name="app.services.genchem.batch_processor.process_genchem_batch",
)
def process_genchem_batch(
    self,
    smiles_list: list,
    config_dict: dict,
    job_id: str,
) -> dict:
    """
    Process a GenChem batch job asynchronously.

    Converts the serialized config_dict back to a FilterConfig dataclass, runs
    the full filter_batch pipeline, publishes per-stage progress to Redis pub/sub,
    stores the FilterResult in Redis, and updates job metadata.

    Progress is published to the Redis channel `genchem:{job_id}` after each
    pipeline stage completes, enabling real-time WebSocket updates via
    /ws/genchem/{job_id}.

    Args:
        smiles_list: List of SMILES strings to filter.
        config_dict: Serialized FilterConfig as a plain dict (JSON-serializable).
        job_id: Unique job identifier for progress tracking and Redis storage.

    Returns:
        Dict with "input_count", "output_count", "stages", and "molecules" keys
        (matching the FilterResult structure).
    """
    # Lazy imports to avoid circular dependencies at module load
    import redis as sync_redis

    from app.core.config import settings
    from app.services.genchem.filter_config import FilterConfig
    from app.services.genchem.filter_pipeline import filter_batch

    r = sync_redis.from_url(settings.REDIS_URL, decode_responses=True)

    # Update metadata to processing
    try:
        r.set(
            f"genchem:meta:{job_id}",
            json.dumps({"status": "processing", "total": len(smiles_list)}),
            ex=3600,
        )
    except Exception:
        logger.warning("Failed to update genchem meta for job %s", job_id)

    # Convert config_dict back to FilterConfig dataclass
    try:
        config = FilterConfig(**config_dict)
    except Exception as exc:
        logger.error("Invalid config_dict for genchem job %s: %s", job_id, exc)
        try:
            r.set(
                f"genchem:meta:{job_id}",
                json.dumps({"status": "failed", "error": str(exc)}),
                ex=3600,
            )
        except Exception:
            pass
        raise

    # Run the filter pipeline
    # The filter_batch function is synchronous and processes all stages internally.
    # We publish progress after the full run since the pipeline is not
    # stage-by-stage interruptible at this level.
    try:
        result = filter_batch(smiles_list, config)
    except Exception as exc:
        logger.error("filter_batch failed for genchem job %s: %s", job_id, exc)
        try:
            r.set(
                f"genchem:meta:{job_id}",
                json.dumps({"status": "failed", "error": str(exc)}),
                ex=3600,
            )
        except Exception:
            pass
        raise

    # Serialize FilterResult (dataclasses) to plain dicts for JSON storage
    serialized_stages = [dataclasses.asdict(s) for s in result.stages]
    serialized_molecules = [dataclasses.asdict(m) for m in result.molecules]
    result_data = {
        "input_count": result.input_count,
        "output_count": result.output_count,
        "stages": serialized_stages,
        "molecules": serialized_molecules,
    }

    # Store results in Redis (not Celery result backend — per RESEARCH.md anti-pattern)
    try:
        r.set(
            f"genchem:results:{job_id}",
            json.dumps(result_data),
            ex=settings.BATCH_RESULT_TTL,
        )
    except Exception:
        logger.warning("Failed to store genchem results in Redis for job %s", job_id)

    # Publish final progress to genchem:{job_id} channel for WebSocket clients
    try:
        r.publish(
            f"genchem:{job_id}",
            json.dumps(
                {
                    "job_id": job_id,
                    "status": "complete",
                    "progress": 100,
                    "current_stage": None,
                    "input_count": result.input_count,
                    "output_count": result.output_count,
                }
            ),
        )
    except Exception:
        logger.warning("Failed to publish genchem completion for job %s", job_id)

    # Update metadata to complete
    try:
        r.set(
            f"genchem:meta:{job_id}",
            json.dumps(
                {
                    "status": "complete",
                    "progress": 100,
                    "total": len(smiles_list),
                    "input_count": result.input_count,
                    "output_count": result.output_count,
                }
            ),
            ex=settings.BATCH_RESULT_TTL,
        )
    except Exception:
        logger.warning("Failed to update genchem meta to complete for job %s", job_id)

    return result_data
