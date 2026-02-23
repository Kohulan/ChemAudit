"""
Analytics Celery Tasks

Implements cheap (auto-dispatched) and expensive (user-triggered) analytics tasks
for batch post-processing. Both tasks are defensive: they check for expired results
and gracefully handle missing analytics service modules.
"""

import logging

from app.celery_app import celery_app
from app.services.analytics.storage import analytics_storage
from app.services.batch.result_aggregator import result_storage

logger = logging.getLogger(__name__)


@celery_app.task(bind=True, queue="default", ignore_result=True)
def run_cheap_analytics(self, job_id: str) -> None:
    """
    Run cheap analytics automatically after batch aggregation completes.

    Computes deduplication and property statistics — analyses that are always
    valuable and have acceptable runtime for any batch size.

    Dispatched by aggregate_batch_results / aggregate_batch_results_priority.

    Args:
        job_id: Batch job identifier.
    """
    results = result_storage.get_all_results(job_id)
    if not results:
        logger.warning(
            "run_cheap_analytics: no results for job %s — "
            "batch results expired; resubmit batch.",
            job_id,
        )
        analytics_storage.update_status(
            job_id,
            "deduplication",
            "failed",
            error="Batch results expired — resubmit batch",
        )
        analytics_storage.update_status(
            job_id,
            "statistics",
            "failed",
            error="Batch results expired — resubmit batch",
        )
        return

    analytics_storage.init_status(
        job_id, auto_analyses=["deduplication", "statistics"]
    )

    # --- Deduplication ---
    dedup_available = False
    try:
        from app.services.analytics.deduplication import compute_all_dedup_levels  # noqa: F401

        dedup_available = True
    except ImportError:
        logger.warning(
            "run_cheap_analytics: deduplication module not yet available for job %s",
            job_id,
        )
        analytics_storage.update_status(job_id, "deduplication", "skipped")

    if dedup_available:
        try:
            from app.services.analytics.deduplication import compute_all_dedup_levels

            dedup_result = compute_all_dedup_levels(results)
            analytics_storage.store_result(
                job_id, "deduplication", dedup_result.model_dump()
            )
            analytics_storage.update_status(job_id, "deduplication", "complete")
        except Exception as exc:
            logger.exception(
                "run_cheap_analytics: deduplication failed for job %s", job_id
            )
            analytics_storage.update_status(
                job_id, "deduplication", "failed", error=str(exc)
            )

    # --- Statistics ---
    stats_available = False
    try:
        from app.services.analytics.statistics import compute_all_statistics  # noqa: F401

        stats_available = True
    except ImportError:
        logger.warning(
            "run_cheap_analytics: statistics module not yet available for job %s",
            job_id,
        )
        analytics_storage.update_status(job_id, "statistics", "skipped")

    if stats_available:
        try:
            from app.services.analytics.statistics import compute_all_statistics

            stats_result = compute_all_statistics(results)
            analytics_storage.store_result(
                job_id, "statistics", stats_result.model_dump()
            )
            analytics_storage.update_status(job_id, "statistics", "complete")
        except Exception as exc:
            logger.exception(
                "run_cheap_analytics: statistics failed for job %s", job_id
            )
            analytics_storage.update_status(
                job_id, "statistics", "failed", error=str(exc)
            )


@celery_app.task(bind=True, queue="default", ignore_result=True)
def run_expensive_analytics(
    self,
    job_id: str,
    analysis_type: str,
    params: dict | None = None,
) -> None:
    """
    Run an expensive analytics computation triggered by the user.

    Dispatched via POST /batch/{job_id}/analytics/{analysis_type}.

    Supported analysis types:
    - scaffold: Murcko scaffold diversity analysis.
    - chemical_space: PCA or t-SNE 2-D embedding.
    - mmp: Matched molecular pair analysis.
    - similarity_search: Nearest-neighbor similarity search.
    - rgroup: R-group decomposition.

    Args:
        job_id: Batch job identifier.
        analysis_type: The type of analytics to compute.
        params: Optional parameters for the analysis (e.g. method, activity_column).
    """
    if params is None:
        params = {}

    results = result_storage.get_all_results(job_id)
    if not results:
        logger.warning(
            "run_expensive_analytics: no results for job %s (%s) — batch results expired.",
            job_id,
            analysis_type,
        )
        analytics_storage.update_status(
            job_id,
            analysis_type,
            "failed",
            error="Batch results expired — resubmit batch",
        )
        return

    analytics_storage.update_status(job_id, analysis_type, "computing")

    try:
        if analysis_type == "scaffold":
            from app.services.analytics.scaffold_analysis import compute_scaffold_analysis

            result = compute_scaffold_analysis(results)
            analytics_storage.store_result(job_id, "scaffold", result)

        elif analysis_type == "chemical_space":
            from app.services.analytics.chemical_space import compute_chemical_space

            method = params.get("method", "pca")
            result = compute_chemical_space(results, method=method)
            data = result.model_dump() if hasattr(result, "model_dump") else result
            analytics_storage.store_result(job_id, "chemical_space", data)

        elif analysis_type == "mmp":
            from app.services.analytics.mmp import compute_mmp

            activity_column = params.get("activity_column")
            result = compute_mmp(results, activity_column=activity_column)
            analytics_storage.store_result(job_id, "mmp", result.model_dump())

        elif analysis_type == "similarity_search":
            from app.services.analytics.chemical_space import find_similar_molecules

            query_smiles = params.get("query_smiles")
            query_index = params.get("query_index")
            top_k = params.get("top_k", 10)
            result = find_similar_molecules(
                results,
                query_smiles=query_smiles,
                query_index=query_index,
                top_k=top_k,
            )
            analytics_storage.store_result(
                job_id, "similarity_search", result
            )

        elif analysis_type == "rgroup":
            from app.services.analytics.scaffold_analysis import compute_rgroup_decomposition

            core_smarts = params.get("core_smarts")
            result = compute_rgroup_decomposition(results, core_smarts=core_smarts)
            analytics_storage.store_result(job_id, "rgroup", result)

        else:
            analytics_storage.update_status(
                job_id, analysis_type, "failed", error="Unknown analysis type"
            )
            return

        analytics_storage.update_status(job_id, analysis_type, "complete")

    except Exception as exc:
        logger.exception(
            "run_expensive_analytics: %s failed for job %s", analysis_type, job_id
        )
        analytics_storage.update_status(
            job_id, analysis_type, "failed", error=str(exc)
        )
