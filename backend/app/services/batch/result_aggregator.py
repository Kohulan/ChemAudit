"""
Result Aggregator Module

Computes statistics from batch processing results.
"""
import json
from dataclasses import dataclass, asdict, field
from typing import List, Dict, Any, Optional
from collections import Counter

import redis

from app.core.config import settings


@dataclass
class BatchStatisticsData:
    """Statistics computed from batch processing results."""

    total: int = 0
    successful: int = 0
    errors: int = 0
    avg_validation_score: Optional[float] = None
    avg_ml_readiness_score: Optional[float] = None
    score_distribution: Dict[str, int] = field(default_factory=dict)
    alert_summary: Dict[str, int] = field(default_factory=dict)
    processing_time_seconds: Optional[float] = None


def compute_statistics(results: List[Dict[str, Any]]) -> BatchStatisticsData:
    """
    Compute aggregate statistics from batch results.

    Args:
        results: List of result dictionaries from molecule processing

    Returns:
        BatchStatisticsData with computed statistics
    """
    stats = BatchStatisticsData()
    stats.total = len(results)

    validation_scores = []
    ml_readiness_scores = []
    alert_counts: Counter = Counter()

    for result in results:
        if result.get("status") == "error" or result.get("error"):
            stats.errors += 1
        else:
            stats.successful += 1

            # Collect validation scores
            if "validation" in result and result["validation"]:
                score = result["validation"].get("overall_score")
                if score is not None:
                    validation_scores.append(score)

            # Collect ML-readiness scores
            if "scoring" in result and result["scoring"]:
                ml_score = result["scoring"].get("ml_readiness", {}).get("score")
                if ml_score is not None:
                    ml_readiness_scores.append(ml_score)

            # Count alerts
            if "alerts" in result and result["alerts"]:
                for alert in result["alerts"].get("alerts", []):
                    alert_type = alert.get("catalog", "Unknown")
                    alert_counts[alert_type] += 1

    # Calculate averages
    if validation_scores:
        stats.avg_validation_score = round(sum(validation_scores) / len(validation_scores), 1)

    if ml_readiness_scores:
        stats.avg_ml_readiness_score = round(
            sum(ml_readiness_scores) / len(ml_readiness_scores), 1
        )

    # Score distribution buckets
    stats.score_distribution = _compute_score_distribution(validation_scores)

    # Alert summary
    stats.alert_summary = dict(alert_counts)

    return stats


def _compute_score_distribution(scores: List[int]) -> Dict[str, int]:
    """
    Compute histogram buckets for validation scores.

    Buckets:
    - excellent: 90-100
    - good: 70-89
    - moderate: 50-69
    - poor: 0-49
    """
    distribution = {"excellent": 0, "good": 0, "moderate": 0, "poor": 0}

    for score in scores:
        if score >= 90:
            distribution["excellent"] += 1
        elif score >= 70:
            distribution["good"] += 1
        elif score >= 50:
            distribution["moderate"] += 1
        else:
            distribution["poor"] += 1

    return distribution


class ResultStorage:
    """
    Stores and retrieves batch results in Redis.

    Results are stored with pagination support and expiration.
    """

    RESULT_EXPIRY = 3600  # 1 hour
    PAGE_SIZE = 50  # Default page size

    def __init__(self, redis_url: str = None):
        self._redis_url = redis_url or settings.REDIS_URL
        self._redis: Optional[redis.Redis] = None

    def _get_redis(self) -> redis.Redis:
        """Get or create Redis connection."""
        if self._redis is None:
            self._redis = redis.from_url(self._redis_url)
        return self._redis

    def store_results(
        self,
        job_id: str,
        results: List[Dict[str, Any]],
        statistics: BatchStatisticsData,
    ) -> None:
        """
        Store batch results and statistics.

        Args:
            job_id: Job identifier
            results: List of result dictionaries
            statistics: Computed statistics
        """
        r = self._get_redis()

        # Store results as a JSON list (for smaller batches, this is fine)
        # For very large batches (10K+), consider chunked storage
        r.set(
            f"batch:results:{job_id}",
            json.dumps(results),
            ex=self.RESULT_EXPIRY,
        )

        # Store statistics separately
        r.set(
            f"batch:stats:{job_id}",
            json.dumps(asdict(statistics)),
            ex=self.RESULT_EXPIRY,
        )

    def get_results(
        self,
        job_id: str,
        page: int = 1,
        page_size: int = 50,
        status_filter: Optional[str] = None,
        min_score: Optional[int] = None,
        max_score: Optional[int] = None,
    ) -> Dict[str, Any]:
        """
        Get paginated results for a job with optional filtering.

        Args:
            job_id: Job identifier
            page: Page number (1-indexed)
            page_size: Results per page
            status_filter: Filter by status ('success', 'error')
            min_score: Minimum validation score
            max_score: Maximum validation score

        Returns:
            Dictionary with results, pagination info, and statistics
        """
        r = self._get_redis()

        # Get all results
        results_data = r.get(f"batch:results:{job_id}")
        if not results_data:
            return {
                "results": [],
                "page": page,
                "page_size": page_size,
                "total_results": 0,
                "total_pages": 0,
            }

        results = json.loads(results_data)

        # Apply filters
        filtered = self._apply_filters(results, status_filter, min_score, max_score)

        # Paginate
        total_results = len(filtered)
        total_pages = (total_results + page_size - 1) // page_size if total_results > 0 else 0
        start_idx = (page - 1) * page_size
        end_idx = start_idx + page_size
        page_results = filtered[start_idx:end_idx]

        return {
            "results": page_results,
            "page": page,
            "page_size": page_size,
            "total_results": total_results,
            "total_pages": total_pages,
        }

    def get_statistics(self, job_id: str) -> Optional[BatchStatisticsData]:
        """Get statistics for a job."""
        r = self._get_redis()
        data = r.get(f"batch:stats:{job_id}")
        if data:
            return BatchStatisticsData(**json.loads(data))
        return None

    def _apply_filters(
        self,
        results: List[Dict[str, Any]],
        status_filter: Optional[str],
        min_score: Optional[int],
        max_score: Optional[int],
    ) -> List[Dict[str, Any]]:
        """Apply filters to results list."""
        filtered = results

        if status_filter:
            if status_filter == "success":
                filtered = [r for r in filtered if r.get("status") == "success"]
            elif status_filter == "error":
                filtered = [r for r in filtered if r.get("status") == "error" or r.get("error")]

        if min_score is not None:
            filtered = [
                r
                for r in filtered
                if r.get("validation", {}).get("overall_score", 0) >= min_score
            ]

        if max_score is not None:
            filtered = [
                r
                for r in filtered
                if r.get("validation", {}).get("overall_score", 100) <= max_score
            ]

        return filtered

    def delete_results(self, job_id: str) -> None:
        """Delete stored results for a job."""
        r = self._get_redis()
        r.delete(f"batch:results:{job_id}")
        r.delete(f"batch:stats:{job_id}")


# Singleton instance
result_storage = ResultStorage()
