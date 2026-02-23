"""
Analytics Storage Module

Stores and retrieves analytics results in Redis with 24h TTL.
Used by analytics Celery tasks to persist computed results.
"""

import json
from typing import Optional

import redis

from app.core.config import settings


class AnalyticsStorage:
    """
    Stores and retrieves analytics results in Redis.

    Key format:
    - Results: batch:analytics:{analysis_type}:{job_id}
    - Status:  batch:analytics:status:{job_id}
    """

    ANALYTICS_TTL = 86400  # 24 hours

    def __init__(self, redis_url: str = None):
        self._redis_url = redis_url or settings.REDIS_URL
        self._redis: Optional[redis.Redis] = None

    def _get_redis(self) -> redis.Redis:
        """Get or create Redis connection."""
        if self._redis is None:
            self._redis = redis.from_url(self._redis_url)
        return self._redis

    def store_result(
        self,
        job_id: str,
        analysis_type: str,
        data: dict,
        ttl: int = None,
    ) -> None:
        """
        Store analytics result for a job and analysis type.

        Args:
            job_id: Batch job identifier.
            analysis_type: Type of analytics (e.g. "deduplication", "scaffold").
            data: Result dictionary to serialize and store.
            ttl: Time-to-live in seconds. Defaults to ANALYTICS_TTL (24h).
        """
        r = self._get_redis()
        key = f"batch:analytics:{analysis_type}:{job_id}"
        r.set(key, json.dumps(data), ex=ttl if ttl is not None else self.ANALYTICS_TTL)

    def get_result(self, job_id: str, analysis_type: str) -> Optional[dict]:
        """
        Retrieve analytics result for a job and analysis type.

        Args:
            job_id: Batch job identifier.
            analysis_type: Type of analytics to retrieve.

        Returns:
            Deserialized result dict, or None if not found / expired.
        """
        r = self._get_redis()
        key = f"batch:analytics:{analysis_type}:{job_id}"
        data = r.get(key)
        if data is None:
            return None
        return json.loads(data)

    def get_status(self, job_id: str) -> Optional[dict]:
        """
        Retrieve analytics status dict for a job.

        Returns:
            Status dict keyed by analysis type, or None if not found.
        """
        r = self._get_redis()
        data = r.get(f"batch:analytics:status:{job_id}")
        if data is None:
            return None
        return json.loads(data)

    def update_status(
        self,
        job_id: str,
        analysis_type: str,
        status: str,
        error: Optional[str] = None,
    ) -> None:
        """
        Update the status of a single analysis type within the job's status dict.

        Reads current status dict, updates the sub-key for analysis_type, writes back.

        Args:
            job_id: Batch job identifier.
            analysis_type: The analysis type whose status to update.
            status: New status value ("pending", "computing", "complete", "failed", "skipped").
            error: Optional error message if status is "failed".
        """
        r = self._get_redis()
        status_key = f"batch:analytics:status:{job_id}"
        existing = r.get(status_key)
        status_dict = json.loads(existing) if existing else {}
        entry = status_dict.setdefault(analysis_type, {"status": "pending"})
        entry["status"] = status
        if error is not None:
            entry["error"] = error
        r.set(status_key, json.dumps(status_dict), ex=self.ANALYTICS_TTL)

    def init_status(self, job_id: str, auto_analyses: list[str]) -> None:
        """
        Initialize the analytics status dict for a job.

        Auto analyses start as "computing"; all others listed here default to "pending".

        Args:
            job_id: Batch job identifier.
            auto_analyses: List of analysis types that will be computed automatically
                           (set to "computing" immediately).
        """
        all_types = [
            "deduplication",
            "statistics",
            "scaffold",
            "chemical_space",
            "mmp",
            "similarity_search",
            "rgroup",
        ]
        status_dict = {}
        for analysis_type in all_types:
            if analysis_type in auto_analyses:
                status_dict[analysis_type] = {"status": "computing"}
            else:
                status_dict[analysis_type] = {"status": "pending"}

        r = self._get_redis()
        r.set(
            f"batch:analytics:status:{job_id}",
            json.dumps(status_dict),
            ex=self.ANALYTICS_TTL,
        )


# Singleton instance
analytics_storage = AnalyticsStorage()
