"""
Tests for analytics infrastructure (Phase 03-01 INFRA-01).

Covers:
- BATCH_RESULT_TTL config setting
- ResultStorage.get_all_results
- AnalyticsStorage CRUD operations
- Analytics Pydantic schemas
- GET /batch/{job_id}/analytics endpoint (404 path)
"""

import json
from unittest.mock import MagicMock, patch

import pytest
from fastapi.testclient import TestClient

from app.core.config import settings
from app.schemas.analytics import (
    AnalysisStatus,
    AnalyticsTriggerResponse,
    BatchAnalyticsResponse,
    DeduplicationGroup,
    DeduplicationResult,
    PropertyStats,
    QualityScore,
    ScaffoldGroup,
    ScaffoldResult,
    StatisticsResult,
)
from app.services.analytics.storage import AnalyticsStorage
from app.services.batch.result_aggregator import ResultStorage


# ---------------------------------------------------------------------------
# Config settings
# ---------------------------------------------------------------------------


def test_batch_result_ttl_is_24h():
    """BATCH_RESULT_TTL must equal 86400 (24 hours in seconds)."""
    assert settings.BATCH_RESULT_TTL == 86400


# ---------------------------------------------------------------------------
# ResultStorage.get_all_results
# ---------------------------------------------------------------------------


def test_get_all_results_returns_empty_for_missing_job():
    """get_all_results returns [] when there are no results in Redis."""
    storage = ResultStorage()
    mock_redis = MagicMock()
    mock_redis.get.return_value = None
    storage._redis = mock_redis

    result = storage.get_all_results("nonexistent-job-id")

    assert result == []
    mock_redis.get.assert_called_once_with("batch:results:nonexistent-job-id")


def test_get_all_results_deserializes_stored_data():
    """get_all_results correctly deserializes stored JSON list."""
    storage = ResultStorage()
    raw = [{"smiles": "CCO", "status": "success"}, {"smiles": "c1ccccc1", "status": "success"}]
    mock_redis = MagicMock()
    mock_redis.get.return_value = json.dumps(raw).encode()
    storage._redis = mock_redis

    result = storage.get_all_results("test-job-123")

    assert result == raw
    assert len(result) == 2


# ---------------------------------------------------------------------------
# AnalyticsStorage
# ---------------------------------------------------------------------------


def _make_storage_with_mock():
    """Return (AnalyticsStorage, mock_redis) for unit tests."""
    storage = AnalyticsStorage()
    mock_redis = MagicMock()
    storage._redis = mock_redis
    return storage, mock_redis


def test_analytics_storage_store_result():
    """store_result serializes data and calls Redis SET with TTL."""
    storage, mock_redis = _make_storage_with_mock()
    data = {"groups": [], "total_unique": {"exact": 5}}

    storage.store_result("job-1", "deduplication", data)

    mock_redis.set.assert_called_once()
    call_args = mock_redis.set.call_args
    key = call_args[0][0] if call_args[0] else call_args[1].get("name", call_args[0][0])
    # Verify key and TTL
    args, kwargs = call_args
    assert "batch:analytics:deduplication:job-1" in args[0]
    assert json.loads(args[1]) == data
    assert kwargs.get("ex") == AnalyticsStorage.ANALYTICS_TTL


def test_analytics_storage_get_result_returns_none_when_missing():
    """get_result returns None when key does not exist in Redis."""
    storage, mock_redis = _make_storage_with_mock()
    mock_redis.get.return_value = None

    result = storage.get_result("job-1", "deduplication")

    assert result is None


def test_analytics_storage_round_trip():
    """store_result followed by get_result returns the original dict."""
    # Use an in-memory dict to simulate Redis
    store: dict = {}

    def fake_set(key, value, ex=None):
        store[key] = value

    def fake_get(key):
        return store.get(key)

    storage = AnalyticsStorage()
    mock_redis = MagicMock()
    mock_redis.set.side_effect = fake_set
    mock_redis.get.side_effect = fake_get
    storage._redis = mock_redis

    data = {"scaffold_smiles": "c1ccccc1", "count": 10}
    storage.store_result("job-99", "scaffold", data)
    result = storage.get_result("job-99", "scaffold")

    assert result == data


def test_analytics_storage_init_status():
    """init_status writes correct status dict with computing/pending states."""
    store: dict = {}

    def fake_set(key, value, ex=None):
        store[key] = value

    storage = AnalyticsStorage()
    mock_redis = MagicMock()
    mock_redis.set.side_effect = fake_set
    storage._redis = mock_redis

    storage.init_status("job-2", auto_analyses=["deduplication", "statistics"])

    status_key = "batch:analytics:status:job-2"
    assert status_key in store
    status_dict = json.loads(store[status_key])
    assert status_dict["deduplication"]["status"] == "computing"
    assert status_dict["statistics"]["status"] == "computing"
    assert status_dict["scaffold"]["status"] == "pending"
    assert status_dict["chemical_space"]["status"] == "pending"


def test_analytics_storage_update_status():
    """update_status correctly updates an entry in the status dict."""
    store: dict = {}

    def fake_set(key, value, ex=None):
        store[key] = value

    def fake_get(key):
        return store.get(key)

    storage = AnalyticsStorage()
    mock_redis = MagicMock()
    mock_redis.set.side_effect = fake_set
    mock_redis.get.side_effect = fake_get
    storage._redis = mock_redis

    # Initialize
    storage.init_status("job-3", auto_analyses=["deduplication"])
    # Update to complete
    storage.update_status("job-3", "deduplication", "complete")

    status_dict = json.loads(store["batch:analytics:status:job-3"])
    assert status_dict["deduplication"]["status"] == "complete"


def test_analytics_storage_update_status_with_error():
    """update_status stores error message when status is 'failed'."""
    store: dict = {}

    def fake_set(key, value, ex=None):
        store[key] = value

    def fake_get(key):
        return store.get(key)

    storage = AnalyticsStorage()
    mock_redis = MagicMock()
    mock_redis.set.side_effect = fake_set
    mock_redis.get.side_effect = fake_get
    storage._redis = mock_redis

    storage.update_status("job-4", "scaffold", "failed", error="Something went wrong")

    status_dict = json.loads(store["batch:analytics:status:job-4"])
    assert status_dict["scaffold"]["status"] == "failed"
    assert status_dict["scaffold"]["error"] == "Something went wrong"


# ---------------------------------------------------------------------------
# Schema instantiation
# ---------------------------------------------------------------------------


def test_analysis_status_schema():
    """AnalysisStatus validates all literal status values."""
    for s in ("pending", "computing", "complete", "failed", "skipped"):
        obj = AnalysisStatus(status=s)
        assert obj.status == s


def test_deduplication_group_schema():
    """DeduplicationGroup instantiates correctly."""
    group = DeduplicationGroup(
        level="exact",
        representative_index=0,
        duplicate_indices=[1, 2],
        group_key="InChI=1/CH4/h1H4",
        count=3,
    )
    assert group.count == 3
    assert group.level == "exact"


def test_batch_analytics_response_schema_minimal():
    """BatchAnalyticsResponse validates with only required fields."""
    resp = BatchAnalyticsResponse(
        job_id="test-job",
        status={"deduplication": AnalysisStatus(status="pending")},
    )
    assert resp.job_id == "test-job"
    assert resp.deduplication is None
    assert resp.scaffold is None


def test_analytics_trigger_response_schema():
    """AnalyticsTriggerResponse instantiates correctly."""
    resp = AnalyticsTriggerResponse(
        job_id="test-job",
        analysis_type="scaffold",
        status="queued",
    )
    assert resp.status == "queued"


# ---------------------------------------------------------------------------
# Analytics GET endpoint — 404 path
# ---------------------------------------------------------------------------


def test_get_batch_analytics_returns_404_when_no_analytics():
    """GET /api/v1/batch/{job_id}/analytics returns 404 when analytics haven't been computed."""
    import os

    os.environ.setdefault("RATE_LIMIT_ENABLED", "false")

    from app.main import app

    with patch(
        "app.api.routes.batch.analytics_storage.get_status", return_value=None
    ):
        client = TestClient(app)
        response = client.get("/api/v1/batch/nonexistent-job/analytics")

    assert response.status_code == 404
    assert "not available" in response.json()["detail"].lower()
