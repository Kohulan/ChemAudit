"""ResultStorage chunked (Redis list) storage + paginated reads.

Batch results are stored as a Redis list (one JSON-encoded result per element)
instead of a single JSON blob, so the default paginated view can be served with
an LRANGE range query without deserializing the entire dataset.
"""

import json

import fakeredis
import pytest

from app.services.batch.result_aggregator import (
    BatchStatisticsData,
    ResultStorage,
    compute_statistics,
)


def _results(n: int, *, statuses=None):
    out = []
    for i in range(n):
        status = "success"
        if statuses and i in statuses:
            status = "error"
        out.append(
            {
                "index": i,
                "smiles": f"C{'C' * (i % 5)}O",
                "status": status,
                "error": "bad" if status == "error" else None,
                "validation": {"overall_score": (i * 7) % 101, "issues": []},
            }
        )
    return out


@pytest.fixture
def storage():
    s = ResultStorage()
    s._redis = fakeredis.FakeRedis()  # sync, byte responses (matches real client)
    return s


def _store(storage, results):
    stats = compute_statistics(results)
    storage.store_results("job1", results, stats)


def test_store_results_uses_redis_list(storage):
    _store(storage, _results(3))
    r = storage._redis
    assert r.type("batch:results:job1") == b"list"
    assert r.llen("batch:results:job1") == 3


def test_get_all_results_round_trip(storage):
    data = _results(5)
    _store(storage, data)
    assert storage.get_all_results("job1") == data


def test_get_all_results_empty_for_missing(storage):
    assert storage.get_all_results("missing") == []


def test_default_view_paginates(storage):
    _store(storage, _results(120))
    page = storage.get_results("job1", page=2, page_size=50)
    assert page["total_results"] == 120
    assert page["total_pages"] == 3
    assert len(page["results"]) == 50
    # Page 2 of stored order = indices 50..99
    assert [r["index"] for r in page["results"]] == list(range(50, 100))


def test_filtered_view_returns_only_matches(storage):
    _store(storage, _results(10, statuses={2, 5, 7}))
    page = storage.get_results("job1", page=1, page_size=50, status_filter="error")
    assert page["total_results"] == 3
    assert {r["index"] for r in page["results"]} == {2, 5, 7}


def test_sorted_view_orders_results(storage):
    _store(storage, _results(10))
    page = storage.get_results(
        "job1", page=1, page_size=50, sort_by="score", sort_dir="desc"
    )
    scores = [r["validation"]["overall_score"] for r in page["results"]]
    assert scores == sorted(scores, reverse=True)


def test_backward_compat_legacy_string_blob(storage):
    """A pre-existing single-JSON-blob value must still be readable."""
    data = _results(4)
    storage._redis.set("batch:results:job1", json.dumps(data))  # legacy format
    assert storage.get_all_results("job1") == data
    page = storage.get_results("job1", page=1, page_size=50)
    assert page["total_results"] == 4


def test_store_overwrites_previous(storage):
    _store(storage, _results(10))
    _store(storage, _results(3))
    assert storage._redis.llen("batch:results:job1") == 3
    assert storage.get_results("job1", page=1, page_size=50)["total_results"] == 3


def test_empty_results(storage):
    _store(storage, [])
    assert storage.get_all_results("job1") == []
    assert storage.get_results("job1", page=1, page_size=50)["total_results"] == 0


def test_statistics_round_trip(storage):
    _store(storage, _results(6, statuses={1}))
    stats = storage.get_statistics("job1")
    assert isinstance(stats, BatchStatisticsData)
    assert stats.total == 6
    assert stats.errors == 1
