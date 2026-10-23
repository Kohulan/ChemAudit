"""ProgressTracker's lazy Redis init must be thread-safe."""

import threading

from app.services.batch.progress_tracker import ProgressTracker, progress_tracker


def test_tracker_has_lock():
    assert hasattr(progress_tracker, "_lock")
    assert hasattr(progress_tracker._lock, "acquire")
    assert hasattr(progress_tracker._lock, "release")


def test_get_redis_returns_single_instance_under_concurrency():
    """Concurrent _get_redis() calls must all receive the same pooled client."""
    tracker = ProgressTracker()
    results = []
    barrier = threading.Barrier(16)

    def worker():
        barrier.wait()  # release all threads simultaneously to maximise the race
        results.append(tracker._get_redis())

    threads = [threading.Thread(target=worker) for _ in range(16)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert len(results) == 16
    assert len({id(r) for r in results}) == 1


def test_mark_complete_reports_all_molecules_processed():
    """Throttled progress writes can lag; a complete job must report processed == total."""
    import json
    from unittest.mock import MagicMock

    tracker = ProgressTracker()
    redis = MagicMock()
    redis.get.return_value = json.dumps(
        {"job_id": "j1", "status": "processing", "progress": 33, "processed": 1, "total": 3}
    )
    tracker._get_redis = lambda: redis

    tracker.mark_complete("j1")

    stored = json.loads(redis.set.call_args.args[1])
    assert stored["status"] == "complete"
    assert stored["processed"] == 3
