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
