"""End-to-end Celery chord aggregation test (eager mode).

Exercises the real group(chunk_tasks) -> aggregate chord wiring so serialization
or aggregation-callback bugs are caught in CI instead of only in production.
DB-audit and analytics dispatch (separate infra with their own tests) are stubbed.
"""

import fakeredis
import pytest

from app.celery_app import celery_app
from app.services.batch import tasks
from app.services.batch.progress_tracker import progress_tracker
from app.services.batch.result_aggregator import result_storage


@pytest.fixture
def eager_chord(monkeypatch):
    conf = celery_app.conf
    prev = (conf.task_always_eager, conf.task_eager_propagates)
    conf.task_always_eager = True
    conf.task_eager_propagates = True

    monkeypatch.setattr(progress_tracker, "_redis", fakeredis.FakeRedis())
    monkeypatch.setattr(result_storage, "_redis", fakeredis.FakeRedis())

    async def _noop_audit(**kwargs):
        return None

    monkeypatch.setattr(tasks, "_log_batch_audit", _noop_audit)
    monkeypatch.setattr(tasks, "_init_analytics_and_dispatch", lambda job_id: None)
    yield
    conf.task_always_eager, conf.task_eager_propagates = prev


def test_chord_aggregation_stores_results_and_completes(eager_chord):
    molecules = [
        {"smiles": "CCO", "name": "ethanol", "index": 0, "properties": {}},
        {"smiles": "c1ccccc1", "name": "benzene", "index": 1, "properties": {}},
        {"smiles": "not_a_valid_smiles_string", "name": "bad", "index": 2, "properties": {}},
    ]
    job_id = "itest-chord-1"

    returned = tasks.process_batch_job(job_id, molecules)
    assert returned == job_id

    stored = result_storage.get_all_results(job_id)
    assert len(stored) == 3

    stats = result_storage.get_statistics(job_id)
    assert stats is not None
    assert stats.total == 3
    assert stats.successful >= 2  # the two valid molecules

    progress = progress_tracker.get_progress(job_id)
    assert progress is not None
    assert progress.status == "complete"
