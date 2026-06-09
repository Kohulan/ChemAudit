"""The /validate/async handler must not block the event loop on task.get()."""

import threading
from unittest.mock import MagicMock

import app.api.routes.validation as v


async def test_async_validate_offloads_blocking_get(client, monkeypatch):
    """The blocking Celery result fetch must run in a worker thread."""
    main_thread = threading.current_thread().name
    captured = {}

    fake_task = MagicMock()

    def fake_get(timeout=None):
        captured["thread"] = threading.current_thread().name
        return {
            "validation": {"overall_score": 99, "issues": []},
            "alerts": None,
            "scoring": None,
        }

    fake_task.get.side_effect = fake_get

    mock_task_fn = MagicMock()
    mock_task_fn.delay.return_value = fake_task
    monkeypatch.setattr(v, "validate_single_molecule", mock_task_fn)

    resp = await client.post(
        "/api/v1/validate/async",
        json={"molecule": "CCO", "format": "smiles"},
    )

    assert resp.status_code == 200, resp.text
    assert resp.json()["overall_score"] == 99
    # The blocking .get must have executed off the event-loop thread.
    assert captured["thread"] != main_thread
