"""Tests for the WebSocket ConnectionManager (Redis pub/sub forwarding).

Covers the previously-untested connect/disconnect lifecycle, pub/sub forwarding,
broadcast fan-out + dead-connection pruning, initial-status delivery, and the
ownership guard that prevents a superseded subscriber from broadcasting.
"""

import asyncio
import json
from unittest.mock import AsyncMock

import fakeredis.aioredis as fakeaioredis
import pytest

from app.websockets.manager import ConnectionManager


def make_ws() -> AsyncMock:
    """A fake WebSocket with awaitable accept()/send_json()."""
    ws = AsyncMock()
    return ws


@pytest.fixture
async def mgr():
    m = ConnectionManager()
    m._redis = fakeaioredis.FakeRedis(decode_responses=True)
    await m._redis.ping()
    yield m
    for t in list(m._subscriber_tasks.values()):
        t.cancel()
    await asyncio.sleep(0)
    try:
        await m._redis.aclose()
    except Exception:
        pass


async def test_connect_registers_and_accepts(mgr):
    ws = make_ws()
    ok = await mgr.connect("job1", ws)
    assert ok is True
    ws.accept.assert_awaited_once()
    assert ws in mgr.active_connections["job1"]
    assert "job1" in mgr._subscriber_tasks


async def test_broadcast_sends_to_all_connections(mgr):
    a, b = make_ws(), make_ws()
    mgr.active_connections["job1"] = [a, b]
    await mgr._broadcast_to_job("job1", {"progress": 50})
    a.send_json.assert_awaited_once_with({"progress": 50})
    b.send_json.assert_awaited_once_with({"progress": 50})


async def test_broadcast_prunes_dead_connections(mgr):
    good, bad = make_ws(), make_ws()
    bad.send_json.side_effect = RuntimeError("closed")
    mgr.active_connections["job1"] = [good, bad]
    await mgr._broadcast_to_job("job1", {"x": 1})
    assert bad not in mgr.active_connections.get("job1", [])
    assert good in mgr.active_connections.get("job1", [])


async def test_pubsub_forwards_published_message(mgr):
    ws = make_ws()
    await mgr.connect("job1", ws)
    await mgr._redis.publish(
        "batch:progress:job1", json.dumps({"status": "processing", "progress": 25})
    )
    for _ in range(20):
        await asyncio.sleep(0.05)
        if ws.send_json.await_count:
            break
    ws.send_json.assert_awaited()


async def test_disconnect_cancels_subscriber_when_last_leaves(mgr):
    ws = make_ws()
    await mgr.connect("job1", ws)
    task = mgr._subscriber_tasks["job1"]
    mgr.disconnect("job1", ws)
    assert "job1" not in mgr.active_connections
    assert "job1" not in mgr._subscriber_tasks
    await asyncio.sleep(0.05)
    assert task.cancelled() or task.done()


async def test_send_initial_status_from_redis(mgr):
    ws = make_ws()
    await mgr._redis.set(
        "batch:job:job1", json.dumps({"status": "processing", "progress": 10})
    )
    ok = await mgr.send_initial_status("job1", ws)
    assert ok is True
    ws.send_json.assert_awaited_once()


async def test_send_initial_status_missing_returns_false(mgr):
    ws = make_ws()
    ok = await mgr.send_initial_status("nope", ws)
    assert ok is False


async def test_superseded_subscriber_does_not_broadcast(mgr):
    """A subscriber that is no longer the registered owner must not broadcast."""
    ws = make_ws()
    mgr.active_connections["job1"] = [ws]
    mgr._subscription_ready["job1"] = asyncio.Event()
    other_owner = asyncio.create_task(asyncio.sleep(10))
    mgr._subscriber_tasks["job1"] = other_owner  # a different task "owns" job1

    sub = asyncio.create_task(mgr._subscribe_to_job("job1"))
    await asyncio.wait_for(mgr._subscription_ready["job1"].wait(), timeout=2)
    await mgr._redis.publish(
        "batch:progress:job1", json.dumps({"status": "processing", "progress": 50})
    )
    await asyncio.sleep(0.2)

    assert sub.done()
    ws.send_json.assert_not_awaited()
    other_owner.cancel()
