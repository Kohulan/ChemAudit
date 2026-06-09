"""Persistent SYBA worker manager: IPC, lifecycle, restart, timeout, concurrency.

SYBA itself is not installed, so these tests drive the manager with fake worker
scripts that speak the same JSON-line protocol — fully exercising the manager's
infrastructure without the GPL dependency.
"""

import sys
import textwrap
import threading

import app.services.profiler.sa_comparison as sac

_PRELUDE = "import json, sys, time\n"

GOOD = """
sys.stdout.write(json.dumps({"ready": True}) + "\\n"); sys.stdout.flush()
for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    json.loads(line)
    sys.stdout.write(json.dumps({"score": 42.5}) + "\\n"); sys.stdout.flush()
"""

FATAL = """
sys.stdout.write(json.dumps({"ready": False, "fatal": True, "error": "no syba"}) + "\\n")
sys.stdout.flush()
"""

HANG = """
sys.stdout.write(json.dumps({"ready": True}) + "\\n"); sys.stdout.flush()
for line in sys.stdin:
    time.sleep(30)
"""

CRASH_ONCE = """
sys.stdout.write(json.dumps({"ready": True}) + "\\n"); sys.stdout.flush()
sys.stdin.readline()
sys.exit(1)
"""


def _command(tmp_path, name, body):
    script = tmp_path / name
    script.write_text(_PRELUDE + textwrap.dedent(body))
    return [sys.executable, "-u", str(script)]


def test_predict_returns_score_and_reuses_worker(tmp_path):
    w = sac._SybaWorker(command=_command(tmp_path, "good.py", GOOD))
    try:
        assert w.predict("CCO") == 42.5
        proc = w._proc
        assert w.predict("c1ccccc1") == 42.5
        assert w._proc is proc, "worker must be reused, not respawned"
    finally:
        w.shutdown()


def test_fatal_handshake_is_sticky(tmp_path):
    w = sac._SybaWorker(command=_command(tmp_path, "fatal.py", FATAL))
    try:
        assert w.predict("CCO") is None
        assert w._fatal is True
        assert w._proc is None  # not left running
    finally:
        w.shutdown()


def test_input_guard_rejects_without_spawning(tmp_path):
    w = sac._SybaWorker(command=_command(tmp_path, "good.py", GOOD))
    try:
        assert w.predict("C" * 5000) is None
        assert w.predict("CC\nO") is None
        assert w.predict("CC\x00O") is None
        assert w._proc is None, "invalid input must be rejected before any spawn"
    finally:
        w.shutdown()


def test_timeout_kills_worker(tmp_path, monkeypatch):
    monkeypatch.setattr(sac, "_SYBA_PREDICT_TIMEOUT", 1.0)
    w = sac._SybaWorker(command=_command(tmp_path, "hang.py", HANG))
    try:
        assert w.predict("CCO") is None
        assert w._proc is None  # hung worker was killed
    finally:
        w.shutdown()


def test_crash_then_transparent_restart(tmp_path):
    w = sac._SybaWorker(command=_command(tmp_path, "crash.py", CRASH_ONCE))
    try:
        assert w.predict("CCO") is None  # worker exits on the request
        # A subsequent call must not raise — the manager respawns a fresh worker.
        assert w.predict("CCO") is None
    finally:
        w.shutdown()


def test_concurrent_predictions_are_serialised(tmp_path):
    w = sac._SybaWorker(command=_command(tmp_path, "good.py", GOOD))
    results = []

    def call():
        results.append(w.predict("CCO"))

    threads = [threading.Thread(target=call) for _ in range(8)]
    try:
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        assert results == [42.5] * 8
    finally:
        w.shutdown()
