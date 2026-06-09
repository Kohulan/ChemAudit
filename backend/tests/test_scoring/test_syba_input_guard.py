"""SYBA subprocess input must be validated before a process is spawned."""

import subprocess

import app.services.profiler.sa_comparison as sac


def test_syba_rejects_oversized_input(monkeypatch):
    spawned = {"n": 0}

    def fake_run(*a, **k):
        spawned["n"] += 1
        raise AssertionError("subprocess must not be spawned for invalid input")

    monkeypatch.setattr(subprocess, "run", fake_run)
    assert sac._syba_via_subprocess("C" * 5000) is None
    assert spawned["n"] == 0


def test_syba_rejects_nonprintable_input(monkeypatch):
    spawned = {"n": 0}

    def fake_run(*a, **k):
        spawned["n"] += 1
        raise AssertionError("subprocess must not be spawned for invalid input")

    monkeypatch.setattr(subprocess, "run", fake_run)
    assert sac._syba_via_subprocess("CC\nO") is None
    assert sac._syba_via_subprocess("CC\x00O") is None
    assert spawned["n"] == 0


def test_syba_runs_for_valid_input(monkeypatch):
    class FakeProc:
        returncode = 0
        stdout = "42.5"
        stderr = ""

    spawned = {"n": 0}

    def fake_run(*a, **k):
        spawned["n"] += 1
        assert k.get("shell") is False  # explicit, never shell=True
        return FakeProc()

    monkeypatch.setattr(subprocess, "run", fake_run)
    assert sac._syba_via_subprocess("CCO") == 42.5
    assert spawned["n"] == 1
