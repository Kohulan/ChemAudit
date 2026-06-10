"""Batch tasks must only auto-retry transient infrastructure errors."""

from app.services.batch import tasks


def test_transient_errors_constant_excludes_generic_exception():
    assert Exception not in tasks.TRANSIENT_ERRORS
    assert ValueError not in tasks.TRANSIENT_ERRORS


def test_transient_errors_constant_includes_connection_errors():
    import redis.exceptions

    assert ConnectionError in tasks.TRANSIENT_ERRORS
    assert TimeoutError in tasks.TRANSIENT_ERRORS
    assert redis.exceptions.ConnectionError in tasks.TRANSIENT_ERRORS


def test_tasks_use_transient_errors_tuple():
    # Assert the decorators reference the shared constant rather than a blanket
    # (Exception,) retry policy.
    import inspect

    src = inspect.getsource(tasks)
    assert "autoretry_for=(Exception,)" not in src
    assert src.count("autoretry_for=TRANSIENT_ERRORS") == 3
