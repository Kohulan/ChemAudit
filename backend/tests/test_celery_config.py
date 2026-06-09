"""Tests for Celery application configuration hardening."""

from app.celery_app import celery_app


def test_task_time_limits_are_configured():
    """Hard and soft time limits must be set to prevent worker starvation."""
    conf = celery_app.conf
    assert conf.task_time_limit == 3600, "hard time limit should be 1 hour"
    assert conf.task_soft_time_limit == 3300, "soft limit should be 55 minutes"
    assert conf.task_soft_time_limit < conf.task_time_limit
