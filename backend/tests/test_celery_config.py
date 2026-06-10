"""Tests for Celery application configuration hardening."""

from app.celery_app import celery_app


def test_task_time_limits_are_configured():
    """Hard and soft time limits must be set to prevent worker starvation."""
    conf = celery_app.conf
    assert conf.task_time_limit == 3600, "hard time limit should be 1 hour"
    assert conf.task_soft_time_limit == 3300, "soft limit should be 55 minutes"
    assert conf.task_soft_time_limit < conf.task_time_limit


def test_analytics_queue_defined_and_routed():
    """Expensive analytics must route to a dedicated 'analytics' queue."""
    queue_names = {q.name for q in celery_app.conf.task_queues}
    assert "analytics" in queue_names

    routes = celery_app.conf.task_routes
    key = "app.services.batch.analytics_tasks.run_expensive_analytics"
    assert routes[key]["queue"] == "analytics"


def test_expensive_analytics_has_own_time_limit():
    """The expensive analytics task must carry a tighter per-task time limit."""
    from app.services.batch.analytics_tasks import run_expensive_analytics

    assert run_expensive_analytics.soft_time_limit == 600
    assert run_expensive_analytics.time_limit == 660
