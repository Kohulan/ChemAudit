"""
Celery Application Configuration

Configures Celery for batch processing of molecules with Redis as broker/backend.
"""
from celery import Celery

from app.core.config import settings

celery_app = Celery(
    "chemstructval",
    broker=settings.REDIS_URL,
    backend=settings.REDIS_URL,
    include=["app.services.batch.tasks"],
)

celery_app.conf.update(
    task_serializer="json",
    result_serializer="json",
    accept_content=["json"],
    result_expires=3600,  # Results expire after 1 hour
    task_track_started=True,
    worker_prefetch_multiplier=1,  # Process one task at a time for accurate progress
    task_acks_late=True,  # Acknowledge tasks after completion for reliability
    task_reject_on_worker_lost=True,  # Requeue tasks if worker dies
)
