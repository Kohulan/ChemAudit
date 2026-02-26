"""
Webhook Notification Service

Sends HMAC-SHA256 signed HTTP POST on batch completion with 3 retries
and exponential backoff via Celery autoretry.
"""

import hashlib
import hmac
import json
import logging
from datetime import datetime, timezone

import httpx

from app.celery_app import celery_app
from app.core.config import settings

logger = logging.getLogger(__name__)


def compute_hmac_signature(secret: str, payload_bytes: bytes) -> str:
    """
    Compute HMAC-SHA256 signature for webhook payload.

    Args:
        secret: Webhook secret key
        payload_bytes: JSON-encoded payload bytes

    Returns:
        Hex digest string prefixed with 'sha256='
    """
    sig = hmac.new(secret.encode("utf-8"), payload_bytes, hashlib.sha256).hexdigest()
    return f"sha256={sig}"


def build_webhook_payload(
    job_id: str,
    status: str,
    molecule_count: int,
    pass_count: int,
    fail_count: int,
    avg_score: float,
) -> dict:
    """
    Build webhook payload dict.

    Args:
        job_id: Batch job ID
        status: Job completion status ('complete' or 'failed')
        molecule_count: Total molecules processed
        pass_count: Number passing validation
        fail_count: Number failing validation
        avg_score: Average validation score

    Returns:
        Payload dict suitable for JSON serialization
    """
    return {
        "event": "batch.complete",
        "job_id": job_id,
        "status": status,
        "molecule_count": molecule_count,
        "pass_count": pass_count,
        "fail_count": fail_count,
        "avg_score": round(avg_score, 2),
        "report_url": f"{settings.BASE_URL}/batch/{job_id}",
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }


@celery_app.task(
    bind=True,
    autoretry_for=(Exception,),
    retry_backoff=True,
    max_retries=3,
    queue="default",
)
def send_webhook(self, url: str, secret: str, payload: dict):
    """
    Send HMAC-SHA256 signed webhook POST request.

    Retries up to 3 times with exponential backoff on failure.

    Args:
        url: Webhook endpoint URL
        secret: HMAC signing secret
        payload: JSON-serializable payload dict
    """
    body = json.dumps(payload, sort_keys=True).encode("utf-8")
    signature = compute_hmac_signature(secret, body)

    headers = {
        "Content-Type": "application/json",
        "X-Signature": signature,
        "X-Hook-Event": "batch.complete",
    }

    with httpx.Client(timeout=10) as client:
        resp = client.post(url, content=body, headers=headers)
        resp.raise_for_status()

    logger.info("Webhook sent to %s for job %s (status=%d)", url, payload.get("job_id"), resp.status_code)
