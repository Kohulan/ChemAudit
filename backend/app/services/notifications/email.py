"""
Email Notification Service

Sends batch completion email via SMTP through a Celery task.
Uses a Jinja2 HTML template for professional email formatting.
"""

import logging
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from pathlib import Path

from jinja2 import Environment, FileSystemLoader

from app.celery_app import celery_app
from app.core.config import settings

logger = logging.getLogger(__name__)

# Template directory
_TEMPLATE_DIR = Path(__file__).parent.parent.parent / "templates" / "emails"


def _render_email_template(job_id: str, stats: dict) -> str:
    """
    Render the batch completion email template.

    Args:
        job_id: Batch job ID
        stats: Stats dict with molecule_count, pass_count, fail_count, avg_score

    Returns:
        Rendered HTML string
    """
    env = Environment(loader=FileSystemLoader(str(_TEMPLATE_DIR)))
    template = env.get_template("batch_complete.html")
    return template.render(
        job_id=job_id,
        base_url=settings.BASE_URL,
        molecule_count=stats.get("molecule_count", 0),
        pass_count=stats.get("pass_count", 0),
        fail_count=stats.get("fail_count", 0),
        avg_score=stats.get("avg_score", 0),
    )


@celery_app.task(queue="default")
def send_batch_complete_email(to: str, job_id: str, stats: dict):
    """
    Send batch completion email via SMTP.

    Fire-and-forget with logging. Does not raise on SMTP errors.

    Args:
        to: Recipient email address
        job_id: Batch job ID
        stats: Stats dict with molecule_count, pass_count, fail_count, avg_score
    """
    try:
        html_content = _render_email_template(job_id, stats)

        msg = MIMEMultipart("alternative")
        msg["Subject"] = f"ChemAudit: Batch {job_id[:8]} complete"
        msg["From"] = settings.SMTP_FROM
        msg["To"] = to

        # Plain text fallback
        plain = (
            f"Your batch validation job {job_id[:8]} is complete.\n"
            f"Molecules: {stats.get('molecule_count', 0)}, "
            f"Passed: {stats.get('pass_count', 0)}, "
            f"Failed: {stats.get('fail_count', 0)}\n"
            f"View results: {settings.BASE_URL}/batch/{job_id}"
        )
        msg.attach(MIMEText(plain, "plain"))
        msg.attach(MIMEText(html_content, "html"))

        with smtplib.SMTP(settings.SMTP_HOST, settings.SMTP_PORT) as server:
            if settings.SMTP_TLS:
                server.starttls()
            if settings.SMTP_USER:
                server.login(settings.SMTP_USER, settings.SMTP_PASSWORD)
            server.sendmail(settings.SMTP_FROM, [to], msg.as_string())

        logger.info("Batch completion email sent to %s for job %s", to, job_id[:8])

    except Exception as e:
        logger.warning("Failed to send email to %s: %s", to, e)
