"""
Celery Beat task for purging expired session data.

Runs daily to remove bookmarks and audit entries tied to anonymous sessions
older than SESSION_MAX_AGE (30 days). API-key-scoped data is not affected.
"""

import logging
from datetime import datetime, timedelta

from sqlalchemy import delete

from app.celery_app import celery_app
from app.core.session import SESSION_MAX_AGE
from app.db.models.audit import ValidationAuditEntry
from app.db.models.bookmark import Bookmark

logger = logging.getLogger(__name__)


@celery_app.task(name="app.services.session.cleanup.purge_expired_sessions")
def purge_expired_sessions():
    """Delete session-scoped data older than 30 days.

    Only affects rows where session_id IS NOT NULL (anonymous browser sessions).
    API-key-scoped rows (api_key_hash IS NOT NULL) are preserved.

    Runs synchronously against the sync engine — Celery tasks have no event loop
    under the default prefork pool, so no asyncio bridging is needed.
    """
    from app.db import sync_session

    cutoff = datetime.utcnow() - timedelta(seconds=SESSION_MAX_AGE)

    with sync_session() as db:
        # Delete expired bookmarks (session-scoped only)
        bm_count = db.execute(
            delete(Bookmark).where(
                Bookmark.session_id.isnot(None),
                Bookmark.api_key_hash.is_(None),
                Bookmark.created_at < cutoff,
            )
        ).rowcount

        # Delete expired audit entries (session-scoped only)
        audit_count = db.execute(
            delete(ValidationAuditEntry).where(
                ValidationAuditEntry.session_id.isnot(None),
                ValidationAuditEntry.api_key_hash.is_(None),
                ValidationAuditEntry.created_at < cutoff,
            )
        ).rowcount

        db.commit()

    logger.info(
        "Session purge complete: %d bookmarks, %d audit entries removed (cutoff: %s)",
        bm_count,
        audit_count,
        cutoff.isoformat(),
    )

    return {"bookmarks_deleted": bm_count, "audit_entries_deleted": audit_count}
