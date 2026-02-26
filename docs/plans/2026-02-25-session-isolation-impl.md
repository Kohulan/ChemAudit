# Session-Based Data Isolation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Isolate bookmarks and history per-browser-session (anonymous) or per-API-key, enforced at both application and database levels, with GDPR-compliant data purge and accurate Privacy page.

**Architecture:** HttpOnly cookie carries UUID v4 session ID. All bookmark/history data scoped by `session_id` (anonymous) or `api_key_hash` (API key users). PostgreSQL Row-Level Security as defense-in-depth. Celery Beat auto-purges stale sessions after 30 days.

**Tech Stack:** FastAPI middleware, SQLAlchemy/asyncpg, PostgreSQL 16 RLS, Alembic migration, Axios `withCredentials`, Celery Beat periodic task.

---

### Task 1: Add `session_id` Column to ORM Models

**Files:**
- Modify: `backend/app/db/models/bookmark.py:28` (add session_id after api_key_hash)
- Modify: `backend/app/db/models/audit.py:29` (add session_id after api_key_hash)

**Step 1: Add `session_id` to Bookmark model**

In `backend/app/db/models/bookmark.py`, add after line 30 (the `api_key_hash` column):

```python
session_id: Mapped[str | None] = mapped_column(
    String(36), nullable=True, index=True
)
```

**Step 2: Add `session_id` to ValidationAuditEntry model**

In `backend/app/db/models/audit.py`, add after line 29 (the `api_key_hash` column):

```python
session_id: Mapped[str | None] = mapped_column(
    String(36), nullable=True, index=True
)
```

Import `String` if not already imported (it is — line 9 in both files).

**Step 3: Commit**

```bash
git add backend/app/db/models/bookmark.py backend/app/db/models/audit.py
git commit -m "feat: add session_id column to bookmark and audit models"
```

---

### Task 2: Alembic Migration — Add Columns, Purge Unscoped Data, Enable RLS

**Files:**
- Create: `backend/alembic/versions/002_session_isolation.py`

**Step 1: Create the migration**

```python
"""Add session_id columns and enable row-level security.

Revision ID: 002
Revises: 001
"""

from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa

revision: str = "002"
down_revision: Union[str, None] = "001"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # --- 1. Add session_id columns ---
    op.add_column("bookmarks", sa.Column("session_id", sa.String(36), nullable=True))
    op.create_index("ix_bookmarks_session_id", "bookmarks", ["session_id"])

    op.add_column("validation_audit", sa.Column("session_id", sa.String(36), nullable=True))
    op.create_index("ix_validation_audit_session_id", "validation_audit", ["session_id"])

    # --- 2. Purge unscoped legacy data ---
    op.execute("DELETE FROM bookmarks WHERE session_id IS NULL AND api_key_hash IS NULL")
    op.execute(
        "DELETE FROM validation_audit WHERE session_id IS NULL AND api_key_hash IS NULL"
    )

    # --- 3. Enable Row-Level Security ---
    # Bookmarks
    op.execute("ALTER TABLE bookmarks ENABLE ROW LEVEL SECURITY")
    op.execute("ALTER TABLE bookmarks FORCE ROW LEVEL SECURITY")
    op.execute("""
        CREATE POLICY bookmark_isolation ON bookmarks
        USING (
            session_id = NULLIF(current_setting('app.session_id', true), '')
            OR api_key_hash = NULLIF(current_setting('app.api_key_hash', true), '')
        )
    """)

    # Audit trail
    op.execute("ALTER TABLE validation_audit ENABLE ROW LEVEL SECURITY")
    op.execute("ALTER TABLE validation_audit FORCE ROW LEVEL SECURITY")
    op.execute("""
        CREATE POLICY audit_isolation ON validation_audit
        USING (
            session_id = NULLIF(current_setting('app.session_id', true), '')
            OR api_key_hash = NULLIF(current_setting('app.api_key_hash', true), '')
        )
    """)


def downgrade() -> None:
    # Remove RLS
    op.execute("DROP POLICY IF EXISTS audit_isolation ON validation_audit")
    op.execute("ALTER TABLE validation_audit DISABLE ROW LEVEL SECURITY")
    op.execute("DROP POLICY IF EXISTS bookmark_isolation ON bookmarks")
    op.execute("ALTER TABLE bookmarks DISABLE ROW LEVEL SECURITY")

    # Remove columns
    op.drop_index("ix_validation_audit_session_id", "validation_audit")
    op.drop_column("validation_audit", "session_id")
    op.drop_index("ix_bookmarks_session_id", "bookmarks")
    op.drop_column("bookmarks", "session_id")
```

**Step 2: Run the migration locally to test**

```bash
cd backend && alembic upgrade head
```

Expected: migration applies cleanly. If DB has existing data, it will be purged.

**Step 3: Commit**

```bash
git add backend/alembic/versions/002_session_isolation.py
git commit -m "feat: add session isolation migration with RLS policies"
```

**Note on RLS + table owner:** PostgreSQL's `FORCE ROW LEVEL SECURITY` makes RLS apply even to the table owner. If the app connects as the table owner (`chemaudit`), `FORCE` ensures policies still apply. For production, consider a separate `chemaudit_app` role, but `FORCE` handles the single-role case.

---

### Task 3: Session Middleware (`backend/app/core/session.py`)

**Files:**
- Create: `backend/app/core/session.py`

**Step 1: Create session management module**

```python
"""
Session Management for Data Isolation

Provides per-browser session scoping via HttpOnly cookie.
API key users are scoped by api_key_hash instead.
"""

import uuid
from typing import Optional

from fastapi import Request, Response
from sqlalchemy import text
from sqlalchemy.ext.asyncio import AsyncSession

from app.core.security import get_api_key, hash_api_key_for_lookup

SESSION_COOKIE = "chemaudit_sid"
SESSION_MAX_AGE = 30 * 24 * 3600  # 30 days


def get_session_id(request: Request) -> Optional[str]:
    """Read session ID from cookie. Returns None if no cookie."""
    return request.cookies.get(SESSION_COOKIE)


def ensure_session_cookie(response: Response, session_id: str) -> None:
    """Set session cookie on response if not already present."""
    response.set_cookie(
        key=SESSION_COOKIE,
        value=session_id,
        max_age=SESSION_MAX_AGE,
        httponly=True,
        secure=True,
        samesite="lax",
    )


def create_session_id() -> str:
    """Generate a new cryptographically random session ID."""
    return str(uuid.uuid4())


async def get_data_scope(request: Request) -> tuple[Optional[str], Optional[str]]:
    """
    Determine data scope from request.

    Returns:
        (session_id, api_key_hash) — API key takes precedence.
    """
    # Check for API key first
    api_key = request.headers.get("X-API-Key")
    if api_key:
        api_key_hash = hash_api_key_for_lookup(api_key)
        return None, api_key_hash

    # Fall back to session cookie
    session_id = get_session_id(request)
    return session_id, None


async def set_rls_context(
    db: AsyncSession, session_id: Optional[str], api_key_hash: Optional[str]
) -> None:
    """Set PostgreSQL session variables for RLS policy evaluation."""
    await db.execute(
        text("SET LOCAL app.session_id = :sid"),
        {"sid": session_id or ""},
    )
    await db.execute(
        text("SET LOCAL app.api_key_hash = :akh"),
        {"akh": api_key_hash or ""},
    )
```

**Step 2: Commit**

```bash
git add backend/app/core/session.py
git commit -m "feat: add session management module for data isolation"
```

---

### Task 4: Integrate Session into Bookmark Routes

**Files:**
- Modify: `backend/app/api/routes/bookmarks.py` (all endpoints)

**Step 1: Add session imports and dependency**

At the top of `bookmarks.py`, add these imports (after existing imports around line 16):

```python
from fastapi import APIRouter, Depends, HTTPException, Query, Request, Response
from app.core.session import (
    create_session_id,
    ensure_session_cookie,
    get_data_scope,
    get_session_id,
    set_rls_context,
)
from app.core.security import get_api_key, hash_api_key_for_lookup
```

Note: `Request` and `Response` need to be added to the `fastapi` import line.

**Step 2: Update `create_bookmark` (line 132-149)**

Replace the endpoint with:

```python
@router.post("/bookmarks", response_model=BookmarkResponse, status_code=201)
async def create_bookmark(
    body: BookmarkCreate,
    request: Request,
    response: Response,
    db: AsyncSession = Depends(get_db),
    api_key: Optional[str] = Depends(get_api_key),
):
    """Create a new bookmark. Auto-computes InChIKey from SMILES."""
    # Determine scope
    session_id, api_key_hash = await get_data_scope(request)
    if api_key:
        api_key_hash = hash_api_key_for_lookup(api_key)
    if not session_id and not api_key_hash:
        session_id = create_session_id()
    if session_id:
        ensure_session_cookie(response, session_id)

    # Set RLS context
    await set_rls_context(db, session_id, api_key_hash)

    inchikey = _compute_inchikey(body.smiles)
    bm = Bookmark(
        smiles=body.smiles,
        name=body.name,
        inchikey=inchikey,
        tags=_tags_to_str(body.tags),
        notes=body.notes,
        source=body.source,
        job_id=body.job_id,
        session_id=session_id,
        api_key_hash=api_key_hash,
    )
    db.add(bm)
    await db.commit()
    await db.refresh(bm)
    return _bookmark_to_response(bm)
```

**Step 3: Update `list_bookmarks` (line 76-119)**

Add `Request` parameter and session scoping. Replace the function signature and add scoping after line 87:

```python
@router.get("/bookmarks", response_model=BookmarkListResponse)
async def list_bookmarks(
    request: Request,
    page: int = Query(default=1, ge=1),
    page_size: int = Query(default=50, ge=1, le=200),
    tag: Optional[str] = Query(None, description="Single tag filter"),
    tags: Optional[str] = Query(None, description="Comma-separated tag filter"),
    source: Optional[str] = Query(None, description="Source filter"),
    search: Optional[str] = Query(None, description="SMILES substring search"),
    db: AsyncSession = Depends(get_db),
    api_key: Optional[str] = Depends(get_api_key),
):
    """List bookmarks with optional filters and pagination."""
    session_id, api_key_hash = await get_data_scope(request)
    if api_key:
        api_key_hash = hash_api_key_for_lookup(api_key)
    await set_rls_context(db, session_id, api_key_hash)

    base_query = select(Bookmark)

    # App-level filtering (defense-in-depth with RLS)
    if api_key_hash:
        base_query = base_query.where(Bookmark.api_key_hash == api_key_hash)
    elif session_id:
        base_query = base_query.where(Bookmark.session_id == session_id)
    else:
        # No session, no API key — return nothing
        return {"bookmarks": [], "total": 0, "page": page, "page_size": page_size}

    # ... rest of filters unchanged (tag, source, search) ...
```

**Step 4: Update `get_bookmark`, `update_bookmark`, `delete_bookmark`, `bulk_delete_bookmarks`**

Add `Request` and `api_key` params to each. Add scope check. The ownership check pattern for single-item endpoints:

```python
@router.get("/bookmarks/{bookmark_id}", response_model=BookmarkResponse)
async def get_bookmark(
    bookmark_id: int,
    request: Request,
    db: AsyncSession = Depends(get_db),
    api_key: Optional[str] = Depends(get_api_key),
):
    """Get a single bookmark by ID."""
    session_id, api_key_hash = await get_data_scope(request)
    if api_key:
        api_key_hash = hash_api_key_for_lookup(api_key)
    await set_rls_context(db, session_id, api_key_hash)

    result = await db.execute(select(Bookmark).where(Bookmark.id == bookmark_id))
    bm = result.scalar_one_or_none()
    if bm is None:
        raise HTTPException(status_code=404, detail="Bookmark not found")
    return _bookmark_to_response(bm)
```

Apply the same pattern to `update_bookmark`, `delete_bookmark`, `bookmark_batch_submit`, and `bulk_delete_bookmarks`. Each gets:
1. `request: Request` parameter
2. `api_key: Optional[str] = Depends(get_api_key)` parameter
3. `get_data_scope()` + `set_rls_context()` at the top

**Step 5: Commit**

```bash
git add backend/app/api/routes/bookmarks.py
git commit -m "feat: scope bookmark endpoints by session/API key"
```

---

### Task 5: Integrate Session into History Routes

**Files:**
- Modify: `backend/app/api/routes/history.py`

**Step 1: Add imports and update `get_history`**

Add session imports. Remove the public `api_key_hash` query parameter. Replace with automatic scoping:

```python
from fastapi import APIRouter, Depends, Query, Request
from app.core.session import get_data_scope, set_rls_context
from app.core.security import get_api_key, hash_api_key_for_lookup
```

Update `get_history`:

```python
@router.get("/history", response_model=AuditHistoryResponse)
async def get_history(
    request: Request,
    page: int = Query(default=1, ge=1),
    page_size: int = Query(default=50, ge=1, le=200),
    date_from: Optional[datetime] = Query(None),
    date_to: Optional[datetime] = Query(None),
    outcome: Optional[str] = Query(None),
    source: Optional[str] = Query(None),
    smiles_search: Optional[str] = Query(None),
    db: AsyncSession = Depends(get_db),
    api_key: Optional[str] = Depends(get_api_key),
):
    """Get paginated audit trail scoped to current session/API key."""
    session_id, api_key_hash = await get_data_scope(request)
    if api_key:
        api_key_hash = hash_api_key_for_lookup(api_key)
    await set_rls_context(db, session_id, api_key_hash)

    query = select(ValidationAuditEntry)
    count_query = select(func.count(ValidationAuditEntry.id))

    # App-level scope filtering (defense-in-depth with RLS)
    if api_key_hash:
        query = query.where(ValidationAuditEntry.api_key_hash == api_key_hash)
        count_query = count_query.where(ValidationAuditEntry.api_key_hash == api_key_hash)
    elif session_id:
        query = query.where(ValidationAuditEntry.session_id == session_id)
        count_query = count_query.where(ValidationAuditEntry.session_id == session_id)
    else:
        return AuditHistoryResponse(entries=[], total=0, page=page, page_size=page_size)

    # ... rest of filters unchanged (date_from, date_to, outcome, source, smiles_search) ...
```

**Step 2: Update `get_history_stats` to be session-scoped**

```python
@router.get("/history/stats")
async def get_history_stats(
    request: Request,
    db: AsyncSession = Depends(get_db),
    api_key: Optional[str] = Depends(get_api_key),
):
    """Get summary statistics scoped to current session/API key."""
    session_id, api_key_hash = await get_data_scope(request)
    if api_key:
        api_key_hash = hash_api_key_for_lookup(api_key)
    await set_rls_context(db, session_id, api_key_hash)

    # Build scope filter
    scope_filter = []
    if api_key_hash:
        scope_filter.append(ValidationAuditEntry.api_key_hash == api_key_hash)
    elif session_id:
        scope_filter.append(ValidationAuditEntry.session_id == session_id)
    else:
        return {"total_validations": 0, "outcome_distribution": {}, "source_distribution": {}}

    total_result = await db.execute(
        select(func.count(ValidationAuditEntry.id)).where(*scope_filter)
    )
    # ... rest uses *scope_filter on each query ...
```

**Step 3: Commit**

```bash
git add backend/app/api/routes/history.py
git commit -m "feat: scope history endpoints by session/API key"
```

---

### Task 6: Wire Session into Validation Audit Logging

**Files:**
- Modify: `backend/app/api/routes/validation.py:205-221`
- Modify: `backend/app/services/audit/service.py:18-65` (add session_id parameter)
- Modify: `backend/app/services/batch/tasks.py:470-499` (pass session_id from Redis job metadata)

**Step 1: Add `session_id` parameter to audit service**

In `backend/app/services/audit/service.py`, add `session_id: Optional[str] = None` parameter to both functions. Update the `ValidationAuditEntry` constructor to include `session_id=session_id`.

```python
async def log_validation_event(
    db: AsyncSession,
    smiles: str,
    inchikey: Optional[str],
    outcome: str,
    score: Optional[int],
    source: str,
    job_id: Optional[str] = None,
    molecule_count: Optional[int] = None,
    pass_count: Optional[int] = None,
    fail_count: Optional[int] = None,
    api_key_hash: Optional[str] = None,
    session_id: Optional[str] = None,
) -> ValidationAuditEntry:
    # ...
    entry = ValidationAuditEntry(
        smiles=smiles,
        inchikey=inchikey,
        outcome=outcome,
        score=score,
        source=source,
        job_id=job_id,
        molecule_count=molecule_count,
        pass_count=pass_count,
        fail_count=fail_count,
        api_key_hash=api_key_hash,
        session_id=session_id,
    )
    # ...
```

Also update `log_batch_event` to accept and forward `session_id`.

**Step 2: Pass session scope from validation route**

In `backend/app/api/routes/validation.py`, the `validate_molecule` endpoint (around line 60) already has `api_key: Optional[str] = Depends(get_api_key)`. Add session imports and pass scope to audit logging.

Add imports:
```python
from app.core.session import create_session_id, ensure_session_cookie, get_data_scope
from app.core.security import hash_api_key_for_lookup
```

Add `request: Request, response: Response` to endpoint params (if not already present).

Before the audit logging block (line 205), compute scope:
```python
    session_id, api_key_hash = await get_data_scope(request)
    if api_key:
        api_key_hash = hash_api_key_for_lookup(api_key)
    if not session_id and not api_key_hash:
        session_id = create_session_id()
    if session_id:
        ensure_session_cookie(response, session_id)
```

Then pass to `log_validation_event`:
```python
    await log_validation_event(
        db,
        smiles=body.molecule,
        inchikey=mol_info.inchikey,
        outcome=audit_outcome,
        score=score,
        source="single",
        api_key_hash=api_key_hash,
        session_id=session_id,
    )
```

**Step 3: Handle batch audit logging**

For batch jobs (`backend/app/services/batch/tasks.py`), the session context is not available inside the Celery worker (no HTTP request). Two options:
- **Option A:** Store `session_id` in the Redis job metadata when the batch is submitted, then read it in the worker. This is the correct approach.
- **Option B:** Skip session scoping for batch audit (batch jobs are already keyed by UUID job_id which is unguessable).

Implement Option A: In the batch submission endpoint, store `session_id` and `api_key_hash` in the job metadata dict that goes to Redis. In the Celery worker, read those from the job dict and pass to `log_batch_event`.

**Step 4: Commit**

```bash
git add backend/app/services/audit/service.py backend/app/api/routes/validation.py backend/app/services/batch/tasks.py
git commit -m "feat: pass session/API key scope through audit logging pipeline"
```

---

### Task 7: Frontend — Enable Cookie Transmission

**Files:**
- Modify: `frontend/src/services/api.ts:167-173`

**Step 1: Add `withCredentials` to axios instance**

In `frontend/src/services/api.ts`, update the `api` instance (line 167-173):

```typescript
export const api = axios.create({
  baseURL: API_BASE_URL,
  timeout: 30000,
  withCredentials: true,  // Send HttpOnly session cookie with every request
  headers: {
    'Content-Type': 'application/json',
  },
});
```

Also check the CSRF token fetch at line 149 — it uses raw `axios.get()` not the `api` instance. Update it to also include credentials:

```typescript
const response = await axios.get(`${API_BASE_URL}/csrf-token`, {
  withCredentials: true,
});
```

**Step 2: Commit**

```bash
git add frontend/src/services/api.ts
git commit -m "feat: enable cookie credentials on API client"
```

---

### Task 8: Backend — Update CORS Configuration

**Files:**
- Modify: `backend/app/core/config.py:46` (verify `allow_credentials` compatible origins)
- Verify: `backend/app/main.py:157-164` (CORS middleware already has `allow_credentials=True`)

**Step 1: Verify CORS config**

Check that `CORS_ORIGINS_STR` does NOT include `*` (it doesn't — it lists specific origins). The existing CORS middleware at `main.py:157-164` already has `allow_credentials=True`. No change needed here.

**Step 2: Verify cookie works in dev**

In development, frontend is at `http://localhost:3002` and backend at `http://localhost:8001`. The `Secure` flag on the cookie means it won't be sent over HTTP. For local dev, we need to conditionally skip `Secure`:

In `backend/app/core/session.py`, update `ensure_session_cookie`:

```python
from app.core.config import settings

def ensure_session_cookie(response: Response, session_id: str) -> None:
    """Set session cookie on response."""
    is_dev = "localhost" in settings.CORS_ORIGINS_STR or "127.0.0.1" in settings.CORS_ORIGINS_STR
    response.set_cookie(
        key=SESSION_COOKIE,
        value=session_id,
        max_age=SESSION_MAX_AGE,
        httponly=True,
        secure=not is_dev,  # Secure=False for localhost dev
        samesite="lax",
    )
```

**Step 3: Commit**

```bash
git add backend/app/core/session.py
git commit -m "fix: allow session cookie on localhost for development"
```

---

### Task 9: Purge My Data Endpoint

**Files:**
- Create: `backend/app/api/routes/session.py`
- Modify: `backend/app/main.py` (register router)

**Step 1: Create session route**

```python
"""
Session Data Management Routes

GDPR-compliant data purge endpoint.
"""

from fastapi import APIRouter, Depends, HTTPException, Request, Response
from sqlalchemy import delete, func, select
from sqlalchemy.ext.asyncio import AsyncSession

from app.core.security import get_api_key, hash_api_key_for_lookup
from app.core.session import SESSION_COOKIE, get_data_scope, set_rls_context
from app.db import get_db
from app.db.models.audit import ValidationAuditEntry
from app.db.models.bookmark import Bookmark

router = APIRouter()


@router.delete("/me/data")
async def purge_my_data(
    request: Request,
    response: Response,
    db: AsyncSession = Depends(get_db),
    api_key: str | None = Depends(get_api_key),
):
    """
    Delete all data associated with the current session or API key.

    GDPR Article 17 — Right to Erasure.
    """
    session_id, api_key_hash = await get_data_scope(request)
    if api_key:
        api_key_hash = hash_api_key_for_lookup(api_key)

    if not session_id and not api_key_hash:
        raise HTTPException(status_code=404, detail="No active session or API key")

    # Count before delete
    bm_filter = (
        Bookmark.api_key_hash == api_key_hash if api_key_hash
        else Bookmark.session_id == session_id
    )
    audit_filter = (
        ValidationAuditEntry.api_key_hash == api_key_hash if api_key_hash
        else ValidationAuditEntry.session_id == session_id
    )

    bm_count = (await db.execute(select(func.count()).where(bm_filter))).scalar() or 0
    audit_count = (await db.execute(select(func.count()).where(audit_filter))).scalar() or 0

    # Delete
    await db.execute(delete(Bookmark).where(bm_filter))
    await db.execute(delete(ValidationAuditEntry).where(audit_filter))
    await db.commit()

    # Clear session cookie
    if session_id:
        response.delete_cookie(SESSION_COOKIE)

    return {
        "status": "purged",
        "deleted": {"bookmarks": bm_count, "history": audit_count},
    }
```

**Step 2: Register the router in `main.py`**

In `backend/app/main.py`, add to the imports block (around line 44-59):

```python
from app.api.routes import session as session_routes
```

And register the router (with the other `app.include_router` calls):

```python
app.include_router(session_routes.router, prefix="/api/v1", tags=["session"])
```

**Step 3: Commit**

```bash
git add backend/app/api/routes/session.py backend/app/main.py
git commit -m "feat: add GDPR data purge endpoint DELETE /me/data"
```

---

### Task 10: Celery Beat Auto-Purge Task

**Files:**
- Create: `backend/app/services/session/cleanup.py`
- Modify: `backend/app/services/batch/celery_app.py` (add beat schedule)

**Step 1: Create cleanup task**

```python
"""
Session Data Cleanup

Periodic task to purge expired anonymous session data.
Runs daily via Celery Beat.
"""

import logging

from sqlalchemy import delete, func, text
from sqlalchemy.ext.asyncio import AsyncSession

from app.db import async_session
from app.db.models.audit import ValidationAuditEntry
from app.db.models.bookmark import Bookmark

logger = logging.getLogger(__name__)

SESSION_TTL_DAYS = 30


async def _purge_expired_sessions() -> dict:
    """Delete session-scoped data older than SESSION_TTL_DAYS."""
    async with async_session() as db:
        # Bypass RLS for cleanup (use superuser SET)
        await db.execute(text("SET LOCAL app.session_id = ''"))
        await db.execute(text("SET LOCAL app.api_key_hash = ''"))

        cutoff = func.now() - text(f"INTERVAL '{SESSION_TTL_DAYS} days'")

        bm_result = await db.execute(
            delete(Bookmark)
            .where(Bookmark.session_id.isnot(None))
            .where(Bookmark.created_at < cutoff)
            .returning(func.count())
        )
        bm_count = bm_result.scalar() or 0

        audit_result = await db.execute(
            delete(ValidationAuditEntry)
            .where(ValidationAuditEntry.session_id.isnot(None))
            .where(ValidationAuditEntry.created_at < cutoff)
            .returning(func.count())
        )
        audit_count = audit_result.scalar() or 0

        await db.commit()

    logger.info(
        "Session cleanup: purged %d bookmarks, %d audit entries (>%d days)",
        bm_count, audit_count, SESSION_TTL_DAYS,
    )
    return {"bookmarks": bm_count, "audit": audit_count}
```

**Step 2: Register Celery Beat schedule**

Find the Celery app configuration and add a beat schedule entry. Check `backend/app/services/batch/celery_app.py` for the Celery instance. Add:

```python
app.conf.beat_schedule = {
    **getattr(app.conf, "beat_schedule", {}),
    "purge-expired-sessions": {
        "task": "app.services.session.cleanup.purge_expired_sessions",
        "schedule": 86400.0,  # daily
    },
}
```

**Note:** The exact wiring depends on how Celery tasks are registered in this project. The cleanup function is async, so it may need a sync wrapper using `asyncio.run()` similar to how batch tasks handle async DB operations.

**Step 3: Commit**

```bash
git add backend/app/services/session/cleanup.py backend/app/services/batch/celery_app.py
git commit -m "feat: add daily Celery Beat task for session data cleanup"
```

---

### Task 11: Frontend — Purge My Data Button

**Files:**
- Modify: `frontend/src/pages/Privacy.tsx`
- Modify: `frontend/src/services/api.ts` (add purge API call)

**Step 1: Add purge API function**

In `frontend/src/services/api.ts`, add to the exports:

```typescript
export const sessionApi = {
  purgeMyData: () => api.delete<{ status: string; deleted: { bookmarks: number; history: number } }>('/me/data'),
};
```

**Step 2: Add Purge button to Privacy page**

In `frontend/src/pages/Privacy.tsx`, add a "Delete All My Data" section after the existing content. Include a confirmation dialog (simple `window.confirm` is sufficient):

```tsx
const [purging, setPurging] = useState(false);
const [purgeResult, setPurgeResult] = useState<{ bookmarks: number; history: number } | null>(null);

const handlePurge = async () => {
  if (!window.confirm('This will permanently delete ALL your bookmarks and validation history. This cannot be undone. Continue?')) return;
  setPurging(true);
  try {
    const { data } = await sessionApi.purgeMyData();
    setPurgeResult(data.deleted);
    // Also clear IndexedDB snapshots
    const { clearAllSnapshots } = await import('../lib/bookmarkStore');
    await clearAllSnapshots();
  } catch {
    alert('Failed to purge data. Please try again.');
  } finally {
    setPurging(false);
  }
};
```

Add a `clearAllSnapshots` function to `frontend/src/lib/bookmarkStore.ts` if it doesn't exist.

**Step 3: Commit**

```bash
git add frontend/src/pages/Privacy.tsx frontend/src/services/api.ts frontend/src/lib/bookmarkStore.ts
git commit -m "feat: add Purge My Data button on Privacy page"
```

---

### Task 12: Fix Privacy Page Content

**Files:**
- Modify: `frontend/src/pages/Privacy.tsx`

**Step 1: Update inaccurate claims**

Replace the TL;DR section's false claims with accurate ones:

| Current (Wrong) | Updated (Correct) |
|---|---|
| "No data retention - structures processed and immediately discarded" | "Session-scoped storage - bookmarks and history tied to your browser session, auto-deleted after 30 days of inactivity" |
| "Local storage only - only store theme preference" | "Minimal storage - theme preference locally, bookmarks and history server-side (scoped to your session)" |
| "No cookies for tracking" | "Session cookie only - a functional cookie identifies your session (not for tracking)" |
| "No storage of chemical structures on our servers" | REMOVE this line entirely (we DO store bookmarks/history) |

Update the "What We Store" section to include:
- Session cookie (`chemaudit_sid`) — functional, not tracking
- Bookmarks (SMILES, name, tags, notes) — scoped to session
- Validation history (SMILES, score, outcome) — scoped to session, auto-purges after 30 days

Update the "Chemical Structure Processing" section to accurately describe that structures ARE logged in the audit trail for history purposes, scoped to the user's session.

Add a "Your Rights" section:
- View your data (Bookmarks page, History page)
- Delete all your data ("Purge My Data" button)
- Data auto-expires after 30 days of inactivity

**Step 2: Commit**

```bash
git add frontend/src/pages/Privacy.tsx
git commit -m "fix: update Privacy page to accurately reflect data practices"
```

---

### Task 13: Write Tests

**Files:**
- Create: `backend/tests/test_session_isolation.py`

**Step 1: Write integration tests**

```python
"""
Session Isolation Tests

Verifies that:
1. Bookmarks are scoped to session cookies
2. History is scoped to session cookies
3. API key users are scoped by api_key_hash
4. Cross-session access is blocked
5. Purge endpoint deletes all session data
"""
import pytest
from httpx import AsyncClient

# Test: Create bookmark sets session cookie
# Test: List bookmarks returns only session's bookmarks
# Test: Different session cookie sees different bookmarks
# Test: API key user sees only their bookmarks
# Test: No session/no API key returns empty list
# Test: Purge endpoint deletes all data for session
# Test: Purge endpoint clears cookie
# Test: History endpoint scoped to session
# Test: History stats scoped to session
```

Key test cases to implement:

1. **Session cookie lifecycle**: POST bookmark → response has `Set-Cookie: chemaudit_sid=...`
2. **Isolation**: Client A creates bookmark → Client B (different cookie) lists bookmarks → empty
3. **API key scoping**: Client with `X-API-Key` creates bookmark → same key sees it, different key doesn't
4. **Purge**: Create data → DELETE /me/data → verify empty
5. **No auth returns nothing**: GET /bookmarks without cookie/API key → `{"bookmarks": [], "total": 0}`

**Step 2: Run tests**

```bash
cd backend && pytest tests/test_session_isolation.py -v
```

**Step 3: Commit**

```bash
git add backend/tests/test_session_isolation.py
git commit -m "test: add session isolation integration tests"
```

---

### Task 14: End-to-End Verification

**Step 1: Start the full stack**

```bash
docker-compose up -d
```

**Step 2: Run the migration**

```bash
docker-compose exec backend alembic upgrade head
```

**Step 3: Manual verification checklist**

- [ ] Open browser A (Chrome): create a bookmark → see it in list
- [ ] Open browser B (Firefox): bookmarks list is empty
- [ ] Browser A: run validation → see it in history
- [ ] Browser B: history is empty
- [ ] Browser A: click "Purge My Data" → bookmarks and history empty
- [ ] Check `chemaudit_sid` cookie exists in browser A dev tools (HttpOnly)
- [ ] Use API key: create bookmark with `X-API-Key` header → visible only with same key

**Step 4: Run full test suite**

```bash
cd backend && pytest -x
cd frontend && npm test
```

**Step 5: Final commit**

```bash
git add -A
git commit -m "feat: complete session-based data isolation with RLS"
```

---

## Summary of All Files Changed

| File | Action | Purpose |
|------|--------|---------|
| `backend/app/db/models/bookmark.py` | Modify | Add `session_id` column |
| `backend/app/db/models/audit.py` | Modify | Add `session_id` column |
| `backend/alembic/versions/002_session_isolation.py` | Create | Migration: columns + RLS |
| `backend/app/core/session.py` | Create | Session management module |
| `backend/app/api/routes/bookmarks.py` | Modify | Scope all endpoints |
| `backend/app/api/routes/history.py` | Modify | Scope all endpoints |
| `backend/app/api/routes/validation.py` | Modify | Pass scope to audit |
| `backend/app/services/audit/service.py` | Modify | Add `session_id` param |
| `backend/app/services/batch/tasks.py` | Modify | Forward session from job metadata |
| `backend/app/api/routes/session.py` | Create | Purge endpoint |
| `backend/app/main.py` | Modify | Register session router |
| `backend/app/services/session/cleanup.py` | Create | Celery Beat cleanup task |
| `frontend/src/services/api.ts` | Modify | `withCredentials: true` + purge API |
| `frontend/src/pages/Privacy.tsx` | Modify | Fix inaccurate claims + Purge button |
| `frontend/src/lib/bookmarkStore.ts` | Modify | Add `clearAllSnapshots` |
| `backend/tests/test_session_isolation.py` | Create | Integration tests |
