# Session-Based Data Isolation Design

**Date:** 2026-02-25
**Status:** Approved
**Branch:** `features_v1`

## Problem

ChemAudit stores bookmarks and validation history in PostgreSQL with **zero data isolation**. All users on a multi-tenant SaaS deployment see each other's bookmarks and history. The `api_key_hash` column exists on both tables but is never populated. The Privacy page incorrectly claims "no data is stored."

## Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Deployment model | Multi-tenant SaaS | Strict isolation required |
| Anonymous user isolation | HttpOnly cookie with UUID v4 session ID | Industry standard; no login required |
| API key user isolation | Scoped by `api_key_hash` (existing column) | Already in schema, just unwired |
| DB-level enforcement | PostgreSQL Row-Level Security (RLS) | AWS/Crunchy Data recommend for multi-tenant SaaS |
| App-level enforcement | WHERE clauses on all queries | Readable, testable, defense-in-depth with RLS |
| Audit trail | Always-on, session-scoped | Visible only to the session that created it |
| Existing data | Purge all unscoped rows | Clean slate for pre-production app |
| Cookie consent | Lazy init (set on first bookmark/history action) | "Strictly necessary" exemption under ePrivacy Directive |
| Data retention | 30 days inactive for anonymous sessions | Auto-purge via Celery Beat |

## Architecture

### Session Token Flow

```
User visits ChemAudit
       |
       v
Any request hits middleware
  -> Read session cookie (if exists)
  -> Read X-API-Key header (if exists)
  -> Determine scope: api_key_hash > session_id > None
  -> SET LOCAL app.session_id / app.api_key_hash on DB connection
       |
       v
First state-changing request (create bookmark, run validation)
  -> If no session cookie: generate UUID v4, set HttpOnly cookie
  -> Store session_id on the created row
       |
       v
All CRUD queries scoped by session_id or api_key_hash
  -> App-level: explicit WHERE clause
  -> DB-level: RLS policy as safety net
       |
       v
"Purge My Data" (DELETE /api/v1/me/data)
  -> Deletes all rows for current session/API key
  -> Clears session cookie
```

### Cookie Configuration

```python
SESSION_COOKIE = "chemaudit_sid"
response.set_cookie(
    key=SESSION_COOKIE,
    value=str(uuid.uuid4()),
    max_age=30 * 24 * 3600,  # 30 days
    httponly=True,             # Not accessible to JS (XSS protection)
    secure=True,               # HTTPS only
    samesite="lax",            # CSRF protection
)
```

### Database Changes

#### New Column

```sql
-- Add session_id to bookmarks and audit tables
ALTER TABLE bookmarks ADD COLUMN session_id VARCHAR(36);
ALTER TABLE validation_audit_entries ADD COLUMN session_id VARCHAR(36);

-- Indexes
CREATE INDEX ix_bookmarks_session_id ON bookmarks(session_id);
CREATE INDEX ix_validation_audit_session_id ON validation_audit_entries(session_id);
```

#### Row-Level Security

```sql
-- Create application role (non-owner, subject to RLS)
CREATE ROLE chemaudit_app LOGIN PASSWORD '...';
GRANT SELECT, INSERT, UPDATE, DELETE ON bookmarks TO chemaudit_app;
GRANT SELECT, INSERT, UPDATE, DELETE ON validation_audit_entries TO chemaudit_app;

-- Enable RLS
ALTER TABLE bookmarks ENABLE ROW LEVEL SECURITY;
ALTER TABLE bookmarks FORCE ROW LEVEL SECURITY;

ALTER TABLE validation_audit_entries ENABLE ROW LEVEL SECURITY;
ALTER TABLE validation_audit_entries FORCE ROW LEVEL SECURITY;

-- Policies: row visible if session_id OR api_key_hash matches
CREATE POLICY bookmark_isolation ON bookmarks
  USING (
    session_id = NULLIF(current_setting('app.session_id', true), '')
    OR api_key_hash = NULLIF(current_setting('app.api_key_hash', true), '')
  );

CREATE POLICY audit_isolation ON validation_audit_entries
  USING (
    session_id = NULLIF(current_setting('app.session_id', true), '')
    OR api_key_hash = NULLIF(current_setting('app.api_key_hash', true), '')
  );
```

#### Data Purge

```sql
-- Migration: purge unscoped data
DELETE FROM bookmarks WHERE session_id IS NULL AND api_key_hash IS NULL;
DELETE FROM validation_audit_entries WHERE session_id IS NULL AND api_key_hash IS NULL;
```

### Backend Changes

#### 1. Session Middleware (`backend/app/core/session.py` — new file)

- `get_session_id(request)`: reads `chemaudit_sid` cookie, returns UUID or None
- `get_data_scope(request)`: returns `(session_id, api_key_hash)` — API key takes precedence
- `set_rls_context(db, session_id, api_key_hash)`: executes `SET LOCAL` statements
- `ensure_session_cookie(response, session_id)`: sets cookie if not already present

#### 2. DB Middleware / Dependency

FastAPI dependency that runs before route handlers:

```python
async def scoped_db(request: Request, response: Response, db: AsyncSession):
    session_id, api_key_hash = get_data_scope(request)
    await set_rls_context(db, session_id, api_key_hash)
    # Lazy cookie: set only if this is a state-changing endpoint
    if request.method in ("POST", "PUT", "DELETE") and session_id:
        ensure_session_cookie(response, session_id)
    return db, session_id, api_key_hash
```

#### 3. Route Updates

**Bookmarks** (`backend/app/api/routes/bookmarks.py`):
- `create_bookmark`: populate `session_id` and `api_key_hash` on the Bookmark row
- `list_bookmarks`: add `WHERE session_id = :sid OR api_key_hash = :akh`
- `get_bookmark`: add ownership check
- `update_bookmark`: add ownership check
- `delete_bookmark`: add ownership check
- `bulk_delete`: scope to session

**History** (`backend/app/api/routes/history.py`):
- `get_history`: filter by session_id/api_key_hash (remove public api_key_hash query param)
- `get_history_stats`: scope to session (not global)

**Validation** (`backend/app/api/routes/validation.py`):
- `log_validation_event`: pass `session_id` and `api_key_hash`

**Batch** (`backend/app/api/routes/batch.py`):
- Batch audit logging: pass `session_id` and `api_key_hash`

#### 4. New Endpoint: Purge My Data

```
DELETE /api/v1/me/data
```

- Requires active session cookie or API key
- Deletes all bookmarks, audit entries for the session/key
- Clears session cookie
- Returns `{"status": "purged", "deleted": {"bookmarks": N, "history": M}}`

#### 5. Celery Beat: Auto-Purge Expired Sessions

Daily task:
- Delete bookmarks where `session_id IS NOT NULL AND updated_at < NOW() - INTERVAL '30 days'`
- Delete audit entries where `session_id IS NOT NULL AND created_at < NOW() - INTERVAL '30 days'`
- Log purge counts

### Frontend Changes

#### 1. Axios Configuration

```typescript
// api.ts — enable credentials for cookie transmission
const api = axios.create({
  baseURL: API_BASE,
  withCredentials: true,  // Send HttpOnly cookie with every request
});
```

#### 2. CORS Update

```python
# main.py — allow credentials
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3002", "https://chemaudit.example.com"],
    allow_credentials=True,  # Required for cookies
    ...
)
```

**Note:** `allow_origins` cannot be `["*"]` when `allow_credentials=True`. Must list specific origins.

#### 3. Purge My Data Button

Add to Privacy page or a Settings dropdown:
- "Delete All My Data" button
- Confirmation dialog: "This will permanently delete all your bookmarks and validation history."
- Calls `DELETE /api/v1/me/data`
- On success: clear IndexedDB snapshots, show confirmation

#### 4. IndexedDB (No Change)

Bookmark snapshots in IndexedDB are already browser-scoped. No changes needed.

### Privacy Page Fix

Update `frontend/src/pages/Privacy.tsx` to accurately state:

1. **What is stored**: SMILES strings, validation scores, timestamps (for bookmarks and history)
2. **How it's scoped**: Data is tied to your browser session (cookie) or API key
3. **Retention**: Anonymous session data auto-deletes after 30 days of inactivity
4. **Your rights**: "Purge My Data" button to delete everything immediately
5. **What is NOT stored**: IP addresses, personal identifiers, browser fingerprints
6. **Shared cache**: Validation results are cached by molecule structure (InChIKey), not by user

## GDPR / Privacy Compliance

| Requirement | How We Meet It |
|-------------|----------------|
| Data minimization (Art. 5(1)(c)) | Only store SMILES, scores, timestamps. No personal identifiers. |
| Right to erasure (Art. 17) | "Purge My Data" endpoint + auto-expiry after 30 days |
| Lawful basis (Art. 6) | Legitimate interest (anonymous sessions) / Contract (API key users) |
| Cookie consent (ePrivacy) | Lazy init — cookie set only when user creates bookmark/runs validation. Strictly necessary exemption. |
| Pseudonymisation (Art. 25) | Session IDs are opaque UUIDs with no link to personal identity |
| Security (Art. 32) | RLS + app-level filtering + HttpOnly/Secure/SameSite cookies |

## Pitfalls to Avoid

1. **CORS + cookies**: Must set `withCredentials: true` in Axios and list specific origins (not `*`)
2. **RLS + connection pool**: Use `SET LOCAL` (transaction-scoped), not `SET` (connection-scoped)
3. **RLS + table owner**: App must connect as non-owner role for RLS to apply
4. **Cookie not on GET**: Session cookie is only set on state-changing requests (POST/PUT/DELETE). GET requests that need the session (list bookmarks) read from existing cookie.
5. **Race condition**: Two rapid requests without a cookie could create two sessions. Middleware must handle this (first request sets cookie, subsequent requests read it).

## Out of Scope

- User registration / login system
- Cross-device sync for anonymous sessions
- Per-tenant encryption keys
- Scoring profiles isolation (future — profiles don't contain sensitive molecule data)
- Batch permalinks isolation (intentionally public/shareable)

## References

- [AWS: Multi-tenant data isolation with PostgreSQL RLS](https://aws.amazon.com/blogs/database/multi-tenant-data-isolation-with-postgresql-row-level-security/)
- [Crunchy Data: RLS for Tenants in Postgres](https://www.crunchydata.com/blog/row-level-security-for-tenants-in-postgres)
- [Nile: Shipping multi-tenant SaaS using Postgres RLS](https://www.thenile.dev/blog/multi-tenant-rls)
- [OWASP: Session Management Cheat Sheet](https://cheatsheetseries.owasp.org/cheatsheets/Session_Management_Cheat_Sheet.html)
- [GDPR Art. 6 — Lawfulness of Processing](https://gdpr-info.eu/art-6-gdpr/)
- [GDPR Recital 26 — Anonymous Data](https://gdpr-info.eu/recitals/no-26/)
- [ComplyDog: Multi-Tenant SaaS Privacy Architecture](https://complydog.com/blog/multi-tenant-saas-privacy-data-isolation-compliance-architecture)
