---
phase: 06-export-api-workflow
plan: "06"
subsystem: backend
tags: [jpype, opsin, iupac, java, celery, audit-trail, webhook, batch, notifications]

# Dependency graph
requires:
  - phase: 06-export-api-workflow
    provides: "Audit trail service (log_batch_event), webhook task (send_webhook), OPSIN converter code, DB session factory"
provides:
  - "jpype1>=1.5.0 dependency in pyproject.toml for OPSIN JAR loading"
  - "Java JRE (default-jre-headless) and opsin-cli-2.8.0-jar-with-dependencies.jar in Docker images"
  - "Batch audit trail wiring via _log_batch_audit() async helper in tasks.py"
  - "Webhook dispatch via send_webhook.delay() in both aggregate task variants"
  - "WEBHOOK_URL and WEBHOOK_SECRET settings in config.py"
affects: [06-export-api-workflow, docker, deployment]

# Tech tracking
tech-stack:
  added: [jpype1>=1.5.0, default-jre-headless, opsin-cli-2.8.0-jar]
  patterns:
    - "asyncio.run() to call async DB functions from sync Celery tasks"
    - "Deferred imports inside Celery tasks for optional features (webhook, analytics)"
    - "fail_count = stats.errors; pass_count = total - fail_count from BatchStatisticsData"

key-files:
  created: []
  modified:
    - backend/pyproject.toml
    - backend/Dockerfile
    - backend/Dockerfile.prod
    - backend/app/core/config.py
    - backend/app/services/batch/tasks.py

key-decisions:
  - "asyncio.run() pattern for calling async log_batch_event from sync Celery task — safe because Celery workers have no existing event loop"
  - "pass_count computed as total - stats.errors (not from validation score buckets) for simplicity and accuracy"
  - "WEBHOOK_URL defaults to empty string (disabled); webhook dispatch guarded by truthiness check"
  - "Java JRE provisioned as default-jre-headless in both dev Dockerfile and prod Dockerfile.prod runtime stage"

patterns-established:
  - "Deferred DB imports inside _log_batch_audit() async helper to avoid circular imports at module load time"
  - "Both aggregate task variants (default + priority queue) get identical audit + webhook wiring"

requirements-completed: [WORK-10, WORK-12, WORK-14]

# Metrics
duration: 7min
completed: 2026-02-24
---

# Phase 06 Plan 06: Backend Integration Gap Closure Summary

**OPSIN runtime deps provisioned (jpype1 + Java JRE + opsin.jar) and batch audit trail + webhook dispatch wired into both Celery aggregator variants**

## Performance

- **Duration:** 7 min
- **Started:** 2026-02-24T20:06:14Z
- **Completed:** 2026-02-24T20:13:00Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments
- Added jpype1>=1.5.0 to pyproject.toml, resolving the missing Python-Java bridge dependency for OPSIN
- Provisioned Java JRE (default-jre-headless) and downloaded opsin-cli-2.8.0-jar to /app/data/opsin.jar in both Dockerfile and Dockerfile.prod
- Wired audit trail logging via `_log_batch_audit()` async helper (uses asyncio.run() from Celery sync context) into `aggregate_batch_results` and `aggregate_batch_results_priority`
- Wired webhook dispatch (send_webhook.delay) in both aggregator variants when WEBHOOK_URL is configured
- Added WEBHOOK_URL and WEBHOOK_SECRET settings to config.py with empty-string defaults (disabled by default)

## Task Commits

Each task was committed atomically:

1. **Task 1: Provision OPSIN Runtime Dependencies** - `42508d1` (feat)
2. **Task 2: Wire Batch Audit Trail and Webhook Dispatch** - `79bb305` (feat)

**Plan metadata:** (docs commit — created next)

## Files Created/Modified
- `backend/pyproject.toml` - Added jpype1>=1.5.0 dependency
- `backend/Dockerfile` - Added default-jre-headless to apt install, RUN to download opsin.jar
- `backend/Dockerfile.prod` - Same Java JRE + opsin.jar provisioning in runtime stage
- `backend/app/core/config.py` - Added WEBHOOK_URL and WEBHOOK_SECRET settings
- `backend/app/services/batch/tasks.py` - Added asyncio/logging imports, module-level logger, settings import, _log_batch_audit() helper, audit trail + webhook wiring in both aggregators

## Decisions Made
- Used `asyncio.run()` to call `log_batch_event` (async) from Celery task (sync): safe because Celery workers run in isolated processes with no pre-existing event loop.
- `pass_count = total - stats.errors`, `fail_count = stats.errors`: uses the `BatchStatisticsData.errors` field which accurately tracks all error-status molecules, avoiding the need to inspect validation scores per-molecule.
- `WEBHOOK_URL` defaults to `""` (empty = disabled); guarded by `if webhook_url:` to skip dispatch entirely when not configured.
- Both `aggregate_batch_results` (default queue) and `aggregate_batch_results_priority` (high_priority queue) receive identical wiring — ensures consistent behavior for small and large batch jobs.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 2 - Missing Critical] Added WEBHOOK_URL and WEBHOOK_SECRET to config.py**
- **Found during:** Task 2 (Wire Batch Audit Trail and Webhook Dispatch)
- **Issue:** Plan's webhook dispatch code referenced `settings.WEBHOOK_URL` and `settings.WEBHOOK_SECRET`, but these fields did not exist in `config.py`. Without them, `tasks.py` would raise `AttributeError` at runtime.
- **Fix:** Added `WEBHOOK_URL: str = ""` and `WEBHOOK_SECRET: str = ""` to `Settings` class in `config.py`.
- **Files modified:** `backend/app/core/config.py`
- **Verification:** `ruff check app/core/config.py` passes; settings instantiates without error.
- **Committed in:** `79bb305` (Task 2 commit)

**2. [Rule 2 - Missing Critical] Replaced plan's warn_count parameter with actual audit service signature**
- **Found during:** Task 2 (Wire Batch Audit Trail and Webhook Dispatch)
- **Issue:** Plan's `log_batch_event()` call snippet included `warn_count=warn_count`, but the actual `log_batch_event` in `audit/service.py` does not accept a `warn_count` parameter. Also, the plan's call omitted the required `db: AsyncSession` parameter (the service requires a DB session, not standalone DB access).
- **Fix:** Created `_log_batch_audit()` async helper that creates its own DB session via `async_session()` and calls `log_batch_event()` with the correct signature (db, job_id, molecule_count, pass_count, fail_count).
- **Files modified:** `backend/app/services/batch/tasks.py`
- **Verification:** `ruff check` passes; function signature matches audit service definition exactly.
- **Committed in:** `79bb305` (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (2 missing critical)
**Impact on plan:** Both auto-fixes essential for runtime correctness. No scope creep.

## Issues Encountered
- The `log_batch_event` audit service function requires a `db: AsyncSession` parameter that was not mentioned in the plan's code snippet. Created the `_log_batch_audit()` async wrapper pattern to handle DB session lifecycle correctly from the Celery sync context.

## User Setup Required
To enable webhook notifications, set these environment variables:
- `WEBHOOK_URL` — HTTPS endpoint to receive batch completion POSTs
- `WEBHOOK_SECRET` — HMAC signing secret (for `X-Signature: sha256=...` header verification)

Leave both empty (defaults) to disable webhook dispatch.

## Next Phase Readiness
- All three UAT blockers/gaps closed: IUPAC 503 blocker fixed (OPSIN provisioned), batch audit trail wired, webhook dispatch wired
- Docker rebuild required for OPSIN/Java changes to take effect (`docker-compose build backend`)
- IUPAC conversion ready end-to-end once container is rebuilt with Java JRE + opsin.jar

---
*Phase: 06-export-api-workflow*
*Completed: 2026-02-24*
