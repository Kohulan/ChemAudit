---
phase: 06-export-api-workflow
plan: 09
subsystem: api
tags: [celery, email, smtp, redis, batch, notifications]

# Dependency graph
requires:
  - phase: 06-export-api-workflow
    provides: send_batch_complete_email Celery task (fully implemented but unwired)
  - phase: 06-export-api-workflow
    provides: webhook dispatch pattern in tasks.py (template for email wiring)
provides:
  - Email dispatch wired in both aggregate_batch_results and aggregate_batch_results_priority
  - NOTIFICATION_EMAIL global config setting (empty = disabled, matching WEBHOOK_URL pattern)
  - Per-job notification_email Form param on /batch/upload (overrides global)
  - Redis-backed notification email storage using batch:email:{job_id} key
  - TestEmailDispatchWiring tests verifying .delay() call and no-op when unconfigured
affects: [06-export-api-workflow, batch-processing, notifications]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - Redis key batch:email:{job_id} stores recipient at upload time, consumed at aggregation
    - Email dispatch mirrors webhook pattern exactly: try/except, fire-and-forget, conditional on setting
    - Per-job override via Form param takes precedence over global NOTIFICATION_EMAIL setting

key-files:
  created:
    - backend/tests/test_api/test_notifications.py (TestEmailDispatchWiring class appended)
  modified:
    - backend/app/core/config.py
    - backend/app/services/batch/tasks.py
    - backend/app/api/routes/batch.py

key-decisions:
  - "Redis key batch:email:{job_id} created at upload time (not job start) to survive async aggregation worker"
  - "Email dispatch in both aggregate_batch_results and aggregate_batch_results_priority for symmetry with webhook"
  - "Tests import BatchStatisticsData directly (not MagicMock) because asdict() requires real dataclass instances"

patterns-established:
  - "Notification dispatch pattern: Redis key read at aggregation, delete after reading to avoid stale data"

requirements-completed: [WORK-13]

# Metrics
duration: 4min
completed: 2026-02-24
---

# Phase 06 Plan 09: Email Dispatch Wiring Summary

**Wired send_batch_complete_email.delay() into both Celery batch aggregation tasks via Redis-backed per-job email storage with global NOTIFICATION_EMAIL fallback**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-24T21:07:07Z
- **Completed:** 2026-02-24T21:11:00Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- `NOTIFICATION_EMAIL` global setting added to config.py (empty string default = disabled, mirrors WEBHOOK_URL pattern)
- `notification_email` optional Form parameter added to `POST /batch/upload` route with global config fallback
- Email dispatch wired in `aggregate_batch_results` (default queue) after webhook block, before analytics dispatch
- Email dispatch wired in `aggregate_batch_results_priority` (priority queue) — identical pattern
- `TestEmailDispatchWiring` class with 2 tests: dispatch when Redis has email, no dispatch when Redis returns None
- WORK-13 marked `[x]` complete in REQUIREMENTS.md

## Task Commits

Each task was committed atomically:

1. **Task 1: Add NOTIFICATION_EMAIL setting and wire email dispatch in both aggregation tasks** - `1226a93` (feat)
2. **Task 2: Add email dispatch wiring test and mark WORK-13 done** - `9349aad` (test)

## Files Created/Modified
- `backend/app/core/config.py` - Added `NOTIFICATION_EMAIL: str = ""` setting after SMTP_TLS
- `backend/app/services/batch/tasks.py` - Added email dispatch blocks in both aggregation functions
- `backend/app/api/routes/batch.py` - Added `notification_email` Form param and Redis storage logic
- `backend/tests/test_api/test_notifications.py` - Appended `TestEmailDispatchWiring` class (2 tests)

## Decisions Made
- Redis key `batch:email:{job_id}` is written at upload time (in the route) rather than looked up from config at aggregation time, enabling per-job override without coupling the Celery task to the request context
- Used real `BatchStatisticsData` dataclass instance in tests (not MagicMock) because `asdict()` in the aggregation function requires a genuine dataclass — discovered during test iteration
- Deleted Redis key after reading to avoid stale notification emails on job retry scenarios

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- Test iteration needed: initial test used `MagicMock()` for stats, but `asdict(stats)` at the end of `aggregate_batch_results` requires a real dataclass instance. Fixed by importing `BatchStatisticsData` and constructing a real instance in the test helper.

## User Setup Required
None - no external service configuration required. NOTIFICATION_EMAIL defaults to empty (disabled). Users who want email notifications set `NOTIFICATION_EMAIL=user@example.com` in their `.env` file along with the existing SMTP settings.

## Next Phase Readiness
- WORK-13 complete — all 14 WORK items are now checked in REQUIREMENTS.md
- Phase 6 gap-closure complete with 9 plans (06-01 through 06-09)
- No blockers

---
*Phase: 06-export-api-workflow*
*Completed: 2026-02-24*
