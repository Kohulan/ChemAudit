---
phase: 03-batch-analytics
plan: 01
subsystem: batch-analytics
tags: [infra, redis, celery, analytics, schemas, pydantic]
dependency_graph:
  requires: []
  provides:
    - BATCH_RESULT_TTL config setting (24h Redis TTL for batch results)
    - ResultStorage.get_all_results() for full dataset access
    - AnalyticsStorage singleton (store/get/update/init_status)
    - Analytics Pydantic schemas (BatchAnalyticsResponse + 15 sub-models)
    - GET /batch/{job_id}/analytics endpoint
    - POST /batch/{job_id}/analytics/{analysis_type} endpoint
    - run_cheap_analytics Celery task (auto-dispatched on aggregation)
    - run_expensive_analytics Celery task (user-triggered per type)
  affects:
    - backend/app/services/batch/tasks.py (dispatch wiring)
    - backend/app/api/routes/batch.py (new endpoints)
tech_stack:
  added:
    - fakeredis (via MagicMock for test isolation)
  patterns:
    - Lazy Redis singleton (AnalyticsStorage mirrors ResultStorage pattern)
    - Import-guarded analytics modules in Celery tasks (graceful skipping)
    - Defensive empty-results check in all analytics tasks
key_files:
  created:
    - backend/app/services/analytics/__init__.py
    - backend/app/services/analytics/storage.py
    - backend/app/services/batch/analytics_tasks.py
    - backend/app/schemas/analytics.py
    - backend/tests/test_infra_analytics.py
  modified:
    - backend/app/core/config.py (BATCH_RESULT_TTL setting)
    - backend/app/services/batch/result_aggregator.py (TTL from config, get_all_results)
    - backend/app/api/routes/batch.py (analytics GET + POST endpoints)
    - backend/app/services/batch/tasks.py (dispatch run_cheap_analytics in both aggregators)
decisions:
  - BATCH_RESULT_TTL=86400 in config so all subsequent analytics plans inherit automatically
  - Import-guarded analytics module calls in run_cheap_analytics allow incremental delivery of plans 03-02 through 03-06 without breaking the dispatcher
  - analytics_storage key format batch:analytics:{type}:{job_id} keeps type and job scoped separately for easy cleanup
  - run_cheap_analytics dispatched via lazy import inside try/except to prevent aggregation failure if analytics module fails to load
metrics:
  duration_seconds: 228
  completed_date: "2026-02-23"
  tasks_completed: 2
  tasks_total: 2
  files_created: 5
  files_modified: 4
  tests_added: 14
  tests_passing: 14
---

# Phase 03 Plan 01: Analytics Infrastructure Skeleton Summary

**One-liner:** 24h TTL policy (INFRA-01) and full analytics infrastructure — storage, schemas, Celery task dispatcher, and GET/POST endpoints — wired and tested as the foundation for plans 03-02 through 03-06.

## What Was Built

### Task 1: INFRA-01 TTL Policy and get_all_results

- Added `BATCH_RESULT_TTL: int = 86400` to `config.py` Caching section with comment referencing INFRA-01.
- Changed `ResultStorage.RESULT_EXPIRY` from hardcoded `3600` to `settings.BATCH_RESULT_TTL`, extending batch result lifetime from 1h to 24h.
- Added `get_all_results(job_id) -> list[dict]` to `ResultStorage` — returns full raw results list for analytics consumption without pagination overhead.

### Task 2: Analytics Infrastructure Skeleton

**Storage (`backend/app/services/analytics/storage.py`):**
- `AnalyticsStorage` class with lazy Redis singleton pattern.
- Key schema: `batch:analytics:{analysis_type}:{job_id}` for results, `batch:analytics:status:{job_id}` for status dict.
- Methods: `store_result`, `get_result`, `get_status`, `update_status`, `init_status`.
- `ANALYTICS_TTL = 86400` (24h). Singleton: `analytics_storage = AnalyticsStorage()`.

**Schemas (`backend/app/schemas/analytics.py`):**
- 18 Pydantic v2 models covering all analytics types.
- `AnalysisStatus`, `DeduplicationGroup/Result`, `ScaffoldGroup/Result`, `RGroupResult`, `ChemSpaceCoordinates`, `SimilarityHit/Result/NearestNeighborResult/SimilarityMatrixResult`, `MMPPair/ActivityCliff/MMPResult`, `PropertyStats/OutlierInfo/PropertyCorrelation/QualityScore/StatisticsResult`.
- Top-level: `BatchAnalyticsResponse` (all analytics optional), `AnalyticsTriggerResponse`.

**Celery Tasks (`backend/app/services/batch/analytics_tasks.py`):**
- `run_cheap_analytics(job_id)` — auto-dispatched; computes deduplication + statistics; import-guarded (skips gracefully if service modules don't exist yet).
- `run_expensive_analytics(job_id, analysis_type, params)` — user-triggered; routes scaffold/chemical_space/mmp/similarity_search/rgroup to service modules; defensive empty-results check.

**Routes (`backend/app/api/routes/batch.py`):**
- `GET /batch/{job_id}/analytics` (30/min rate limit) — returns `BatchAnalyticsResponse` with status dict + available results; 404 if analytics not initialized.
- `POST /batch/{job_id}/analytics/{analysis_type}` (10/min rate limit) — validates allowed set, dispatches `run_expensive_analytics.delay()`, returns `AnalyticsTriggerResponse(status="queued")`.

**Task Dispatch Wiring (`backend/app/services/batch/tasks.py`):**
- Both `aggregate_batch_results` and `aggregate_batch_results_priority` now dispatch `run_cheap_analytics.delay(job_id)` after `progress_tracker.mark_complete(job_id)`.
- Dispatch wrapped in `try/except` so analytics failure cannot block aggregation.

**Tests (`backend/tests/test_infra_analytics.py`):**
- 14 tests covering: TTL config, `get_all_results` (missing + deserialize), `AnalyticsStorage` store/get/round-trip/init_status/update_status/error, schema instantiation, 404 endpoint path.
- All 14 tests pass.

## Deviations from Plan

None — plan executed exactly as written.

## Verification Results

All plan verification criteria confirmed:

1. `settings.BATCH_RESULT_TTL` returns `86400` — OK
2. `ResultStorage.RESULT_EXPIRY` reads from `settings.BATCH_RESULT_TTL` (not hardcoded) — OK
3. `get_all_results()` method exists and returns `list[dict]` — OK
4. `analytics_storage` singleton can store/retrieve JSON data — OK
5. `BatchAnalyticsResponse` schema validates with all optional fields — OK
6. `GET /batch/{job_id}/analytics` returns 404 for unknown jobs — OK
7. `POST /batch/{job_id}/analytics/scaffold` registered (returns 200 queued in production) — OK
8. Both `aggregate_batch_results` and `aggregate_batch_results_priority` dispatch cheap analytics — OK

## Self-Check: PASSED
