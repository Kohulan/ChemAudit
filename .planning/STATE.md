# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-23)

**Core value:** Every chemical structure submitted gets a thorough, transparent, and reproducible quality assessment — from basic validity through ML-readiness — so scientists can trust their molecular data.
**Current focus:** Phase 1 — Deep Validation (in progress)

## Current Position

Phase: 1 of 6 (Deep Validation)
Plan: 1 of 3 in current phase
Status: In progress
Last activity: 2026-02-23 — Plan 01-01 complete: 4 M1.1 deep validation checks + CHECKS_VERSION=v2 cache key versioning

Progress: [█░░░░░░░░░] 5%

## Performance Metrics

**Velocity:**
- Total plans completed: 1
- Average duration: 5 min
- Total execution time: 0.1 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 1. Deep Validation | 1/3 | 5 min | 5 min |
| 2. Standardization Intelligence | 0/2 | — | — |
| 3. Batch Analytics | 0/6 | — | — |
| 4. Scoring Expansion | 0/3 | — | — |
| 5. Visualizations | 0/2 | — | — |
| 6. Export, API & Workflow | 0/3 | — | — |

**Recent Trend:**
- Last 5 plans: 01-01 (5 min)
- Trend: —

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Roadmap: Follow recommended build order from Features/ROADMAP.md (M1.1 → M1.2 → M1.3 → M2.1 → M3.1 → M4.3 → M3.2 → M3.3 → M3.5 → M4.1 → M4.2 → M5.1 → M5.2 → M3.4 → M2.2 → M6.1 → M6.2 → M6.3)
- Phase 3: INFRA-01 (batch pagination + Redis TTL) is plan 03-01 and must ship before any analytics endpoint
- Tautomer canonicalization: Document stereo stripping in all provenance and dedup responses (stereo_stripped flag)
- Validation cache: CHECKS_VERSION=v2 added to cache key in 01-01; key format is now validation:v2:{inchikey}:{checks_hash}
- Batch analytics: Run as post-aggregation Celery chord; frontend polls /batch/{job_id}/analytics separately from batch completion
- 01-01 (M1.1): TautomerDetection uses enumerator.Canonicalize() (GetCanonicalTautomer does not exist); stereoisomer cap=128; multi-fragment checks run on largest fragment only

### Pending Todos

None yet.

### Blockers/Concerns

- Phase 1: Startup assertion for CheckRegistry deferred — will address in a future plan (M1.1 complete without it)
- Phase 3: Redis memory ceiling under combined analytics load — profile during M3.1 setup
- Phase 6: Java availability in Docker for py2opsin (WORK-10) — validate with spike before committing full implementation

## Session Continuity

Last session: 2026-02-23
Stopped at: Completed 01-01-PLAN.md — 4 M1.1 deep stereo/tautomer checks implemented and tested; CHECKS_VERSION=v2 cache key versioning complete
Resume file: None
