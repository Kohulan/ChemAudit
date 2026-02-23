# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-23)

**Core value:** Every chemical structure submitted gets a thorough, transparent, and reproducible quality assessment — from basic validity through ML-readiness — so scientists can trust their molecular data.
**Current focus:** Phase 1 — Deep Validation (in progress)

## Current Position

Phase: 1 of 6 (Deep Validation)
Plan: 3 of 4 in current phase
Status: In progress
Last activity: 2026-02-23 — Plan 01-03 complete: 6 M1.3 structural complexity checks (hypervalent, polymer, ring strain, macrocycle, charged species, explicit H audit)

Progress: [███░░░░░░░] 15%

## Performance Metrics

**Velocity:**
- Total plans completed: 3
- Average duration: 6 min
- Total execution time: 0.3 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 1. Deep Validation | 3/4 | 18 min | 6 min |
| 2. Standardization Intelligence | 0/2 | — | — |
| 3. Batch Analytics | 0/6 | — | — |
| 4. Scoring Expansion | 0/3 | — | — |
| 5. Visualizations | 0/2 | — | — |
| 6. Export, API & Workflow | 0/3 | — | — |

**Recent Trend:**
- Last 5 plans: 01-01 (5 min), 01-02 (6 min), 01-03 (7 min)
- Trend: Stable

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
- 01-02 (M1.2): MolVS REMOVE_FRAGMENTS broad SMARTS ([#7], [#8]) need heavy-atom-count guard to avoid false positives; largest carbon fragment always = drug overriding MW < 50 heuristic
- 01-03 (M1.3): Ring strain heuristic only (size 3/4); macrocycle threshold >12 not >=12; startup check uses logger.warning not hard assert; polymer severity INFO; zwitterion = net 0 with both + and - atoms

### Pending Todos

None yet.

### Blockers/Concerns

- Phase 1: Startup assertion implemented (logger.warning) in 01-03; CI hard assert in test_deep_complexity_checks.py
- Phase 3: Redis memory ceiling under combined analytics load — profile during M3.1 setup
- Phase 6: Java availability in Docker for py2opsin (WORK-10) — validate with spike before committing full implementation

## Session Continuity

Last session: 2026-02-23
Stopped at: Completed 01-03-PLAN.md — 6 M1.3 structural complexity checks implemented and tested; 199 total validation tests passing; all 16 deep validation checks registered
Resume file: None
