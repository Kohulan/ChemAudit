---
phase: 03-batch-analytics
plan: 03
subsystem: batch-analytics
tags: [scaffold-analysis, murcko-scaffold, r-group-decomposition, shannon-entropy, rdkit, analytics]
dependency_graph:
  requires:
    - 03-01 (AnalyticsStorage, analytics schemas, run_expensive_analytics Celery task)
  provides:
    - compute_scaffold_analysis: groups molecules by Murcko scaffold SMILES, generic scaffold (double-GetScaffoldForMol), Shannon entropy, frequency distribution capped at top 50 + Other
    - compute_rgroup_decomposition: validates SMARTS, runs RGroupDecompose, returns per-molecule R-groups with unmatched count
  affects:
    - backend/app/services/batch/analytics_tasks.py (scaffold and rgroup branches now resolve import guard)
tech_stack:
  added: []
  patterns:
    - Double-GetScaffoldForMol pattern for generic scaffold (MakeScaffoldGeneric converts exocyclic doubles to singles; second call removes them)
    - Acyclic molecules grouped under empty string scaffold key
    - Shannon entropy from scaffold frequency distribution with log2 and edge-case guards
    - Frequency distribution capped at top 50 scaffolds + "Other" bucket for remainder
    - SMARTS validation via Chem.MolFromSmarts before calling RGroupDecompose
    - Both scaffold and rgroup return plain dicts (not Pydantic models) — analytics_tasks stores them directly
key_files:
  created:
    - backend/app/services/analytics/scaffold_analysis.py
    - backend/tests/test_analytics_scaffold.py
  modified:
    - backend/app/services/batch/analytics_tasks.py
decisions:
  - "Double-GetScaffoldForMol for generic scaffold: MakeScaffoldGeneric converts exocyclic =O/=N to single-bonded substituents; second GetScaffoldForMol removes these. Required to match existing app/services/scoring/scaffold.py pattern."
  - "Acyclic molecules use empty string scaffold key '' — no ring system, no scaffold, grouped together for consistent schema structure"
  - "Shannon entropy handles single-scaffold edge case (entropy=0) and empty-batch edge case separately before log2 computation"
  - "Frequency distribution capped at 50 + Other bucket: prevents unbounded dict size for very diverse batches"
  - "Both functions return plain dicts instead of Pydantic models: analytics_tasks.py stores result directly via analytics_storage.store_result which accepts dict"
  - "analytics_tasks.py updated from app.services.analytics.scaffold to scaffold_analysis and from app.services.analytics.rgroup to scaffold_analysis — both functions consolidated into one module"
patterns_established:
  - "Scaffold extraction: check ring count first, then MurckoScaffold.GetScaffoldForMol, then generic via MakeScaffoldGeneric + second GetScaffoldForMol"
  - "R-group decomposition: validate SMARTS first (Chem.MolFromSmarts), then parse all mols with indices, then RGroupDecompose, then map matched rows back to original indices"
requirements_completed: [BATCH-05, BATCH-06, BATCH-07, BATCH-08]
metrics:
  duration_seconds: 415
  completed_date: "2026-02-23"
  tasks_completed: 2
  tasks_total: 2
  files_created: 2
  files_modified: 1
  tests_added: 12
  tests_passing: 12
---

# Phase 03 Plan 03: Scaffold Analysis Service Summary

**Murcko scaffold grouping with double-GetScaffoldForMol generic scaffold, Shannon entropy diversity, frequency distribution capped at top 50, and R-group decomposition around user-specified SMARTS cores — fully tested.**

## Performance

- **Duration:** ~7 min
- **Started:** 2026-02-23T16:14:15Z
- **Completed:** 2026-02-23T16:21:04Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments

- Implemented `compute_scaffold_analysis`: extracts Murcko and generic scaffolds for every molecule in a batch, groups by scaffold SMILES, computes Shannon entropy diversity, builds frequency distribution capped at 50 + "Other"
- Implemented `compute_rgroup_decomposition`: validates SMARTS before decomposition, maps matched R-groups back to original batch indices, counts unmatched molecules
- Created 12-test suite covering scaffold grouping, generic scaffold (benzene/pyridine collapse to C1CCCCC1), entropy (0 for single scaffold, log2(N) for equal-frequency N scaffolds), frequency cap (51+ scaffolds produce "Other" bucket), R-group decomposition, invalid SMARTS, no-matches case, and empty batch

## Task Commits

Each task was committed atomically (bundled with related analytics services in the same session):

1. **Task 1: Scaffold Analysis Service** - `a987b09` (feat — scaffold_analysis.py created alongside statistics.py in a single session commit)
2. **Task 2: Scaffold Analysis Tests** - `c4726dc` (docs — test_analytics_scaffold.py added in plan 04 docs commit)

## Files Created/Modified

- `backend/app/services/analytics/scaffold_analysis.py` — Murcko scaffold decomposition, generic scaffold extraction, Shannon entropy, frequency distribution, R-group decomposition (248 lines)
- `backend/tests/test_analytics_scaffold.py` — 12 tests covering all plan requirements (248 lines)
- `backend/app/services/batch/analytics_tasks.py` — Updated scaffold and rgroup import paths from `app.services.analytics.scaffold` / `rgroup` to `scaffold_analysis`; changed `result.model_dump()` to `result` for dict-returning functions

## Decisions Made

- **Double-GetScaffoldForMol pattern**: `MakeScaffoldGeneric` converts exocyclic =O/=N to single-bonded substituents. A second `GetScaffoldForMol` call removes them and gives the true Bemis-Murcko framework. Pattern matched existing `app/services/scoring/scaffold.py` implementation.
- **Acyclic molecules under `""` key**: Molecules with no ring system have no Murcko scaffold. Grouping them under `""` keeps the data structure consistent — they appear in the scaffolds list rather than being silently dropped.
- **Plain dict returns, not Pydantic models**: `analytics_tasks.py` stores results via `analytics_storage.store_result(job_id, type, result)` which serializes to Redis. Returning dicts avoids a `.model_dump()` call and is consistent with how the route endpoint assembles the full `BatchAnalyticsResponse`.
- **analytics_tasks.py import path fix**: Plan specified `from app.services.analytics.scaffold_analysis import` but existing `analytics_tasks.py` had `from app.services.analytics.scaffold import` and `from app.services.analytics.rgroup import`. Updated to match the single consolidated module.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Updated analytics_tasks.py import paths to match scaffold_analysis module**
- **Found during:** Task 1 (verifying the key_link requirement)
- **Issue:** `analytics_tasks.py` referenced `app.services.analytics.scaffold` and `app.services.analytics.rgroup` as separate modules, but the plan creates a single `scaffold_analysis.py` with both functions. The `run_expensive_analytics` scaffold/rgroup branches would fail ImportError at runtime.
- **Fix:** Updated both import statements in `analytics_tasks.py` to `from app.services.analytics.scaffold_analysis import compute_scaffold_analysis` / `compute_rgroup_decomposition`. Also changed `result.model_dump()` to `result` since functions return plain dicts.
- **Files modified:** `backend/app/services/batch/analytics_tasks.py`
- **Verification:** Syntax check passed; import verified in isolation; all 12 scaffold tests pass.
- **Committed in:** `a987b09` (bundled with Task 1 implementation)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Required for `run_expensive_analytics` scaffold/rgroup branches to function. No scope creep — only updated pre-existing import placeholder paths to match delivered module.

## Issues Encountered

None — implementation matched plan specification exactly. The double-GetScaffoldForMol pattern was pre-documented in the existing `scoring/scaffold.py` with a GitHub reference, making implementation straightforward.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- `compute_scaffold_analysis` importable by `run_expensive_analytics` scaffold branch
- `compute_rgroup_decomposition` importable by `run_expensive_analytics` rgroup branch
- ScaffoldResult-compatible dict returned matching `analytics.py` schema
- All BATCH-05 through BATCH-08 requirements satisfied
- Ready for plan 03-05 (MMP / activity cliff analytics)

## Self-Check: PASSED

- scaffold_analysis.py: FOUND
- test_analytics_scaffold.py: FOUND
- 03-03-SUMMARY.md: FOUND
- Commit a987b09 (scaffold_analysis.py): FOUND
- Commit c4726dc (test_analytics_scaffold.py): FOUND

---
*Phase: 03-batch-analytics*
*Completed: 2026-02-23*
