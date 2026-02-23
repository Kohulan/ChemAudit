---
phase: 03-batch-analytics
plan: 06
subsystem: analytics
tags: [numpy, rdkit, statistics, batch-analytics, outlier-detection, quality-score, pearson-correlation]

# Dependency graph
requires:
  - phase: 03-01
    provides: "analytics infrastructure skeleton: AnalyticsStorage, StatisticsResult schema, run_cheap_analytics task stub"

provides:
  - "compute_all_statistics: property stats + correlations + IQR outliers + quality score"
  - "compute_property_stats: numpy-based mean/median/std/quartiles/IQR/min/max per property"
  - "compute_correlations: pairwise Pearson with NaN guard and <10 data point exclusion"
  - "compute_outliers: IQR fence method [Q1-1.5*IQR, Q3+1.5*IQR] per property"
  - "compute_quality_score: composite 0-100 metric (40% validity, 35% scaffold entropy, 25% Lipinski)"
  - "21 tests covering all functions, edge cases, and empty batch behavior"

affects:
  - "analytics endpoint (GET /batch/{job_id}/analytics) — statistics field now populated"
  - "run_cheap_analytics task — statistics branch no longer skipped"

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "PROPERTY_EXTRACTORS dict maps name -> lambda for clean property extraction from nested dicts"
    - "compute_all_statistics returns Pydantic StatisticsResult model (not plain dict) because analytics_tasks.py calls .model_dump()"
    - "Near-zero threshold abs < 1e-10 filters floating-point noise from constant arrays in corrcoef"
    - "Scaffold diversity computed inline via MurckoScaffold (no dependency on scaffold_analysis service)"
    - "Lipinski fallback to 50% neutral when no druglikeness data available"

key-files:
  created:
    - backend/app/services/analytics/statistics.py
    - backend/tests/test_analytics_statistics.py
  modified: []

key-decisions:
  - "compute_all_statistics returns StatisticsResult Pydantic model (not plain dict) because analytics_tasks.py calls .model_dump() on the result — plan said dict but tasks code was authoritative"
  - "abs(pearson_r) < 1e-10 threshold to filter floating-point noise from constant arrays — numpy corrcoef returns -3.4e-16 (not NaN) for constant vs varying arrays on this numpy version"
  - "Scaffold diversity computed inline with MurckoScaffold — plan explicitly required no dependency on scaffold_analysis service"

patterns-established:
  - "Property extractor lambdas: safe nested .get() chains handle missing/None scoring blocks gracefully"
  - "IQR outlier method: skip < 4 values (statistical minimum), correlation: skip < 10 values"
  - "Quality score components rounded to 1 decimal place before composite calculation"

requirements-completed: [BATCH-16, BATCH-17, BATCH-18, BATCH-19]

# Metrics
duration: 4min
completed: 2026-02-23
---

# Phase 3 Plan 06: Batch Statistics Summary

**numpy-based property statistics (mean/median/std/IQR), Pearson correlation matrix, IQR outlier detection, and composite quality score (40% validity + 35% scaffold entropy + 25% Lipinski) — all wired into run_cheap_analytics**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-23T16:14:30Z
- **Completed:** 2026-02-23T16:18:30Z
- **Tasks:** 2
- **Files created:** 2

## Accomplishments

- Implemented `statistics.py` with 5 public functions covering property stats, correlations, outlier detection, and composite quality scoring
- 21 tests pass covering all 4 computation paths, edge cases (insufficient data, constant arrays, empty batch), and integration via `compute_all_statistics`
- `run_cheap_analytics` statistics branch now fully operational — no longer skipped on import

## Task Commits

1. **Task 1: Statistics, Quality Score, and Outlier Service** - `a987b09` (feat)
2. **Task 2: Statistics and Quality Score Tests** - `b4a6b62` (test)

**Plan metadata:** (docs commit follows)

## Files Created/Modified

- `backend/app/services/analytics/statistics.py` - Property stats, correlations, outlier detection, quality score service; exports `compute_all_statistics`
- `backend/tests/test_analytics_statistics.py` - 21 tests covering all functions and edge cases

## Decisions Made

1. `compute_all_statistics` returns a `StatisticsResult` Pydantic model instead of a plain dict — `analytics_tasks.py` calls `.model_dump()` on the result, making the Pydantic model the correct return type (the plan specified "dict matching StatisticsResult schema" but the task code was authoritative).

2. Near-zero threshold `abs(pearson_r) < 1e-10` used instead of strict `== 0.0` check — numpy's `corrcoef` for constant arrays returns floating-point noise (`-3.4e-16`, not `0.0` or `NaN`) on this numpy version. The threshold correctly filters these artifacts while preserving genuine near-zero correlations.

3. Scaffold diversity computed inline with `MurckoScaffold.GetScaffoldForMol` — plan explicitly prohibited dependency on the scaffold_analysis service. Shannon entropy of scaffold frequency distribution normalized by log2(unique_scaffolds).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Added abs < 1e-10 threshold for floating-point noise in correlations**
- **Found during:** Task 2 (test verification)
- **Issue:** numpy `corrcoef` returns `-3.4e-16` (not `NaN` or exactly `0.0`) for constant arrays, causing the `== 0.0` skip check to fail — constant-value properties appeared in correlation results with spurious non-zero values
- **Fix:** Changed `if pearson_r == 0.0: continue` to `if abs(pearson_r) < 1e-10: continue` in `compute_correlations`; also updated test assertion to verify this behavior
- **Files modified:** `backend/app/services/analytics/statistics.py`, `backend/tests/test_analytics_statistics.py`
- **Verification:** Test `test_constant_qed_excluded_from_correlations` passes; all 21 tests pass
- **Committed in:** `b4a6b62` (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (Rule 1 - Bug)
**Impact on plan:** Fix was necessary for correct behavior — constant-value properties should not produce correlation pairs. No scope creep.

## Issues Encountered

- numpy `corrcoef` floating-point noise: for constant arrays, returns ~1e-16 (not NaN or 0.0). Fixed by using absolute threshold. Standard numpy warning about `invalid value encountered in divide` is expected and harmless for this edge case.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- `compute_all_statistics` ready for import by `run_cheap_analytics` — statistics branch now active
- BATCH-16 through BATCH-19 fully implemented and tested
- Phase 03 plans 02-05 (deduplication, scaffold, chemical space, MMP) can proceed independently
- Quality score can be displayed in frontend batch analytics panel alongside deduplication results

## Self-Check: PASSED

- FOUND: `backend/app/services/analytics/statistics.py`
- FOUND: `backend/tests/test_analytics_statistics.py`
- FOUND: `.planning/phases/03-batch-analytics/03-06-SUMMARY.md`
- FOUND: commit `a987b09` (Task 1 — statistics service)
- FOUND: commit `b4a6b62` (Task 2 — 21 tests)
- 21 tests collected and passing

---
*Phase: 03-batch-analytics*
*Completed: 2026-02-23*
