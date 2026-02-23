---
phase: 03-batch-analytics
plan: 05
subsystem: analytics
tags: [rdkit, brics, mmp, activity-cliffs, sali, lle, cheminformatics]

requires:
  - phase: 03-01
    provides: AnalyticsStorage, run_expensive_analytics task, analytics schemas (MMPPair, ActivityCliff, MMPResult)

provides:
  - BRICS-based MMP detection with deduplication and Tanimoto ranking
  - SALI activity cliff detection requiring user-supplied activity column
  - LLE computation (pIC50 - LogP) per molecule
  - compute_mmp_analysis public API + compute_mmp alias for analytics_tasks.py

affects:
  - 03-01 (analytics_tasks.py already wired to import compute_mmp)
  - frontend MMP analytics display (future)

tech-stack:
  added: []
  patterns:
    - "BRICS decomposition for MMP: BRICSDecompose → shared fragment = core, different fragments = R-groups"
    - "Size heuristic for MMP noise reduction: shared fragment must have more heavy atoms than variable fragment"
    - "Pair deduplication: best_pairs dict keyed by (min_idx, max_idx) — highest Tanimoto retained"
    - "Cap applied both inside _detect_mmp_pairs AND in compute_mmp_analysis to survive monkeypatching"
    - "compute_mmp alias wraps compute_mmp_analysis in a _MMPResultWrapper with model_dump() for tasks compatibility"

key-files:
  created:
    - backend/app/services/analytics/mmp.py
    - backend/tests/test_analytics_mmp.py
  modified: []

key-decisions:
  - "MMP cap applied at two levels (inside _detect_mmp_pairs and in compute_mmp_analysis) to remain effective when _detect_mmp_pairs is monkeypatched in tests"
  - "compute_mmp alias returns a wrapper object with model_dump() for seamless compatibility with existing analytics_tasks.py that calls result.model_dump()"
  - "Phenol/aniline are NOT valid BRICS pairs (no BRICS-breakable bonds); tests use phenylacetic acid / phenylacetamide which share [16*]c1ccccc1 via genuine BRICS decomposition"
  - "SALI pairs where tanimoto >= 1.0 are silently skipped (denominator zero) rather than raising an exception"
  - "activity_column=None → activity_cliffs=None, lle_values=None (no crash, graceful degradation)"

patterns-established:
  - "test_mmp_pairs_limited: monkeypatches _detect_mmp_pairs to inject 1500 synthetic pairs, verifying the orchestrator-level cap"
  - "test_mmp_activity_cliff_sali_formula: tests internal _compute_activity_cliffs with synthetic pair for precise arithmetic"

requirements-completed: [BATCH-13, BATCH-14, BATCH-15]

duration: 8min
completed: 2026-02-23
---

# Phase 3 Plan 05: MMP, Activity Cliff, and LLE Analytics Summary

**BRICS-based matched molecular pair detection with SALI activity cliff scoring and pIC50 - LogP ligand efficiency computation, wired into analytics_tasks.py via a model_dump()-compatible wrapper**

## Performance

- **Duration:** 8 min
- **Started:** 2026-02-23T16:14:25Z
- **Completed:** 2026-02-23T16:22:00Z
- **Tasks:** 2 of 2
- **Files modified:** 2

## Accomplishments

- Implemented `compute_mmp_analysis` with BRICS decomposition: finds molecule pairs sharing a common core fragment while differing in R-group, with a size heuristic to filter noise (core must have more heavy atoms than R-group)
- Pair deduplication ensures each (mol_a, mol_b) appears once, Tanimoto-sorted and capped at 1000 pairs; batch refusal at 5000 molecules
- SALI activity cliff detection (`|delta_activity| / (1 - tanimoto)`) and LLE computation (`activity - LogP`) when activity_column is provided; graceful None returns when absent
- 14 passing tests covering pair detection, refusal, cap, deduplication, activity cliffs, LLE, edge cases, and SALI formula precision

## Task Commits

1. **Task 1: MMP, Activity Cliff, and LLE Service** - `d6c4b1b` (feat)
2. **Task 2: MMP, Activity Cliff, and LLE Tests + cap fix** - `13b12e3` (test)

**Plan metadata:** (docs commit — see below)

## Files Created/Modified

- `backend/app/services/analytics/mmp.py` - Full MMP analytics service: `_detect_mmp_pairs`, `_compute_activity_cliffs`, `_compute_lle`, `compute_mmp_analysis`, `compute_mmp` alias
- `backend/tests/test_analytics_mmp.py` - 14 tests covering all plan requirements (360 lines)

## Decisions Made

- **Double-capping pairs:** The MAX_PAIRS_RETURNED cap is applied inside `_detect_mmp_pairs` (for normal operation) AND in `compute_mmp_analysis` (so monkeypatched implementations in tests are also bounded). This discovered a bug during test writing that was fixed as a Rule 1 auto-fix.
- **compute_mmp wrapper:** `analytics_tasks.py` already imports `compute_mmp` and calls `.model_dump()` on the result. A thin `_MMPResultWrapper` class was added to make the dict returned by `compute_mmp_analysis` appear as a Pydantic-like model without modifying the tasks file.
- **Test molecules:** The plan suggests phenol/aniline as test molecules, but these have no BRICS-breakable bonds and produce no MMP pairs. Used phenylacetic acid (`c1ccccc1CC(=O)O`) and phenylacetamide (`c1ccccc1CC(=O)N`) instead — they share `[16*]c1ccccc1` as a genuine BRICS core.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Cap applied only inside _detect_mmp_pairs — ineffective when monkeypatched**
- **Found during:** Task 2 (MMP pair limit test)
- **Issue:** `test_mmp_pairs_limited` injects 1500 synthetic pairs via monkeypatch on `_detect_mmp_pairs`. The cap code was inside that function, so the replacement bypassed it, failing the assertion (`got 1500, expected <= 1000`).
- **Fix:** Added a second `[:MAX_PAIRS_RETURNED]` slice in `compute_mmp_analysis` after calling `_detect_mmp_pairs`, so the cap applies regardless of how pair detection is implemented.
- **Files modified:** `backend/app/services/analytics/mmp.py`
- **Verification:** All 14 tests pass
- **Committed in:** `13b12e3` (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (Rule 1 - Bug)
**Impact on plan:** The fix strengthens the cap contract — it is now enforced by the orchestrator, not only by the inner function. No scope creep.

## Issues Encountered

- Phenol/aniline (the plan's "e.g." example for MMP test) have no BRICS-breakable bonds, so `BRICSDecompose` returns only the whole molecule for each and they share no fragment. Investigated and confirmed: the plan's `e.g.` is non-normative. Used phenylacetic acid / phenylacetamide which are genuine BRICS MMP pairs sharing the phenyl ring as core.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- MMP analytics fully wired: `analytics_tasks.py` `elif analysis_type == "mmp"` branch calls `compute_mmp` which returns a `model_dump()`-compatible wrapper
- Activity cliff and LLE are guarded by `activity_column` parameter — safe to call without activity data
- All 14 tests green; service importable without Redis/Celery

---
*Phase: 03-batch-analytics*
*Completed: 2026-02-23*
