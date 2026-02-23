---
phase: 01-deep-validation
plan: 03
subsystem: validation
tags: [rdkit, structural-complexity, hypervalent, polymer, ring-strain, macrocycle, zwitterion, explicit-hydrogen]

# Dependency graph
requires:
  - phase: 01-deep-validation/01-01
    provides: "BaseCheck, CheckRegistry decorator pattern, deep_stereo_tautomer checks"
  - phase: 01-deep-validation/01-02
    provides: "deep_composition checks, __init__.py export pattern"
provides:
  - "6 new structural complexity checks registered in CheckRegistry"
  - "hypervalent_atoms check (DVAL-12): flags atoms exceeding max allowed valence"
  - "polymer_detection check (DVAL-13): SGroup markers, MW > 1500 Da, dummy atoms"
  - "ring_strain check (DVAL-14): 3/4-membered ring heuristic detection"
  - "macrocycle_detection check (DVAL-15): rings > 12 atoms with SSSR note"
  - "charged_species check (DVAL-16): net charge, zwitterion detection"
  - "explicit_hydrogen_audit check (DVAL-17): explicit H count and H atom objects"
  - "EXPECTED_DEEP_VALIDATION_CHECKS constant in main.py with logger.warning startup check"
  - "65 tests covering all 6 M1.3 checks plus CI hard-assert for all 16 deep validation checks"
affects: [scoring, batch-analytics, ml-readiness, export]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Structural complexity checks use ring size heuristic (not force field) for ring strain"
    - "Macrocycle threshold is >12 (not >=12) atoms per SSSR rings"
    - "Zwitterion = net charge 0 AND both positive and negative atoms present"
    - "Startup warning via logger.warning for gradual deployment; CI hard assert in test file"
    - "Polymer detection uses three heuristics: SGroup markers, MW threshold, dummy atoms"

key-files:
  created:
    - backend/app/services/validation/checks/deep_complexity.py
    - backend/tests/test_validation/test_deep_complexity_checks.py
  modified:
    - backend/app/services/validation/checks/__init__.py
    - backend/app/services/validation/engine.py
    - backend/app/main.py

key-decisions:
  - "Ring strain uses heuristic (ring size 3/4) only — no force field calculation per research decision"
  - "Macrocycle threshold is >12 atoms (not >=12) per plan spec"
  - "Startup check uses logger.warning not hard assert to allow partial deployments; CI test has hard assert"
  - "Polymer severity is INFO (macrolides, natural products valid); ring strain is WARNING (stability concern)"
  - "ExplicitHydrogenAudit uses GetNumExplicitHs() for H on heavy atoms + counts H atom objects from AddHs()"

patterns-established:
  - "All structural_complexity checks return INFO or WARNING severity (never CRITICAL/ERROR for passing molecules)"
  - "None molecule check at top of run() returns Severity.ERROR"
  - "details dict always includes all keys even when check passes (empty lists, zero counts)"
  - "affected_atoms is union of all flagged atom indices (sorted for determinism)"

requirements-completed: [DVAL-12, DVAL-13, DVAL-14, DVAL-15, DVAL-16, DVAL-17]

# Metrics
duration: 7min
completed: 2026-02-23
---

# Phase 01 Plan 03: Deep Validation M1.3 Structural Complexity Summary

**6 structural complexity checks (hypervalent atoms, polymer detection, ring strain, macrocycle, charged species/zwitterion, explicit hydrogen audit) registered in CheckRegistry with 65 passing tests and a startup assertion covering all 16 deep validation checks.**

## Performance

- **Duration:** 7 minutes
- **Started:** 2026-02-23T10:08:58Z
- **Completed:** 2026-02-23T10:16:23Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments

- Implemented 6 new check classes in `deep_complexity.py` registered via `@CheckRegistry.register()` decorator
- `engine.py` now imports all three deep check modules (`deep_stereo_tautomer`, `deep_composition`, `deep_complexity`)
- `main.py` lifespan adds startup check with `EXPECTED_DEEP_VALIDATION_CHECKS` covering all 16 DVAL checks
- 65 comprehensive tests covering all 6 M1.3 checks: normal molecules pass, flag conditions detected, details structure verified, None handling
- CI cross-check (`TestAllDeepValidationChecks`) hard-asserts all 16 deep validation checks registered

## Task Commits

Each task was committed atomically:

1. **Task 1: Implement 6 structural complexity checks + startup assertion** - `ab0eb4c` (feat)
2. **Task 2: Write comprehensive tests for M1.3 checks + CI registration test** - (committed as part of docs commit `6da19e7` from prior plan execution)

## Files Created/Modified

- `backend/app/services/validation/checks/deep_complexity.py` - 6 new check classes: HypervalentAtomCheck, PolymerDetectionCheck, RingStrainCheck, MacrocycleDetectionCheck, ChargedSpeciesCheck, ExplicitHydrogenAuditCheck
- `backend/tests/test_validation/test_deep_complexity_checks.py` - 65 tests across 8 test classes including CI cross-check
- `backend/app/services/validation/checks/__init__.py` - Added deep_complexity imports and __all__ exports
- `backend/app/services/validation/engine.py` - Added deep_complexity and deep_composition imports
- `backend/app/main.py` - Added EXPECTED_DEEP_VALIDATION_CHECKS constant and logger.warning startup check in lifespan

## Decisions Made

- **Ring strain heuristic only**: 3/4-membered ring size check, no force field energy calculation (per research decision in 01-RESEARCH.md)
- **Macrocycle threshold >12**: rings with more than 12 atoms are flagged, exactly 12-atom rings pass (boundary condition)
- **Startup check uses logger.warning**: allows partial deployment when not all plans are yet executed; CI test file has the hard assert
- **Polymer severity INFO**: high-MW natural products and macrolide drugs are valid compounds; polymer flag is informational
- **Ring strain severity WARNING**: strain can affect stability and is more actionable than INFO

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Added missing deep_composition import to engine.py**
- **Found during:** Task 1 (engine.py update)
- **Issue:** Plan 02's engine.py import for `deep_composition` was missing — composition checks would not register without it
- **Fix:** Added `import app.services.validation.checks.deep_composition  # noqa: F401` to engine.py alongside the deep_complexity import
- **Files modified:** backend/app/services/validation/engine.py
- **Verification:** Verified all 16 deep validation checks register when engine.py is imported
- **Committed in:** ab0eb4c (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (Rule 3 - blocking import missing from prior plan)
**Impact on plan:** Auto-fix necessary for all 16 checks to be discoverable. No scope creep.

## Issues Encountered

- The test file `test_deep_complexity_checks.py` was pre-committed in a prior docs commit (`6da19e7`) from plan 01's execution — the test file content was already there when Task 2 began. No re-commit was needed.
- The venv required `pip install -e ".[dev]"` before tests could run (slowapi not installed in system Python).

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- All 16 deep validation checks (DVAL-01 through DVAL-17) from Phase 1 are now registered and tested
- Startup assertion in `main.py` confirms all checks present at runtime
- Phase 2 (Standardization Intelligence) can proceed; the validation engine is complete
- ValidationEngine singleton already includes all M1.3 checks for immediate use in API routes

## Self-Check: PASSED

- deep_complexity.py: FOUND
- test_deep_complexity_checks.py: FOUND
- 01-03-SUMMARY.md: FOUND
- Commit ab0eb4c: FOUND
- Commit ba1e3ef (docs/final): FOUND
- All 199 validation tests: PASSING

---
*Phase: 01-deep-validation*
*Completed: 2026-02-23*
