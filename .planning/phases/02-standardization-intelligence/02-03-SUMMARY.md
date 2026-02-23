---
phase: 02-standardization-intelligence
plan: "03"
subsystem: api
tags: [rdkit, fragment-dict, provenance, dval, stereo, standardization]

# Dependency graph
requires:
  - phase: 02-standardization-intelligence
    provides: ProvenancePipeline, ProvStageRecord, StereoProvenance schemas and stereo per-center tracking (STD-05, STD-06)
provides:
  - Deduplicated COUNTERION_NAMES: O=C(O)O=carbonic acid, O=CO=formic acid, 55 total unique entries
  - DVAL cross-reference population in stereo provenance (dval_cross_refs on StereoProvenance)
  - DVAL cross-reference population on tautomer stage ProvStageRecord
  - dval_results optional field on StandardizationOptions Pydantic schema
affects:
  - 03-batch-analytics
  - frontend provenance UI (dval_cross_refs now populated in API responses)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - Fragment dict deduplication guard via regex key-count test prevents future silent overwrite
    - DVAL cross-ref population: optional dict param flows route->service, populated post-stage-build
    - Backward-compatible optional field pattern: None default preserves existing API behavior

key-files:
  created: []
  modified:
    - backend/app/services/standardization/fragment_dict.py
    - backend/app/services/standardization/provenance.py
    - backend/app/schemas/standardization.py
    - backend/app/api/routes/standardization.py
    - backend/tests/test_standardization/test_provenance.py

key-decisions:
  - "Fragment dict: remove duplicate O=C(O)O key (formic acid carbonic) and add O=CO for real formic acid — surgical single-entry fix"
  - "DVAL cross-refs on StereoProvenance: add dval_cross_refs List[str] field to StereoProvenance schema (not just ProvStageRecord) to expose DVAL-01 linkage"
  - "DVAL cross-refs on tautomer stage: mutate taut_stage.dval_cross_refs after _capture_tautomer_provenance() returns, keeping capture logic pure"
  - "dval_results accepted as kwarg on standardize_with_provenance(), not inside StandardizationOptions dataclass — keeps internal dataclass clean; Pydantic schema carries it"

patterns-established:
  - "Regression test pattern: inspect module source with regex to count key occurrences, preventing future duplicate-key overwrites"
  - "DVAL cross-ref injection: post-stage annotation pattern — build stage first, then annotate with DVAL refs, keeping stage builder methods pure"

requirements-completed: [STD-01, STD-02, STD-03, STD-04, STD-05, STD-06]

# Metrics
duration: 4min
completed: 2026-02-23
---

# Phase 02 Plan 03: Gap Closure Summary

**Deduplicated COUNTERION_NAMES (O=CO=formic acid, 55 unique entries) and DVAL-01/DVAL-03 cross-reference population in stereo and tautomer provenance, lifting Phase 02 verification from 10/11 to 11/11**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-23T13:59:40Z
- **Completed:** 2026-02-23T14:03:13Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments

- Fixed silent data correctness bug: `COUNTERION_NAMES` had a duplicate `"O=C(O)O"` key where "formic acid (carbonic)" silently overwrote "carbonic acid" at runtime — corrected by removing the duplicate and adding `"O=CO"` for real formic acid (HCOOH, MW~46)
- Implemented DVAL-01 cross-referencing: when callers pass `dval_results={"undefined_stereo": {"count": N}}`, the stereo provenance's `dval_cross_refs` is populated with `"DVAL-01: N undefined stereocenters detected"`
- Implemented DVAL-03 cross-referencing: when `dval_results={"tautomer_detection": {"count": N}}`, the tautomer stage's `dval_cross_refs` is populated with `"DVAL-03: N tautomers enumerated"`
- Added `dval_results: Optional[dict[str, Any]]` field to Pydantic `StandardizationOptions` schema, route passes it through; fully backward compatible (None default, empty lists when omitted)
- Added regression-preventing test: `test_no_duplicate_keys_in_fragment_dict` uses source inspection to count key occurrences, failing immediately if a duplicate is introduced
- All 82 standardization tests pass (77 prior + 5 new)

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix duplicate fragment dictionary key and add formic acid** - `361fa5e` (fix)
2. **Task 2: Populate DVAL cross-references in stereo provenance** - `47f13c9` (feat)

**Plan metadata:** see final docs commit below

## Files Created/Modified

- `/Volumes/Data_Drive/Project/2026/chemstructval/backend/app/services/standardization/fragment_dict.py` - Removed duplicate `"O=C(O)O"` entry; added `"O=CO": {"name": "formic acid", "role": "salt"}`
- `/Volumes/Data_Drive/Project/2026/chemstructval/backend/app/services/standardization/provenance.py` - Added `dval_results` kwarg; populate `stereo_dval_cross_refs` and `taut_dval_cross_refs`; removed TODO comment
- `/Volumes/Data_Drive/Project/2026/chemstructval/backend/app/schemas/standardization.py` - Added `dval_cross_refs` field to `StereoProvenance`; added `dval_results` field to `StandardizationOptions`; added `Any` import
- `/Volumes/Data_Drive/Project/2026/chemstructval/backend/app/api/routes/standardization.py` - Pass `dval_results=body.options.dval_results` to `standardize_with_provenance()`
- `/Volumes/Data_Drive/Project/2026/chemstructval/backend/tests/test_standardization/test_provenance.py` - Added 5 tests: `test_no_duplicate_keys_in_fragment_dict`, `test_formic_acid_classified`, `test_dval_cross_refs_populated_with_stereo`, `test_dval_cross_refs_empty_without_dval_results`, `test_dval_cross_refs_tautomer`

## Decisions Made

- Fragment dict fix: surgical removal of one duplicate line + one new entry; no other entries touched
- DVAL cross-refs placed on `StereoProvenance` (added `dval_cross_refs` field) rather than only on `ProvStageRecord` stages — this puts the DVAL-01 ref where it belongs conceptually (stereo summary context, not standardizer stage)
- Tautomer DVAL-03 ref placed on `ProvStageRecord.dval_cross_refs` (already existing field) for the `tautomer_canonicalization` stage
- `dval_results` stays as a kwarg to `standardize_with_provenance()` (not embedded in the internal `StandardizationOptions` dataclass) — keeps the dataclass internal-only; only the Pydantic schema exposes it to callers

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

- Missing `prometheus_fastapi_instrumentator` package caused conftest import error during initial test run; installed (`pip install prometheus-fastapi-instrumentator`) and tests passed immediately.

## Next Phase Readiness

- Phase 02 verification score: 11/11 must-haves verified
- Fragment dict has 55 unique entries, no duplicates
- DVAL cross-reference infrastructure complete end-to-end (schema -> service -> route)
- Ready to advance to Phase 03 (Batch Analytics), beginning with INFRA-01 (batch pagination + Redis TTL)

---
*Phase: 02-standardization-intelligence*
*Completed: 2026-02-23*
