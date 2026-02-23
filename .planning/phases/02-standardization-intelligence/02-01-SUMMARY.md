---
phase: 02-standardization-intelligence
plan: "01"
subsystem: api
tags: [rdkit, chembl, standardization, provenance, pydantic, fastapi]

# Dependency graph
requires:
  - phase: 01-deep-validation
    provides: validation infrastructure and RDKit setup

provides:
  - ProvenancePipeline wrapping ChEMBL pipeline with per-stage atom-level diffs
  - COUNTERION_NAMES fragment dictionary with 57+ curated entries
  - classify_fragment() for named fragment role/MW classification
  - ChargeChange, BondChange, RadicalChange, FragmentRemoval, TautomerProvenance Pydantic models
  - StandardizationProvenance schema extending StandardizationResult
  - include_provenance option in /api/v1/standardize endpoint (backward compatible)

affects:
  - 02-02 (stereo provenance plan — STD-06 builds on StereoProvenance placeholder)
  - frontend standardization UI (provenance display, charge visualization)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - ProvenancePipeline wrapper pattern (wraps existing pipeline, adds provenance externally)
    - SMILES fragment diffing for removed fragment detection (avoids atom-idx issues on count change)
    - Atom-by-atom diff only when atom count preserved (standardizer stage)
    - NORMALIZE_RULE_NAMES lookup for normalization rule identification
    - Separate typed provenance fields per stage (no generic changes: List[dict] to avoid TS any[])

key-files:
  created:
    - backend/app/services/standardization/fragment_dict.py
    - backend/app/services/standardization/provenance.py
    - backend/tests/test_standardization/test_provenance.py
  modified:
    - backend/app/schemas/standardization.py
    - backend/app/services/standardization/chembl_pipeline.py
    - backend/app/api/routes/standardization.py

key-decisions:
  - "ProvenancePipeline wraps (not inherits from) StandardizationPipeline — zero modifications to existing pipeline internals, provenance captured by re-calling same ChEMBL functions stage-by-stage"
  - "Use SMILES fragment diffing (set(before.split('.')) - set(after.split('.'))) for get_parent provenance — avoids atom-idx comparison pitfall when atom count changes after fragment removal"
  - "Atom-level diff only runs when before.GetNumAtoms() == after.GetNumAtoms() — prevents index mismatches for standardizer stage"
  - "Dictionary keys must be RDKit-canonical SMILES — fixed 15 entries during Task 2 (benzoic acid, sulfuric acid, p-TSA, phosphoric acid, citric acid, fumaric acid, etc.)"
  - "Separate typed fields per stage (charge_changes, bond_changes, radical_changes, fragment_removals) instead of generic changes: List[dict] — prevents TypeScript any[] on frontend"
  - "NORMALIZE_RULE_NAMES dict covers 9 common patterns; unmatched changes get rule_name='unknown_normalization'"
  - "Tautomer provenance uses result.tautomers (Python-iterable Mol objects) not result.smiles (C++ vector) — avoids crash"
  - "complexity_flag set when > 100 tautomers enumerated to flag high-complexity molecules without bloating response"

patterns-established:
  - "Rule 1 auto-fix: Fragment dict keys non-canonical — fixed to RDKit-canonical SMILES during test verification"

requirements-completed: [STD-01, STD-02, STD-03, STD-04]

# Metrics
duration: 7min
completed: 2026-02-23
---

# Phase 02 Plan 01: Standardization Provenance Pipeline Summary

**Per-stage atom-level provenance for ChEMBL standardization pipeline: charge diffs (STD-02), fragment classification (STD-04), tautomer modified atoms/bonds (STD-01), and normalization rule identification (STD-03)**

## Performance

- **Duration:** 7 min
- **Started:** 2026-02-23T13:10:39Z
- **Completed:** 2026-02-23T13:18:18Z
- **Tasks:** 2
- **Files modified:** 6

## Accomplishments

- ProvenancePipeline wraps ChEMBL pipeline stage-by-stage, capturing atom-level diffs without modifying pipeline internals
- Fragment dictionary with 57 curated COUNTERION_NAMES entries (chloride, sodium, water, DMSO, benzoic acid, TFA, etc.) and classify_fragment() with graceful fallback
- All four STD provenance types implemented: charge changes (STD-02), normalization rules (STD-03), fragment removals with name/role/MW (STD-04), tautomer modified atoms/bonds/stereo (STD-01)
- 23 new tests passing covering all four STD requirements plus backward compatibility and endpoint integration
- All 65 standardization tests pass — no existing behavior changed

## Task Commits

Each task was committed atomically:

1. **Task 1: ProvenancePipeline, fragment dictionary, and schema extension** - `858345c` (feat)
2. **Task 2: Comprehensive provenance tests** - `ac0e517` (test)

**Plan metadata:** (docs commit — this summary)

## Files Created/Modified

- `backend/app/services/standardization/fragment_dict.py` - COUNTERION_NAMES dict (57 entries) and classify_fragment()
- `backend/app/services/standardization/provenance.py` - ProvenancePipeline with per-stage atom-level diff capture
- `backend/app/schemas/standardization.py` - 7 new Pydantic models (ChargeChange, BondChange, RadicalChange, FragmentRemoval, TautomerProvenance, ProvStageRecord, StereoProvenance, StandardizationProvenance) + provenance field on StandardizationResult + include_provenance on StandardizationOptions
- `backend/app/services/standardization/chembl_pipeline.py` - Added include_provenance: bool = False to internal StandardizationOptions dataclass
- `backend/app/api/routes/standardization.py` - Added ProvenancePipeline singleton, conditional provenance path when include_provenance=True, updated _convert_pipeline_result to accept provenance param
- `backend/tests/test_standardization/test_provenance.py` - 23 tests: 6 fragment dict, 13 pipeline, 4 endpoint integration

## Decisions Made

- ProvenancePipeline wraps (not inherits from) StandardizationPipeline to avoid touching internals
- SMILES fragment diffing used for get_parent stage (avoids atom-idx issues when atom count changes)
- Atom-level diff only runs when before/after atom counts match (standardizer stage)
- Separate typed fields per stage (not generic List[dict]) to prevent TypeScript any[] on frontend
- TautomerEnumerator.Enumerate().tautomers (Python iterable) used, not .smiles (C++ vector — crashes)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fragment dictionary had non-canonical SMILES keys**
- **Found during:** Task 2 (test_common_salts_present and test_fragment_removal_unknown failures)
- **Issue:** 15 dictionary entries used non-canonical SMILES as keys (e.g., `OC(=O)c1ccccc1` instead of `O=C(O)c1ccccc1` for benzoic acid; similar for sulfuric acid, p-TSA, phosphoric acid, oxalic acid, etc.). classify_fragment() canonicalizes input before lookup, so non-canonical keys would never match.
- **Fix:** Corrected all 15 keys to RDKit-canonical SMILES using Chem.MolToSmiles(Chem.MolFromSmiles()) verification. Also removed duplicate `OC(O)=O` (formic acid) that collided with `OC(=O)O` (carbonic acid) after canonicalization — both become `O=C(O)O`.
- **Files modified:** backend/app/services/standardization/fragment_dict.py
- **Verification:** test_common_salts_present and classify_fragment('[Cl-]') tests both pass. classify_fragment('OC(=O)c1ccccc1') now correctly returns name='benzoic acid'.
- **Committed in:** ac0e517 (Task 2 commit, also updated fragment_dict.py)

---

**Total deviations:** 1 auto-fixed (Rule 1 — bug)
**Impact on plan:** Fix was necessary for correctness — non-canonical keys would cause all dictionary lookups to silently fail. No scope creep.

## Issues Encountered

- Pentane (`CCCCC`) is in COUNTERION_NAMES as a solvent — initial test used it as "unknown fragment" and failed. Fixed test to use isobutyric acid (`CC(C)C(=O)O`) which is genuinely not in the dictionary.

## Next Phase Readiness

- Provenance schemas are ready for Plan 02-02 (stereo provenance: STD-06 builds on StereoProvenance placeholder)
- Frontend can now display per-stage provenance when include_provenance=True
- No blockers

## Self-Check: PASSED

- fragment_dict.py: FOUND
- provenance.py: FOUND
- test_provenance.py: FOUND
- 02-01-SUMMARY.md: FOUND
- 858345c (Task 1 commit): FOUND
- ac0e517 (Task 2 commit): FOUND

---
*Phase: 02-standardization-intelligence*
*Completed: 2026-02-23*
