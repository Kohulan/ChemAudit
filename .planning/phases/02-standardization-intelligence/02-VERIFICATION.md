---
phase: 02-standardization-intelligence
verified: 2026-02-23T15:00:00Z
status: passed
score: 11/11 must-haves verified
re_verification:
  previous_status: gaps_found
  previous_score: 10/11
  gaps_closed:
    - "Fragment dictionary has no duplicate canonical SMILES keys — O=C(O)O appears exactly once (carbonic acid), O=CO added for formic acid"
    - "DVAL cross-references are populated in stereo provenance when dval_results is provided — DVAL-01 on StereoProvenance, DVAL-03 on tautomer ProvStageRecord"
  gaps_remaining: []
  regressions: []
human_verification:
  - test: "Open standardization UI, run a molecule with include_provenance=true, observe ProvenanceTimeline renders"
    expected: "Vertical timeline with expandable stage cards showing before/after SMILES, change badges, and 'Show structure' button"
    why_human: "Visual rendering and interactive expand/collapse cannot be verified programmatically"
  - test: "Click 'Show structure' button on a stage card with charge changes"
    expected: "MoleculeViewer renders with blue-highlighted atoms corresponding to changed atom indices"
    why_human: "RDKit.js on-demand rendering and atom color highlighting requires browser execution"
  - test: "Submit a molecule without include_provenance, verify StepsList still renders"
    expected: "Standard StepsList with steps_applied is shown instead of ProvenanceTimeline"
    why_human: "UI fallback path requires visual confirmation"
---

# Phase 02: Standardization Intelligence Verification Report

**Phase Goal:** Users running standardization see exactly what changed and why — every stage of the ChEMBL pipeline (Checker, Standardizer, GetParent, Tautomer) produces a provenance record with atom-level before/after diffs, removed fragment names, and stereo change tracking

**Verified:** 2026-02-23T15:00:00Z
**Status:** passed
**Re-verification:** Yes — after gap closure (Plan 03 closed both gaps from initial verification)

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Standardization request with include_provenance=true returns provenance with per-stage records (checker, standardizer, get_parent, tautomer) | VERIFIED | 28/28 provenance tests pass; `test_basic_provenance_structure` confirms 4 stage names present |
| 2 | Each provenance stage record has input_smiles, output_smiles, applied flag, and typed change arrays | VERIFIED | `ProvStageRecord` Pydantic model has all required fields; 82/82 standardization tests pass |
| 3 | Neutralization provenance (STD-02) reports per-atom charge changes with atom_idx, element, before_charge, after_charge, rule_name, and smarts | VERIFIED | `ChargeChange` schema fully typed; `test_neutralization_charge_tracking` passes for alanine zwitterion |
| 4 | Parent extraction provenance (STD-04) reports removed fragments with smiles, name, role (salt/solvent/counterion/unknown), and mw | VERIFIED | `FragmentRemoval` schema with all fields; `test_parent_extraction_provenance` confirms HCl named "hydrochloric acid" |
| 5 | Tautomer provenance (STD-01) reports input_smiles, canonical_smiles, num_tautomers_enumerated, modified_atoms, modified_bonds, stereo_stripped | VERIFIED | `TautomerProvenance` model has all 7 required fields; `test_tautomer_provenance` passes |
| 6 | Functional group audit (STD-03) reports atom-level changes with rule_name identification | VERIFIED | `NORMALIZE_RULE_NAMES` dict with 9 patterns; `test_functional_group_audit` passes for DMSO |
| 7 | Existing API consumers without include_provenance receive identical responses (provenance is null) | VERIFIED | `StandardizationResult.provenance` defaults to None; backward compat test passes; all 82 tests pass without regression |
| 8 | Kekulization/aromaticity changes tracked as ring-level changes with ring_atoms, ring_size, before_type, after_type (STD-05) | VERIFIED | `_capture_ring_aromaticity_changes()` integrated into standardizer stage; 4 ring tracking tests pass |
| 9 | Stereo provenance includes per-center breakdown with atom_idx, type, before_config, after_config, reason (STD-06) | VERIFIED | `StereoTracker.compare()` returns `per_center_detail`; `StereoProvenance.dval_cross_refs` field populated when `dval_results` provided; 4 stereo tracker tests + 4 stereo normalization tests pass |
| 10 | Frontend displays provenance as vertical timeline with expandable stage cards when provenance is present | VERIFIED | `ProvenanceTimeline.tsx` (144 lines) and `ProvenanceStageCard.tsx` (398 lines) exist; TypeScript compiles without errors; conditional render in `StandardizationResults.tsx` confirmed |
| 11 | DVAL cross-references are populated in stereo provenance when dval_results are provided | VERIFIED | `StereoProvenance.dval_cross_refs` populated with "DVAL-01: N undefined stereocenters detected" when `dval_results={"undefined_stereo": {"count": N}}`; tautomer stage `ProvStageRecord.dval_cross_refs` populated with "DVAL-03: N tautomers enumerated"; backward compatible (empty list when dval_results is None); 3 new DVAL cross-ref tests pass |

**Score:** 11/11 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `backend/app/services/standardization/provenance.py` | ProvenancePipeline with per-stage diff capture and DVAL cross-ref injection | VERIFIED | `standardize_with_provenance(mol, options, dval_results=None)` — dval_results kwarg added; stereo cross-refs at lines 271-287; tautomer cross-refs injected post-stage-build; TODO comment removed |
| `backend/app/services/standardization/fragment_dict.py` | COUNTERION_NAMES with 55 unique entries, no duplicate keys | VERIFIED | 55 entries confirmed; `O=C(O)O` (carbonic acid) appears exactly once at line 52; `O=CO` (formic acid) added at line 53; live dict has 55 unique keys |
| `backend/app/schemas/standardization.py` | StereoProvenance.dval_cross_refs and StandardizationOptions.dval_results | VERIFIED | `StereoProvenance.dval_cross_refs: List[str]` at line 150; `StandardizationOptions.dval_results: Optional[dict[str, Any]]` at line 190; backward compatible with None default |
| `backend/app/api/routes/standardization.py` | Route passes dval_results to standardize_with_provenance | VERIFIED | Line 104: `dval_results=body.options.dval_results` passed in provenance call |
| `backend/tests/test_standardization/test_provenance.py` | 28 tests including 5 gap-closure tests | VERIFIED | 28 tests, all pass in 0.44s; 5 new: `test_no_duplicate_keys_in_fragment_dict`, `test_formic_acid_classified`, `test_dval_cross_refs_populated_with_stereo`, `test_dval_cross_refs_empty_without_dval_results`, `test_dval_cross_refs_tautomer` |
| `backend/tests/test_standardization/test_provenance_m22.py` | 12 tests for STD-05 and STD-06 | VERIFIED | 12 tests pass; no regressions from gap closure changes |
| `frontend/src/types/standardization.ts` | TypeScript interfaces for all provenance types including dval_cross_refs on StereoProvenance | VERIFIED | All 9 provenance interfaces present; `StereoProvenance` has `dval_cross_refs: string[]` |
| `frontend/src/components/standardization/ProvenanceTimeline.tsx` | Vertical timeline rendering provenance stages | VERIFIED | 144 lines; full implementation confirmed in initial verification; no changes made in Plan 03 |
| `frontend/src/components/standardization/ProvenanceStageCard.tsx` | Expandable stage card with typed change tables and on-demand MoleculeViewer | VERIFIED | 398 lines; full implementation confirmed in initial verification; no changes made in Plan 03 |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|----|--------|---------|
| `backend/app/api/routes/standardization.py` | `backend/app/services/standardization/provenance.py` | dval_results passed from request options to ProvenancePipeline | VERIFIED | Line 104: `dval_results=body.options.dval_results` — route reads from Pydantic schema, passes as kwarg |
| `backend/app/services/standardization/fragment_dict.py` | `backend/tests/test_standardization/test_provenance.py` | test_no_duplicate_keys_in_fragment_dict validates all keys are unique via source inspection | VERIFIED | Test uses `inspect.getsource()` + regex to count key occurrences; will fail immediately on any future duplicate |
| `backend/app/services/standardization/provenance.py` | `backend/app/services/standardization/fragment_dict.py` | classify_fragment() for removed fragment naming | VERIFIED (regression check) | Unchanged from initial verification — `from app.services.standardization.fragment_dict import classify_fragment` |
| `backend/app/schemas/standardization.py` | `backend/app/api/routes/standardization.py` | StandardizationProvenance and dval_results in response conversion | VERIFIED (regression check) | Unchanged wiring from initial verification; `dval_results` field flows route -> service |
| `frontend/src/components/standardization/StandardizationResults.tsx` | `frontend/src/components/standardization/ProvenanceTimeline.tsx` | Conditional render when result.provenance is truthy | VERIFIED (regression check) | Frontend unchanged in Plan 03; conditional render confirmed in initial verification |

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|------------|------------|-------------|--------|---------|
| STD-01 | 02-01 | Canonical tautomer generation with provenance (input to canonical mapping) | SATISFIED | `TautomerProvenance` with 7 required fields; `test_tautomer_provenance` passes; DVAL-03 cross-ref now injectable via `dval_results` |
| STD-02 | 02-01 | Neutralization report (atom-level charge changes with before/after) | SATISFIED | `ChargeChange` records with atom_idx, element, before_charge, after_charge, rule_name, smarts; live-verified for alanine zwitterion |
| STD-03 | 02-01 | Functional group standardization audit (which groups were normalized and how) | SATISFIED | `NORMALIZE_RULE_NAMES` with 9 patterns; `test_functional_group_audit` passes for DMSO normalization |
| STD-04 | 02-01 | Parent extraction provenance (removed fragments with names, MW, structure) | SATISFIED | `FragmentRemoval` schema; COUNTERION_NAMES with 55 unique curated entries (no duplicates); formic acid correctly keyed as `O=CO` |
| STD-05 | 02-02 | Kekulization/aromaticity normalization report (ring system changes) | SATISFIED | `RingChange` schema; `_capture_ring_aromaticity_changes()` integrated; 4 ring tracking tests pass |
| STD-06 | 02-02, 02-03 | Stereochemistry normalization tracking (stereo changes during standardization) with DVAL cross-references | SATISFIED | `StereoCenterDetail` with atom_idx, type, before_config, after_config, reason; `StereoProvenance.dval_cross_refs` populated from `dval_results`; 4 stereo tracker tests + 3 DVAL cross-ref tests pass |

### Anti-Patterns Found

No anti-patterns found in gap-closure files. The TODO comment in `provenance.py` has been replaced with a documented implementation. No new TODOs, placeholders, or empty returns introduced.

### Human Verification Required

#### 1. ProvenanceTimeline renders in browser

**Test:** Submit a molecule (e.g., `[NH3+]CC([O-])=O`) with `include_provenance: true` checked in the UI
**Expected:** Vertical timeline with expandable stage cards — Checker, Standardizer, Parent Extraction, Tautomer Canonicalization. Standardizer card auto-expands showing 2 charge change rows.
**Why human:** Visual rendering, Framer Motion animations, and interactive expand/collapse require browser execution

#### 2. On-demand structure rendering with atom highlighting

**Test:** Expand a stage card that has charge changes, click "Show structure" button
**Expected:** MoleculeViewer renders the output SMILES with atom indices highlighted (blue for charge-changed atoms)
**Why human:** RDKit.js on-demand rendering and canvas-based atom color highlighting cannot be verified programmatically

#### 3. StepsList fallback when provenance is null

**Test:** Submit a molecule without the include_provenance option checked
**Expected:** The standard StepsList (not ProvenanceTimeline) renders in the Standardization Pipeline section
**Why human:** Visual confirmation of conditional rendering path in browser required

---

## Re-verification Summary

**Both gaps from the initial verification were fully resolved by Plan 03.**

**Gap 1 resolved — Duplicate fragment dictionary key fixed:**

`fragment_dict.py` previously had two entries using key `"O=C(O)O"` (lines 52 and 54), causing Python to silently overwrite "carbonic acid" with "formic acid (carbonic)". Plan 03 removed the duplicate entry and added a correctly-keyed entry `"O=CO"` for real formic acid (HCOOH, MW 46.025). The dictionary now has 55 unique entries. A regression-preventing test (`test_no_duplicate_keys_in_fragment_dict`) uses source inspection to count key occurrences and will fail immediately if a duplicate is introduced in future. `test_formic_acid_classified` confirms `classify_fragment("OC=O")` returns `name="formic acid"`, `role="salt"`, `mw~46.0`.

**Gap 2 resolved — DVAL cross-references now populated:**

`provenance.py` previously contained a TODO comment deferring DVAL cross-ref population. Plan 03 implemented full cross-ref injection:
- `StereoProvenance.dval_cross_refs` is populated with `"DVAL-01: N undefined stereocenters detected"` when `dval_results={"undefined_stereo": {"count": N}}` is provided
- `ProvStageRecord.dval_cross_refs` on the tautomer stage is populated with `"DVAL-03: N tautomers enumerated"` when `dval_results={"tautomer_detection": {"count": N}}` is provided
- Backward compatible: when `dval_results` is None (default), all `dval_cross_refs` remain empty lists — all 77 pre-existing tests still pass unchanged
- `StandardizationOptions.dval_results: Optional[dict[str, Any]]` added to the Pydantic schema; route passes it through at line 104

**Zero regressions.** All 82 standardization tests pass (77 prior + 5 new). Commits `361fa5e` (fragment dict fix) and `47f13c9` (DVAL cross-refs) verified in git log.

---

_Verified: 2026-02-23T15:00:00Z_
_Verifier: Claude (gsd-verifier)_
