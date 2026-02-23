---
phase: 01-deep-validation
plan: "01"
subsystem: backend-validation
tags:
  - validation
  - stereochemistry
  - tautomers
  - rdkit
  - cache
dependency_graph:
  requires:
    - backend/app/services/validation/checks/base.py
    - backend/app/services/validation/registry.py
    - backend/app/services/validation/engine.py
    - backend/app/core/cache.py
  provides:
    - stereoisomer_enumeration check (DVAL-01/02)
    - tautomer_detection check (DVAL-03)
    - aromatic_system_validation check (DVAL-04)
    - coordinate_dimension check (DVAL-05)
    - CHECKS_VERSION cache key versioning
  affects:
    - backend/app/services/validation/engine.py
    - backend/app/core/cache.py
    - backend/app/services/validation/checks/__init__.py
tech_stack:
  added: []
  patterns:
    - "@CheckRegistry.register() decorator for self-registering checks"
    - "rdkit.Chem.EnumerateStereoisomers with StereoEnumerationOptions(unique=True)"
    - "rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumerator.Canonicalize()"
    - "CHECKS_VERSION prefix in cache keys for invalidation on new check addition"
key_files:
  created:
    - backend/app/services/validation/checks/deep_stereo_tautomer.py
    - backend/tests/test_validation/test_deep_stereo_tautomer_checks.py
  modified:
    - backend/app/services/validation/checks/__init__.py
    - backend/app/services/validation/engine.py
    - backend/app/core/cache.py
decisions:
  - "CHECKS_VERSION set to v2; invalidation scan pattern updated from validation:{inchikey}:* to validation:*:{inchikey}:* to match versioned keys"
  - "Stereoisomer enumeration cap set to 128 (2^7) — returns count-only with cap_exceeded=True above cap, empty stereoisomer_smiles list"
  - "TautomerDetection uses enumerator.Canonicalize() not GetCanonicalTautomer() (does not exist in RDKit)"
  - "Multi-fragment molecules: both StereoisomerEnumeration and TautomerDetection operate on largest fragment only"
  - "AromaticSystemValidation extends existing AromaticityCheck (which handles Kekulize failures) with ring-size and charge analysis only"
  - "CoordinateDimensionCheck always passes (informational) — reports 2d/3d/no_coordinates/degenerate"
metrics:
  duration_seconds: 294
  completed_date: "2026-02-23"
  tasks_completed: 2
  tasks_total: 2
  files_created: 2
  files_modified: 3
---

# Phase 1 Plan 1: M1.1 Deep Stereo/Tautomer Checks Summary

**One-liner:** 4 deep validation checks (stereoisomer enumeration, tautomer detection, aromatic system validation, coordinate dimension) registered via RDKit APIs with CHECKS_VERSION=v2 cache key invalidation.

## What Was Built

### Task 1: 4 M1.1 Check Classes + Cache Key Versioning

**`backend/app/services/validation/checks/deep_stereo_tautomer.py`** — new module containing:

1. **`StereoisomerEnumerationCheck`** (`stereoisomer_enumeration`, covers DVAL-01/02):
   - Finds undefined tetrahedral stereocenters via `FindMolChiralCenters(includeUnassigned=True)`
   - Enumerates stereoisomers with `EnumerateStereoisomers(StereoEnumerationOptions(unique=True))`
   - Cap: 128 (2^7) — above cap returns empty `stereoisomer_smiles` and `cap_exceeded=True`
   - Returns: `undefined_count`, `total_centers`, `atom_indices`, `stereoisomer_smiles`, `enumeration_cap`, `cap_exceeded`
   - Severity: WARNING (undefined stereo) / INFO (all defined)

2. **`TautomerDetectionCheck`** (`tautomer_detection`, covers DVAL-03):
   - Uses `rdMolStandardize.TautomerEnumerator()` — not MolVS (per research mandate)
   - Uses `enumerator.Canonicalize(mol)` — not `GetCanonicalTautomer` (does not exist)
   - Operates on largest fragment for multi-fragment molecules
   - Returns: `tautomer_count`, `canonical_smiles`, `is_canonical_form`, `tautomer_smiles`
   - Severity: INFO (always passes — purely informational)

3. **`AromaticSystemValidationCheck`** (`aromatic_system_validation`, covers DVAL-04):
   - Extends existing `AromaticityCheck` (Kekulize failures) with ring-size and charge analysis
   - Detects aromatic rings not in size {5, 6} via `RingInfo.AtomRings()`
   - Detects aromatic atoms with `GetFormalCharge() != 0`
   - Returns: `unusual_ring_sizes`, `charged_aromatics` (each with atom indices and details)
   - Severity: WARNING / INFO

4. **`CoordinateDimensionCheck`** (`coordinate_dimension`, covers DVAL-05):
   - Checks `mol.GetNumConformers()` — 0 → `no_coordinates`
   - Checks z-coordinates: all 0.0 → `2d`, any non-zero → `3d`
   - Detects degenerate coordinates (all atoms at identical position)
   - Returns: `dimension` ("2d"|"3d"|"no_coordinates"|"degenerate"), `num_conformers`
   - Severity: INFO (always passes — purely informational)

**Cache key versioning (`backend/app/core/cache.py`)**:
- Added `CHECKS_VERSION = "v2"` constant (satisfying STATE.md mandate)
- Updated `validation_cache_key()` format: `validation:{CHECKS_VERSION}:{inchikey}:{checks_hash}`
- Updated `invalidate_cached_validation()` scan pattern: `validation:*:{inchikey}:*` (matches across all versions)

**Engine and init updates**:
- `engine.py`: added `import app.services.validation.checks.deep_stereo_tautomer # noqa: F401`
- `__init__.py`: added all 4 new check classes to imports and `__all__`

### Task 2: Comprehensive Tests

**`backend/tests/test_validation/test_deep_stereo_tautomer_checks.py`** — 31 tests:

- `TestStereoisomerEnumeration` (8 tests): no stereocenters pass, undefined centers detected + enumerated, fully defined pass, cap exceeded, None guard, details structure, affected_atoms, message content
- `TestTautomerDetection` (7 tests): phenol tautomers, ethane no tautomers, canonical True, non-canonical False (cyclohexadienone keto form), None guard, details structure, message content
- `TestAromaticSystemValidation` (7 tests): benzene pass, pyrrole pass, charged aromatic detection, None guard, affected_atoms populated, azulene 7-membered ring flagged, details structure
- `TestCoordinateDimension` (7 tests): no conformer, 2D coords via Compute2DCoords, 3D coords via EmbedMolecule, None guard, details structure, always-passes contract, message content
- `TestRegistration` (2 tests): all 4 checks in registry, all have correct category `stereo_tautomer`

## Verification Results

All verifications pass:

1. 4 M1.1 checks registered in `stereo_tautomer` category: confirmed
2. All 31 new tests pass: confirmed
3. Cache key format includes `:v2:`: confirmed (`validation:v2:{inchikey}:{hash}`)
4. All 79 existing validation tests unbroken: confirmed

## Deviations from Plan

### Auto-fixed Issues

None — plan executed exactly as written with one minor test correction:

**[Rule 1 - Bug] Fixed incorrect naphthalene SMILES in test**
- **Found during:** Task 2 test execution
- **Issue:** Test `test_affected_atoms_populated` used `c1ccc2cccccc2c1` (malformed — 11 carbons, can't kekulize), which returned `None` from RDKit instead of naphthalene
- **Fix:** Replaced with `c1ccc2cccccc12` (azulene, correct 10-atom structure with 5+7 rings) and also updated `test_azulene_unusual_ring` to use this correct SMILES
- **Files modified:** `backend/tests/test_validation/test_deep_stereo_tautomer_checks.py`
- **Commit:** 58e4b11

### Background Discovery

During execution, the system-triggered linter/formatter updated `__init__.py` and `engine.py` to also include `deep_composition.py` and `deep_complexity.py` check imports. These files were pre-existing in the codebase (from prior work). This is a beneficial side-effect — all 27 checks are now properly registered. No action taken as these files were already present and correct.

## Commits

| Task | Commit | Message |
|------|--------|---------|
| Task 1 | `51bf66e` | `feat(01-01): implement 4 M1.1 deep validation checks + cache key versioning` |
| Task 2 | `58e4b11` | `test(01-01): add comprehensive tests for M1.1 deep stereo/tautomer checks` |

## Self-Check: PASSED
