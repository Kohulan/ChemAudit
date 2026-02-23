---
phase: 01-deep-validation
verified: 2026-02-23T13:00:00Z
status: passed
score: 12/12 must-haves verified
re_verification: false
human_verification:
  - test: "Deep Validation tab visual appearance and interaction"
    expected: "Tab appears after 'Validate & Score' tab, 3 collapsible domain sections visible, segmented control switches view mode, gear icon opens severity config panel"
    why_human: "Visual rendering and interactive UX elements cannot be verified programmatically without a running browser"
  - test: "Atom cross-linking from DeepCheckCard to molecule viewer"
    expected: "Clicking an atom index badge in a check result (e.g., stereo center at atom 1) highlights that atom in the RDKit.js molecule viewer"
    why_human: "Requires rendered DOM interaction in a live browser session"
  - test: "Severity override localStorage persistence"
    expected: "After changing a check from WARNING to ERROR and refreshing the page, the override is still applied"
    why_human: "Requires browser localStorage + page refresh cycle"
---

# Phase 1: Deep Validation Verification Report

**Phase Goal:** Implement 16 deep validation checks (stereo/tautomer, chemical composition, structural complexity) with backend checks, tests, startup assertion, and a frontend Deep Validation tab.
**Verified:** 2026-02-23T13:00:00Z
**Status:** PASSED
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| #  | Truth | Status | Evidence |
|----|-------|--------|----------|
| 1  | A molecule with undefined stereocenters returns a check result listing each undefined center's atom index and all enumerated stereoisomer SMILES (up to 128 cap) | VERIFIED | `StereoisomerEnumerationCheck` in `deep_stereo_tautomer.py` uses `FindMolChiralCenters(includeUnassigned=True)` and `EnumerateStereoisomers` with ENUMERATION_CAP=128; returns `atom_indices`, `stereoisomer_smiles`, `cap_exceeded` in details |
| 2  | A molecule with tautomeric forms returns the canonical tautomer SMILES, a count of enumerable tautomers, and a flag indicating whether the input is the canonical form | VERIFIED | `TautomerDetectionCheck` uses `rdMolStandardize.TautomerEnumerator().Canonicalize()` (not the non-existent `GetCanonicalTautomer`); returns `tautomer_count`, `canonical_smiles`, `is_canonical_form` |
| 3  | A molecule with aromatic system issues (unusual ring sizes, charged aromatics) returns a structured warning with affected atom indices | VERIFIED | `AromaticSystemValidationCheck` checks ring sizes not in {5,6} and `GetFormalCharge() != 0` on aromatic atoms; returns `unusual_ring_sizes`, `charged_aromatics`, `affected_atoms` |
| 4  | A molecule with 3D coordinates is identified as 3D; one with 2D as 2D; one with no conformer as 'no_coordinates' | VERIFIED | `CoordinateDimensionCheck` checks `GetNumConformers()` and z-coordinate values; correctly classifies all four states |
| 5  | Validation cache keys include a CHECKS_VERSION prefix so adding new checks invalidates stale cached results | VERIFIED | `CHECKS_VERSION = "v2"` in `cache.py`; key format `validation:v2:{inchikey}:{hash}` confirmed by live test |
| 6  | A mixture input (SMILES with dot separator) is flagged with each fragment classified as drug, salt, solvent, or unknown | VERIFIED | `MixtureDetectionCheck` uses `GetMolFrags` + MolVS `REMOVE_FRAGMENTS` patterns + heuristic fallback; returns `fragments` list with `classification`, `smiles`, `molecular_weight` |
| 7  | A molecule containing a known solvent fragment returns a solvent contamination warning | VERIFIED | `SolventContaminationCheck` with 15-solvent curated list + canonical SMILES + bidirectional substructure matching |
| 8  | An inorganic molecule (no carbon) or organometallic (contains metals) is flagged with atom details | VERIFIED | `InorganicFilterCheck` checks `GetAtomicNum() == 6` for carbon, METAL_ATOMIC_NUMS set for metals; ERROR for inorganic, WARNING for organometallic |
| 9  | A molecule with radical electrons is flagged with affected atom indices | VERIFIED | `RadicalDetectionCheck` iterates `atom.GetNumRadicalElectrons()` |
| 10 | A molecule with isotope labels is flagged with isotope details per atom | VERIFIED | `IsotopeLabelDetectionCheck` iterates `atom.GetIsotope()` with common_name mapping (deuterium, 13C, etc.) |
| 11 | A trivially small molecule (heavy atom count <= 3) is flagged as too small | VERIFIED | `TrivialMoleculeCheck` uses `mol.GetNumHeavyAtoms()` with threshold=3; ERROR severity |
| 12 | A molecule with hypervalent atoms returns a warning with affected atom indices, actual vs allowed valences | VERIFIED | `HypervalentAtomCheck` uses `GetPeriodicTable().GetValenceList()` against `GetTotalValence()` |
| 13 | A molecule with SGroup markers or MW > 1500 is flagged as a possible polymer | VERIFIED | `PolymerDetectionCheck` checks SGroups via `GetMolSubstanceGroups`, MW threshold, and dummy atoms |
| 14 | A molecule with 3 or 4-membered rings returns a ring strain warning with atom indices | VERIFIED | `RingStrainCheck` uses `AtomRings()` and checks `len(ring) in (3, 4)` |
| 15 | A molecule with rings > 12 atoms returns a macrocycle warning | VERIFIED | `MacrocycleDetectionCheck` uses `AtomRings()` and checks `len(ring) > 12` (not >=12 — boundary condition confirmed) |
| 16 | A charged molecule returns charge details; zwitterions (net charge 0 with both + and -) are specifically flagged | VERIFIED | `ChargedSpeciesCheck` computes net_charge, is_zwitterion (net==0 AND has both charges) |
| 17 | A molecule with unusual explicit hydrogen patterns is flagged | VERIFIED | `ExplicitHydrogenAuditCheck` uses `GetNumExplicitHs()` for H on heavy atoms + counts H atom objects |
| 18 | Startup assertion confirms all 16 deep validation check names are registered | VERIFIED | `EXPECTED_DEEP_VALIDATION_CHECKS` set in `main.py` lifespan; `logger.warning` if missing; CI hard-assert in `TestAllDeepValidationChecks` |
| 19 | Users see a 'Deep Validation' tab in single validation results | VERIFIED | `TabType` includes `'deep-validation'`; tab config with Microscope icon added to TABS array in `SingleValidation.tsx`; `DeepValidationTab` rendered when `activeTab === 'deep-validation'` |
| 20 | Deep Validation tab shows collapsible sections grouped by domain | VERIFIED | `DomainSection` component renders 3 sections using `DEEP_CHECK_DOMAINS`; Framer Motion AnimatePresence for collapse |
| 21 | Users can toggle between Category view and Severity view | VERIFIED | `viewMode` state with segmented control buttons; severity view groups checks by effective severity |
| 22 | Users can configure severity per check via gear icon; overrides persist in localStorage | VERIFIED | `SeverityConfigPanel` opened by gear/Settings icon; `useDeepValidationConfig` hook uses `useLocalStorage` under key `chemaudit:deep-validation-config` |
| 23 | Clicking atom indices in check results highlights corresponding atoms on molecule viewer | VERIFIED | `onHighlightAtoms={setHighlightedAtoms}` wired from `SingleValidation` → `DeepValidationTab` → `DeepCheckCard`; `highlightAtoms={highlightedAtoms}` passed to `MoleculeViewer` |
| 24 | Stereoisomer enumeration results show a collapsible list of enumerated SMILES | VERIFIED | `StereoisomerList` component in `DeepCheckCard` for `stereoisomer_enumeration` check name |
| 25 | Mixture detection results show a fragment classification table | VERIFIED | `FragmentClassificationTable` component in `DeepCheckCard` for `mixture_detection` check name |

**Score:** 12/12 must-haves (from PLAN frontmatter) verified. All 25 observable truths derived from the 4 plan must_haves sections also VERIFIED.

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `backend/app/services/validation/checks/deep_stereo_tautomer.py` | 4 checks: stereoisomer_enumeration, tautomer_detection, aromatic_system_validation, coordinate_dimension; contains `@CheckRegistry.register` | VERIFIED | 493 lines; 4 `@CheckRegistry.register` decorators confirmed |
| `backend/tests/test_validation/test_deep_stereo_tautomer_checks.py` | Tests for DVAL-01..05; contains `class TestStereoisomerEnumeration` | VERIFIED | 419 lines; 31 tests; all pass |
| `backend/app/core/cache.py` | CHECKS_VERSION constant in cache key | VERIFIED | `CHECKS_VERSION = "v2"`; key format `validation:v2:...` confirmed live |
| `backend/app/services/validation/checks/deep_composition.py` | 6 checks: mixture_detection, solvent_contamination, inorganic_filter, radical_detection, isotope_label_detection, trivial_molecule | VERIFIED | 859 lines; 6 `@CheckRegistry.register` decorators |
| `backend/tests/test_validation/test_deep_composition_checks.py` | Tests for DVAL-06..11; contains `class TestMixtureDetection` | VERIFIED | 656 lines; 55 tests; all pass |
| `backend/app/services/validation/checks/deep_complexity.py` | 6 checks: hypervalent_atoms, polymer_detection, ring_strain, macrocycle_detection, charged_species, explicit_hydrogen_audit | VERIFIED | 613 lines; 6 `@CheckRegistry.register` decorators |
| `backend/tests/test_validation/test_deep_complexity_checks.py` | Tests for DVAL-12..17; contains `class TestHypervalentAtoms` | VERIFIED | 906 lines; 65 tests; all pass |
| `backend/app/main.py` | Startup assertion with `EXPECTED_DEEP_VALIDATION_CHECKS` | VERIFIED | `EXPECTED_DEEP_VALIDATION_CHECKS` set declared at module level; used in lifespan to `logger.warning` if missing |
| `frontend/src/components/validation/DeepValidationTab.tsx` | Main tab container with segmented control and gear icon; contains `DeepValidationTab` | VERIFIED | 369 lines; segmented control, Settings icon, DomainSection, SeverityConfigPanel all present |
| `frontend/src/components/validation/DeepCheckCard.tsx` | Individual check result card; contains `DeepCheckCard` | VERIFIED | 746 lines; per-check-name detail rendering, `onHighlightAtoms` callback |
| `frontend/src/components/validation/SeverityConfigPanel.tsx` | Severity configuration modal; contains `SeverityConfigPanel` | VERIFIED | 263 lines; 3-button toggle (ERROR/WARNING/INFO), reset all, close button |
| `frontend/src/hooks/useDeepValidationConfig.ts` | localStorage-backed severity config hook; contains `useDeepValidationConfig` | VERIFIED | 48 lines; uses `useLocalStorage` hook; `setSeverityOverride`, `removeSeverityOverride`, `resetAllOverrides`, `getEffectiveSeverity` exported |
| `frontend/src/pages/SingleValidation.tsx` | Deep Validation tab entry; contains `deep-validation` | VERIFIED | `TabType` includes `'deep-validation'`; `DeepValidationTab` rendered with `checks={result.all_checks}` and `onHighlightAtoms={setHighlightedAtoms}` |
| `frontend/src/components/validation/StereoisomerList.tsx` | Collapsible stereoisomer SMILES list | VERIFIED | 107 lines; expand/collapse with count |
| `frontend/src/components/validation/FragmentClassificationTable.tsx` | Fragment classification mini-table | VERIFIED | 93 lines; colored classification badges |
| `frontend/src/tests/components/DeepValidationTab.test.tsx` | Frontend tests for all deep validation components | VERIFIED | 764 lines; 41 tests; all pass |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `deep_stereo_tautomer.py` | `registry.py` | `@CheckRegistry.register()` decorator | WIRED | 4 registrations found at lines 38, 166, 251, 369 |
| `engine.py` | `deep_stereo_tautomer.py` | `import app.services.validation.checks.deep_stereo_tautomer` | WIRED | Line 16 of `engine.py` |
| `engine.py` | `deep_composition.py` | `import app.services.validation.checks.deep_composition` | WIRED | Line 17 of `engine.py` |
| `engine.py` | `deep_complexity.py` | `import app.services.validation.checks.deep_complexity` | WIRED | Line 15 of `engine.py` |
| `cache.py` | validation cache keys | `CHECKS_VERSION` prefix | WIRED | Key format `validation:v2:{inchikey}:{hash}` confirmed live |
| `deep_composition.py` | `registry.py` | `@CheckRegistry.register()` | WIRED | 6 registrations: lines 238, 335, 478, 604, 689, 781 |
| `deep_complexity.py` | `registry.py` | `@CheckRegistry.register()` | WIRED | 6 registrations: lines 24, 123, 227, 317, 411, 514 |
| `main.py` | `registry.py` | `EXPECTED_DEEP_VALIDATION_CHECKS` startup check | WIRED | Lines 22-100 of `main.py`; imports `CheckRegistry` in lifespan |
| `DeepValidationTab.tsx` | `useDeepValidationConfig.ts` | `useDeepValidationConfig` hook import | WIRED | Line 20 of `DeepValidationTab.tsx` |
| `DeepCheckCard.tsx` | `SingleValidation.tsx` | `onHighlightAtoms` callback | WIRED | `onHighlightAtoms={setHighlightedAtoms}` at line 785; `highlightAtoms={highlightedAtoms}` to MoleculeViewer at line 1126 |
| `SingleValidation.tsx` | `DeepValidationTab.tsx` | Tab rendering | WIRED | Lines 780-792; `checks={result.all_checks}` — `all_checks` field confirmed in both backend schema and frontend type |

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|------------|-------------|--------|----------|
| DVAL-01 | 01-01, 01-04 | Undefined stereocenter detection with count and enumeration | SATISFIED | `StereoisomerEnumerationCheck`; `FindMolChiralCenters(includeUnassigned=True)`; 8 tests pass |
| DVAL-02 | 01-01, 01-04 | Stereoisomer enumeration returning actual SMILES | SATISFIED | Same check; `EnumerateStereoisomers`; cap=128; SMILES returned in `stereoisomer_smiles`; `StereoisomerList` in frontend |
| DVAL-03 | 01-01, 01-04 | Tautomer flagging with canonical form detection | SATISFIED | `TautomerDetectionCheck`; `Canonicalize()` used; `is_canonical_form` flag; 7 tests pass |
| DVAL-04 | 01-01, 01-04 | Aromatic system validation | SATISFIED | `AromaticSystemValidationCheck`; ring size + charge analysis; 7 tests pass |
| DVAL-05 | 01-01, 01-04 | Coordinate dimension check | SATISFIED | `CoordinateDimensionCheck`; 2d/3d/no_coordinates/degenerate; 7 tests pass |
| DVAL-06 | 01-02, 01-04 | Mixture detection with fragment classification | SATISFIED | `MixtureDetectionCheck`; drug/salt/solvent/unknown; `FragmentClassificationTable` in frontend; 7 tests pass |
| DVAL-07 | 01-02, 01-04 | Solvent contamination detection | SATISFIED | `SolventContaminationCheck`; 15-solvent list; false-positive guard for broad MolVS SMARTS; tests pass |
| DVAL-08 | 01-02, 01-04 | Inorganic/organometallic filter | SATISFIED | `InorganicFilterCheck`; carbon check + METAL_ATOMIC_NUMS; ERROR/WARNING severity distinction; tests pass |
| DVAL-09 | 01-02, 01-04 | Radical electron detection | SATISFIED | `RadicalDetectionCheck`; `GetNumRadicalElectrons()`; affected_atoms populated; tests pass |
| DVAL-10 | 01-02, 01-04 | Isotope label detection | SATISFIED | `IsotopeLabelDetectionCheck`; `GetIsotope()`; common_name mapping; INFO severity; tests pass |
| DVAL-11 | 01-02, 01-04 | Trivial molecule check | SATISFIED | `TrivialMoleculeCheck`; `GetNumHeavyAtoms() <= 3`; ERROR severity; tests pass |
| DVAL-12 | 01-03, 01-04 | Hypervalent atom detection | SATISFIED | `HypervalentAtomCheck`; `GetValenceList()` vs `GetTotalValence()`; tests pass |
| DVAL-13 | 01-03, 01-04 | Polymer/repeating unit detection | SATISFIED | `PolymerDetectionCheck`; SGroup markers, MW > 1500, dummy atoms; INFO severity; tests pass |
| DVAL-14 | 01-03, 01-04 | Ring strain detection | SATISFIED | `RingStrainCheck`; 3/4-membered rings; heuristic only (no force field); WARNING severity; tests pass |
| DVAL-15 | 01-03, 01-04 | Macrocycle detection | SATISFIED | `MacrocycleDetectionCheck`; > 12 atoms (not >= 12); SSSR note included; INFO severity; tests pass |
| DVAL-16 | 01-03, 01-04 | Charged species with zwitterion identification | SATISFIED | `ChargedSpeciesCheck`; net_charge + is_zwitterion (net==0 AND both charges); tests pass |
| DVAL-17 | 01-03, 01-04 | Explicit hydrogen audit | SATISFIED | `ExplicitHydrogenAuditCheck`; `GetNumExplicitHs()` + H atom object detection; INFO severity; tests pass |

**Note:** The REQUIREMENTS.md traceability table shows "Pending" status for all DVAL requirements — this is stale documentation, not a gap. The requirement checkboxes at the top of REQUIREMENTS.md correctly show `[x]` for all DVAL-01 through DVAL-17. Implementation was verified programmatically against the actual codebase.

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `DeepValidationTab.tsx` | 315 | `return null` | Info | Legitimate: returns null only when a domain has zero matching checks (empty domain guard, not a stub) |
| `DeepCheckCard.tsx` | 114 | `return null` | Info | Legitimate: returns null when details entries list is empty (generic fallback with no entries) |

No blockers or warnings found. All `return null` instances are guarded conditionals on legitimate empty-state data, not implementation stubs.

### Human Verification Required

#### 1. Deep Validation Tab Visual Appearance

**Test:** Start the dev server (`npm run dev` in frontend), navigate to the single validation page, enter a SMILES string (e.g., `CC(O)C(N)CC` for stereo test or `CCN.Cl` for mixture test), click Validate, then click the "Deep Validation" tab.
**Expected:** Tab appears in the tab bar after "Validate & Score" with a Microscope icon. Three collapsible domain sections are visible: "Stereo & Tautomers", "Chemical Composition", "Structural Complexity". Each section shows check cards with severity badges, messages, and atom index badges.
**Why human:** Visual rendering, icon display, and layout cannot be verified without a live browser.

#### 2. Atom Cross-Linking to Molecule Viewer

**Test:** With a validation result open on the Deep Validation tab, find a check card with atom indices (e.g., stereoisomer_enumeration for a molecule with undefined stereocenters). Click one of the atom index badges.
**Expected:** The 2D molecule viewer highlights the clicked atom (changes color or shows selection indicator on the correct atom).
**Why human:** Requires RDKit.js rendering in a live browser with DOM interaction to verify the `highlightAtoms` prop is consumed by the MoleculeViewer.

#### 3. Severity Override Persistence

**Test:** Open the severity config panel (gear icon), change one check from its default severity to ERROR. Refresh the page. Navigate back to the Deep Validation tab.
**Expected:** The severity override is still applied (check shows ERROR severity), demonstrating localStorage persistence.
**Why human:** Requires browser localStorage read/write cycle across page refresh.

### Test Results Summary

| Test Suite | Tests | Result |
|-----------|-------|--------|
| `test_deep_stereo_tautomer_checks.py` | 31 | All pass |
| `test_deep_composition_checks.py` | 55 | All pass |
| `test_deep_complexity_checks.py` | 65 | All pass |
| **Backend total** | **151** | **All pass** |
| `DeepValidationTab.test.tsx` | 41 | All pass |
| **Frontend total** | **41** | **All pass** |

### Registration Verification

Live Python verification (conda `cheminformatics` environment):

```
Missing: set()
Extra (other registered checks): {'inchi_roundtrip', 'connectivity', 'valence', ...}
All 16 deep checks present: True
CHECKS_VERSION: v2
Cache key sample: validation:v2:BSYNRYMUTXBXSQ-UHFFFAOYSA-N:5ef5ef0364b6939c
v2 in key: True
```

All 16 deep validation checks (`stereoisomer_enumeration`, `tautomer_detection`, `aromatic_system_validation`, `coordinate_dimension`, `mixture_detection`, `solvent_contamination`, `inorganic_filter`, `radical_detection`, `isotope_label_detection`, `trivial_molecule`, `hypervalent_atoms`, `polymer_detection`, `ring_strain`, `macrocycle_detection`, `charged_species`, `explicit_hydrogen_audit`) are registered and discoverable via `CheckRegistry.get_all()`.

### TypeScript Compilation

```
npx tsc --noEmit --pretty
(no output — zero errors)
```

---

_Verified: 2026-02-23T13:00:00Z_
_Verifier: Claude (gsd-verifier)_
