# Design: Fix Regex-Based Molecule Info Bugs

**Date**: 2026-02-24
**Branch**: `fix/molecule-info-regex-bugs` from `development`
**Scope**: Validation endpoint + frontend `useMoleculeInfo` hook + `SingleValidation` page

## Problem

The frontend `useMoleculeInfo.ts` hook uses regex to estimate molecular properties (atom count, bond count, ring count, stereocenters) from SMILES strings. These regexes break on molecules with isotope labels, charges, or bracket atoms:

- `/\d/g` counts isotope digits (`[63Ni]`, `[13C]`) as ring closures → phantom rings
- `/[A-Z][a-z]?/g` misses aromatic atoms (`c`, `n`) → benzene shows 0 atoms
- `replace(/[^A-Z]/gi, '').length` counts bracket contents → inflated atom counts
- `kekulizeSmiles()` regex heuristic is fragile and unreliable

The backend `MoleculeInfo` schema only returns `num_atoms`, forcing the frontend to compute `num_bonds`, `num_rings`, and stereocenter info client-side.

## Solution: Approach B — Backend Authority + RDKit.js Preview

### Backend Changes

**`backend/app/schemas/validation.py` — Enrich `MoleculeInfo`:**

Add fields:
- `num_bonds: Optional[int]`
- `num_rings: Optional[int]`
- `num_aromatic_rings: Optional[int]`
- `num_stereocenters: Optional[int]`
- `has_stereochemistry: Optional[bool]`

**`backend/app/api/routes/validation.py` — Update `extract_molecule_info()`:**

Compute new fields using RDKit:
- `mol.GetNumBonds()`
- `rdMolDescriptors.CalcNumRings(mol)`
- `Descriptors.NumAromaticRings(mol)`
- `len(FindMolChiralCenters(mol, includeUnassigned=True))`
- Derive `has_stereochemistry` from stereocenters + E/Z bond detection

### Frontend Changes

**`frontend/src/types/validation.ts` — Mirror schema:**

Add matching optional fields to `MoleculeInfo` interface.

**`frontend/src/hooks/useMoleculeInfo.ts` — Replace regex with `get_descriptors()`:**

- Delete `kekulizeSmiles()` function entirely
- Replace all regex fallbacks with `mol.get_descriptors()` which returns `NumRings`, `NumAromaticRings`, `NumHeavyAtoms`, `NumAtoms`, `NumRotatableBonds`, etc.
- Keep kekulized SMILES via RDKit.js `get_smiles()` with kekulize option; return `null` if unavailable (no regex fallback)

**`frontend/src/pages/SingleValidation.tsx` — Single source post-validation:**

- Post-validation tiles: use ALL fields from `result.molecule_info` (backend)
- Pre-validation tiles: use `useMoleculeInfo` (now powered by `get_descriptors()`)
- Remove redundant `/[/\\]/.test()` regex for E/Z detection

### What Gets Removed

- `kekulizeSmiles()` regex function (60 lines)
- All SMILES regex fallbacks (`/\d/g`, `/[A-Z][a-z]?/g`, `/@+/g`, `/[/\\]/`)
- Mixed data source pattern in validated view

### What Stays Unchanged

- `useRDKit.ts`, `useMolecule.ts` — unrelated
- Alerts/scoring `MoleculeInfoSchema` — separate PR
- Batch processing — unaffected

## Testing

- Verify with isotope molecules: `[63Ni]`, `[13C]`, `[2H]`, `[99Tc]`
- Verify aromatic molecules: `c1ccccc1`, ferrocene
- Verify normal molecules: aspirin, imatinib
- Verify pre-validation preview matches post-validation backend values
