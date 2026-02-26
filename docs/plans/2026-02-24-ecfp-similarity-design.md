# ECFP Similarity in Molecule Comparison

**Date**: 2026-02-24
**Branch**: `features_v1`

## Goal

When comparing 2 molecules in the MoleculeComparisonPanel, show their Tanimoto similarity score computed from ECFP4 (Morgan) fingerprints.

## Approach: New Lightweight API Endpoint

### Backend: `POST /validate/similarity`

**Request**:
```json
{
  "smiles_a": "CCO",
  "smiles_b": "CCCO"
}
```

**Response**:
```json
{
  "tanimoto_similarity": 0.78,
  "fingerprint_type": "ECFP4",
  "radius": 2,
  "n_bits": 2048,
  "common_bits": 45,
  "bits_a": 52,
  "bits_b": 58
}
```

**Implementation**:
- File: `backend/app/api/routes/validation.py` (add new route)
- Schema: `backend/app/schemas/validation.py` (add SimilarityRequest/Response)
- Logic: Use `AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)` + `DataStructs.TanimotoSimilarity()`
- Existing reference: `backend/app/services/analytics/chemical_space.py` already uses Morgan fingerprints

### Frontend: MoleculeComparisonPanel Enhancement

**When 2 molecules are loaded**:
1. Call `POST /validate/similarity` with both SMILES
2. Display Tanimoto score in the panel header between the two structures
3. Visual: circular gauge or progress bar with color coding:
   - >= 0.85: "Very Similar" (green)
   - >= 0.70: "Similar" (amber)
   - >= 0.50: "Moderate" (orange)
   - < 0.50: "Dissimilar" (red)
4. Show fingerprint details on hover (bits_a, bits_b, common_bits)

### Files to Modify

| File | Change |
|------|--------|
| `backend/app/api/routes/validation.py` | Add `POST /validate/similarity` route |
| `backend/app/schemas/validation.py` | Add `SimilarityRequest`, `SimilarityResponse` |
| `frontend/src/services/api.ts` | Add `getSimilarity()` API call |
| `frontend/src/components/batch/MoleculeComparisonPanel.tsx` | Fetch + display similarity score |

### Existing Backend Code Reference

`backend/app/services/analytics/chemical_space.py` line ~350:
- `compute_similarity_matrix()` already uses Morgan fingerprints + Tanimoto
- We can extract a simple pairwise function from this pattern

### Edge Cases
- One or both SMILES invalid → return error, panel shows "N/A"
- Identical molecules → 1.0 similarity
- Loading state while fetching
