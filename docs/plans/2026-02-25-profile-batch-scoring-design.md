# Profile-Integrated Batch Scoring

**Date**: 2026-02-25
**Branch**: `features_v1`
**Status**: Approved

## Problem

The "Apply" button on the Profiles page is cosmetic — it sets local React state (`activeProfileId`) that is lost on navigation and has zero effect on scoring. The profile system has full CRUD but no functional integration with batch or single-molecule scoring.

## Solution

Move the profile picker into the batch validation flow as a sidebar, and wire profile-based scoring into batch processing using desirability-based MPO scoring.

## Architecture

### 1. Backend — Add Missing Properties to Batch Results

**File**: `backend/app/services/batch/tasks.py`

The batch result dict currently stores `mw`, `logp`, `hbd`, `hba`, `tpsa`, `fsp3` but is missing `rotatable_bonds` and `aromatic_rings`. These are needed for profile threshold evaluation.

**Change**: Add `rotatable_bonds` (from Veber result) and `aromatic_rings` (compute via `Lipinski.NumAromaticRings(mol)`) to the druglikeness section of batch results.

### 2. Backend — Profile-Aware Scoring in Batch Tasks

**File**: `backend/app/services/batch/tasks.py`, `backend/app/schemas/batch.py`

- Add optional `profile_id` to `BatchUploadRequest` schema
- Pass `profile_id` through to `_process_single_molecule()` via `safety_options`
- If `profile_id` is present:
  1. Fetch profile thresholds/weights from DB (cache for the job duration)
  2. For each property in the profile, compute a **desirability score** (0-1):
     - 1.0 if value is within [min, max] range
     - Smooth linear falloff outside the range (e.g., 10% beyond max → 0.5)
     - 0.0 if value is >50% beyond the range boundary
  3. Combine per-property desirability scores using **weighted geometric mean** (QED-style)
  4. Scale to 0-100 → "Profile Score"
  5. Store in `result["scoring"]["profile"]`:
     ```json
     {
       "profile_id": 1,
       "profile_name": "Drug-like (Lipinski)",
       "score": 85.2,
       "properties": {
         "mw": {"value": 342.1, "min": 150, "max": 500, "in_range": true, "desirability": 1.0},
         "logp": {"value": 5.3, "min": -2, "max": 5, "in_range": false, "desirability": 0.7}
       }
     }
     ```

### 3. Frontend — Profile Sidebar in Batch Upload

**Files**: `frontend/src/pages/BatchValidation.tsx`, new `frontend/src/components/batch/ProfileSidebar.tsx`

- Remove `/profiles` as a standalone route (or keep as redirect)
- Add a collapsible left sidebar/drawer on the batch validation page
- Toggled by a "Scoring Profile" toggle in the upload section
- Shows preset + custom profiles (reuse `PresetPicker` grid)
- "Customize" button opens `ProfileBuilder` inline or as a modal
- Selected `profile_id` sent with batch job submission
- Sidebar state persisted in `BatchCacheContext`

### 4. Frontend — Profile Results in Batch Table

**Files**: `frontend/src/components/batch/BatchResultsTable.tsx`, `frontend/src/components/batch/BatchAnalyticsPanel.tsx`

- New "Profile Score" column in results table (conditional — only when job used a profile)
- Color-coded badge: green (>=80), amber (50-79), red (<50)
- Row detail expansion shows per-property pass/fail breakdown
- Analytics panel: histogram of profile score distribution
- Filter: "Profile compliant" / "Non-compliant" toggle

## Data Flow

```
Upload page → user toggles profile sidebar → picks a profile
  → clicks "Start Processing"
  → POST /batch/upload with profile_id in request body
  → Celery task receives profile_id
  → _process_single_molecule() fetches profile thresholds (cached)
  → computes per-property desirability + weighted MPO score
  → stores in result["scoring"]["profile"]
  → Frontend receives results with profile scores
  → BatchResultsTable shows Profile Score column
  → BatchAnalyticsPanel shows profile score distribution
```

## Properties Supported

| Key | Source in Batch Results | Profile Threshold |
|-----|----------------------|-------------------|
| `mw` | `druglikeness.mw` | min/max Da |
| `logp` | `druglikeness.logp` | min/max |
| `tpsa` | `druglikeness.tpsa` | min/max Å² |
| `hbd` | `druglikeness.hbd` | min/max |
| `hba` | `druglikeness.hba` | min/max |
| `rotatable_bonds` | **NEW** — add to druglikeness | min/max |
| `aromatic_rings` | **NEW** — add to druglikeness | min/max |
| `fsp3` | `admet.fsp3` | min/max |

## Desirability Function

Based on QED methodology (Bickerton et al., 2012):

```python
def desirability(value, min_val, max_val):
    """Compute 0-1 desirability for a property value given range."""
    if min_val <= value <= max_val:
        return 1.0
    range_width = max(max_val - min_val, 1e-6)
    if value < min_val:
        distance = (min_val - value) / range_width
    else:
        distance = (value - max_val) / range_width
    return max(0.0, 1.0 - distance)

def profile_score(properties, thresholds, weights):
    """Weighted geometric mean of desirabilities → 0-100."""
    scores = []
    total_weight = 0
    for key, thresh in thresholds.items():
        if key in properties and properties[key] is not None:
            d = desirability(properties[key], thresh["min"], thresh["max"])
            w = weights.get(key, 1.0)
            scores.append(d ** w)
            total_weight += w
    if not scores or total_weight == 0:
        return None
    geometric_mean = product(scores) ** (1.0 / total_weight)
    return round(geometric_mean * 100, 1)
```

## What Stays Unchanged

- Profile CRUD API (already functional)
- Existing scoring (ML-readiness, QED, safety filters, ADMET)
- Single validation page (no profile integration)
- Batch analytics (scaffold, chemical space) — unaffected

## Implementation Order

1. Backend: add `rotatable_bonds` + `aromatic_rings` to batch results
2. Backend: add `profile_id` to batch upload schema + task processing
3. Backend: implement desirability scoring function
4. Frontend: create `ProfileSidebar` component
5. Frontend: wire sidebar into `BatchValidation.tsx` upload flow
6. Frontend: add Profile Score column to `BatchResultsTable`
7. Frontend: add profile score distribution to analytics
8. Frontend: remove standalone `/profiles` route (or redirect)
