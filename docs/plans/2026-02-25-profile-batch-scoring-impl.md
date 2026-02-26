# Profile-Integrated Batch Scoring — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Wire scoring profiles into batch processing so each molecule gets a desirability-based "Profile Score" (0-100) when a profile is selected at upload time.

**Architecture:** The upload endpoint fetches profile thresholds/weights from DB and passes them through `safety_options` to Celery tasks (no async-in-Celery). A new `profile_scoring.py` module computes per-property desirability scores and a weighted geometric mean. Frontend adds a profile picker sidebar to the batch upload flow, a Profile Score column to results, and a score distribution chart.

**Tech Stack:** Python/FastAPI/Celery/RDKit (backend), React/TypeScript/Recharts/Tailwind (frontend)

---

### Task 1: Backend — Add `rotatable_bonds` + `aromatic_rings` to batch results

**Files:**
- Modify: `backend/app/services/batch/tasks.py:344-356`
- Test: `backend/tests/test_api/test_batch_properties.py` (create)

**Step 1: Write the failing test**

Create `backend/tests/test_api/test_batch_properties.py`:

```python
"""Test that batch results include rotatable_bonds and aromatic_rings."""
from app.services.batch.tasks import _process_single_molecule


def test_druglikeness_includes_rotatable_bonds_and_aromatic_rings():
    """rotatable_bonds and aromatic_rings should appear in druglikeness dict."""
    mol_data = {"smiles": "c1ccc(CC(=O)O)cc1", "name": "test", "index": 0}
    result = _process_single_molecule(mol_data)
    dl = result["scoring"]["druglikeness"]
    assert "rotatable_bonds" in dl
    assert "aromatic_rings" in dl
    assert isinstance(dl["rotatable_bonds"], int)
    assert isinstance(dl["aromatic_rings"], int)
    # Phenylacetic acid: 1 aromatic ring, 2 rotatable bonds
    assert dl["aromatic_rings"] == 1
    assert dl["rotatable_bonds"] >= 1
```

**Step 2: Run test to verify it fails**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -m pytest backend/tests/test_api/test_batch_properties.py -v`
Expected: FAIL — `rotatable_bonds` not in dl

**Step 3: Add properties to `_process_single_molecule`**

In `backend/app/services/batch/tasks.py`, modify lines 347-356. After `"tpsa"` (line 355), add `rotatable_bonds` and `aromatic_rings`:

```python
        try:
            dl_result = calculate_druglikeness(mol, include_extended=False)
            result["scoring"]["druglikeness"] = {
                "qed_score": dl_result.qed.score,
                "lipinski_passed": dl_result.lipinski.passed,
                "lipinski_violations": dl_result.lipinski.violations,
                "mw": dl_result.lipinski.mw,
                "logp": dl_result.lipinski.logp,
                "hbd": dl_result.lipinski.hbd,
                "hba": dl_result.lipinski.hba,
                "tpsa": dl_result.veber.tpsa if dl_result.veber else None,
                "rotatable_bonds": dl_result.veber.rotatable_bonds if dl_result.veber else None,
                "aromatic_rings": dl_result.qed.properties.get("arom", None) if dl_result.qed and dl_result.qed.properties else None,
            }
```

> **Note:** If `dl_result.qed.properties` doesn't have `"arom"`, fall back to computing directly:
> ```python
> from rdkit.Chem import Lipinski
> "aromatic_rings": Lipinski.NumAromaticRings(mol),
> ```

**Step 4: Run test to verify it passes**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -m pytest backend/tests/test_api/test_batch_properties.py -v`
Expected: PASS

**Step 5: Update frontend type**

In `frontend/src/types/batch.ts`, add to the `druglikeness` interface (after line 64):

```typescript
    druglikeness?: {
      qed_score: number;
      lipinski_passed: boolean;
      lipinski_violations: number;
      mw?: number;
      logp?: number;
      hbd?: number;
      hba?: number;
      tpsa?: number;
      rotatable_bonds?: number;
      aromatic_rings?: number;
    };
```

**Step 6: Commit**

```bash
git add backend/app/services/batch/tasks.py backend/tests/test_api/test_batch_properties.py frontend/src/types/batch.ts
git commit -m "feat: add rotatable_bonds and aromatic_rings to batch results"
```

---

### Task 2: Backend — Create profile scoring module

**Files:**
- Create: `backend/app/services/scoring/profile_scoring.py`
- Test: `backend/tests/test_services/test_profile_scoring.py` (create)

**Step 1: Write the failing test**

Create `backend/tests/test_services/test_profile_scoring.py`:

```python
"""Test desirability and profile scoring functions."""
import pytest
from app.services.scoring.profile_scoring import desirability, profile_score


class TestDesirability:
    def test_value_in_range_returns_one(self):
        assert desirability(300, 0, 500) == 1.0

    def test_value_at_boundary_returns_one(self):
        assert desirability(500, 0, 500) == 1.0
        assert desirability(0, 0, 500) == 1.0

    def test_value_above_max_linear_falloff(self):
        # 50% beyond range width → d = 0.5
        result = desirability(750, 0, 500)
        assert result == pytest.approx(0.5, abs=0.01)

    def test_value_below_min_linear_falloff(self):
        # 100 below min of 200, range width 300 → distance = 100/300 ≈ 0.33
        result = desirability(100, 200, 500)
        assert result == pytest.approx(0.667, abs=0.01)

    def test_value_far_outside_returns_zero(self):
        # > 100% beyond range → clamped to 0
        assert desirability(1100, 0, 500) == 0.0

    def test_none_value_returns_none(self):
        assert desirability(None, 0, 500) is None


class TestProfileScore:
    def test_all_in_range_returns_100(self):
        properties = {"mw": 300, "logp": 2.5}
        thresholds = {"mw": {"min": 0, "max": 500}, "logp": {"min": -5, "max": 5}}
        weights = {"mw": 1.0, "logp": 1.0}
        assert profile_score(properties, thresholds, weights) == 100.0

    def test_one_out_of_range_reduces_score(self):
        properties = {"mw": 300, "logp": 7.5}
        thresholds = {"mw": {"min": 0, "max": 500}, "logp": {"min": -5, "max": 5}}
        weights = {"mw": 1.0, "logp": 1.0}
        score = profile_score(properties, thresholds, weights)
        assert 0 < score < 100

    def test_missing_property_skipped(self):
        properties = {"mw": 300}
        thresholds = {"mw": {"min": 0, "max": 500}, "logp": {"min": -5, "max": 5}}
        weights = {"mw": 1.0, "logp": 1.0}
        # Only mw evaluated, in range → 100
        assert profile_score(properties, thresholds, weights) == 100.0

    def test_empty_thresholds_returns_none(self):
        assert profile_score({"mw": 300}, {}, {}) is None

    def test_weights_affect_score(self):
        properties = {"mw": 600, "logp": 2.5}
        thresholds = {"mw": {"min": 0, "max": 500}, "logp": {"min": -5, "max": 5}}
        # Heavy weight on mw (which is out of range) → lower score
        score_heavy = profile_score(properties, thresholds, {"mw": 3.0, "logp": 1.0})
        score_light = profile_score(properties, thresholds, {"mw": 1.0, "logp": 1.0})
        assert score_heavy < score_light
```

**Step 2: Run test to verify it fails**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -m pytest backend/tests/test_services/test_profile_scoring.py -v`
Expected: FAIL — ModuleNotFoundError

**Step 3: Implement profile scoring module**

Create `backend/app/services/scoring/profile_scoring.py`:

```python
"""
Profile-based desirability scoring.

Computes per-property desirability (0-1) and a weighted geometric mean
profile score (0-100) following QED-style methodology (Bickerton et al., 2012).
"""

from math import prod
from typing import Any, Dict, Optional


def desirability(
    value: Optional[float],
    min_val: float,
    max_val: float,
) -> Optional[float]:
    """Compute 0-1 desirability for a property value given an acceptable range.

    Returns 1.0 if value is within [min_val, max_val].
    Linear falloff outside the range, reaching 0.0 at 100% beyond the range width.
    """
    if value is None:
        return None
    if min_val <= value <= max_val:
        return 1.0
    range_width = max(max_val - min_val, 1e-6)
    if value < min_val:
        distance = (min_val - value) / range_width
    else:
        distance = (value - max_val) / range_width
    return max(0.0, 1.0 - distance)


def profile_score(
    properties: Dict[str, Any],
    thresholds: Dict[str, Dict[str, Any]],
    weights: Dict[str, float],
) -> Optional[float]:
    """Compute weighted geometric mean of desirabilities scaled to 0-100.

    Args:
        properties: Property name → value (e.g. {"mw": 342.1, "logp": 2.3})
        thresholds: Property name → {"min": float, "max": float}
        weights: Property name → weight (default 1.0)

    Returns:
        Score 0-100, or None if no evaluable properties.
    """
    scores = []
    total_weight = 0.0
    for key, thresh in thresholds.items():
        val = properties.get(key)
        if val is None:
            continue
        t_min = thresh.get("min", float("-inf"))
        t_max = thresh.get("max", float("inf"))
        d = desirability(val, t_min, t_max)
        if d is None:
            continue
        w = weights.get(key, 1.0)
        scores.append(d ** w)
        total_weight += w

    if not scores or total_weight == 0:
        return None

    geometric_mean = prod(scores) ** (1.0 / total_weight)
    return round(geometric_mean * 100, 1)


def compute_profile_result(
    properties: Dict[str, Any],
    profile_id: int,
    profile_name: str,
    thresholds: Dict[str, Dict[str, Any]],
    weights: Dict[str, float],
) -> Dict[str, Any]:
    """Build the full profile scoring result dict for a single molecule.

    Returns dict with score, profile metadata, and per-property breakdown.
    """
    score = profile_score(properties, thresholds, weights)
    prop_details = {}
    for key, thresh in thresholds.items():
        val = properties.get(key)
        t_min = thresh.get("min", float("-inf"))
        t_max = thresh.get("max", float("inf"))
        d = desirability(val, t_min, t_max) if val is not None else None
        in_range = (t_min <= val <= t_max) if val is not None else None
        prop_details[key] = {
            "value": val,
            "min": t_min if t_min != float("-inf") else None,
            "max": t_max if t_max != float("inf") else None,
            "in_range": in_range,
            "desirability": round(d, 3) if d is not None else None,
        }

    return {
        "profile_id": profile_id,
        "profile_name": profile_name,
        "score": score,
        "properties": prop_details,
    }
```

**Step 4: Run test to verify it passes**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -m pytest backend/tests/test_services/test_profile_scoring.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
git add backend/app/services/scoring/profile_scoring.py backend/tests/test_services/test_profile_scoring.py
git commit -m "feat: add desirability-based profile scoring module"
```

---

### Task 3: Backend — Wire profile into batch upload + Celery tasks

**Files:**
- Modify: `backend/app/api/routes/batch.py:59-198`
- Modify: `backend/app/services/batch/tasks.py:198-408`
- Test: `backend/tests/test_api/test_batch_profile_scoring.py` (create)

**Step 1: Write the failing test**

Create `backend/tests/test_api/test_batch_profile_scoring.py`:

```python
"""Test that profile scoring is applied when profile data is in safety_options."""
from app.services.batch.tasks import _process_single_molecule


def test_profile_score_in_result_when_profile_provided():
    """When safety_options includes profile data, result should have scoring.profile."""
    mol_data = {"smiles": "c1ccc(CC(=O)O)cc1", "name": "test", "index": 0}
    safety_options = {
        "profile_id": 1,
        "profile_name": "Drug-like (Lipinski)",
        "profile_thresholds": {
            "mw": {"min": 0, "max": 500},
            "logp": {"min": -5, "max": 5},
            "hbd": {"min": 0, "max": 5},
            "hba": {"min": 0, "max": 10},
        },
        "profile_weights": {"mw": 1.0, "logp": 1.0, "hbd": 1.0, "hba": 1.0},
    }
    result = _process_single_molecule(mol_data, safety_options=safety_options)
    profile = result["scoring"]["profile"]
    assert profile["profile_id"] == 1
    assert profile["profile_name"] == "Drug-like (Lipinski)"
    assert 0 <= profile["score"] <= 100
    assert "mw" in profile["properties"]
    assert profile["properties"]["mw"]["in_range"] is True


def test_no_profile_score_when_profile_not_provided():
    """When no profile data in safety_options, result should not have scoring.profile."""
    mol_data = {"smiles": "c1ccc(CC(=O)O)cc1", "name": "test", "index": 0}
    result = _process_single_molecule(mol_data)
    assert "profile" not in result.get("scoring", {})
```

**Step 2: Run test to verify it fails**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -m pytest backend/tests/test_api/test_batch_profile_scoring.py -v`
Expected: FAIL — `profile` not in scoring

**Step 3: Add profile scoring to `_process_single_molecule`**

In `backend/app/services/batch/tasks.py`:

1. Add import at top (after line 26):
```python
from app.services.scoring.profile_scoring import compute_profile_result
```

2. After the ADMET block (after line 400, before `result["status"] = "success"`), add profile scoring:

```python
        # Calculate profile score (if profile was selected at upload)
        opts = safety_options or {}
        if opts.get("profile_id") is not None:
            try:
                # Gather all computed properties for profile evaluation
                dl = result.get("scoring", {}).get("druglikeness", {})
                admet = result.get("scoring", {}).get("admet", {})
                mol_properties = {
                    "mw": dl.get("mw"),
                    "logp": dl.get("logp"),
                    "hbd": dl.get("hbd"),
                    "hba": dl.get("hba"),
                    "tpsa": dl.get("tpsa"),
                    "rotatable_bonds": dl.get("rotatable_bonds"),
                    "aromatic_rings": dl.get("aromatic_rings"),
                    "fsp3": admet.get("fsp3"),
                }
                result["scoring"]["profile"] = compute_profile_result(
                    properties=mol_properties,
                    profile_id=opts["profile_id"],
                    profile_name=opts.get("profile_name", ""),
                    thresholds=opts.get("profile_thresholds", {}),
                    weights=opts.get("profile_weights", {}),
                )
            except Exception as e:
                result["scoring"]["profile"] = {"error": str(e)}
```

**Step 4: Run test to verify it passes**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -m pytest backend/tests/test_api/test_batch_profile_scoring.py -v`
Expected: PASS

**Step 5: Add `profile_id` form param to upload endpoint**

In `backend/app/api/routes/batch.py`:

1. Add form parameter after `notification_email` (line 86-87):
```python
    profile_id: Optional[int] = Form(
        default=None,
        description="Scoring profile ID for profile-based desirability scoring",
    ),
```

2. After building `safety_options` dict (line 195), add profile fetching:
```python
    # If profile selected, fetch thresholds and attach to safety_options
    if profile_id is not None:
        from app.db.session import async_session_factory
        from app.services.profiles.service import ProfileService
        import json

        async with async_session_factory() as db:
            profile = await ProfileService().get(db, profile_id)
        if profile is None:
            raise HTTPException(status_code=404, detail=f"Profile {profile_id} not found")
        safety_options["profile_id"] = profile.id
        safety_options["profile_name"] = profile.name
        safety_options["profile_thresholds"] = json.loads(profile.thresholds) if isinstance(profile.thresholds, str) else profile.thresholds
        safety_options["profile_weights"] = json.loads(profile.weights) if isinstance(profile.weights, str) else profile.weights
```

**Step 6: Run all batch tests**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -m pytest backend/tests/test_api/test_batch_profile_scoring.py backend/tests/test_api/test_batch_properties.py backend/tests/test_services/test_profile_scoring.py -v`
Expected: All PASS

**Step 7: Run ruff check**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval/backend && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -m ruff check .`
Expected: No errors (fix any that appear)

**Step 8: Commit**

```bash
git add backend/app/api/routes/batch.py backend/app/services/batch/tasks.py backend/tests/test_api/test_batch_profile_scoring.py
git commit -m "feat: wire profile scoring into batch upload and Celery processing"
```

---

### Task 4: Backend — Add profile score to sort extractors + statistics

**Files:**
- Modify: `backend/app/services/batch/result_aggregator.py:237-254`

**Step 1: Add `profile_score` sort extractor**

In `backend/app/services/batch/result_aggregator.py`, add to `SORT_EXTRACTORS` dict (after the `issues` entry at line 253):

```python
        "profile_score": lambda r: ((r.get("scoring") or {}).get("profile") or {}).get(
            "score", -1
        ),
```

**Step 2: Add profile stats to `compute_statistics`**

In `compute_statistics()`, after extracting other averages (around line 115), add:

```python
        # Profile score stats
        profile_scores = []
        for r in successful:
            ps = ((r.get("scoring") or {}).get("profile") or {}).get("score")
            if ps is not None:
                profile_scores.append(ps)

        avg_profile_score = (
            round(sum(profile_scores) / len(profile_scores), 1)
            if profile_scores
            else None
        )
        profile_compliant_count = sum(1 for s in profile_scores if s >= 80)
        profile_compliance_rate = (
            round(profile_compliant_count / len(profile_scores) * 100, 1)
            if profile_scores
            else None
        )
```

Add these to the `BatchStatisticsData` dataclass and `BatchStatistics` Pydantic schema if needed, OR simply store them in the statistics dict as extra fields. Since `BatchStatistics` uses `Dict[str, Any]` patterns extensively, the simplest approach is to add them to the returned dataclass.

**Step 3: Run ruff check + existing tests**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval/backend && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -m ruff check .`

**Step 4: Commit**

```bash
git add backend/app/services/batch/result_aggregator.py
git commit -m "feat: add profile_score sort extractor and profile stats to aggregator"
```

---

### Task 5: Frontend — Add `profile_id` to upload API + types

**Files:**
- Modify: `frontend/src/services/api.ts:418-456`
- Modify: `frontend/src/types/batch.ts:51-82` and `101-119`
- Modify: `frontend/src/components/batch/BatchUpload.tsx:8-11,143-199`

**Step 1: Add profile types to `batch.ts`**

In `frontend/src/types/batch.ts`, add `profile` to the `scoring` interface (after `admet`, around line 80):

```typescript
    profile?: {
      profile_id: number;
      profile_name: string;
      score: number | null;
      properties: Record<string, {
        value: number | null;
        min: number | null;
        max: number | null;
        in_range: boolean | null;
        desirability: number | null;
      }>;
      error?: string;
    };
```

Add `'profile_score'` to `SortField` (line 161):
```typescript
export type SortField = 'index' | 'name' | 'smiles' | 'score' | 'qed' | 'safety' | 'status' | 'issues' | 'profile_score';
```

Add profile stats to `BatchStatistics` (after `processing_time_seconds`, line 119):
```typescript
  avg_profile_score?: number | null;
  profile_compliance_rate?: number | null;
```

**Step 2: Update `uploadBatch` API function**

In `frontend/src/services/api.ts`, modify the `uploadBatch` signature (line 418-423) to accept `profileId`:

```typescript
  uploadBatch: async (
    file: File,
    smilesColumn?: string,
    nameColumn?: string,
    onUploadProgress?: (progress: number) => void,
    safetyOptions?: { includeExtended?: boolean; includeChembl?: boolean; includeStandardization?: boolean },
    profileId?: number | null
  ): Promise<BatchUploadResponse> => {
```

Add to FormData construction (after line 441):
```typescript
    if (profileId != null) {
      formData.append('profile_id', String(profileId));
    }
```

**Step 3: Update `BatchUpload` component props**

In `frontend/src/components/batch/BatchUpload.tsx`, modify:

1. Props interface (line 8-12):
```typescript
interface BatchUploadProps {
  onUploadSuccess: (jobId: string, totalMolecules: number, options?: { includeAnalytics: boolean }) => void;
  onUploadError: (error: string) => void;
  disabled?: boolean;
  profileId?: number | null;
}
```

2. Destructure prop (line 33-37):
```typescript
export function BatchUpload({
  onUploadSuccess,
  onUploadError,
  disabled = false,
  profileId,
}: BatchUploadProps) {
```

3. Pass to `uploadBatch` call (line 188-198):
```typescript
      const response = await batchApi.uploadBatch(
        fileToUpload,
        csvColumns ? smilesColForUpload : undefined,
        csvColumns && nameColForUpload ? nameColForUpload : undefined,
        (progress) => setUploadProgress(progress),
        {
          includeExtended: includeExtendedSafety,
          includeChembl: includeChemblAlerts,
          includeStandardization: includeStandardization,
        },
        profileId
      );
```

**Step 4: Run TypeScript check**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval/frontend && npx tsc --noEmit`
Expected: No errors

**Step 5: Commit**

```bash
git add frontend/src/types/batch.ts frontend/src/services/api.ts frontend/src/components/batch/BatchUpload.tsx
git commit -m "feat: add profile_id to batch upload API and types"
```

---

### Task 6: Frontend — Create ProfileSidebar component

**Files:**
- Create: `frontend/src/components/batch/ProfileSidebar.tsx`

**Step 1: Create the ProfileSidebar component**

This is a collapsible panel that appears in the batch upload section when toggled. It fetches profiles from the API and shows a compact grid of presets with radio-style selection.

Create `frontend/src/components/batch/ProfileSidebar.tsx`:

```tsx
import { useState, useEffect } from 'react';
import { ChevronDown, ChevronRight, Beaker, Check } from 'lucide-react';
import { profilesApi } from '../../services/api';
import { ClayButton } from '../ui/ClayButton';
import { cn } from '../../lib/utils';
import type { ScoringProfile, ThresholdRange } from '../../types/workflow';

interface ProfileSidebarProps {
  selectedProfileId: number | null;
  onProfileChange: (profileId: number | null) => void;
  disabled?: boolean;
}

function formatRange(range: ThresholdRange | undefined): string {
  if (!range) return '--';
  if (range.min !== undefined && range.max !== undefined) return `${range.min}-${range.max}`;
  if (range.max !== undefined) return `≤${range.max}`;
  if (range.min !== undefined) return `≥${range.min}`;
  return '--';
}

export function ProfileSidebar({ selectedProfileId, onProfileChange, disabled }: ProfileSidebarProps) {
  const [isExpanded, setIsExpanded] = useState(selectedProfileId != null);
  const [profiles, setProfiles] = useState<ScoringProfile[]>([]);
  const [isLoading, setIsLoading] = useState(false);

  useEffect(() => {
    if (isExpanded && profiles.length === 0) {
      setIsLoading(true);
      profilesApi.getProfiles()
        .then(res => setProfiles(res.data))
        .catch(() => {})
        .finally(() => setIsLoading(false));
    }
  }, [isExpanded, profiles.length]);

  const handleToggle = () => {
    const next = !isExpanded;
    setIsExpanded(next);
    if (!next) onProfileChange(null);
  };

  const selected = profiles.find(p => p.id === selectedProfileId);

  return (
    <div className="bg-[var(--color-surface-elevated)] rounded-xl border border-[var(--color-border)] overflow-hidden">
      <button
        onClick={handleToggle}
        disabled={disabled}
        className="w-full p-4 flex items-center justify-between bg-[var(--color-surface-sunken)]/50 hover:bg-[var(--color-surface-sunken)] transition-colors"
      >
        <div className="flex items-center gap-2">
          <Beaker className="w-4 h-4 text-[var(--color-primary)]" />
          <span className="text-sm font-medium text-[var(--color-text-primary)]">Scoring Profile</span>
          {selected && (
            <span className="text-xs px-2 py-0.5 rounded-full bg-[var(--color-primary)]/10 text-[var(--color-primary)] font-medium">
              {selected.name}
            </span>
          )}
        </div>
        {isExpanded ? (
          <ChevronDown className="w-4 h-4 text-[var(--color-text-muted)]" />
        ) : (
          <ChevronRight className="w-4 h-4 text-[var(--color-text-muted)]" />
        )}
      </button>

      {isExpanded && (
        <div className="p-4 space-y-3">
          <p className="text-xs text-[var(--color-text-muted)]">
            Select a profile to score each molecule against property thresholds.
          </p>

          {isLoading ? (
            <div className="flex justify-center py-4">
              <div className="animate-spin w-6 h-6 border-2 border-[var(--color-primary)] border-t-transparent rounded-full" />
            </div>
          ) : (
            <div className="grid grid-cols-1 gap-2">
              {profiles.map(profile => {
                const isSelected = profile.id === selectedProfileId;
                const keys = Object.keys(profile.thresholds).slice(0, 4);
                return (
                  <button
                    key={profile.id}
                    onClick={() => onProfileChange(isSelected ? null : profile.id)}
                    disabled={disabled}
                    className={cn(
                      'w-full text-left p-3 rounded-lg border transition-all',
                      isSelected
                        ? 'border-[var(--color-primary)] bg-[var(--color-primary)]/5 ring-1 ring-[var(--color-primary)]'
                        : 'border-[var(--color-border)] hover:border-[var(--color-text-muted)] bg-[var(--color-surface-elevated)]'
                    )}
                  >
                    <div className="flex items-center justify-between">
                      <span className="text-sm font-medium text-[var(--color-text-primary)]">{profile.name}</span>
                      {isSelected && <Check className="w-4 h-4 text-[var(--color-primary)]" />}
                    </div>
                    <div className="flex flex-wrap gap-1.5 mt-1.5">
                      {keys.map(key => (
                        <span key={key} className="text-[10px] px-1.5 py-0.5 rounded bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)]">
                          {key}: {formatRange(profile.thresholds[key])}
                        </span>
                      ))}
                    </div>
                  </button>
                );
              })}
            </div>
          )}

          {selectedProfileId != null && (
            <ClayButton
              variant="ghost"
              size="sm"
              onClick={() => onProfileChange(null)}
              className="w-full"
            >
              Clear selection
            </ClayButton>
          )}
        </div>
      )}
    </div>
  );
}
```

**Step 2: Run TypeScript check**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval/frontend && npx tsc --noEmit`
Expected: No errors

**Step 3: Commit**

```bash
git add frontend/src/components/batch/ProfileSidebar.tsx
git commit -m "feat: create ProfileSidebar component for batch upload flow"
```

---

### Task 7: Frontend — Wire ProfileSidebar into BatchValidation upload flow

**Files:**
- Modify: `frontend/src/pages/BatchValidation.tsx`
- Modify: `frontend/src/contexts/BatchCacheContext.tsx`

**Step 1: Add `selectedProfileId` to BatchCacheContext**

In `frontend/src/contexts/BatchCacheContext.tsx`, add to `BatchCache` interface (after `includeAnalytics`, line 18):

```typescript
  selectedProfileId: number | null;
```

**Step 2: Wire ProfileSidebar into BatchValidation page**

In `frontend/src/pages/BatchValidation.tsx`:

1. Add import:
```typescript
import { ProfileSidebar } from '../components/batch/ProfileSidebar';
```

2. Add state (near other state declarations):
```typescript
const [selectedProfileId, setSelectedProfileId] = useState<number | null>(null);
```

3. Include `selectedProfileId` when persisting to cache and restoring from cache.

4. Add `<ProfileSidebar>` in the upload section (before the `<BatchUpload>` component, inside the upload pageState block):
```tsx
<ProfileSidebar
  selectedProfileId={selectedProfileId}
  onProfileChange={setSelectedProfileId}
  disabled={pageState !== 'upload'}
/>
```

5. Pass `profileId` to `<BatchUpload>`:
```tsx
<BatchUpload
  onUploadSuccess={handleUploadSuccess}
  onUploadError={handleError}
  profileId={selectedProfileId}
/>
```

**Step 3: Run TypeScript check**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval/frontend && npx tsc --noEmit`
Expected: No errors

**Step 4: Commit**

```bash
git add frontend/src/pages/BatchValidation.tsx frontend/src/contexts/BatchCacheContext.tsx
git commit -m "feat: wire ProfileSidebar into batch validation upload flow"
```

---

### Task 8: Frontend — Add Profile Score column to BatchResultsTable

**Files:**
- Modify: `frontend/src/components/batch/BatchResultsTable.tsx:170-239`

**Step 1: Add Profile Score column header**

In `BatchResultsTable.tsx`, after the Safety column header (line 219) and before the Status column header (line 221), add:

```tsx
              {results.some(r => r.scoring?.profile) && (
                <th
                  className="px-4 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase cursor-pointer hover:bg-[var(--color-surface-elevated)]"
                  onClick={() => handleSort('profile_score')}
                >
                  Profile {sortBy === 'profile_score' && (sortDir === 'asc' ? '↑' : '↓')}
                </th>
              )}
```

**Step 2: Add Profile Score cell in each row**

In the row rendering section (where QED, Safety, etc. cells are rendered), add a cell for the profile score. Use color coding: green (>=80), amber (50-79), red (<50):

```tsx
{results.some(r => r.scoring?.profile) && (
  <td className="px-4 py-3 text-center">
    {result.scoring?.profile?.score != null ? (
      <span className={cn(
        'inline-flex items-center px-2 py-0.5 rounded-full text-xs font-medium',
        result.scoring.profile.score >= 80
          ? 'bg-green-500/10 text-green-600 dark:text-green-400'
          : result.scoring.profile.score >= 50
          ? 'bg-amber-500/10 text-amber-600 dark:text-amber-400'
          : 'bg-red-500/10 text-red-600 dark:text-red-400'
      )}>
        {result.scoring.profile.score.toFixed(1)}
      </span>
    ) : (
      <span className="text-[var(--color-text-muted)]">--</span>
    )}
  </td>
)}
```

**Step 3: Update colspan for loading/empty states**

Update `colSpan={9}` references to account for the conditional extra column. The simplest approach: compute column count based on whether profile data exists.

**Step 4: Run TypeScript check**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval/frontend && npx tsc --noEmit`
Expected: No errors

**Step 5: Commit**

```bash
git add frontend/src/components/batch/BatchResultsTable.tsx
git commit -m "feat: add conditional Profile Score column to batch results table"
```

---

### Task 9: Frontend — Add profile score distribution chart to analytics

**Files:**
- Modify: `frontend/src/components/batch/BatchAnalyticsPanel.tsx`

**Step 1: Add profile score histogram**

In `BatchAnalyticsPanel.tsx`, add a new chart card in the Distributions tab (after the Alert Frequency card, around line 481). Use Recharts `BarChart` to show profile score distribution in buckets: Excellent (80-100), Good (50-79), Moderate (20-49), Poor (0-19).

The chart only renders if `results.some(r => r.scoring?.profile?.score != null)`.

```tsx
{/* VIZ: Profile Score Distribution (conditional) */}
{results.some(r => r.scoring?.profile?.score != null) && (
  <ChartCard
    title="Profile Score Distribution"
    icon={<Beaker className="w-4 h-4" />}
    ariaLabel="Profile score distribution histogram"
    caption="Distribution of profile desirability scores"
  >
    <ProfileScoreHistogram results={results} />
  </ChartCard>
)}
```

Create a simple inline `ProfileScoreHistogram` component that:
1. Buckets results by score range
2. Renders a `BarChart` with 4 bars
3. Color-codes: green/amber/orange/red

**Step 2: Import Beaker icon**

Add `Beaker` to lucide-react imports.

**Step 3: Run TypeScript check**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval/frontend && npx tsc --noEmit`
Expected: No errors

**Step 4: Commit**

```bash
git add frontend/src/components/batch/BatchAnalyticsPanel.tsx
git commit -m "feat: add profile score distribution chart to batch analytics"
```

---

### Task 10: Frontend — Remove standalone `/profiles` route

**Files:**
- Modify: `frontend/src/main.tsx`

**Step 1: Replace `/profiles` route with redirect**

In `frontend/src/main.tsx`, replace the `/profiles` route with a redirect to `/batch`:

```tsx
<Route path="/profiles" element={<Navigate to="/batch" replace />} />
```

Import `Navigate` from `react-router-dom` if not already imported. Remove the lazy import of `ProfilesPage` if it's no longer used elsewhere.

**Step 2: Remove `/profiles` from navigation**

If there's a nav link to `/profiles` in the header/sidebar, remove it or redirect it.

**Step 3: Run TypeScript check**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval/frontend && npx tsc --noEmit`
Expected: No errors

**Step 4: Commit**

```bash
git add frontend/src/main.tsx
git commit -m "refactor: redirect /profiles route to /batch"
```

---

### Task 11: Final verification

**Step 1: Run all backend tests**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -m pytest backend/tests/ -v --tb=short`

**Step 2: Run frontend type check**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval/frontend && npx tsc --noEmit`

**Step 3: Run ruff check**

Run: `cd /Volumes/Data_Drive/Project/2026/chemstructval/backend && /Users/kohulanrajan/anaconda3/envs/cheminformatics/bin/python -m ruff check .`

**Step 4: Manual smoke test**

1. Start backend + frontend dev servers
2. Navigate to `/batch`
3. Toggle "Scoring Profile" → select "Drug-like (Lipinski)"
4. Upload a small CSV (5-10 molecules)
5. Verify results show "Profile" column with scores
6. Verify analytics shows profile score distribution chart
7. Navigate to `/profiles` → should redirect to `/batch`
