---
phase: 03-batch-analytics
plan: 04
subsystem: backend-analytics
tags: [chemical-space, pca, tsne, similarity, morgan-fingerprints, opentsne]
requirements: [BATCH-09, BATCH-10, BATCH-11, BATCH-12]

dependency_graph:
  requires: ["03-01"]
  provides: ["chemical_space.py service", "chemical space test coverage"]
  affects: ["analytics_tasks.py (chemical_space branch)", "batch analytics GET response"]

tech_stack:
  added:
    - openTSNE>=1.0.0 (t-SNE dimensionality reduction)
  patterns:
    - Randomized SVD PCA: numpy-only, seeded RNG, no sklearn dependency
    - Size-gated matrix representation: dense (<=500), sparse (500-2000), refused (>2000)
    - t-SNE perplexity auto-adjustment: min(perplexity, max(5, n//3))
    - Non-deprecated Morgan fingerprints: rdFingerprintGenerator.GetMorganGenerator

key_files:
  created:
    - backend/app/services/analytics/chemical_space.py
    - backend/tests/test_analytics_chemspace.py
  modified:
    - backend/pyproject.toml

decisions:
  - "PCA uses randomized numpy SVD (no sklearn) — seeded with np.random.default_rng(42) for reproducibility"
  - "compute_chemical_space dispatcher added for analytics_tasks.py integration (method param routes to pca/tsne)"
  - "pytest 'slow' mark registered in pyproject.toml to suppress unknown-mark warnings from t-SNE tests"

metrics:
  duration: "3m 27s"
  completed: "2026-02-23"
  tasks_completed: 2
  files_changed: 3
---

# Phase 03 Plan 04: Chemical Space Analysis Summary

Chemical space analysis service: randomized-SVD PCA, openTSNE t-SNE (capped at 2000 molecules), Tanimoto similarity search by SMILES or batch index, nearest-neighbor isolation scores, and size-gated dense/sparse pairwise similarity matrix.

## What Was Built

### chemical_space.py

Five public functions plus a Celery dispatcher:

- **compute_pca**: Randomized SVD PCA via pure numpy (no sklearn). Centers data, performs random projection, QR decomposition, then full SVD on the reduced matrix. Returns 2D coordinates + variance explained. Seeded with `np.random.default_rng(42)` for reproducibility. Requires >= 3 molecules.

- **compute_tsne**: openTSNE TSNE with hard 2000-molecule cap. Auto-adjusts perplexity via `min(perplexity, max(5, n//3))` to prevent crashes on small batches. `n_jobs=-1` for parallel execution.

- **find_similar_molecules**: Tanimoto similarity search via `DataStructs.BulkTanimotoSimilarity`. Accepts external `query_smiles` or `query_index` (excludes query itself from results). Returns top-k neighbors sorted descending.

- **compute_nearest_neighbors**: Per-molecule nearest neighbor with `isolation_score = 1.0 - max_neighbor_similarity`. O(n^2) BulkTanimotoSimilarity calls.

- **compute_similarity_matrix**: Pairwise upper-triangle matrix. Dense list-of-lists for ≤500 molecules; sparse `{i, j, similarity}` dicts for 500-2000; error dict refusal for >2000.

- **compute_chemical_space**: Dispatcher for `analytics_tasks.py` — routes to `compute_pca` or `compute_tsne` based on `method` param, returns `ChemSpaceCoordinates` schema object.

### Shared fingerprint utility

Module-level `_morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)` — non-deprecated API, ECFP4 equivalent.

### Test coverage (test_analytics_chemspace.py)

13 tests in total:
- 3 PCA tests (shape, error guard, reproducibility)
- 3 t-SNE tests (2 marked `@pytest.mark.slow`, 1 fast refusal test)
- 2 similarity search tests (by SMILES, by index)
- 1 nearest neighbor test (isolation score math)
- 3 similarity matrix tests (dense, sparse, refusal)
- 1 API correctness test (verifies non-deprecated Morgan fingerprint API)

**Non-slow test results: 11/11 passed in 3.18s**

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 2 - Missing Config] Registered pytest 'slow' mark in pyproject.toml**
- **Found during:** Task 2 test run
- **Issue:** `@pytest.mark.slow` produced `PytestUnknownMarkWarning` since the mark was not registered
- **Fix:** Added `markers = ["slow: marks tests as slow..."]` to `[tool.pytest.ini_options]` in pyproject.toml
- **Files modified:** `backend/pyproject.toml`
- **Commit:** 5fa254f

**2. [Rule 2 - Integration] Added compute_chemical_space dispatcher**
- **Found during:** Task 1 implementation
- **Issue:** `analytics_tasks.py` calls `compute_chemical_space(results, method=method)` but the plan only specified `compute_pca`/`compute_tsne` as exports
- **Fix:** Added `compute_chemical_space` function that routes to `compute_pca` or `compute_tsne` and returns `ChemSpaceCoordinates` schema object
- **Files modified:** `backend/app/services/analytics/chemical_space.py`
- **Commit:** ab46c1a

## Self-Check: PASSED

| Item | Status |
|------|--------|
| `backend/app/services/analytics/chemical_space.py` | FOUND |
| `backend/tests/test_analytics_chemspace.py` | FOUND (288 lines) |
| `.planning/phases/03-batch-analytics/03-04-SUMMARY.md` | FOUND |
| Commit ab46c1a (Task 1) | FOUND |
| Commit 5fa254f (Task 2) | FOUND |
| 11 non-slow tests passing | PASSED |
