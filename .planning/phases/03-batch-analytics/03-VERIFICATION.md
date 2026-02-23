---
phase: 03-batch-analytics
verified: 2026-02-23T17:30:00Z
status: gaps_found
score: 19/20 must-haves verified
re_verification: false
gaps:
  - truth: "Similarity search (analytics type 'similarity_search') is callable via POST /batch/{job_id}/analytics/similarity_search and executes correctly"
    status: failed
    reason: "run_expensive_analytics imports from app.services.analytics.similarity which does not exist. The actual function (find_similar_molecules) lives in chemical_space.py under a different name. At runtime the ImportError is caught by try/except and the job status is marked 'failed' — but the analysis type is silently broken."
    artifacts:
      - path: "backend/app/services/batch/analytics_tasks.py"
        issue: "Line 183 imports 'from app.services.analytics.similarity import compute_similarity_search' — this module does not exist"
      - path: "backend/app/services/analytics/similarity.py"
        issue: "File does not exist"
    missing:
      - "Create backend/app/services/analytics/similarity.py that exports compute_similarity_search, OR rename find_similar_molecules in chemical_space.py to compute_similarity_search and re-export it, OR update analytics_tasks.py to import find_similar_molecules from chemical_space"
human_verification:
  - test: "Submit a batch job, wait for it to complete, then GET /batch/{job_id}/analytics"
    expected: "HTTP 200 with status dict showing deduplication and statistics as 'complete', remaining types as 'pending'"
    why_human: "Requires running Redis + Celery workers; cannot verify task dispatch and storage round-trip with grep alone"
  - test: "POST /batch/{job_id}/analytics/similarity_search with a valid query_smiles param"
    expected: "HTTP 200 with status='queued'; subsequent GET shows status='failed' with ImportError message"
    why_human: "Documents the current broken behavior so it is observed before the gap fix"
---

# Phase 03: Batch Analytics Verification Report

**Phase Goal:** A completed batch job exposes a second analytics layer — multi-level deduplication groups, scaffold families, chemical space projections, MMP pairs, and statistical summaries — all computed asynchronously by a post-aggregation Celery chord, accessible via a dedicated analytics endpoint.

**Verified:** 2026-02-23T17:30:00Z
**Status:** gaps_found
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Batch results stored with 24h TTL (INFRA-01) | VERIFIED | `settings.BATCH_RESULT_TTL == 86400`; `ResultStorage.RESULT_EXPIRY = settings.BATCH_RESULT_TTL` (not hardcoded) |
| 2 | `get_all_results(job_id)` returns all raw results without pagination | VERIFIED | Method exists on `ResultStorage`, reads `batch:results:{job_id}`, returns `json.loads(data)` |
| 3 | `GET /batch/{job_id}/analytics` returns analytics status and cached results | VERIFIED | Endpoint registered at line 398 of `batch.py`, 30/min rate limit, 404 for unknown jobs |
| 4 | `POST /batch/{job_id}/analytics/{analysis_type}` triggers expensive analytics | VERIFIED | Endpoint registered at line 461 of `batch.py`, 10/min rate limit, dispatches `run_expensive_analytics.delay()` |
| 5 | `aggregate_batch_results` dispatches `run_cheap_analytics.delay(job_id)` after storing | VERIFIED | Lines 432-442 of `tasks.py`: try/except dispatch in both `aggregate_batch_results` and `aggregate_batch_results_priority` |
| 6 | Analytics results cached under `batch:analytics:{type}:{job_id}` keys with 24h TTL | VERIFIED | `AnalyticsStorage.ANALYTICS_TTL = 86400`; key format confirmed in `storage.py` |
| 7 | Multi-level deduplication returns exact, tautomeric, stereo_insensitive, salt_form groups | VERIFIED | `compute_all_dedup_levels` functional; all 4 dedup functions import and execute; 20 tests pass |
| 8 | Scaffold families grouped by Murcko scaffold with Shannon entropy metric | VERIFIED | `compute_scaffold_analysis` uses double-GetScaffoldForMol, entropy computed inline; 12 tests pass |
| 9 | Chemical space PCA 2D projection uses randomized numpy SVD | VERIFIED | `compute_pca` in `chemical_space.py` uses `np.linalg.svd` (no sklearn); seeded `np.random.default_rng(42)` |
| 10 | t-SNE projection via openTSNE, limited to 2000 molecules | VERIFIED | `compute_tsne` enforces `> 2000` refusal; `openTSNE>=1.0.0` in `pyproject.toml` |
| 11 | Similarity matrix: dense ≤500, sparse 500-2000, refused >2000 | VERIFIED | Logic confirmed in `compute_similarity_matrix`; tests pass for all three tiers |
| 12 | Similarity search executable via the analytics endpoint | FAILED | `run_expensive_analytics` imports `from app.services.analytics.similarity import compute_similarity_search` — this module does not exist |
| 13 | MMP pairs detected via BRICS fragmentation with 5000-molecule limit | VERIFIED | `MAX_MMP_BATCH_SIZE = 5000`; BRICS in `mmp.py`; 14 tests pass |
| 14 | Activity cliff SALI computed when activity_column provided; skipped otherwise | VERIFIED | `_compute_activity_cliffs` guards on activity dict; `activity_cliffs=None` when no column |
| 15 | LLE = pIC50 - LogP per molecule when activity_column provided | VERIFIED | `_compute_lle` calls `Descriptors.MolLogP`; formula confirmed in tests |
| 16 | Property distribution stats: mean, median, std, quartiles, IQR, min, max | VERIFIED | `compute_property_stats` uses numpy; all 5 properties extracted via `PROPERTY_EXTRACTORS` lambdas |
| 17 | Pearson correlation matrix with NaN guard | VERIFIED | `abs(pearson_r) < 1e-10` filter handles numpy floating-point noise from constant arrays |
| 18 | Batch quality score: 40% validity + 35% diversity + 25% drug-likeness | VERIFIED | Formula at line 352: `score = validity_pct * 0.40 + diversity_pct * 0.35 + druglikeness_pct * 0.25` |
| 19 | IQR-based outlier detection flags molecules outside [Q1-1.5*IQR, Q3+1.5*IQR] | VERIFIED | `compute_outliers` uses numpy Q1/Q3 and IQR fence; 21 statistics tests pass |
| 20 | All analytics are accessible via `BatchAnalyticsResponse` schema with all fields optional | VERIFIED | 18 Pydantic v2 schemas verified importable; `BatchAnalyticsResponse` has all optional fields |

**Score:** 19/20 truths verified

---

## Required Artifacts

| Artifact | Plan | Status | Details |
|----------|------|--------|---------|
| `backend/app/core/config.py` | 03-01 | VERIFIED | `BATCH_RESULT_TTL: int = 86400` present |
| `backend/app/services/batch/result_aggregator.py` | 03-01 | VERIFIED | `RESULT_EXPIRY = settings.BATCH_RESULT_TTL`; `get_all_results` method exists |
| `backend/app/services/batch/analytics_tasks.py` | 03-01 | VERIFIED | `run_cheap_analytics` and `run_expensive_analytics` Celery tasks registered |
| `backend/app/services/analytics/storage.py` | 03-01 | VERIFIED | `AnalyticsStorage` with all 5 methods; `analytics_storage` singleton |
| `backend/app/schemas/analytics.py` | 03-01 | VERIFIED | 18 Pydantic v2 schemas; `BatchAnalyticsResponse` with optional analytics fields |
| `backend/app/api/routes/batch.py` | 03-01 | VERIFIED | GET + POST analytics endpoints at `/batch/{job_id}/analytics` |
| `backend/app/services/analytics/deduplication.py` | 03-02 | VERIFIED | All 5 exports present; `compute_all_dedup_levels` returns `DeduplicationResult` |
| `backend/tests/test_analytics_dedup.py` | 03-02 | VERIFIED | 374 lines; 20 tests pass |
| `backend/app/services/analytics/scaffold_analysis.py` | 03-03 | VERIFIED | `compute_scaffold_analysis` and `compute_rgroup_decomposition` exported |
| `backend/tests/test_analytics_scaffold.py` | 03-03 | VERIFIED | 248 lines; 12 tests pass |
| `backend/pyproject.toml` | 03-04 | VERIFIED | `openTSNE>=1.0.0` in dependencies |
| `backend/app/services/analytics/chemical_space.py` | 03-04 | VERIFIED | 5 public functions + `compute_chemical_space` dispatcher; `_morgan_gen` uses non-deprecated API |
| `backend/tests/test_analytics_chemspace.py` | 03-04 | VERIFIED | 288 lines; 11/11 non-slow tests pass |
| `backend/app/services/analytics/mmp.py` | 03-05 | VERIFIED | `compute_mmp_analysis`, `compute_mmp` alias, `_MMPResultWrapper` for `.model_dump()` compatibility |
| `backend/tests/test_analytics_mmp.py` | 03-05 | VERIFIED | 360 lines; 14 tests pass |
| `backend/app/services/analytics/statistics.py` | 03-06 | VERIFIED | `compute_all_statistics` returns `StatisticsResult` Pydantic model |
| `backend/tests/test_analytics_statistics.py` | 03-06 | VERIFIED | 585 lines; 21 tests pass |
| `backend/app/services/analytics/similarity.py` | (none) | MISSING | Referenced by `analytics_tasks.py` but file was never created |

---

## Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `tasks.py` | `analytics_tasks.py` | `run_cheap_analytics.delay(job_id)` in both aggregators | WIRED | Lines 434-442 and 504-514 of `tasks.py` |
| `batch.py` | `storage.py` | `analytics_storage.get_status()` and `get_result()` | WIRED | Lines 413, 441 of `batch.py` |
| `result_aggregator.py` | `config.py` | `settings.BATCH_RESULT_TTL` for TTL value | WIRED | Line 191: `RESULT_EXPIRY = settings.BATCH_RESULT_TTL` |
| `analytics_tasks.py` | `deduplication.py` | `from app.services.analytics.deduplication import` in `run_cheap_analytics` | WIRED | Lines 59, 71 of `analytics_tasks.py` |
| `analytics_tasks.py` | `statistics.py` | `from app.services.analytics.statistics import` in `run_cheap_analytics` | WIRED | Lines 89, 101 of `analytics_tasks.py` |
| `analytics_tasks.py` | `scaffold_analysis.py` | `from app.services.analytics.scaffold_analysis import` in scaffold/rgroup branches | WIRED | Lines 163, 199 of `analytics_tasks.py` |
| `analytics_tasks.py` | `chemical_space.py` | `from app.services.analytics.chemical_space import compute_chemical_space` | WIRED | Line 169 of `analytics_tasks.py` |
| `analytics_tasks.py` | `mmp.py` | `from app.services.analytics.mmp import compute_mmp` | WIRED | Line 176 of `analytics_tasks.py` |
| `analytics_tasks.py` | `similarity.py` | `from app.services.analytics.similarity import compute_similarity_search` | NOT_WIRED | `similarity.py` does not exist; runtime ImportError caught, status set to "failed" |

---

## Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|------------|-------------|--------|----------|
| INFRA-01 | 03-01 | Batch result pagination and Redis TTL policy | SATISFIED | `BATCH_RESULT_TTL=86400` in config; `RESULT_EXPIRY` reads from config; `get_all_results` method functional |
| BATCH-01 | 03-02 | Exact duplicate detection (canonical SMILES matching) | SATISFIED | `compute_exact_dedup` groups by `Chem.MolToSmiles(mol)`; min-index representative; verified by test suite |
| BATCH-02 | 03-02 | Tautomeric duplicate detection | SATISFIED | `compute_tautomer_dedup` with per-call `TautomerEnumerator()` (thread-safe); canonical tautomer SMILES key |
| BATCH-03 | 03-02 | Stereoisomer-insensitive deduplication | SATISFIED | `compute_stereo_dedup` strips `/t /m /s` InChI layers via `frozenset` filter |
| BATCH-04 | 03-02 | Salt-form duplicate detection | SATISFIED | `compute_saltform_dedup` uses `get_parent_mol` from `chembl_structure_pipeline.standardizer`; fast-path for pre-standardized SMILES |
| BATCH-05 | 03-03 | Murcko scaffold decomposition and grouping | SATISFIED | `compute_scaffold_analysis` groups by Murcko scaffold SMILES; acyclic under `""` key |
| BATCH-06 | 03-03 | Generic scaffold grouping | SATISFIED | Double-`GetScaffoldForMol` pattern after `MakeScaffoldGeneric`; matches existing `scaffold.py` pattern |
| BATCH-07 | 03-03 | Scaffold diversity metrics | SATISFIED | `unique_scaffold_count`, `shannon_entropy`, `frequency_distribution` (capped at 50 + "Other") |
| BATCH-08 | 03-03 | R-group decomposition around common core | SATISFIED | `compute_rgroup_decomposition` validates SMARTS, calls `rdRGroupDecomposition.RGroupDecompose` |
| BATCH-09 | 03-04 | Chemical space PCA/t-SNE visualization | SATISFIED | `compute_pca` (randomized SVD), `compute_tsne` (openTSNE); both return `ChemSpaceCoordinates` schema |
| BATCH-10 | 03-04 | Similarity search within batch | PARTIAL | `find_similar_molecules` exists in `chemical_space.py` but the `similarity_search` analytics type in `analytics_tasks.py` imports from non-existent `similarity.py` — function unreachable via the API endpoint |
| BATCH-11 | 03-04 | Nearest neighbor analysis | SATISFIED | `compute_nearest_neighbors` returns isolation_score per molecule; tests pass |
| BATCH-12 | 03-04 | Similarity matrix/heatmap data | SATISFIED | `compute_similarity_matrix` with dense/sparse/refused tiers; tests pass |
| BATCH-13 | 03-05 | Matched molecular pair detection | SATISFIED | BRICS fragmentation; size heuristic; deduplication; 1000-pair cap |
| BATCH-14 | 03-05 | Activity cliff detection (SALI index) | SATISFIED | `_compute_activity_cliffs` computes `|delta_activity| / (1 - tanimoto)`; skips pairs without activity |
| BATCH-15 | 03-05 | Lipophilic ligand efficiency | SATISFIED | `_compute_lle` computes `activity_value - MolLogP(mol)`; None when no activity column |
| BATCH-16 | 03-06 | Property distribution statistics | SATISFIED | numpy mean/median/std/Q1/Q3/IQR/min/max for 5 properties via `PROPERTY_EXTRACTORS` |
| BATCH-17 | 03-06 | Property correlation matrix | SATISFIED | Pairwise Pearson with `abs < 1e-10` NaN guard; min 10 values per property |
| BATCH-18 | 03-06 | Batch quality score | SATISFIED | `score = validity_pct*0.40 + diversity_pct*0.35 + druglikeness_pct*0.25`; inline scaffold entropy |
| BATCH-19 | 03-06 | Outlier detection (IQR-based) | SATISFIED | `compute_outliers` uses `[Q1-1.5*IQR, Q3+1.5*IQR]` fences; min 4 values required |

**Note:** REQUIREMENTS.md shows BATCH-01 through BATCH-19 and INFRA-01 as "Pending" (checkbox and tracking table not updated post-execution). This is a documentation gap, not an implementation gap — all algorithms are implemented and tested.

---

## Anti-Patterns Found

| File | Pattern | Severity | Impact |
|------|---------|----------|--------|
| `analytics_tasks.py:183` | `from app.services.analytics.similarity import compute_similarity_search` — module does not exist | Blocker | `similarity_search` analysis type silently fails at runtime with ImportError; status set to "failed" |

No other stub, placeholder, or TODO anti-patterns detected. All analytics files contain substantive implementations with proper docstrings, type hints, and error handling.

---

## Test Results Summary

| Test File | Tests | Result |
|-----------|-------|--------|
| `test_infra_analytics.py` | 14 | All pass |
| `test_analytics_dedup.py` | 20 | All pass |
| `test_analytics_scaffold.py` | 12 | All pass |
| `test_analytics_chemspace.py` | 11 (non-slow) + 2 (slow/t-SNE) | 11/11 non-slow pass; slow tests skipped in CI |
| `test_analytics_mmp.py` | 14 | All pass |
| `test_analytics_statistics.py` | 21 | All pass |
| **Total** | **94** | **92/92 automated; 2 skipped (@slow)** |

---

## Human Verification Required

### 1. End-to-End Analytics Pipeline

**Test:** Submit a batch of 20+ molecules via `/api/v1/batch/upload`, wait for processing to complete (poll status until "complete"), then call `GET /api/v1/batch/{job_id}/analytics`.

**Expected:** HTTP 200 response with `status.deduplication == "complete"` and `status.statistics == "complete"`, results populated in the `deduplication` and `statistics` fields of `BatchAnalyticsResponse`.

**Why human:** Requires running Redis + Celery workers. The Celery chord dispatch (`run_cheap_analytics.delay`) and the analytics storage round-trip (`AnalyticsStorage.store_result` → Redis → `get_result`) cannot be verified without live infrastructure.

### 2. Similarity Search Broken Behavior

**Test:** POST to `/api/v1/batch/{job_id}/analytics/similarity_search` with `{"query_smiles": "CCO"}`. Then poll `GET /api/v1/batch/{job_id}/analytics`.

**Expected (current broken behavior):** HTTP 200 with `status="queued"` from POST; subsequent GET shows `similarity_search` status as `"failed"` with an error message about the missing module.

**Why human:** Documents the observable impact of the gap before it is fixed; confirms the try/except guard works as designed.

---

## Gaps Summary

One gap found: the `similarity_search` analytics type in `run_expensive_analytics` imports from `app.services.analytics.similarity` which was never created. The actual implementation (`find_similar_molecules`) lives in `chemical_space.py` under a different name.

**Root cause:** Plan 03-04 specified `find_similar_molecules` as the exported function name, but the analytics_tasks.py scaffold written in plan 03-01 used `compute_similarity_search` from a separate `similarity.py` module. Plan 03-04 did not detect the mismatch and did not create the missing module.

**Impact:** `POST /batch/{job_id}/analytics/similarity_search` queues successfully, but the Celery task fails immediately on import with `ModuleNotFoundError`. The job status is set to `"failed"` — no crash, but the feature is non-functional.

**Fix scope:** Small — either create `backend/app/services/analytics/similarity.py` with a `compute_similarity_search` wrapper around `find_similar_molecules`, or update the import in `analytics_tasks.py` to use `find_similar_molecules` from `chemical_space.py` directly.

All other analytics — deduplication, scaffold, chemical space (PCA/t-SNE), nearest neighbors, similarity matrix, MMP, statistics, and quality scoring — are fully implemented, wired, and tested.

---

_Verified: 2026-02-23T17:30:00Z_
_Verifier: Claude (gsd-verifier)_
