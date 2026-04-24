---
name: chemaudit-batch-validation
description: >
  Process CSV, TSV, TXT, or SDF files of molecules through ChemAudit's batch
  pipeline with Redis-backed progress tracking, on-demand analytics (scaffold,
  chemical space, clustering, MMP, taxonomy, R-group), and nine export formats.
  Use when user says "validate this file", "batch validation", "process these
  molecules", "upload CSV of SMILES", "export results", "scaffold analysis",
  "cluster compounds", "Butina clustering", or "PDF report for these molecules".
  Deployment-profile-gated file size and molecule-count limits — query /config
  first.
license: MIT
allowed-tools: "Bash(curl:*) Bash(python:*) Bash(jq:*)"
compatibility: >
  Requires a running ChemAudit instance with Redis and Celery workers
  (docker compose up). WebSocket support requires the browser or a WS client.
metadata:
  author: Kohulan Rajan
  version: 1.0.0
  mcp-server: chemaudit
  category: cheminformatics
  tags: [batch-processing, chemistry, CSV, SDF, clustering, scaffold, analytics]
---

# ChemAudit Batch Validation

## Overview

File-in, file-out bulk validation with three phases:

1. **Upload → Celery enqueue** — instant `job_id` response.
2. **Processing** — per-molecule validation + optional scoring, safety, profiling, standardization. Progress via WebSocket or `/status` polling.
3. **Analytics (opt-in)** — scaffold, chemical space, clustering, MMP, taxonomy, similarity search, R-group decomposition. Computed on demand after the base job completes.

All state lives in Redis keyed by `job_id`. Session cookies scope access (one owner per job).

## Deployment limits — check first

```bash
curl -sS http://localhost:8000/api/v1/config
```

Returns `{"limits": {"max_batch_size": N, "max_file_size_mb": M}, "deployment_profile": "medium"}`. Profiles (small/medium/large/xl/coconut) set these via `config/*.yml`.

## Workflow

### 1. CSV column detection (optional, CSV/TSV/TXT only)

Before upload, probe a CSV for SMILES and Name columns:

```bash
curl -sS -X POST http://localhost:8000/api/v1/batch/detect-columns \
  -F "file=@compounds.csv"
```

Returns `{columns, suggested_smiles, suggested_name, column_samples, row_count_estimate, file_size_mb}`.

### 2. Upload

```bash
curl -sS -X POST http://localhost:8000/api/v1/batch/upload \
  -F "file=@compounds.csv" \
  -F "smiles_column=SMILES" \
  -F "name_column=Name" \
  -F "include_extended_safety=false" \
  -F "include_chembl_alerts=false" \
  -F "include_standardization=true" \
  -F "include_profiling=false" \
  -F "include_safety_assessment=false" \
  -F "profile_id=" \
  -F "notification_email="
```

Accepts `.sdf`, `.csv`, `.tsv`, `.txt`. Response:

```json
{
  "job_id": "9e3b2c1e-....",
  "status": "pending",
  "total_molecules": 450,
  "message": "Job submitted. Processing 450 molecules."
}
```

Form fields:
- `smiles_column` — CSV only; defaults to `SMILES`.
- `name_column` — CSV only; optional.
- `include_extended_safety` — add NIH and ZINC catalogs.
- `include_chembl_alerts` — add 7 ChEMBL pharma filter sets (BMS, Dundee, Glaxo, Inpharmatica, LINT, MLSMR, SureChEMBL).
- `include_standardization` — run ChEMBL standardization pipeline per molecule.
- `include_profiling` — PFI, #stars, Abbott bioavailability, CNS MPO.
- `include_safety_assessment` — CYP, hERG, bRo5, REOS, complexity panel.
- `profile_id` — scoring profile ID; adds profile desirability score per molecule.
- `notification_email` — sends a completion email (validated server-side).

### 3. Track progress

**WebSocket (preferred):**

```javascript
const ws = new WebSocket(`ws://localhost:8000/ws/batch/${job_id}`);
ws.onmessage = (e) => console.log(JSON.parse(e.data));
// Send "ping" to keep alive; server replies "pong".
```

Messages: `{job_id, status, progress, processed, total, eta_seconds}`.

**Polling:**

```bash
curl -sS http://localhost:8000/api/v1/batch/<job_id>/status
```

### 4. Paginated results

```bash
curl -sS "http://localhost:8000/api/v1/batch/<job_id>?page=1&page_size=50&status_filter=error&min_score=0&max_score=40"
```

Filters: `status_filter` (`success`/`error`), `min_score`/`max_score` (0–100), `sort_by` (index/name/smiles/score/qed/safety/status/issues), `sort_dir` (asc/desc), `issue_filter` (failed check name), `alert_filter` (catalog name).

Stats-only endpoint (cheap, no results array):

```bash
curl -sS http://localhost:8000/api/v1/batch/<job_id>/stats
```

### 5. On-demand analytics

Get current analytics status and any cached results:

```bash
curl -sS http://localhost:8000/api/v1/batch/<job_id>/analytics
```

Trigger a specific computation:

```bash
curl -sS -X POST http://localhost:8000/api/v1/batch/<job_id>/analytics/scaffold
curl -sS -X POST http://localhost:8000/api/v1/batch/<job_id>/analytics/chemical_space \
  -H 'Content-Type: application/json' -d '{"method": "tsne"}'
curl -sS -X POST http://localhost:8000/api/v1/batch/<job_id>/analytics/clustering \
  -H 'Content-Type: application/json' -d '{"cutoff": 0.35}'
curl -sS -X POST http://localhost:8000/api/v1/batch/<job_id>/analytics/mmp \
  -H 'Content-Type: application/json' -d '{"activity_column": "pIC50"}'
curl -sS -X POST http://localhost:8000/api/v1/batch/<job_id>/analytics/similarity_search \
  -H 'Content-Type: application/json' -d '{"query_smiles": "CCO", "top_k": 20}'
curl -sS -X POST http://localhost:8000/api/v1/batch/<job_id>/analytics/rgroup \
  -H 'Content-Type: application/json' -d '{"core_smarts": "c1ccccc1"}'
curl -sS -X POST http://localhost:8000/api/v1/batch/<job_id>/analytics/taxonomy
```

Supported analysis types: `scaffold`, `chemical_space`, `mmp`, `similarity_search`, `rgroup`, `clustering`, `taxonomy`.

**Clustering cap**: Butina clustering is hard-capped at 1,000 molecules (D-06). Subsample or filter first on larger jobs.

Triggering an analysis that's already cached returns `already_complete` unless params differ (e.g. different `method` for chemical_space or different `cutoff` for clustering — those trigger recompute).

### 6. Pairwise MCS comparison

Between any two molecules in the batch by index:

```bash
curl -sS -X POST http://localhost:8000/api/v1/batch/<job_id>/mcs \
  -H 'Content-Type: application/json' \
  -d '{"index_a": 0, "index_b": 17}'
```

Returns MCS SMARTS, Tanimoto similarity, matched atom/bond counts, property deltas.

### 7. Subset actions

Apply an operation to a hand-picked subset from the UI or another filter pass:

```bash
# Revalidate — creates new batch with the chosen indices
curl -sS -X POST http://localhost:8000/api/v1/batch/<job_id>/subset/revalidate \
  -H 'Content-Type: application/json' \
  -d '{"indices": [0, 3, 17, 42]}'

# Rescore with a profile — creates new batch or scores inline
curl -sS -X POST http://localhost:8000/api/v1/batch/<job_id>/subset/rescore \
  -H 'Content-Type: application/json' \
  -d '{"indices": [0, 3, 17], "profile_id": 4}'

# Inline rescore (no Celery, sync response)
curl -sS -X POST http://localhost:8000/api/v1/batch/<job_id>/subset/score-inline \
  -H 'Content-Type: application/json' \
  -d '{"indices": [0, 3, 17], "profile_id": 4}'

# Export just the subset
curl -sS -X POST http://localhost:8000/api/v1/batch/<job_id>/subset/export \
  -H 'Content-Type: application/json' \
  -d '{"indices": [0, 3, 17], "format": "sdf"}' -o subset.sdf
```

### 8. Export — 9 formats

```bash
curl -sS "http://localhost:8000/api/v1/batch/<job_id>/export?format=excel&include_images=true&sheet_layout=multi" \
  -o results.xlsx

curl -sS "http://localhost:8000/api/v1/batch/<job_id>/export?format=pdf&include_audit=true&sections=validation_summary,score_distribution" \
  -o report.pdf
```

Formats (details in `references/export-formats.md`):

| Format | Content |
|---|---|
| `csv` | All validation, scoring, and alert columns |
| `excel` | Multi-sheet XLSX, optional 2D images, conditional coloring |
| `sdf` | MOL blocks with properties; `include_audit=true` for full audit |
| `json` | Full result objects with nested metadata |
| `pdf` | Professional report; charts, stats, images; section-selectable |
| `fingerprint` | ZIP with Morgan/MACCS/RDKit fingerprints as CSV, .npy, .npz |
| `dedup` | ZIP with deduplication summary and per-group annotated CSVs |
| `scaffold` | CSV with Murcko scaffold SMILES and scaffold-group assignment |
| `property_matrix` | ZIP with flat CSV + multi-sheet Excel of all properties |

Filters: `score_min`, `score_max`, `status` (`success`/`error`/`warning`), `indices=0,1,5,23` (GET) or POST body `{"indices": [...]}` for large selections.

### 9. Cancel

```bash
curl -sS -X DELETE http://localhost:8000/api/v1/batch/<job_id>
```

Already-completed results are retained.

## Examples

### Example 1 — "Upload compounds.csv, poll, export failures as SDF"

1. `POST /api/v1/batch/upload` with the file → capture `job_id`.
2. Poll `GET /api/v1/batch/<job_id>/status` until `status == "complete"`.
3. `GET /api/v1/batch/<job_id>/export?format=sdf&status=error` → save to `failures.sdf`.

### Example 2 — "Scaffold distribution for this library"

1. Upload, wait for complete.
2. `POST /api/v1/batch/<job_id>/analytics/scaffold`.
3. Poll `GET /api/v1/batch/<job_id>/analytics` until the `scaffold` key moves to `status: complete`.
4. Read `scaffold.groups[]` and render a Murcko-scaffold frequency chart.

### Example 3 — "Cluster at Tanimoto 0.7 and export cluster centroids"

1. Ensure the batch has ≤ 1,000 molecules.
2. `POST /api/v1/batch/<job_id>/analytics/clustering` with `{"cutoff": 0.30}` (distance, not similarity — 0.30 distance ≈ 0.70 similarity).
3. Read `clustering.clusters[]`; each contains `centroid_index`.
4. Subset-export those indices as CSV.

## Rate limits

- `/batch/upload`: 3/minute.
- `/batch/detect-columns`: 10/minute.
- `/batch/<job_id>` (paginated results): 60/minute.
- `/batch/<job_id>/status`, `/stats`, DELETE `/batch/<job_id>`: 10/minute.
- `/batch/<job_id>/analytics` (status poll): 120/minute.
- `/batch/<job_id>/analytics/<type>` (trigger): 10/minute.
- `/batch/<job_id>/mcs`: 30/minute.
- `/batch/<job_id>/subset/revalidate`, `/subset/rescore`, `/subset/export`: 10/minute.
- `/batch/<job_id>/subset/score-inline`: 30/minute.
- `/batch/<job_id>/export` (GET and POST): 30/minute.

API key via `X-API-Key` raises the per-minute tier.

## Troubleshooting

### 400 "File too large" / 400 "File contains N molecules. Maximum allowed is M."

Check `/api/v1/config` — deployment limits are profile-driven. Split the file or redeploy with a larger profile (`./deploy.sh large`).

### 400 "No valid molecules found in file"

CSV column wrong, or every row has a parse error. Run `/batch/detect-columns` first, or pre-validate with `/diagnostics/file-prevalidate` (see `chemaudit-diagnostics`).

### 400 "Invalid file content. Please check the file format and column names."

Malicious-content scan matched a pattern (script tags, macros). Legitimate hits are rare — audit-logged server-side. Sanitize the source file.

### 400 "Invalid notification email address"

Pattern-validated; must be a normal `user@domain.tld` under 254 chars.

### 404 on `/batch/<job_id>` a few minutes after completion

Results TTL in Redis is `BATCH_RESULT_TTL` (default 1 hour). Download exports before expiry or increase the TTL.

### 400 "Clustering is limited to 1,000 molecules"

D-06 hard cap. Filter (by score, alert, or scaffold) before triggering clustering.

### WebSocket closes with code 4003 / 4029

- 4003: session cookie doesn't own the job (uploaded anonymously in a different browser).
- 4029: too many simultaneous WS connections for one job (`_WS_MAX_PER_JOB` = 10).

### Analytics status stuck at `computing`

The Celery worker crashed mid-task. Trigger the same analysis again — the idempotency guard resets `computing` to `queued` when params change; otherwise flush the status key manually.

## Further reading

- `references/export-formats.md` — exhaustive spec of all nine export formats.
- `chemaudit-qsar-ready` — for ML-dataset curation pipelines on the same file input.
- `chemaudit-structure-filter` — for generative-chemistry funnel filtering on a SMILES list.
- `chemaudit-dataset-intelligence` — for dataset health auditing with activity columns.
