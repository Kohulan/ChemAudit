---
name: chemaudit-dataset-intelligence
description: >
  Audit dataset health with a Fourches-style 5-component score, detect
  contradictory bioactivity labels across duplicates, compare two dataset
  versions (added / removed / modified / unchanged), and generate curation
  reports. Use when user says "audit this dataset", "dataset health score",
  "contradictory labels", "duplicate activity conflict", "dataset diff",
  "compare dataset versions", "curation report", "data quality check", or
  "clean this dataset for ML". Accepts CSV and SDF with optional activity
  columns.
license: MIT
allowed-tools: "Bash(curl:*) Bash(python:*)"
compatibility: >
  Requires a running ChemAudit instance with Celery workers for async audits.
metadata:
  author: Kohulan Rajan
  version: 1.0.0
  mcp-server: chemaudit
  category: cheminformatics
  tags: [dataset-curation, quality-control, ML-datasets, bioactivity, data-cleaning]
---

# ChemAudit Dataset Intelligence

## Overview

Whole-dataset auditing — operates on the file, not individual molecules. Four outputs:

1. **Health audit** — composite 0-100 score with 5-component breakdown (Fourches 2010).
2. **Contradictory labels** — duplicate InChIKeys with inconsistent activity values.
3. **Dataset diff** — structural and property-level comparison between two uploads.
4. **Curation report + curated CSV** — JSON report and downloadable cleaned CSV with appended diagnostic columns.

## Health audit — 5 sub-scores (Fourches 2010)

| Sub-score | Default weight | What it measures |
|---|---|---|
| `parsability` | 0.25 | Fraction of rows with a valid, parseable SMILES |
| `stereo` | 0.15 | Fraction with defined stereocenters where chirality is expected |
| `uniqueness` | 0.20 | 1 - (duplicates / total) by canonical InChIKey |
| `alerts` | 0.20 | Fraction of molecules alert-free against the default catalogs |
| `std_consistency` | 0.20 | Fraction agreeing across three standardization pipelines (RDKit vs ChEMBL vs minimal) on a sampled subset |

Weights are customizable per audit; source of literature-backed defaults is `DEFAULT_WEIGHTS` in `frontend/src/types/dataset_intelligence.ts`.

Overall score is the weighted sum × 100.

## Contradictory labels — fold-difference thresholds

For every InChIKey appearing ≥ 2 times with a numeric activity column, reports the max/min fold-difference. Threshold options: **3, 5, 10, 50, 100** (default: **10**). Entries above the threshold are flagged as contradictory.

## Workflow

### 1. Upload dataset

```bash
curl -sS -X POST http://localhost:8000/api/v1/dataset/upload \
  -F "file=@bioactivity.csv" \
  -F "smiles_column=canonical_smiles" \
  -F "activity_column=pIC50"
```

Both column names are optional — auto-detection runs if omitted. Accepts `.csv` and `.sdf`. Response:

```json
{
  "job_id": "3f7a9b...",
  "filename": "bioactivity.csv",
  "file_type": "csv",
  "status": "pending",
  "message": "Dataset audit job submitted"
}
```

### 2. Track progress

WebSocket:

```javascript
const ws = new WebSocket(`ws://localhost:8000/ws/dataset/${job_id}`);
```

Poll:

```bash
curl -sS http://localhost:8000/api/v1/dataset/<job_id>/status
```

Response: `{job_id, status, progress, current_stage, eta_seconds}`. Stages: parse → health_audit → contradictions → standardization_sample → diff_ready.

### 3. Fetch audit results

```bash
curl -sS http://localhost:8000/api/v1/dataset/<job_id>/results
```

While processing, returns HTTP 202 with `{job_id, status, message}`. When `complete`, returns:

```json
{
  "job_id": "...",
  "status": "complete",
  "health_audit": {
    "overall_score": 78.3,
    "sub_scores": [
      {"name": "parsability", "score": 0.985, "weight": 0.25, "count": 15, "total": 1000},
      {"name": "stereo", "score": 0.62, "weight": 0.15, "count": 380, "total": 1000},
      {"name": "uniqueness", "score": 0.93, "weight": 0.20, "count": 70, "total": 1000},
      {"name": "alerts", "score": 0.81, "weight": 0.20, "count": 190, "total": 1000},
      {"name": "std_consistency", "score": 0.88, "weight": 0.20, "count": 12, "total": 100}
    ],
    "weights": {"parsability": 0.25, "stereo": 0.15, "uniqueness": 0.20, "alerts": 0.20, "std_consistency": 0.20},
    "molecule_count": 1000,
    "issues": [
      {"row_index": 42, "smiles": "C(C)(C)(C)(C)C", "issue_type": "parse_error", "severity": "error", "description": "..."}
    ],
    "property_distributions": {
      "mw": {"bins": [...], "counts": [...]},
      "logp": {"bins": [...], "counts": [...]},
      "tpsa": {"bins": [...], "counts": [...]}
    },
    "std_pipeline_comparison": {...},
    "std_sample_size": 100,
    "dedup_groups": [{"inchikey": "...", "rows": [12, 88]}]
  },
  "contradictions": [
    {
      "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
      "entries": [
        {"row_index": 5, "smiles": "...", "activity": 7.2},
        {"row_index": 203, "smiles": "...", "activity": 4.1}
      ],
      "fold_difference": 1259.0,
      "entry_count": 2,
      "smiles": "CC(=O)Oc1ccccc1C(=O)O"
    }
  ],
  "numeric_columns": [
    {"name": "pIC50", "priority": 1},
    {"name": "MW", "priority": 2}
  ],
  "curation_report": {...},
  "curated_csv_available": true
}
```

### 4. Compare two dataset versions (diff)

Upload the comparison file against a completed primary job:

```bash
curl -sS -X POST http://localhost:8000/api/v1/dataset/<job_id>/diff \
  -F "file=@bioactivity_v2.csv"
```

Response:

```json
{
  "added": [{"inchikey": "...", "smiles": "...", "row_index": 12, "properties": {...}, "changes": []}],
  "removed": [...],
  "modified": [{"inchikey": "...", "changes": [{"column": "pIC50", "old_value": 6.2, "new_value": 7.1}]}],
  "added_count": 45,
  "removed_count": 12,
  "modified_count": 89,
  "unchanged_count": 854,
  "unique_columns_primary": 2,
  "unique_columns_comparison": 3
}
```

Diff uses InChIKey as the primary key, so it's robust to SMILES notation differences.

### 5. Downloads

```bash
# Full curation report JSON
curl -sS http://localhost:8000/api/v1/dataset/<job_id>/download/report \
  -o curation_report.json

# Curated CSV — original rows + appended diagnostic columns
curl -sS http://localhost:8000/api/v1/dataset/<job_id>/download/csv \
  -o curated.csv
```

Curated CSV appends: `_health_issues` (joined issue types), `_is_duplicate` (bool), `_standardized_smiles`, `_alert_flags` (joined catalog names).

## Examples

### Example 1 — "Audit this ChEMBL bioactivity CSV for quality issues"

1. `POST /dataset/upload` with `smiles_column` and `activity_column` pointing at your CSV.
2. Poll `/status` or WS until `complete`.
3. `GET /results` → read `health_audit.overall_score` (aim ≥ 70) and `health_audit.sub_scores[]` for the weakest component.
4. Read `contradictions[]` (keeps only groups above the default 10-fold threshold — check `entry_count` and `fold_difference` per group).
5. If any `issues` have `severity="error"`, those rows need manual review.

### Example 2 — "Compare v1 and v2 of a training set"

1. Upload v1 and wait for `complete` → capture `job_id_v1`.
2. `POST /dataset/<job_id_v1>/diff` with the v2 file.
3. Read `added_count`, `removed_count`, `modified_count` for an overview.
4. For details on changed labels: `modified[]` contains `changes[]` arrays — look for activity-column changes.

### Example 3 — "Generate a curation report for my supervisor"

1. Audit as in example 1.
2. `GET /download/report` → JSON you can paste into an appendix or render as markdown.
3. `GET /download/csv` → drop-in replacement dataset with diagnostic columns.

## Rate limits

- `/dataset/upload`: 3/min.
- `/dataset/<job_id>/status`: 60/min.
- `/dataset/<job_id>/results`: 30/min.
- `/dataset/<job_id>/diff`: 10/min.
- `/dataset/<job_id>/download/*`: 10/min.

## Troubleshooting

### 400 "Invalid file type. Please upload a CSV or SDF file."

Only `.csv` and `.sdf` extensions accepted — `.tsv`, `.txt`, etc. are rejected by dataset-intelligence (unlike batch validation which accepts all four).

### 413 "File exceeds maximum size limit"

Check `/api/v1/config` for the deployment limit. Split the file or redeploy with a larger profile.

### 202 response on `/results`

Job is still running. Keep polling `/status` (lighter endpoint) until `status` flips to `complete` or `error`.

### 400 "Primary dataset job is not yet complete"

Diff endpoint needs the primary job to finish first. Wait until the primary shows `status=complete` in `/status`.

### `activity_column` auto-detected wrong column

Specify it explicitly. The auto-detector ranks by name heuristics and numeric-column detection; numeric columns are returned in `numeric_columns[]` sorted by `priority`.

### `contradictions` empty when you expect conflicts

The default fold threshold is 10. Lower-fold conflicts aren't reported unless you customize the threshold (currently set server-side; contact admin to change). Check `entries[]` length — duplicates below threshold are listed but won't appear as contradictions.

### `std_consistency` score very low

All three standardization pipelines (RDKit vs ChEMBL vs minimal) disagreed for many molecules. Usually means the dataset has unusual edge cases (metal complexes, large fragments, unusual charges). Look at `std_pipeline_comparison` for per-pipeline statistics.

### No `curated_csv_available` in results

The Celery worker failed mid-curation-report step. Re-upload. If it persists, check worker logs.

## Further reading

- `chemaudit-qsar-ready` — once the audit reveals issues, curate per molecule with provenance.
- `chemaudit-standardization` — per-molecule investigation of why two SMILES share an InChIKey.
- `chemaudit-diagnostics` — pre-flight file integrity checks before upload.
