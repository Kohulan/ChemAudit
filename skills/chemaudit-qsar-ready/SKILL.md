---
name: chemaudit-qsar-ready
description: >
  Run molecules through ChemAudit's configurable 10-step QSAR-Ready curation
  pipeline for ML-ready SMILES with full per-step provenance and InChIKey-change
  detection. Use when user says "QSAR ready", "prepare for QSAR", "curate for
  ML", "standardize dataset for modeling", "strip salts and desalt", "remove
  duplicates by InChIKey", "canonicalize tautomers for ML", "prepare SMILES for
  machine learning", "QSAR-2D curation", or "QSAR-3D curation". Supports single
  molecules and batch CSV/SDF files with pasted SMILES fallback.
license: MIT
allowed-tools: "Bash(curl:*) Bash(python:*)"
compatibility: >
  Requires a running ChemAudit instance. Batch mode requires Celery workers.
metadata:
  author: Kohulan Rajan
  version: 1.0.0
  mcp-server: chemaudit
  category: cheminformatics
  tags: [QSAR, machine-learning, curation, standardization, dataset-preparation]
---

# ChemAudit QSAR-Ready Pipeline

## Overview

Configurable 10-step curation pipeline that produces ML-ready canonical SMILES with duplicate detection. Per-step provenance records what changed, enabling audit trails.

Pipeline steps (all toggleable except `parse` and `canonical`):

| # | Step | Config flag | Default |
|---|---|---|---|
| 1 | `parse` — parse input SMILES | always on | on |
| 2 | `metals` — `MetalDisconnector` | `enable_metals` | on |
| 3 | `desalt` — `LargestFragmentChooser` | `enable_desalt` | on |
| 4 | `normalize` — functional-group `Normalizer` | `enable_normalize` | on |
| 5 | `neutralize` — `Uncharger` | `enable_neutralize` | on |
| 6 | `tautomer` — `TautomerEnumerator` canonicalization | `enable_tautomer` | on |
| 7 | `stereo` — strip stereochemistry | `enable_stereo_strip` | **off** |
| 8 | `isotope` — strip isotope labels | `enable_isotope_strip` | on |
| 9 | `filter` — heavy-atom / MW / inorganic composition filter | `min_heavy_atoms`, `max_heavy_atoms`, `max_mw`, `remove_inorganics` | 3–100 atoms, ≤1500 Da, no inorganics |
| 10 | `canonical` — canonical SMILES + InChIKey | always on | on |

**QSAR-2D vs QSAR-3D**: set `enable_stereo_strip: true` for 2D descriptors (all stereo removed). Keep it `false` for 3D workflows that need R/S/E/Z.

## InChIKey change detection (D-14)

Every run reports:
- `original_inchikey` — from parsed input, before any steps.
- `standardized_inchikey` — from final SMILES, after all steps.
- `inchikey_changed` — `true` if the two differ.

`inchikey_changed=true` is the canonical signal that the pipeline modified the structure in a way that matters for deduplication and downstream ML.

## Status values

- `ok` — successfully curated.
- `rejected` — failed the composition filter (too few atoms, too heavy, all-inorganic).
- `duplicate` — InChIKey matches an earlier molecule in the same batch.
- `error` — pipeline raised an exception (unparseable input, toolkit failure).

## Workflow

### 1. Single molecule

```bash
curl -sS -X POST http://localhost:8000/api/v1/qsar-ready/single \
  -H 'Content-Type: application/json' \
  -d '{
    "smiles": "CC(=O)Oc1ccccc1C(=O)O.[Na+].[Cl-]",
    "config": {
      "enable_metals": true,
      "enable_desalt": true,
      "enable_normalize": true,
      "enable_neutralize": true,
      "enable_tautomer": true,
      "enable_stereo_strip": false,
      "enable_isotope_strip": true,
      "min_heavy_atoms": 3,
      "max_heavy_atoms": 100,
      "max_mw": 1500.0,
      "remove_inorganics": true
    }
  }'
```

Response shape:

```json
{
  "original_smiles": "CC(=O)Oc1ccccc1C(=O)O.[Na+].[Cl-]",
  "original_inchikey": "XXXXXXXX-XXXXXXXXXX-X",
  "curated_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "standardized_inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
  "inchikey_changed": true,
  "status": "ok",
  "rejection_reason": null,
  "steps": [
    {"step_name": "parse", "step_index": 1, "enabled": true, "status": "applied",
     "before_smiles": null, "after_smiles": "...", "detail": null},
    {"step_name": "metals", "step_index": 2, "enabled": true, "status": "no_change", ...},
    {"step_name": "desalt", "step_index": 3, "enabled": true, "status": "applied",
     "detail": "Removed 2 fragments: Na+, Cl-"},
    ...
  ]
}
```

Per-step `status` values: `applied`, `no_change`, `skipped`, `error`.

### 2. Batch — file upload

CSV/SDF file:

```bash
curl -sS -X POST http://localhost:8000/api/v1/qsar-ready/batch/upload \
  -F 'config={"enable_metals":true,"enable_desalt":true,"enable_normalize":true,"enable_neutralize":true,"enable_tautomer":true,"enable_stereo_strip":false,"enable_isotope_strip":true,"min_heavy_atoms":3,"max_heavy_atoms":100,"max_mw":1500.0,"remove_inorganics":true}' \
  -F "file=@compounds.csv"
```

`config` is form-encoded JSON. File is auto-detected (`.csv`, `.sdf`, `.sd`). Phase 9 file pre-validation runs first (M END terminators, encoding, empty rows) — critical issues return HTTP 422 with the pre-validation report.

### 3. Batch — pasted SMILES (D-03)

When you don't have a file:

```bash
curl -sS -X POST http://localhost:8000/api/v1/qsar-ready/batch/upload \
  -F 'config={"enable_metals":true,"enable_desalt":true,...}' \
  -F "smiles_text=CCO
CC(=O)O
c1ccccc1
CC(=O)Oc1ccccc1C(=O)O"
```

One SMILES per line.

### 4. Track progress

WebSocket:

```javascript
const ws = new WebSocket(`ws://localhost:8000/ws/qsar/${job_id}`);
ws.onmessage = (e) => console.log(JSON.parse(e.data));
```

Or poll:

```bash
curl -sS http://localhost:8000/api/v1/qsar-ready/batch/<job_id>/status
```

### 5. Paginated results

```bash
curl -sS "http://localhost:8000/api/v1/qsar-ready/batch/<job_id>/results?page=1&per_page=50"
```

Response includes `summary` with counts: `total`, `ok`, `rejected`, `duplicate`, `error`, and `steps_applied_counts` (per-step count of `applied` outcomes).

### 6. Download — CSV / SDF / JSON

```bash
curl -sS "http://localhost:8000/api/v1/qsar-ready/batch/<job_id>/download/csv" -o curated.csv
curl -sS "http://localhost:8000/api/v1/qsar-ready/batch/<job_id>/download/sdf" -o curated.sdf
curl -sS "http://localhost:8000/api/v1/qsar-ready/batch/<job_id>/download/json" -o full_provenance.json
```

CSV columns: `original_smiles`, `curated_smiles`, `original_inchikey`, `curated_inchikey`, `inchikey_changed`, `status`, `rejection_reason`, `steps_applied` (comma-joined).

SDF uses `curated_smiles` as structure, `original_inchikey` as the molecule title, and attaches `original_smiles`, `curated_smiles`, `status` as properties. Rejected and duplicate entries are skipped.

JSON is the full per-molecule provenance dump plus `summary`, `config`, and a `duplicates` list.

## Examples

### Example 1 — "Curate this SMILES and show what changed"

1. `POST /api/v1/qsar-ready/single` with `config=default`.
2. Print `original_smiles`, `curated_smiles`, `inchikey_changed`.
3. Walk `steps[]`; for each step where `status == "applied"`, print `step_name → detail`.

### Example 2 — "QSAR-2D curation on this SDF, download cleaned CSV"

1. `POST /api/v1/qsar-ready/batch/upload` with `config.enable_stereo_strip=true` and the file.
2. WS or poll until `status=complete`.
3. `GET /api/v1/qsar-ready/batch/<job_id>/results` → check `summary` (counts per status).
4. `GET /api/v1/qsar-ready/batch/<job_id>/download/csv` → 2D-stripped, deduplicated SMILES.

### Example 3 — "How many molecules changed InChIKey after standardization?"

From the results summary isn't enough — you need per-molecule detail:

1. `GET /api/v1/qsar-ready/batch/<job_id>/results?page=1&per_page=500` (or iterate pages).
2. Count `[r for r in results if r["inchikey_changed"]]`.
3. Slice: `[r["original_smiles"] for r in results if r["inchikey_changed"] and r["status"] == "ok"]`.

## Rate limits

- `/qsar-ready/single`: 30/min.
- `/qsar-ready/batch/upload`: 3/min.
- `/qsar-ready/batch/<job_id>/status`: 60/min.
- `/qsar-ready/batch/<job_id>/results`: 30/min.
- `/qsar-ready/batch/<job_id>/download/<format>`: 10/min.

## Troubleshooting

### 400 "Invalid config JSON"

`config` must be a JSON-encoded string in the multipart form. Quote it once, don't double-encode.

### 400 "Either a file or smiles_text must be provided"

Batch upload needs one of the two. Pasted SMILES goes under the `smiles_text` form field, not `text` or `smiles`.

### 422 "File has critical pre-validation issues"

Response body has a `prevalidation` field with the `FilePreValidationResponse` payload. Most common causes: missing `M  END` in SDF, encoding mismatch in CSV, empty file. Fix the file and re-upload.

### 400 "No valid SMILES found in input"

File parsed OK but every row had a `parse_error`. Ensure the `SMILES` column exists (CSV auto-detects only that exact name) or pre-validate with `/diagnostics/file-prevalidate`.

### 400 "Input contains N molecules. Maximum allowed is M."

Hit the deployment-profile batch cap. `GET /api/v1/config` to read the limit, split the file, or redeploy with a larger profile.

### status = "rejected" with `rejection_reason = "Molecule has no carbon atoms"`

The composition filter (step 9) rejected an inorganic. Set `remove_inorganics: false` in config to keep them.

### status = "duplicate"

Another earlier molecule in the same batch canonicalized to the same InChIKey. Intentional — keeps only the first occurrence. Duplicates are written to the JSON download under `duplicates[]`.

### Tautomer canonicalization stripped E/Z double-bond stereo

Known behavior of `TautomerEnumerator`. If that matters, set `enable_tautomer: false` — but expect that molecules with different tautomeric forms of the same compound will produce distinct InChIKeys.

## Further reading

- `chemaudit-standardization` — single-molecule standardization with richer provenance.
- `chemaudit-dataset-intelligence` — health audit across a dataset, complementary to curation.
- `chemaudit-diagnostics` — pre-flight file validation.
