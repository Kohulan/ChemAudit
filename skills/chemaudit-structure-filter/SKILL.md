---
name: chemaudit-structure-filter
description: >
  Filter generative-chemistry outputs through ChemAudit's 6-stage multi-stage
  funnel with a composite 0-1 scorer and REINVENT 4-compatible scoring endpoint.
  Use when user says "filter molecules", "drug-like filter", "structure filter",
  "filter generative output", "REINVENT scoring", "REINVENT Component", "lead-like
  filtering", "fragment filter", "funnel these SMILES", "novelty filter", or
  "score generative batch". Four built-in presets (drug_like, lead_like,
  fragment_like, permissive). Sync for ≤1000 SMILES, async for larger.
license: MIT
allowed-tools: "Bash(curl:*) Bash(python:*)"
compatibility: >
  Requires a running ChemAudit instance. Async batch mode requires Celery workers.
metadata:
  author: Kohulan Rajan
  version: 1.0.0
  mcp-server: chemaudit
  category: cheminformatics
  tags: [generative-chemistry, filtering, drug-likeness, REINVENT, virtual-screening]
---

# ChemAudit Structure Filter

## Overview

Multi-stage funnel for generative-chemistry outputs and virtual screening hits. Six sequential stages, each of which can reject a molecule independently:

1. `parse` — RDKit parse + sanitization.
2. `valence` — valence validity check.
3. `alerts` — structural alert catalogs (PAINS, Brenk, Kazius, NIBR — configurable).
4. `property` — MW, LogP, TPSA, rotatable bonds, rings within bounds.
5. `sa_score` — synthetic accessibility score ≤ threshold.
6. `dedup` — InChIKey-based deduplication.

Optional 7th stage:
7. `novelty` — Tanimoto similarity vs ChEMBL (gated by `enable_novelty`).

Counts how many molecules enter and leave each stage (funnel diagram data).

Paired with a **composite 0-1 scorer** that blends validity / QED / alert-free / SA into a single reward, useful for generative-model reinforcement learning loops. Invalid SMILES score `null`, never `0.0` (D-14 / Pitfall 4) so generative agents can distinguish "unparseable" from "parseable but terrible".

REINVENT 4 Component API is supported out of the box.

## Four built-in presets

See `references/preset-configs.md` for full threshold details.

| Preset | MW | LogP | TPSA | RotB | Rings | SA | Alerts |
|---|---|---|---|---|---|---|---|
| `drug_like` | 200–500 | -1 to 5 | ≤140 | ≤10 | — | ≤5.0 | PAINS + Brenk + Kazius |
| `lead_like` | 200–350 | -1 to 3.5 | ≤140 | ≤7 | — | ≤4.0 | PAINS + Brenk + Kazius + NIBR |
| `fragment_like` | 100–300 | -1 to 3 | ≤100 | ≤3 | ≤3 | ≤3.0 | PAINS only |
| `permissive` | 100–800 | -5 to 8 | ≤200 | ≤15 | — | ≤7.0 | none |

Weight vectors for the composite score differ per preset (D-15):

| Preset | validity | qed | alert_free | sa |
|---|---|---|---|---|
| `drug_like` | 0.3 | 0.3 | 0.2 | 0.2 |
| `lead_like` | 0.2 | **0.4** | 0.2 | 0.2 |
| `fragment_like` | 0.2 | 0.2 | **0.3** | **0.3** |
| `permissive` | **0.4** | 0.3 | 0.1 | 0.2 |

## Sync vs async split (D-22)

- `POST /structure-filter/filter` with ≤ 1000 SMILES → sync `FilterResponse`.
- `POST /structure-filter/filter` with > 1000 SMILES → returns `StructureFilterBatchUploadResponse` with `job_id`, task runs on Celery.
- `POST /structure-filter/batch/upload` (file upload) → always async.

## Workflow

### 1. Filter with a preset (sync, ≤1000 SMILES)

```bash
curl -sS -X POST http://localhost:8000/api/v1/structure-filter/filter \
  -H 'Content-Type: application/json' \
  -d '{
    "smiles_list": ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", "[Ni]"],
    "preset": "drug_like"
  }'
```

Response (synchronous):

```json
{
  "input_count": 4,
  "output_count": 1,
  "stages": [
    {"stage_name": "parse", "stage_index": 1, "input_count": 4, "passed_count": 4, "rejected_count": 0, "enabled": true},
    {"stage_name": "valence", "stage_index": 2, "input_count": 4, "passed_count": 4, "rejected_count": 0, "enabled": true},
    {"stage_name": "alerts", "stage_index": 3, "input_count": 4, "passed_count": 3, "rejected_count": 1, "enabled": true},
    ...
  ],
  "molecules": [
    {"smiles": "CCO", "status": "rejected", "failed_at": "property", "rejection_reason": "MW 46.07 below minimum 200"},
    {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "status": "passed", "failed_at": null, "rejection_reason": null},
    ...
  ]
}
```

### 2. Filter with explicit config (overrides preset)

```bash
curl -sS -X POST http://localhost:8000/api/v1/structure-filter/filter \
  -H 'Content-Type: application/json' \
  -d '{
    "smiles_list": ["CCO", "c1ccccc1C(=O)O"],
    "config": {
      "min_mw": 100, "max_mw": 400,
      "min_logp": -1, "max_logp": 4,
      "max_tpsa": 120, "max_rot_bonds": 8,
      "max_rings": null, "max_sa_score": 4.5,
      "use_pains": true, "use_brenk": true, "use_kazius": false, "use_nibr": false,
      "enable_novelty": false, "novelty_threshold": 0.85,
      "weight_validity": 0.25, "weight_qed": 0.35,
      "weight_alert_free": 0.2, "weight_sa": 0.2
    }
  }'
```

If both `preset` and `config` are given, `preset` wins.

### 3. Composite score per molecule (0-1)

```bash
curl -sS -X POST http://localhost:8000/api/v1/structure-filter/score \
  -H 'Content-Type: application/json' \
  -d '{
    "smiles_list": ["CCO", "c1ccccc1C(=O)O", "invalid_smiles"],
    "preset": "drug_like"
  }'
```

Response: `{"scores": [0.412, 0.871, null]}`. `null` means unparseable — preserve that in downstream ML pipelines rather than coercing to 0.

### 4. REINVENT 4 Component API

The endpoint accepts the REINVENT 4 input contract verbatim (raw list, not wrapped in a model):

```bash
curl -sS -X POST "http://localhost:8000/api/v1/structure-filter/reinvent-score?preset=drug_like" \
  -H 'Content-Type: application/json' \
  -d '[
    {"input_string": "CCO", "query_id": "q1"},
    {"input_string": "c1ccccc1C(=O)O", "query_id": "q2"},
    {"input_string": "invalid", "query_id": "q3"}
  ]'
```

Response:

```json
{
  "output": {
    "successes_list": [
      {"query_id": "q1", "output_value": 0.412},
      {"query_id": "q2", "output_value": 0.871}
    ]
  }
}
```

**Critical invariant**: invalid SMILES are **omitted**, not scored 0.0 (REINVENT expects this — see Pitfall 4 / D-14). Don't insert zeros yourself.

REINVENT 4 `config.toml` snippet:

```toml
[[parameters.component]]
type = "ExternalProcess"
command = "curl -X POST http://chemaudit:8000/api/v1/structure-filter/reinvent-score?preset=drug_like -H 'Content-Type: application/json' -d @-"
```

### 5. Async batch (file upload, always async)

```bash
curl -sS -X POST http://localhost:8000/api/v1/structure-filter/batch/upload \
  -F "file=@generated_library.csv" \
  -F "preset=lead_like"
```

Or with explicit JSON config:

```bash
curl -sS -X POST http://localhost:8000/api/v1/structure-filter/batch/upload \
  -F "file=@library.sdf" \
  -F 'config={"min_mw":200,"max_mw":450,...}'
```

### 6. Track async progress

WebSocket:

```javascript
const ws = new WebSocket(`ws://localhost:8000/ws/structure-filter/${job_id}`);
```

Poll:

```bash
curl -sS http://localhost:8000/api/v1/structure-filter/batch/<job_id>/status
```

Response: `{job_id, status, progress, current_stage}`.

### 7. Async results

```bash
curl -sS http://localhost:8000/api/v1/structure-filter/batch/<job_id>/results
```

Returns the same `FilterResponse` shape wrapped in `{job_id, status, result}`.

### 8. Download

```bash
# Plain text — one passed SMILES per line
curl -sS http://localhost:8000/api/v1/structure-filter/batch/<job_id>/download/passed_txt \
  -o passed.txt

# Full CSV — every molecule with status, failed_at, rejection_reason
curl -sS http://localhost:8000/api/v1/structure-filter/batch/<job_id>/download/full_csv \
  -o results.csv
```

## Examples

### Example 1 — "Filter 10K generated SMILES and save the passing ones"

1. `POST /structure-filter/batch/upload` with the file and `preset=drug_like` → get `job_id`.
2. Track via WS or poll `/status`.
3. `GET /structure-filter/batch/<job_id>/download/passed_txt` → plain-text file of passed SMILES.

### Example 2 — "Integrate with REINVENT"

Wire the `/reinvent-score` endpoint as an `ExternalProcess` component. Choose the preset via query string. REINVENT's `query_id` round-trip is preserved; agents skip invalid SMILES naturally.

### Example 3 — "Fragment-like filter on a small library"

1. `POST /structure-filter/filter` with `smiles_list` (≤1000) and `preset="fragment_like"`.
2. Read `molecules[]`, keep those with `status="passed"`.
3. Render the funnel: `stages[]` has `input_count` / `passed_count` / `rejected_count` per stage.

### Example 4 — "Turn on novelty filtering vs ChEMBL"

Set `enable_novelty: true` and `novelty_threshold: 0.85` in an explicit config (no preset enables novelty by default):

```bash
curl -sS -X POST http://localhost:8000/api/v1/structure-filter/filter \
  -H 'Content-Type: application/json' \
  -d '{"smiles_list": [...], "config": {"min_mw": 200, "max_mw": 500, ..., "enable_novelty": true, "novelty_threshold": 0.85}}'
```

Molecules with Tanimoto > 0.85 to any ChEMBL compound are rejected at the novelty stage.

## Rate limits

- `/structure-filter/filter`: 20/min.
- `/structure-filter/score`, `/structure-filter/reinvent-score`: 30/min.
- `/structure-filter/batch/upload`: 3/min.
- `/structure-filter/batch/<job_id>/status`: 60/min.
- `/structure-filter/batch/<job_id>/results`: 30/min.
- `/structure-filter/batch/<job_id>/download/<format>`: 10/min.

## Troubleshooting

### 400 "Unknown preset '<name>'"

Valid presets: `drug_like`, `lead_like`, `fragment_like`, `permissive`. Case-sensitive.

### 400 "No valid SMILES found in input"

File parsed but every row errored. For CSV, ensure the column is named exactly `SMILES`. For SDF, ensure molecules have valid MOL blocks.

### Score 0.0 vs null — which is it?

`null` = parse failure or SMILES rejected before scoring. `0.0` = parsed and scored but every component was 0. Critical distinction for generative models; never collapse them.

### REINVENT returns fewer items than I sent

Expected — invalid SMILES are omitted from `successes_list`, not zero-scored. Match by `query_id`, not list position.

### Async job_id found in `/status` but no results

Celery worker crashed. Check `docker compose logs worker` and re-upload. Redis TTL is 1 hour on metadata.

### Funnel shows 0 passed but all stages have 100% input

All molecules pass individual stages but none make it through `dedup`. Means the input contains only duplicates (identical InChIKeys). Deduplicate before filtering, or disable dedup by adjusting the pipeline.

### High rejection count at `alerts` stage for a fragment-screen

`fragment_like` preset has only PAINS enabled. If you want Brenk filtering for fragments, use an explicit config rather than the preset.

## Further reading

- `references/preset-configs.md` — exhaustive threshold tables and weight vectors.
- `chemaudit-qsar-ready` — different goal: ML-ready curation with provenance, not pass/fail filtering.
- `chemaudit-molecule-validation` — single-molecule scoring when you need per-check detail.
