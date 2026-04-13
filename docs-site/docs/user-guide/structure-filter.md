---
sidebar_position: 13
title: Structure Filter
description: Multi-stage funnel filtering for generative chemistry outputs and screening libraries
---

# Structure Filter

The Structure Filter provides multi-stage funnel filtering for chemical libraries, generative chemistry outputs, or any SMILES collection. Molecules pass through a sequence of configurable filter stages, with a visual funnel showing pass/fail counts at each stage.

## Filter Presets

Four built-in presets provide ready-to-use filter configurations:

| Preset | Description | Key Criteria |
|--------|-------------|-------------|
| **drug_like** | Lipinski-based drug-likeness | MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10 |
| **lead_like** | Lead optimization criteria | MW 200–450, LogP −1 to 4, RotB ≤ 7 |
| **fragment_like** | Rule of Three for fragments | MW ≤ 300, LogP ≤ 3, HBD ≤ 3, HBA ≤ 3 |
| **permissive** | Minimal filtering | Basic validity checks only |

You can also define custom filter configurations with arbitrary property ranges and SMARTS-based substructure inclusion/exclusion patterns.

## Using the Web Interface

1. Navigate to **Structure Filter** under the **Library** dropdown in the header
2. Choose your input method:
   - **Paste SMILES**: Enter SMILES strings, one per line
   - **Upload file**: Drag and drop a CSV or SDF file
3. Select a preset or configure custom filter stages
4. Add optional SMARTS patterns for substructure filtering
5. Click **Filter**
6. View the **funnel visualization** showing how many molecules pass or fail at each stage
7. Explore stage-by-stage results with per-molecule pass/fail details
8. Download passing molecules

:::info Async Processing
For datasets with more than 1,000 molecules, filtering runs asynchronously with a WebSocket progress feed. Smaller datasets return results immediately.
:::

## Funnel Visualization

The funnel chart shows the filtering pipeline as a series of stages. Each stage displays:

- **Stage name** and index
- **Input count** — molecules entering the stage
- **Passed count** — molecules passing the stage
- **Rejected count** — molecules filtered out
- Whether the stage is **enabled** or skipped

This makes it easy to identify which filter stage removes the most molecules and adjust your criteria accordingly.

## Scoring Mode

The Structure Filter also provides a continuous 0–1 score for each molecule, useful for integration with generative models:

```bash
curl -X POST http://localhost:8001/api/v1/structure-filter/score \
  -H "Content-Type: application/json" \
  -d '{
    "smiles_list": ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"],
    "preset": "drug_like"
  }'
```

Response:
```json
{
  "scores": [0.85, 0.92, 0.78]
}
```

A score of `null` indicates an invalid or unparseable SMILES.

### REINVENT Integration

A REINVENT-compatible scoring endpoint is available for direct integration with the REINVENT generative model:

```bash
curl -X POST "http://localhost:8001/api/v1/structure-filter/reinvent-score?preset=drug_like" \
  -H "Content-Type: application/json" \
  -d '[{"input_string": "CCO", "query_id": 0}, {"input_string": "c1ccccc1", "query_id": 1}]'
```

## API Reference

### Filter (Synchronous ≤ 1,000 molecules)

```bash
curl -X POST http://localhost:8001/api/v1/structure-filter/filter \
  -H "Content-Type: application/json" \
  -d '{
    "smiles_list": ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"],
    "preset": "drug_like"
  }'
```

Response:
```json
{
  "input_count": 3,
  "output_count": 2,
  "stages": [
    {
      "stage_name": "validity",
      "stage_index": 0,
      "input_count": 3,
      "passed_count": 3,
      "rejected_count": 0,
      "enabled": true
    },
    {
      "stage_name": "property_filter",
      "stage_index": 1,
      "input_count": 3,
      "passed_count": 2,
      "rejected_count": 1,
      "enabled": true
    }
  ],
  "molecules": [
    { "smiles": "CCO", "status": "passed", "failed_at": null, "rejection_reason": null },
    { "smiles": "c1ccccc1", "status": "passed", "failed_at": null, "rejection_reason": null },
    { "smiles": "CC(=O)Oc1ccccc1C(=O)O", "status": "rejected", "failed_at": "property_filter", "rejection_reason": "MW > 500" }
  ]
}
```

### Batch Upload (> 1,000 molecules)

```bash
curl -X POST http://localhost:8001/api/v1/structure-filter/batch/upload \
  -F "file=@molecules.csv" \
  -F "preset=drug_like"
```

### Check Status

```bash
curl http://localhost:8001/api/v1/structure-filter/batch/{job_id}/status
```

### Get Results

```bash
curl http://localhost:8001/api/v1/structure-filter/batch/{job_id}/results
```

### Download

```bash
# Passing molecules only (one SMILES per line)
curl http://localhost:8001/api/v1/structure-filter/batch/{job_id}/download/passed_txt -o passed.txt

# Full results with status
curl http://localhost:8001/api/v1/structure-filter/batch/{job_id}/download/full_csv -o results.csv
```

```python
import requests

# Synchronous filter
response = requests.post(
    "http://localhost:8001/api/v1/structure-filter/filter",
    json={
        "smiles_list": ["CCO", "c1ccccc1", "INVALID"],
        "preset": "drug_like"
    }
)
result = response.json()
print(f"Input: {result['input_count']}, Output: {result['output_count']}")
for mol in result["molecules"]:
    print(f"  {mol['smiles']}: {mol['status']}")
```

## Rate Limits

| Endpoint | Limit |
|----------|-------|
| `POST /structure-filter/filter` | 20 req/min |
| `POST /structure-filter/score` | 30 req/min |
| `POST /structure-filter/reinvent-score` | 30 req/min |
| `POST /structure-filter/batch/upload` | 3 req/min |
| `GET /structure-filter/batch/*/status` | 60 req/min |
| `GET /structure-filter/batch/*/results` | 30 req/min |
| `GET /structure-filter/batch/*/download/*` | 10 req/min |

## WebSocket

For batch jobs, connect to the WebSocket for real-time progress:

```
ws://localhost:8001/ws/structure-filter/{job_id}
```

Same message format and keep-alive protocol as the [batch processing WebSocket](/docs/api/websocket).

## Next Steps

- **[QSAR-Ready Pipeline](/docs/user-guide/qsar-ready)** — Curate molecules for ML before filtering
- **[Batch Processing](/docs/user-guide/batch-processing)** — Full validation and scoring
- **[Scoring Profiles](/docs/user-guide/scoring/profiles)** — Custom scoring criteria
