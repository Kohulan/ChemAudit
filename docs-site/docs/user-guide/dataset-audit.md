---
sidebar_position: 14
title: Dataset Audit
description: Comprehensive dataset health auditing with contradiction detection and curation reports
---

# Dataset Audit

The Dataset Audit tool provides comprehensive health auditing for chemical datasets. Upload a CSV or SDF file and receive a detailed health score, contradictory label detection, property distribution analysis, and a curation report with actionable recommendations.

## Features

| Feature | Description |
|---------|-------------|
| **Health Score** | Overall 0–100 dataset quality score with weighted sub-scores |
| **Contradictory Labels** | Detects molecules with the same structure but different activity labels |
| **Property Distributions** | Statistical analysis of molecular properties across the dataset |
| **Curation Report** | Actionable recommendations for dataset improvement |
| **Dataset Diff** | Compare two dataset versions to find added, removed, and modified molecules |
| **Curated CSV** | Download a cleaned version of the dataset with audit columns appended |

## Using the Web Interface

1. Navigate to **Dataset Audit** under the **Data Preparation** dropdown
2. Upload a CSV or SDF file
   - Optionally specify the SMILES column (auto-detected if not provided)
   - Optionally specify the activity column for contradiction detection
3. Monitor audit progress through the stages:
   - **Health audit** — Structural validity, property distributions, standardization coverage
   - **Contradiction detection** — Same-structure/different-label analysis
   - **Curation** — Generating recommendations and curated output
4. Review results:
   - **Health Score** — Overall score with sub-score breakdown
   - **Issues** — Categorized problems found in the dataset
   - **Contradictions** — Molecules with conflicting labels
   - **Property Distributions** — Histograms and statistics
   - **Treemap** — Interactive drill-down visualization of issues by category
5. Download the curation report or curated CSV

## Health Score

The health score is a weighted combination of sub-scores:

| Sub-Score | What It Measures |
|-----------|-----------------|
| **Validity** | Percentage of molecules that parse and pass basic validation |
| **Diversity** | Scaffold diversity via Shannon entropy |
| **Property Coverage** | Whether molecular properties span reasonable drug-like ranges |
| **Standardization** | Consistency of molecular representations |

The overall score (0–100) provides a quick indicator of dataset quality for ML modeling.

:::warning Low Health Scores
A health score below 50 indicates significant data quality issues. Review the issues list and consider curating the dataset before using it for modeling.
:::

## Contradictory Label Detection

When an activity column is provided, the audit identifies molecules that have the same structure (by InChIKey) but different activity labels. These contradictions can introduce noise into ML models and should be reviewed.

Each contradiction includes:
- The shared structure (SMILES and InChIKey)
- All conflicting labels found
- The number of occurrences of each label

## Dataset Diff

Compare two versions of a dataset to understand what changed:

1. Complete the initial audit
2. Click **Compare** and upload a second file
3. View results categorized as:
   - **Added** — Molecules in the new version but not the original
   - **Removed** — Molecules in the original but not the new version
   - **Modified** — Same molecule with changed properties or labels
   - **Unchanged** — Identical across both versions

## API Reference

### Upload Dataset

```bash
curl -X POST http://localhost:8001/api/v1/dataset/upload \
  -F "file=@dataset.csv" \
  -F "smiles_column=SMILES" \
  -F "activity_column=Activity"
```

Response:
```json
{
  "job_id": "550e8400-...",
  "filename": "dataset.csv",
  "file_type": "csv",
  "status": "pending",
  "message": "Dataset audit submitted"
}
```

### Check Status

```bash
curl http://localhost:8001/api/v1/dataset/{job_id}/status
```

Response:
```json
{
  "job_id": "550e8400-...",
  "status": "processing",
  "progress": 60.0,
  "current_stage": "health_audit",
  "eta_seconds": 15.0
}
```

### Get Results

```bash
curl http://localhost:8001/api/v1/dataset/{job_id}/results
```

Returns 202 if the audit is still processing.

### Dataset Diff

```bash
curl -X POST http://localhost:8001/api/v1/dataset/{job_id}/diff \
  -F "file=@dataset_v2.csv"
```

Response:
```json
{
  "added": [...],
  "removed": [...],
  "modified": [...],
  "added_count": 50,
  "removed_count": 10,
  "modified_count": 25,
  "unchanged_count": 915
}
```

### Download Outputs

```bash
# Curation report (JSON)
curl http://localhost:8001/api/v1/dataset/{job_id}/download/report -o report.json

# Curated CSV (with appended audit columns)
curl http://localhost:8001/api/v1/dataset/{job_id}/download/csv -o curated.csv
```

```python
import requests
import time

# Upload dataset
with open("dataset.csv", "rb") as f:
    response = requests.post(
        "http://localhost:8001/api/v1/dataset/upload",
        files={"file": f},
        data={"smiles_column": "SMILES", "activity_column": "Activity"}
    )
job_id = response.json()["job_id"]

# Poll for completion
while True:
    status = requests.get(f"http://localhost:8001/api/v1/dataset/{job_id}/status").json()
    if status["status"] in ["complete", "error"]:
        break
    time.sleep(2)

# Get results
results = requests.get(f"http://localhost:8001/api/v1/dataset/{job_id}/results").json()
print(f"Health Score: {results['health_audit']['overall_score']}")
print(f"Contradictions: {len(results['contradictions'])}")
```

## Rate Limits

| Endpoint | Limit |
|----------|-------|
| `POST /dataset/upload` | 3 req/min |
| `GET /dataset/*/status` | 60 req/min |
| `GET /dataset/*/results` | 30 req/min |
| `POST /dataset/*/diff` | 10 req/min |
| `GET /dataset/*/download/*` | 10 req/min |

## Next Steps

- **[QSAR-Ready Pipeline](/docs/user-guide/qsar-ready)** — Curate datasets for ML after auditing
- **[Batch Processing](/docs/user-guide/batch-processing)** — Full validation and scoring
- **[Batch Analytics](/docs/user-guide/batch-analytics)** — Clustering, taxonomy, and chemical space analysis
