---
sidebar_position: 12
title: QSAR-Ready Pipeline
description: Prepare chemical datasets for machine learning with multi-step curation
---

# QSAR-Ready Pipeline

The QSAR-Ready Pipeline prepares chemical datasets for machine learning by applying a standardized multi-step curation pipeline. It ensures that molecules are clean, canonical, and deduplicated before use in QSAR/QSPR modeling.

## Pipeline Steps

The pipeline applies the following steps in order:

| Step | Description | Details |
|------|-------------|---------|
| **Standardization** | ChEMBL-compatible structure normalization | Fixes nitro groups, metal disconnection, charge normalization |
| **Salt Stripping** | Remove counterions and salts | Extracts the parent molecule using MolVS fragment patterns |
| **Neutralization** | Neutralize charged species | Converts charged forms to neutral where chemically appropriate |
| **Tautomer Canonicalization** | Canonicalize tautomeric forms | Uses RDKit's tautomer enumerator for a canonical representation |
| **Duplicate Removal** | Remove duplicates by InChIKey | Identifies and flags duplicate structures after standardization |

Each molecule receives a result status indicating its outcome.

## Result Status

| Status | Meaning |
|--------|---------|
| **ok** | Successfully curated — molecule passed all pipeline steps |
| **rejected** | Failed a pipeline step (see `rejection_reason` for details) |
| **duplicate** | Duplicate of another molecule by standardized InChIKey |
| **error** | Processing error (e.g., unparseable SMILES) |

## Using the Web Interface

1. Navigate to **QSAR-Ready** under the **Data Preparation** dropdown in the header
2. Choose your input method:
   - **Paste SMILES**: Enter SMILES strings, one per line
   - **Upload file**: Drag and drop a CSV or SDF file
3. Configure pipeline options if needed
4. Click **Process**
5. Monitor real-time progress
6. Review results:
   - Summary statistics (ok / rejected / duplicate / error counts)
   - Per-molecule details with step-by-step provenance
   - InChIKey change tracking (original vs. standardized)
7. Download the curated dataset in CSV, SDF, or JSON format

:::tip Batch Export
Use the download buttons to export only the successfully curated molecules (`ok` status) for direct use in ML pipelines.
:::

## API Reference

### Single Molecule

Process one molecule through the pipeline:

```bash
curl -X POST http://localhost:8001/api/v1/qsar-ready/single \
  -H "Content-Type: application/json" \
  -d '{
    "smiles": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",
    "config": {}
  }'
```

**Response:**

```json
{
  "original_smiles": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",
  "original_inchikey": "XAKUHBCMIFQRLG-UHFFFAOYSA-M",
  "curated_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "standardized_inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
  "inchikey_changed": true,
  "status": "ok",
  "rejection_reason": null,
  "steps": [
    {
      "step_name": "standardization",
      "applied": true,
      "changes": ["Applied standardization rules"]
    },
    {
      "step_name": "salt_stripping",
      "applied": true,
      "changes": ["Removed fragment: [Na+]"]
    }
  ]
}
```

### Batch Upload

```bash
# Upload a file
curl -X POST http://localhost:8001/api/v1/qsar-ready/batch/upload \
  -F "file=@molecules.csv" \
  -F 'config={}'

# Or paste SMILES text
curl -X POST http://localhost:8001/api/v1/qsar-ready/batch/upload \
  -F "smiles_text=CCO
c1ccccc1
CC(=O)Oc1ccccc1C(=O)O"
```

### Check Status

```bash
curl http://localhost:8001/api/v1/qsar-ready/batch/{job_id}/status
```

### Get Results

```bash
curl "http://localhost:8001/api/v1/qsar-ready/batch/{job_id}/results?page=1&per_page=50"
```

### Download Curated Dataset

```bash
# CSV format
curl http://localhost:8001/api/v1/qsar-ready/batch/{job_id}/download/csv -o curated.csv

# SDF format
curl http://localhost:8001/api/v1/qsar-ready/batch/{job_id}/download/sdf -o curated.sdf

# JSON format
curl http://localhost:8001/api/v1/qsar-ready/batch/{job_id}/download/json -o curated.json
```

```python
import requests

# Single molecule
response = requests.post(
    "http://localhost:8001/api/v1/qsar-ready/single",
    json={"smiles": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]", "config": {}}
)
result = response.json()
print(f"Status: {result['status']}, Curated: {result['curated_smiles']}")

# Batch upload
with open("molecules.csv", "rb") as f:
    response = requests.post(
        "http://localhost:8001/api/v1/qsar-ready/batch/upload",
        files={"file": f},
        data={"config": "{}"}
    )
job_id = response.json()["job_id"]
```

## Rate Limits

| Endpoint | Limit |
|----------|-------|
| `POST /qsar-ready/single` | 30 req/min |
| `POST /qsar-ready/batch/upload` | 3 req/min |
| `GET /qsar-ready/batch/*/status` | 60 req/min |
| `GET /qsar-ready/batch/*/results` | 30 req/min |
| `GET /qsar-ready/batch/*/download/*` | 10 req/min |

## Use Cases

### ML Dataset Preparation

1. Upload your raw compound collection
2. The pipeline standardizes, deduplicates, and cleans all structures
3. Download the curated set and use it directly for model training

### Compound Registration

1. Process incoming compounds through the pipeline
2. Use the standardized InChIKey for unique registration
3. Flag duplicates against your existing collection

## Next Steps

- **[Structure Filter](/docs/user-guide/structure-filter)** — Apply property and substructure filters to curated molecules
- **[Batch Processing](/docs/user-guide/batch-processing)** — Full validation and scoring for large datasets
- **[Standardization](/docs/user-guide/standardization)** — Learn about the underlying standardization pipeline
