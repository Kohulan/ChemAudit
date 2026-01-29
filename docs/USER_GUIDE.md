<div align="center">

<img src="assets/logo.png" alt="ChemVault" width="80" />

# User Guide

### Complete Guide to ChemVault Features

</div>

---

## 📋 Table of Contents

- [Single Molecule Validation](#-single-molecule-validation)
- [Batch Processing](#-batch-processing)
- [Structural Alerts](#-structural-alerts)
- [ML-Readiness Scoring](#-ml-readiness-scoring)
- [Standardization](#-standardization)
- [Database Integrations](#-database-integrations)
- [Exporting Results](#-exporting-results)

---

## 🔬 Single Molecule Validation

### Supported Input Formats

| Format | Example | Auto-Detected |
|--------|---------|---------------|
| **SMILES** | `CCO` | ✅ Yes |
| **InChI** | `InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3` | ✅ Yes |
| **MOL Block** | V2000/V3000 format | ✅ Yes |

### How to Validate

**Web Interface:**
1. Navigate to **Single Validation** page
2. Enter or paste your molecule
3. Click **Validate**
4. Review results in the dashboard

**API:**
```bash
curl -X POST http://localhost:8000/api/v1/validate \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "CCO",
    "format": "auto"
  }'
```

### Validation Checks Explained

<details>
<summary><b>🔴 Critical Checks</b></summary>

| Check | Description | Common Causes |
|-------|-------------|---------------|
| **Valence** | Atoms must have valid bond counts | Typos in SMILES, incorrect charges |
| **Kekulization** | Aromatic rings must be resolvable | Invalid aromatic systems |
| **Connectivity** | Molecule must be valid graph | Disconnected fragments |

</details>

<details>
<summary><b>🟡 Warning Checks</b></summary>

| Check | Description | Recommendation |
|-------|-------------|----------------|
| **Aromaticity** | Validates aromatic perception | Review ring systems |
| **Stereochemistry** | Checks stereo specifications | Verify intended stereochemistry |
| **Charges** | Unusual charge states | Verify if intentional |

</details>

<details>
<summary><b>🔵 Info Checks</b></summary>

| Check | Description |
|-------|-------------|
| **Molecular Weight** | Reports MW for reference |
| **Atom Count** | Total heavy atoms |
| **Ring Count** | Number of ring systems |

</details>

---

## 📦 Batch Processing

Process large datasets efficiently with real-time progress tracking.

### Specifications

| Feature | Limit |
|---------|-------|
| **Max File Size** | 1 GB |
| **Max Molecules** | 1,000,000 |
| **Supported Formats** | SDF, CSV |
| **Progress Updates** | Real-time via WebSocket |

### CSV Format Requirements

Your CSV must have a column containing SMILES strings:

```csv
Name,SMILES,Activity
Aspirin,CC(=O)Oc1ccccc1C(=O)O,Active
Caffeine,Cn1cnc2c1c(=O)n(c(=O)n2C)C,Active
Ethanol,CCO,Inactive
```

### SDF Format

Standard MDL SDF format with multiple molecules:

```
MoleculeName
  RDKit  2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
    ...
M  END
> <Activity>
Active

$$$$
```

### How to Process Batch Files

**Web Interface:**
1. Navigate to **Batch Processing** page
2. Drag & drop your file or click to browse
3. For CSV files, select the SMILES column
4. Click **Upload and Process**
5. Monitor real-time progress
6. Download results when complete

**API:**
```bash
# Upload file
curl -X POST http://localhost:8000/api/v1/batch/upload \
  -F "file=@molecules.sdf"

# Response includes job_id
# {"job_id": "abc123", "status": "pending", "total_molecules": 1000}

# Check progress
curl http://localhost:8000/api/v1/batch/abc123/status

# Get results
curl http://localhost:8000/api/v1/batch/abc123
```

### Filtering Results

Filter batch results by:

| Filter | Description |
|--------|-------------|
| **Status** | `success` or `error` |
| **Min Score** | Minimum validation score (0-100) |
| **Max Score** | Maximum validation score (0-100) |

---

## ⚠️ Structural Alerts

Screen molecules against known problematic substructures.

### Available Alert Catalogs

| Catalog | Description | Alerts |
|---------|-------------|--------|
| **PAINS** | Pan-Assay Interference Compounds | ~480 patterns |
| **BRENK** | Unwanted chemical groups | ~105 patterns |
| **NIH** | NIH MLSMR excluded structures | ~180 patterns |
| **Dundee** | University of Dundee filters | ~95 patterns |

### How to Screen

**Web Interface:**
1. Enter your molecule
2. Navigate to **Structural Alerts** tab
3. Select catalogs to screen against
4. View matched alerts with severity

**API:**
```bash
curl -X POST http://localhost:8000/api/v1/alerts \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "c1ccc2c(c1)nc(n2)Sc3nnnn3C",
    "catalogs": ["PAINS", "BRENK"]
  }'
```

### Understanding Alert Results

```json
{
  "has_alerts": true,
  "total_alerts": 2,
  "alerts": [
    {
      "catalog": "PAINS",
      "name": "thiazole_amine_A",
      "severity": "high",
      "description": "Potential assay interference",
      "matched_atoms": [4, 5, 6, 7]
    }
  ]
}
```

---

## 📊 ML-Readiness Scoring

Evaluate how suitable a molecule is for machine learning applications.

### Score Components

| Component | Weight | What It Measures |
|-----------|--------|------------------|
| **Descriptor Calculability** | 30% | Can standard descriptors be computed? |
| **Fingerprint Generation** | 25% | Can molecular fingerprints be generated? |
| **Structural Complexity** | 20% | Is complexity within trainable range? |
| **Data Availability** | 15% | Similar molecules in training sets? |
| **Standardization** | 10% | Can molecule be standardized consistently? |

### Score Interpretation

| Score | Rating | Meaning |
|-------|--------|---------|
| **80-100** | ⭐⭐⭐⭐⭐ Excellent | Ideal for ML applications |
| **60-79** | ⭐⭐⭐⭐ Good | Suitable with minor caveats |
| **40-59** | ⭐⭐⭐ Moderate | May need preprocessing |
| **20-39** | ⭐⭐ Poor | Significant challenges expected |
| **0-19** | ⭐ Very Poor | Not recommended for ML |

### How to Get Score

**API:**
```bash
curl -X POST http://localhost:8000/api/v1/score \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "CCO",
    "include": ["ml_readiness", "np_likeness", "scaffold"]
  }'
```

---

## 🧹 Standardization

Standardize molecules using a ChEMBL-compatible pipeline.

### Pipeline Steps

| Step | Description | Default |
|------|-------------|---------|
| **Normalize** | Fix common drawing errors | ✅ On |
| **Reionize** | Standardize ionization states | ✅ On |
| **Uncharge** | Neutralize charges where appropriate | ✅ On |
| **Remove Stereo** | Strip stereochemistry | ❌ Off |
| **Remove Salts** | Strip counterions | ✅ On |
| **Canonicalize** | Generate canonical SMILES | ✅ On |
| **Tautomerize** | Enumerate/canonicalize tautomers | ✅ On |

### How to Standardize

**API:**
```bash
curl -X POST http://localhost:8000/api/v1/standardize \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",
    "options": {
      "remove_salts": true,
      "neutralize": true
    }
  }'
```

**Response:**
```json
{
  "original_smiles": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",
  "standardized_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "changes": [
    "Removed salt: [Na+]",
    "Neutralized carboxylate"
  ]
}
```

---

## 🗄️ Database Integrations

Look up molecules in major chemical databases.

### Supported Databases

| Database | Data Available | Rate Limit |
|----------|----------------|------------|
| **PubChem** | Properties, synonyms, safety | 5 req/sec |
| **ChEMBL** | Bioactivity, targets, assays | 10 req/sec |
| **COCONUT** | Natural product data | 5 req/sec |

### How to Search

**Web Interface:**
1. Enter your molecule
2. Navigate to **Database Lookup** tab
3. Select databases to search
4. View cross-references and data

**API:**
```bash
# Search PubChem
curl -X POST http://localhost:8000/api/v1/integrations/pubchem/lookup \
  -H "Content-Type: application/json" \
  -d '{"molecule": "CCO", "format": "smiles"}'

# Search all databases
curl -X POST http://localhost:8000/api/v1/integrations/lookup-all \
  -H "Content-Type: application/json" \
  -d '{"molecule": "CCO", "format": "smiles"}'
```

---

## 📤 Exporting Results

### Available Export Formats

| Format | Use Case |
|--------|----------|
| **CSV** | Spreadsheet analysis |
| **JSON** | Programmatic processing |
| **PDF** | Reports and documentation |

### Batch Export

```bash
# Export as CSV
curl "http://localhost:8000/api/v1/batch/{job_id}/export?format=csv" \
  -o results.csv

# Export as PDF report
curl "http://localhost:8000/api/v1/batch/{job_id}/export?format=pdf" \
  -o report.pdf
```

### PDF Report Contents

- Executive summary with statistics
- Score distribution charts
- Failed molecules with error details
- Alert summary by catalog
- Processing metadata

---

<div align="center">

**Need more help?** Check the [API Reference](API_REFERENCE.md) or [Troubleshooting Guide](TROUBLESHOOTING.md)

</div>
