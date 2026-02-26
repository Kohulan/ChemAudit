---
sidebar_position: 7
title: Exporting Results
description: Export batch processing results in CSV, Excel, SDF, JSON, and PDF formats
---

# Exporting Results

Export your batch processing results in multiple formats for analysis, reporting, and data exchange. ChemAudit supports CSV, Excel, SDF, JSON, and PDF exports with flexible filtering options.

## Available Export Formats

| Format | Extension | Use Case | Features |
|--------|-----------|----------|----------|
| **CSV** | `.csv` | Spreadsheet analysis, data pipelines | Simple tabular data |
| **Excel** | `.xlsx` | Formatted reports, presentations | Conditional coloring, summary sheet |
| **SDF** | `.sdf` | Structure-data exchange | Compatible with chemistry tools |
| **JSON** | `.json` | Programmatic processing | Full data fidelity, nested structures |
| **PDF** | `.pdf` | Professional reports, archival | Charts, molecule images, statistics |
| **Fingerprint Matrix** | `.zip` | ML pipeline input | Morgan/MACCS/RDKit FPs in CSV + numpy |
| **Deduplicated Set** | `.zip` | Dataset curation | Summary + annotated CSVs with group IDs |
| **Scaffold-Grouped** | `.csv` | SAR analysis | Murcko scaffolds + group assignments |
| **Property Matrix** | `.zip` | Comprehensive analysis | Flat CSV + 4-sheet Excel |

## How to Export

### Web Interface

1. Complete a batch processing job
2. Click the **Export** button
3. Select export format
4. (Optional) Apply filters:
   - Minimum/maximum validation score
   - Success or error status
   - Specific molecule indices
5. Click **Download**

The file downloads automatically with a descriptive filename including job ID and timestamp.

### API - Export All Results

```bash
# Export all results as CSV
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=csv" -o results.csv

# Export as Excel
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=excel" -o results.xlsx

# Export as SDF
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=sdf" -o results.sdf

# Export as JSON
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=json" -o results.json

# Export as PDF report
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=pdf" -o report.pdf
```

### API - Export with Filters

```bash
# Export molecules with validation score >= 80
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=csv&score_min=80" -o high_quality.csv

# Export only successful validations
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=excel&status=success" -o successful.xlsx

# Export specific score range
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=sdf&score_min=70&score_max=90" -o medium_scores.sdf
```

### API - Export Selected Molecules

For small selections (up to 200 indices), use GET with comma-separated indices:

```bash
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=csv&indices=0,1,5,23,42" -o selected.csv
```

For large selections (> 200 indices), use POST with JSON body:

```bash
curl -X POST "http://localhost:8001/api/v1/batch/{job_id}/export?format=json" \
  -H "Content-Type: application/json" \
  -d '{"indices": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]}' \
  -o selected.json
```

## Format Details

### CSV Export

Simple tabular format with these columns:

- Index
- SMILES
- Name (if provided)
- Validation score
- Status (success/error)
- Number of issues
- ML-readiness score
- QED score
- Safety pass/fail
- Alert counts (PAINS, BRENK, etc.)
- Molecular formula
- Molecular weight

**Pros:**
- Universal compatibility
- Easy to import into spreadsheets or databases
- Small file size

**Cons:**
- No nested data structures
- No formatting or charts
- Limited molecule information

### Excel Export

Formatted spreadsheet with two sheets:

**Results Sheet:**
- All molecule data in table format
- Conditional formatting (red for errors, green for high scores)
- Frozen header row
- Auto-sized columns

**Summary Sheet:**
- Job statistics
- Score distribution
- Alert summary
- Processing metadata

**Pros:**
- Professional appearance
- Ready for presentations
- Includes summary statistics
- Conditional coloring highlights issues

**Cons:**
- Larger file size than CSV
- Requires Excel or compatible software

### SDF Export

Structure-Data File format with molecular structures and properties:

```
Molecule_0
  RDKit          2D

 13 13  0  0  0  0  0  0  0  0999 V2000
    ... coordinates ...
M  END
> <validation_score>
95

> <qed>
0.65

> <ml_readiness_score>
85

$$$$
```

**Pros:**
- Preserves 2D/3D structure coordinates
- Compatible with chemistry software (RDKit, MOE, Schrodinger)
- Can include custom properties

**Cons:**
- Larger file size
- Not human-readable for data analysis

:::tip Chemistry Tools
SDF is the best format for importing into molecular modeling, docking, or cheminformatics tools.
:::

### JSON Export

Complete data in JSON format with nested structures:

```json
[
  {
    "index": 0,
    "smiles": "CCO",
    "name": "Ethanol",
    "status": "success",
    "validation": {
      "overall_score": 95,
      "issues": [],
      "all_checks": [...]
    },
    "alerts": {
      "total_alerts": 0,
      "alerts": []
    },
    "scoring": {
      "ml_readiness": {...},
      "druglikeness": {...},
      "admet": {...}
    }
  }
]
```

**Pros:**
- Complete data preservation
- Nested structures maintained
- Easy to parse programmatically
- Standard format for APIs

**Cons:**
- Larger file size
- Not suitable for spreadsheet analysis

### PDF Report

Professional report document with:

1. **Executive Summary**
   - Job metadata (date, molecule count, processing time)
   - Overall statistics
   - Key findings

2. **Score Distribution Charts**
   - Validation score histogram
   - ML-readiness distribution
   - QED distribution

3. **Molecule Gallery**
   - Structure images (up to 100 molecules)
   - Key properties
   - Validation scores
   - Alert indicators

4. **Failed Molecules Section**
   - Error messages
   - Failure reasons
   - Recommendations

5. **Alert Summary**
   - Alert counts by catalog
   - Most common patterns
   - Risk assessment

6. **Processing Metadata**
   - Configuration used
   - Processing time
   - Worker statistics

**Pros:**
- Professional presentation
- Self-contained (no external dependencies)
- Suitable for archival and sharing
- Includes visualizations

**Cons:**
- Largest file size
- Not suitable for further processing
- Limited to summary information for large batches

:::warning PDF Size Limits
For batches larger than 1,000 molecules, PDF exports include only summary statistics and the first 100 molecules to keep file size manageable.
:::

### Fingerprint Matrix Export

Generates a ZIP containing 9 files — 3 fingerprint types in 3 formats each:

| Fingerprint | Bits | File Prefix |
|------------|------|-------------|
| Morgan ECFP4 (radius 2) | 2,048 | `morgan` |
| MACCS keys | 167 | `maccs` |
| RDKit path-based | 2,048 | `rdkit` |

Each fingerprint is exported as `.csv`, `.npy` (numpy array), and `.npz` (compressed numpy). Molecules that failed parsing are skipped. Ideal for feeding directly into ML pipelines.

```bash
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=fingerprint" -o fingerprints.zip
```

### Deduplicated Set Export

Generates a ZIP with two CSV files:

- **dedup_summary.csv** — One row per unique canonical SMILES group: group_id, representative SMILES, group size, duplicate indices, dedup level
- **dedup_annotated.csv** — All molecules with group_id, is_representative flag, dedup_level, validation score, and InChIKey

Ideal for dataset curation — quickly identify and remove duplicates while preserving provenance.

```bash
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=dedup" -o deduplicated.zip
```

### Scaffold-Grouped Export

Single CSV with Murcko scaffold assignments:

- `scaffold_smiles` — the Murcko scaffold SMILES (empty string for acyclic molecules)
- `scaffold_group` — integer group ID (0 for acyclic, 1+ for scaffold groups)
- Plus: index, SMILES, name, status, overall score, ML-readiness score, InChIKey

Ideal for SAR analysis and scaffold-based dataset organization.

```bash
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=scaffold" -o scaffolds.csv
```

### Property Matrix Export

Generates a ZIP with two files:

- **properties_flat.csv** — All properties as columns in a single flat table
- **properties.xlsx** — Multi-sheet Excel workbook:
  - **Descriptors**: MW, LogP, TPSA, HBD, HBA, Rotatable Bonds, Aromatic Rings, Fsp3, Heavy Atom Count, Ring Count
  - **Scores**: Validation, ML-readiness, QED, SA, NP-likeness, Lipinski violations
  - **Alerts**: Per-catalog alert counts (PAINS, BRENK, NIH, Glaxo) and alert names
  - **Identifiers**: SMILES, name, InChIKey, status

```bash
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=property_matrix" -o properties.zip
```

### PDF Section Selection

PDF exports support selecting specific sections via the `sections` parameter:

| Section Key | Content |
|-------------|---------|
| `validation_summary` | Overall statistics and pass rates |
| `score_distribution` | Score histogram charts |
| `alert_frequency` | Alert count bar charts |
| `chemical_space` | Chemical space scatter plot |
| `scaffold_treemap` | Scaffold frequency visualization |
| `statistics` | Property statistics tables |
| `correlation_matrix` | Property correlation heatmap |
| `mmp_pairs` | Matched molecular pair summary |

All sections are included by default. To select specific sections:

```bash
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=pdf&sections=validation_summary,score_distribution" -o summary.pdf
```

### Excel Structure Images

Add `include_images=true` to embed 2D structure PNG images in the Excel export:

```bash
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=excel&include_images=true" -o results_with_images.xlsx
```

## Filter Options

### Score Filters

| Parameter | Type | Description |
|-----------|------|-------------|
| `score_min` | integer (0-100) | Minimum validation score |
| `score_max` | integer (0-100) | Maximum validation score |

**Example:**
```bash
# Export only high-quality molecules (score >= 90)
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=csv&score_min=90"
```

### Status Filters

| Parameter | Values | Description |
|-----------|--------|-------------|
| `status` | `success`, `error`, `warning` | Processing outcome |

**Example:**
```bash
# Export only failed molecules for debugging
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=csv&status=error"
```

### Molecule Selection

| Parameter | Type | Description |
|-----------|------|-------------|
| `indices` | string or array | Specific molecule indices |

**GET Example (small selections):**
```bash
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=sdf&indices=0,1,2,3"
```

**POST Example (large selections):**
```bash
curl -X POST "http://localhost:8001/api/v1/batch/{job_id}/export?format=json" \
  -d '{"indices": [0,1,2,3,4,5,6,7,8,9]}'
```

## Combining Filters

Filters can be combined for precise exports:

```bash
# High-quality successful molecules with specific indices
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=excel&status=success&score_min=80&indices=0,5,10,15"
```

## Best Practices

1. **Choose the right format**:
   - CSV for data analysis
   - Excel for reports and presentations
   - SDF for chemistry software
   - JSON for programmatic processing
   - PDF for archival and sharing

2. **Use filters to reduce file size**:
   - Export only relevant molecules
   - Filter by score to focus on quality
   - Separate successes and errors

3. **Export multiple times for different purposes**:
   - CSV for your database
   - Excel for stakeholder reports
   - SDF for molecular modeling
   - PDF for documentation

4. **Preserve raw data**:
   - Keep original batch results for reproducibility
   - Export JSON for complete data preservation

## Next Steps

- **[Batch Processing](/docs/user-guide/batch-processing)** — Process large datasets
- **[Batch Analytics](/docs/user-guide/batch-analytics)** — Analytics that power advanced exports
- **[Subset Actions](/docs/user-guide/subset-actions)** — Export selected molecules only
- **[API Reference](/docs/api/endpoints)** — Full export API documentation
