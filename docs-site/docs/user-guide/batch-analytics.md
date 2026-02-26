---
sidebar_position: 3
title: Batch Analytics
description: Interactive analytics, visualizations, and chemical space exploration for batch results
---

# Batch Analytics & Visualizations

After batch processing completes, ChemAudit automatically computes analytics and provides interactive visualizations for exploring your dataset. Additional on-demand analyses are available for deeper investigation.

## Automatic Analytics

These analyses run immediately after batch completion at no additional cost:

### Deduplication

Identifies duplicate molecules across four comparison levels:

| Level | Method | What It Catches |
|-------|--------|----------------|
| **Exact** | Canonical SMILES | Identical structures |
| **Tautomeric** | Canonical tautomer SMILES | Tautomeric forms (e.g., keto/enol) |
| **Stereo-insensitive** | InChI without stereo layers | Enantiomers and diastereomers |
| **Salt-form** | ChEMBL parent compound | Different salt forms of the same drug |

Each level reports the number of unique compounds and groups duplicate molecules together with a representative structure.

### Statistics

Comprehensive property statistics computed across all successful molecules:

- **Per-property stats**: mean, median, standard deviation, Q1/Q3, IQR, min/max for validation score, QED, SA score, ML-readiness, and Fsp3
- **Pairwise correlations**: Pearson correlation between all property pairs (requires ≥10 values)
- **Outlier detection**: IQR-fence method (below Q1 − 1.5×IQR or above Q3 + 1.5×IQR)
- **Composite quality score** (0–100): weighted combination of validity rate (40%), scaffold diversity via Shannon entropy (35%), and Lipinski pass rate (25%)

## On-Demand Analytics

These computationally expensive analyses are triggered manually and run in the background:

### Scaffold Analysis

Extracts Murcko scaffolds (ring systems + linkers) and generic scaffolds (all atoms replaced with carbon):

- **Unique scaffold count** and **Shannon entropy** as a diversity metric
- **Frequency distribution**: top 50 scaffolds plus an "Other" bucket
- Acyclic molecules share an empty scaffold
- Click any scaffold in the treemap to select all molecules containing it

```bash
# Trigger scaffold analysis
curl -X POST http://localhost:8001/api/v1/batch/{job_id}/analytics/scaffold
```

### Chemical Space Mapping

Projects molecules into 2D space using molecular fingerprints (Morgan ECFP4, 2048-bit):

| Method | Algorithm | Limit | Best For |
|--------|-----------|-------|----------|
| **PCA** | Randomized SVD (pure numpy) | No limit | Fast overview, variance analysis |
| **t-SNE** | openTSNE | ≤ 2,000 molecules | Cluster visualization |

Both methods require at least 3 molecules. PCA also reports variance explained per component.

```bash
# Trigger PCA chemical space mapping
curl -X POST http://localhost:8001/api/v1/batch/{job_id}/analytics/chemical_space \
  -H "Content-Type: application/json" \
  -d '{"method": "pca"}'

# Trigger t-SNE (for smaller datasets)
curl -X POST http://localhost:8001/api/v1/batch/{job_id}/analytics/chemical_space \
  -H "Content-Type: application/json" \
  -d '{"method": "tsne"}'
```

### Matched Molecular Pairs (MMP)

Detects structurally similar molecule pairs that differ by a single chemical transformation:

- **BRICS single-cut fragmentation** identifies core + R-group pairs
- **Core heuristic**: core must have more heavy atoms than the R-group
- **Activity cliffs** (optional): SALI index = |Δactivity| / (1 − Tanimoto) for highly similar pairs with large activity differences
- **Lipophilic Ligand Efficiency** (optional): LLE = pActivity − LogP

```bash
# Trigger MMP analysis with activity cliff detection
curl -X POST http://localhost:8001/api/v1/batch/{job_id}/analytics/mmp \
  -H "Content-Type: application/json" \
  -d '{"activity_column": "validation_score"}'
```

:::warning Batch Size Limit
MMP analysis is limited to 5,000 molecules. For larger datasets, select a subset first.
:::

### Similarity Search

Find molecules similar to a query structure using ECFP4 Tanimoto similarity:

```bash
# Search by SMILES
curl -X POST http://localhost:8001/api/v1/batch/{job_id}/analytics/similarity_search \
  -H "Content-Type: application/json" \
  -d '{"query_smiles": "c1ccccc1", "top_k": 10}'

# Search by molecule index in the batch
curl -X POST http://localhost:8001/api/v1/batch/{job_id}/analytics/similarity_search \
  -H "Content-Type: application/json" \
  -d '{"query_index": 42, "top_k": 10}'
```

### R-Group Decomposition

Decompose molecules around a user-supplied core SMARTS pattern:

```bash
curl -X POST http://localhost:8001/api/v1/batch/{job_id}/analytics/rgroup \
  -H "Content-Type: application/json" \
  -d '{"core_smarts": "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"}'
```

## Interactive Visualizations

The analytics panel below the batch results table provides two tabs of interactive charts with linked selection.

### Distributions Tab

**Score Histogram** — Color-coded bars showing the distribution of validation scores across four quality bands:

| Band | Score Range | Color |
|------|-------------|-------|
| Excellent | 80–100 | Amber |
| Good | 60–80 | Dark amber |
| Moderate | 40–60 | Orange |
| Poor | 0–40 | Red |

Click any bar to select all molecules in that score range.

**Property Scatter Plot** — Configurable X/Y scatter plot for exploring property relationships. Choose axes from: MW, LogP, TPSA, QED, Overall Score, SA Score, Fsp3. Points are colored by a selectable property. Click individual points to toggle selection.

**Alert Frequency Chart** — Horizontal bar chart showing structural alert counts by type, sorted by frequency. Height scales dynamically with the number of alert types.

**Validation Treemap** — Hierarchical treemap grouping validation issues by category and count. Provides a quick visual overview of the most common problems in your dataset.

### Chemical Space Tab

**Scaffold Treemap** — Scaffold frequency treemap showing the most common ring systems. Click any scaffold cell to select all molecules sharing that scaffold. Uses 16 distinct colors for scaffold groups.

**Chemical Space Scatter** — Canvas-rendered 2D scatter plot for PCA or t-SNE projections. Optimized for performance with large datasets (800+ points). Features:

- **Brush selection**: Click and drag to draw a rectangle selecting all enclosed points
- **Click selection**: Click individual points to toggle them
- **Color-by dropdown**: Color points by overall score, QED, SA score, or Fsp3
- **PCA/t-SNE toggle**: Switch between projection methods (t-SNE triggers a background computation)
- **PNG export**: Download the current view as an image

### Linked Brushing

Selections are synchronized across all charts. When you select molecules in one visualization, the selection is reflected everywhere:

- Click a score histogram bar → those molecules highlight in the scatter plot
- Brush-select in chemical space → those molecules highlight in score histogram
- Click a scaffold treemap cell → those molecules appear selected in all charts

A selection toolbar appears showing the count of selected molecules with options to **Compare** (1–2 molecules), open **Actions** (subset operations), or **Clear** the selection.

## Molecule Comparison

Select 1 or 2 molecules from the results table and click the **Compare** button to open the comparison panel.

**Side-by-side view includes:**

- 2D structure images for both molecules
- Property comparison table with color-coded highlighting (green for better, red for worse):
  - Overall Score, QED, SA Score, Fsp3, Lipinski Violations, Alert Count
- Overlaid 6-axis radar chart showing normalized property profiles
- ECFP4 Tanimoto similarity score between the pair (via `POST /validate/similarity`)

### ECFP Similarity API

```bash
curl -X POST http://localhost:8001/api/v1/validate/similarity \
  -H "Content-Type: application/json" \
  -d '{
    "smiles_a": "CC(=O)Oc1ccccc1C(=O)O",
    "smiles_b": "CC(=O)Nc1ccc(O)cc1"
  }'
```

## Batch Timeline

A 4-phase horizontal timeline at the top of the batch view tracks processing progress:

1. **Upload** — File received and parsed
2. **Validation** — Molecules being validated
3. **Analytics** — Post-processing analytics running
4. **Complete** — All processing finished

Each phase shows a status indicator with live color updates.

## API Reference

### Get Analytics Status and Results

```bash
curl http://localhost:8001/api/v1/batch/{job_id}/analytics
```

**Response:**

```json
{
  "job_id": "550e8400-...",
  "status": {
    "deduplication": {"status": "complete", "computed_at": "2026-02-26T10:35:00Z"},
    "statistics": {"status": "complete", "computed_at": "2026-02-26T10:35:01Z"},
    "scaffold": {"status": "pending"},
    "chemical_space": {"status": "pending"},
    "mmp": {"status": "pending"}
  },
  "deduplication": {
    "exact": {"total_unique": 950, "groups": [...]},
    "tautomeric": {"total_unique": 920, "groups": [...]},
    "stereo_insensitive": {"total_unique": 900, "groups": [...]},
    "salt_form": {"total_unique": 880, "groups": [...]}
  },
  "statistics": {
    "property_stats": [...],
    "correlations": [...],
    "outliers": [...],
    "quality_score": {"score": 78.5, "validity": 0.95, "diversity": 0.72, "druglikeness": 0.65}
  }
}
```

### Trigger On-Demand Analytics

```bash
curl -X POST http://localhost:8001/api/v1/batch/{job_id}/analytics/{type}
```

**Available types:** `scaffold`, `chemical_space`, `mmp`, `similarity_search`, `rgroup`

**Response:**

```json
{
  "job_id": "550e8400-...",
  "analysis_type": "scaffold",
  "status": "queued"
}
```

Poll `GET /batch/{job_id}/analytics` to check when the analysis completes (status changes from `computing` to `complete`).

:::tip Analytics Caching
Analytics results are cached in Redis for 24 hours. Subsequent requests return cached results instantly without recomputation.
:::

## Next Steps

- **[Subset Actions](/docs/user-guide/subset-actions)** — Work with selected molecules from analytics
- **[Exporting Results](/docs/user-guide/exporting-results)** — Export analytics data and visualizations
- **[Scoring Profiles](/docs/user-guide/scoring/profiles)** — Score molecules against custom criteria
- **[Structural Alerts](/docs/user-guide/structural-alerts)** — Understand the alert frequency chart
