<div align="center">

<img src="assets/logo.png" alt="ChemAudit" width="80" />

# User Guide

### Complete Guide to ChemAudit Features

</div>

---

## Table of Contents

- [Single Molecule Validation](#single-molecule-validation)
  - [Validation Checks](#validation-checks-explained)
- [Batch Processing](#batch-processing)
- [Structural Alerts](#structural-alerts)
- [Scoring](#scoring)
  - [ML-Readiness](#ml-readiness-scoring)
  - [Drug-Likeness](#drug-likeness)
  - [Safety Filters](#safety-filters)
  - [ADMET Predictions](#admet-predictions)
  - [NP-Likeness](#np-likeness)
  - [Scaffold Analysis](#scaffold-analysis)
  - [Aggregator Likelihood](#aggregator-likelihood)
- [Standardization](#standardization)
- [Database Integrations](#database-integrations)
- [Exporting Results](#exporting-results)
- [Scoring Methodology Reference](SCORING_METHODOLOGY.md)

---

## Single Molecule Validation

### Supported Input Formats

| Format | Example | Auto-Detected |
|--------|---------|---------------|
| **SMILES** | `CCO` | Yes |
| **InChI** | `InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3` | Yes |
| **MOL Block** | V2000/V3000 format | Yes |

### How to Validate

**Web Interface:**
1. Navigate to the **Single Validation** page (home)
2. Enter or paste your molecule
3. Click **Validate**
4. Review results across all tabs: Validation, Alerts, Scoring, Standardization, Database Lookup

**API:**
```bash
curl -X POST http://localhost:8001/api/v1/validate \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "CCO",
    "format": "auto"
  }'
```

### Validation Checks Explained

ChemAudit runs 5 basic checks on every molecule and 17 deep validation checks organized into three domains.

<details>
<summary><b>Basic Checks (always run)</b></summary>

| Check | Severity | Description |
|-------|----------|-------------|
| **Parsability** | Critical | Can the input be parsed into a valid molecule? |
| **Sanitization** | Error | Does the molecule pass RDKit sanitization? |
| **Valence** | Critical | Do all atoms have chemically valid bond counts? |
| **Aromaticity** | Error | Can aromatic systems be kekulized? |
| **Connectivity** | Warning | Is the molecule a single connected component? |

</details>

<details>
<summary><b>Deep Validation — Chemical Composition (6 checks)</b></summary>

| Check | Severity | Description |
|-------|----------|-------------|
| **Mixture Detection** | Warning | Identifies disconnected fragments (drug, salt, solvent, unknown) |
| **Solvent Contamination** | Warning | Matches against 15+ known solvents (water, DMSO, DMF, etc.) |
| **Inorganic Filter** | Warning/Error | Detects inorganic or organometallic compounds |
| **Radical Detection** | Warning | Flags atoms with unpaired electrons |
| **Isotope Labels** | Info | Detects isotope-labeled atoms (deuterium, ¹³C, tritium, etc.) |
| **Trivial Molecule** | Error | Flags molecules with ≤ 3 heavy atoms |

</details>

<details>
<summary><b>Deep Validation — Structural Complexity (6 checks)</b></summary>

| Check | Severity | Description |
|-------|----------|-------------|
| **Hypervalent Atoms** | Warning | Atoms exceeding maximum allowed valence |
| **Polymer Detection** | Info | SGroup markers, MW > 1500 Da, or dummy atoms |
| **Ring Strain** | Warning | 3- or 4-membered rings with significant angle strain |
| **Macrocycle Detection** | Info | Rings with > 12 atoms |
| **Charged Species** | Info | Net charge, charged atoms, zwitterion detection |
| **Explicit Hydrogen Audit** | Info | Mixed explicit/implicit hydrogen representation |

</details>

<details>
<summary><b>Deep Validation — Stereo & Tautomers (5 checks)</b></summary>

| Check | Severity | Description |
|-------|----------|-------------|
| **Stereoisomer Enumeration** | Warning | Enumerates possible stereoisomers from undefined centers (cap: 128) |
| **Undefined Stereocenters** | Warning | Counts undefined chiral centers |
| **Tautomer Detection** | Info | Enumerates tautomers and checks if input is canonical |
| **Aromatic System Validation** | Warning | Unusual aromatic ring sizes or charged aromatic atoms |
| **Coordinate Dimension** | Info | Detects 2D, 3D, or no coordinate data |

</details>

All deep validation check severities can be customized through the severity configuration panel. For detailed logic and thresholds, see the [Scoring Methodology](SCORING_METHODOLOGY.md#validation-checks).

### Options

| Option | Description |
|--------|-------------|
| **preserve_aromatic** | Output SMILES with aromatic notation (lowercase atoms like `c1ccccc1`) instead of kekulized form (`C1=CC=CC=C1`) |

---

## Batch Processing

Process large datasets with real-time progress tracking.

### Specifications

| Feature | Details |
|---------|---------|
| **Max Molecules** | Configurable per deployment profile (1K - 1M) |
| **Max File Size** | Configurable per deployment profile (100MB - 1GB) |
| **Supported Formats** | SDF, CSV, TSV, TXT |
| **Progress Updates** | Real-time via WebSocket |
| **Worker Queues** | Separate default and priority queues |

### CSV Format Requirements

Your CSV must have a column containing SMILES strings:

```csv
Name,SMILES,Activity
Aspirin,CC(=O)Oc1ccccc1C(=O)O,Active
Caffeine,Cn1cnc2c1c(=O)n(c(=O)n2C)C,Active
Ethanol,CCO,Inactive
```

TSV and TXT files with tab-separated columns are also supported.

### How to Process Batch Files

**Web Interface:**
1. Navigate to **Batch Processing** page
2. Drag & drop your file or click to browse
3. For CSV/TSV files, select the SMILES column and optional Name column
4. Toggle extended safety filters, ChEMBL alerts, or standardization if desired
5. Click **Upload and Process**
6. Monitor real-time progress via WebSocket
7. View results with sorting, filtering, and pagination
8. Export results in your preferred format

**API:**
```bash
# Upload file
curl -X POST http://localhost:8001/api/v1/batch/upload \
  -F "file=@molecules.sdf"

# Response: {"job_id": "abc123", "status": "pending", "total_molecules": 1000}

# Check progress
curl http://localhost:8001/api/v1/batch/abc123/status

# Get results (paginated)
curl "http://localhost:8001/api/v1/batch/abc123?page=1&page_size=50"

# Get statistics only
curl http://localhost:8001/api/v1/batch/abc123/stats
```

### Filtering and Sorting Results

| Filter | Description |
|--------|-------------|
| **Status** | `success` or `error` |
| **Min/Max Score** | Validation score range (0-100) |
| **Sort By** | index, name, smiles, score, qed, safety, status, issues |
| **Sort Direction** | ascending or descending |

---

## Structural Alerts

Screen molecules against known problematic substructures.

### Available Alert Catalogs

| Catalog | Description | Patterns |
|---------|-------------|----------|
| **PAINS** | Pan-Assay Interference Compounds (A/B/C) | ~480 |
| **BRENK** | Unwanted chemical moieties | ~105 |
| **NIH** | NIH MLSMR excluded structures | ~180 |
| **ZINC** | ZINC database filters | ~95 |
| **ChEMBL** | Pharma company filters (BMS, Dundee, Glaxo, Inpharmatica, Lint, MLSMR, SureChEMBL) | ~700+ |

**Note:** 87 FDA-approved drugs contain PAINS patterns. Alerts are warnings for investigation, not automatic rejections.

### How to Screen

**Web Interface:**
1. Enter your molecule
2. Navigate to **Structural Alerts** tab
3. View matched alerts with severity and matched atoms

**API:**
```bash
# Full screening
curl -X POST http://localhost:8001/api/v1/alerts \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "c1ccc2c(c1)nc(n2)Sc3nnnn3C",
    "catalogs": ["PAINS", "BRENK"]
  }'

# Quick check (faster, yes/no only)
curl -X POST http://localhost:8001/api/v1/alerts/quick-check \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "c1ccc2c(c1)nc(n2)Sc3nnnn3C",
    "catalogs": ["PAINS"]
  }'
```

---

## Scoring

ChemAudit provides comprehensive molecular scoring across 7 dimensions. For complete calculation formulas, thresholds, and academic references, see the [Scoring Methodology](SCORING_METHODOLOGY.md).

### ML-Readiness Scoring

Evaluate how suitable a molecule is for machine learning applications. The score (0–100) is computed across four dimensions:

| Dimension | Max Points | What It Measures |
|-----------|-----------|------------------|
| **Structural Quality** | 20 | Single component, organic elements, no radicals, reasonable charge, no dummy atoms |
| **Property Profile** | 35 | Desirability-scored MW, LogP, TPSA, HBD, HBA, rotatable bonds, aromatic rings |
| **Complexity & Feasibility** | 25 | QED, SA Score, Fsp3, stereocenter complexity |
| **Representation Quality** | 20 | 451 descriptor completeness, 7 fingerprint types, bit density, conformer generation |

| Score | Tier | Interpretation |
|-------|------|----------------|
| **85–100** | Excellent | Suitable for most ML workflows without modification |
| **70–84** | Good | Minor limitations; generally suitable with standard preprocessing |
| **50–69** | Moderate | Usable but may need careful feature selection or preprocessing |
| **30–49** | Limited | Significant challenges; consider alternatives or specialized models |
| **0–29** | Poor | Not recommended for standard ML pipelines |

Each dimension card is expandable in the UI to show per-item scoring breakdowns with tooltips explaining the scoring logic. See [ML-Readiness Methodology](SCORING_METHODOLOGY.md#ml-readiness-scoring) for all formulas.

### Drug-Likeness

Evaluate compliance with 7 established drug-likeness rules plus a consensus score.

| Filter | Key Criteria | Reference |
|--------|-------------|-----------|
| **QED** | Composite 0–1 score (MW, LogP, HBA, HBD, PSA, RotB, aromatics, alerts) | Bickerton et al. (2012) |
| **Lipinski (Ro5)** | MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10 (≤ 1 violation OK) | Lipinski et al. (2001) |
| **Veber** | Rotatable bonds ≤ 10, TPSA ≤ 140 A² | Veber et al. (2002) |
| **Rule of Three** | MW < 300, LogP ≤ 3, HBD ≤ 3, HBA ≤ 3, RotB ≤ 3, TPSA ≤ 60 | Congreve et al. (2003) |
| **Ghose** | MW 160–480, LogP −0.4–5.6, atoms 20–70, MR 40–130 | Ghose et al. (1999) |
| **Egan** | LogP ≤ 5.88, TPSA ≤ 131.6 A² | Egan et al. (2000) |
| **Muegge** | 9 parameters (MW, LogP, TPSA, rings, C, heteroatoms, RotB, HBD, HBA) | Muegge et al. (2001) |
| **Consensus** | 0–5 score: count of Lipinski + Veber + Egan + Ghose + Muegge passes | — |
| **Lead-Likeness** | MW 200–350, LogP −1 to 3, RotB ≤ 7 | — |

See [Drug-Likeness Methodology](SCORING_METHODOLOGY.md#drug-likeness) for exact thresholds per filter.

### Safety Filters

Screen against structural alert databases used in drug discovery.

| Catalog | Source | Patterns | Purpose |
|---------|--------|----------|---------|
| **PAINS** | Baell & Holloway (2010) | ~480 | Pan-Assay Interference Compounds |
| **Brenk** | Brenk et al. (2008) | ~105 | Unfavorable chemical moieties |
| **NIH** | NIH MLSMR Program | ~180 | Screening exclusion filters |
| **ZINC** | ZINC Database | ~95 | Drug-likeness and reactivity |
| **ChEMBL** | Pharma companies | ~700+ | BMS, Dundee, Glaxo, Inpharmatica, LINT, MLSMR, SureChEMBL |

**Note:** 87 FDA-approved drugs contain PAINS patterns. Alerts are warnings for investigation, not automatic rejections.

### ADMET Predictions

Calculated molecular properties predictive of Absorption, Distribution, Metabolism, Excretion, and Toxicity.

| Property | Method | Output | Key Thresholds |
|----------|--------|--------|----------------|
| **Synthetic Accessibility** | SA Score (Ertl 2009) | 1–10 scale | < 4 easy, 4–6 moderate, > 6 difficult |
| **Aqueous Solubility** | ESOL (Delaney 2004) | LogS + mg/mL | ≥ −1 highly soluble → < −5 insoluble |
| **Complexity** | Fsp3, stereocenters, Bertz CT | Classification | Fsp3 > 0.42 = good 3D character |
| **CNS MPO** | Wager et al. (2010) | 0–6 score | ≥ 5 excellent CNS penetration |
| **Bioavailability** | Lipinski + Veber combined | Oral + CNS flags | Pass both = oral absorption likely |
| **Pfizer 3/75 Rule** | LogP + TPSA | Risk flag | LogP > 3 AND TPSA < 75 = toxicity risk |
| **GSK 4/400 Rule** | MW + LogP | Favorable flag | MW ≤ 400 AND LogP ≤ 4 = favorable |
| **Golden Triangle** | MW + LogP | In/out triangle | MW 200–450 AND LogP −0.5–5 = balanced |

See [ADMET Methodology](SCORING_METHODOLOGY.md#admet-predictions) for calculation formulas including the ESOL equation and CNS MPO component scoring.

### NP-Likeness

Natural product likeness scoring based on fragment analysis.

| Score | Category | Description |
|-------|----------|-------------|
| **≥ 2.0** | Strong NP-like | Natural product features clearly evident |
| **1.0–2.0** | NP-like | Suggestive of natural origin |
| **0.3–1.0** | Moderate NP-like | Some NP features present |
| **−0.3 to 0.3** | Mixed | Both NP and synthetic characteristics |
| **−1.0 to −0.3** | Moderate synthetic | More synthetic than NP features |
| **< −1.0** | Synthetic-like | Typical synthetic compound profile |

Color scale in the UI: green (NP-like) → slate (mixed) → red (synthetic).

### Scaffold Analysis

Murcko scaffold extraction for structure-activity analysis.

Returns:
- **Murcko scaffold SMILES** — Ring systems with linkers
- **Generic scaffold SMILES** — Simplified ring framework (all atoms replaced with carbon, all bonds with single)
- **Has scaffold** — Whether the molecule contains ring systems

### Aggregator Likelihood

Predicts whether a molecule is likely to form colloidal aggregates in biological assays, causing false-positive readouts.

| Risk Score | Classification |
|------------|---------------|
| **≥ 0.6** | High aggregation risk |
| **0.3–0.59** | Moderate risk |
| **< 0.3** | Low risk |

Based on 6 indicators: LogP, TPSA, MW, aromatic ring count, conjugation fraction, and 10 known aggregator SMARTS patterns (rhodanines, quinones, catechols, curcumin-like, etc.). See [Aggregator Methodology](SCORING_METHODOLOGY.md#aggregator-likelihood) for full details.

**How to Get Scores (API):**
```bash
curl -X POST http://localhost:8001/api/v1/score \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "CCO",
    "include": ["ml_readiness", "druglikeness", "safety_filters", "admet", "np_likeness", "scaffold", "aggregator"]
  }'
```

---

## Standardization

Standardize molecules using a ChEMBL-compatible pipeline.

### Pipeline Steps

| Step | Description | Always Runs |
|------|-------------|-------------|
| **Checker** | Detect structural issues before standardization | Yes |
| **Standardizer** | Fix common issues (nitro groups, metals, sulphoxides) | Yes |
| **Get Parent** | Extract parent molecule, remove salts and solvents | Yes |
| **Tautomer** | Canonicalize tautomers | No (opt-in) |

### Options

| Option | Default | Description |
|--------|---------|-------------|
| **include_tautomer** | `false` | Enable tautomer canonicalization (may lose E/Z stereo) |
| **preserve_stereo** | `true` | Attempt to preserve stereochemistry |

### How to Standardize

**API:**
```bash
curl -X POST http://localhost:8001/api/v1/standardize \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",
    "options": {
      "include_tautomer": false,
      "preserve_stereo": true
    }
  }'
```

The response includes:
- Original and standardized SMILES
- Steps applied with changes
- Checker issues found
- Excluded fragments (salts, solvents)
- Stereo comparison (if stereochemistry changed)
- Structure comparison (atom count, formula, mass change)

---

## Database Integrations

Look up molecules in major chemical databases.

### Supported Databases

| Database | Data Available | Rate Limit |
|----------|----------------|------------|
| **PubChem** | Properties, synonyms, IUPAC name | 30 req/min |
| **ChEMBL** | Bioactivity, targets, clinical phase | 30 req/min |
| **COCONUT** | Natural product data, organism source | 30 req/min |

### How to Search

**Web Interface:**
1. Enter your molecule on the Single Validation page
2. Navigate to **Database Lookup** tab
3. View cross-references from PubChem, ChEMBL, and COCONUT

**API:**
```bash
# Search PubChem
curl -X POST http://localhost:8001/api/v1/integrations/pubchem/lookup \
  -H "Content-Type: application/json" \
  -d '{"molecule": "CCO", "format": "smiles"}'

# Search ChEMBL
curl -X POST http://localhost:8001/api/v1/integrations/chembl/bioactivity \
  -H "Content-Type: application/json" \
  -d '{"molecule": "CCO", "format": "smiles"}'

# Search COCONUT
curl -X POST http://localhost:8001/api/v1/integrations/coconut/lookup \
  -H "Content-Type: application/json" \
  -d '{"molecule": "CCO", "format": "smiles"}'
```

---

## Exporting Results

### Available Export Formats

| Format | Extension | Use Case |
|--------|-----------|----------|
| **CSV** | `.csv` | Spreadsheet analysis, data pipelines |
| **Excel** | `.xlsx` | Formatted report with conditional coloring and summary sheet |
| **SDF** | `.sdf` | Structure-data exchange with other chemistry tools |
| **JSON** | `.json` | Programmatic processing, full data fidelity |
| **PDF** | `.pdf` | Professional reports with charts and molecule images |

### How to Export

**Web Interface:**
1. Complete a batch processing job
2. Click the **Export** button
3. Select format and optional filters (score range, status, specific molecules)
4. Download the file

**API:**
```bash
# Export all results as CSV
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=csv" -o results.csv

# Export filtered results as Excel
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=excel&score_min=80" -o results.xlsx

# Export specific molecules as SDF
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=sdf&indices=0,1,5,23" -o results.sdf

# Export as PDF report
curl "http://localhost:8001/api/v1/batch/{job_id}/export?format=pdf" -o report.pdf

# Export large selections via POST
curl -X POST "http://localhost:8001/api/v1/batch/{job_id}/export?format=json" \
  -H "Content-Type: application/json" \
  -d '{"indices": [0, 1, 2, 3, 4, 5]}' -o results.json
```

### PDF Report Contents

- Executive summary with statistics
- Score distribution charts
- Molecule images with annotations
- Failed molecules with error details
- Alert summary by catalog
- Processing metadata

---

<div align="center">

**Need more help?** Check the [API Reference](API_REFERENCE.md) or [Troubleshooting Guide](TROUBLESHOOTING.md)

</div>
