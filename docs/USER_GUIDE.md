<div align="center">

<img src="assets/logo.png" alt="ChemAudit" width="80" />

# User Guide

### Complete Guide to ChemAudit Features

</div>

---

## Table of Contents

- [Single Molecule Validation](#single-molecule-validation)
  - [Validation Checks](#validation-checks-explained)
  - [Compound Profiler & Safety Tabs](#compound-profiler--safety-tabs)
- [Batch Processing](#batch-processing)
  - [Batch Analytics](#batch-analytics)
- [QSAR-Ready Pipeline](#qsar-ready-pipeline)
- [Structure Filter](#structure-filter)
- [Dataset Audit](#dataset-audit)
- [Diagnostics](#diagnostics)
- [Structural Alerts](#structural-alerts)
- [Scoring](#scoring)
  - [ML-Readiness](#ml-readiness-scoring)
  - [Drug-Likeness](#drug-likeness)
  - [Safety Filters](#safety-filters)
  - [ADMET Predictions](#admet-predictions)
  - [NP-Likeness](#np-likeness)
  - [Scaffold Analysis](#scaffold-analysis)
  - [Aggregator Likelihood](#aggregator-likelihood)
  - [Ligand Efficiency](#ligand-efficiency)
  - [Bioavailability Radar](#bioavailability-radar)
  - [Property Breakdown](#property-breakdown)
  - [Salt Inventory](#salt-inventory)
- [Standardization](#standardization)
- [Database Integrations](#database-integrations)
  - [Identifier Resolution](#identifier-resolution)
  - [Database Comparison](#database-comparison)
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
4. Review results across all tabs: Validation, Alerts, Scoring, Scoring Profiles, Profiler, Safety, Standardization, Database Lookup

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

### Compound Profiler & Safety Tabs

The **Profiler** and **Safety** panels are consolidated as tabs within the Single Validation page:

**Profiler tab** provides:
- **3D Shape** — Molecular shape descriptors and visualization
- **Ligand Efficiency** — LE and LLE metrics (using proxy pIC50 from QED when no activity data is available)
- **Property Breakdown** — Per-atom TPSA and LogP contributions with functional group summaries
- **Bioavailability Radar** — 6-axis radar chart (LIPO, SIZE, POLAR, INSOLU, INSATU, FLEX) showing whether the molecule falls within optimal oral bioavailability ranges

**Safety tab** provides:
- Screening against all alert catalogs (PAINS, BRENK, NIH, ZINC, ChEMBL)
- Matched atoms highlighting
- Alert severity and details

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

### Batch Analytics

After batch processing completes, additional analytics tabs become available:

| Analytics | Description | Limit |
|-----------|-------------|-------|
| **Butina Clustering** | Sphere-exclusion clustering using Morgan fingerprints (radius=2, 2048 bits) with configurable Tanimoto distance cutoff (0.2–0.6, default 0.35) | Max 1,000 molecules |
| **Chemical Taxonomy** | Classifies molecules using ~50 curated SMARTS rules across 3 groups: Ring Systems, Functional Groups, and Pharmacophoric/Drug-class features. Multi-category matching per molecule. | No limit |
| **Chemical Space** | t-SNE (openTSNE) and PCA visualization of chemical space using Morgan ECFP4 fingerprints | t-SNE: max 2,000 molecules |
| **Registration Hashing** | Canonical deduplication using RDKit's registration hash with tautomer hash v2, identifying exact duplicates | No limit |
| **Scaffold Analysis** | Murcko scaffold mining with Shannon entropy diversity metric and frequency distribution | No limit |

For detailed analytics documentation, see the [Scoring Methodology](SCORING_METHODOLOGY.md#batch-analytics).

---

## QSAR-Ready Pipeline

The QSAR-Ready Pipeline (`/qsar-ready`) prepares chemical datasets for machine learning by applying a multi-step curation pipeline.

### Pipeline Steps

| Step | Description |
|------|-------------|
| **Standardization** | ChEMBL-compatible structure normalization |
| **Salt Stripping** | Remove counterions and salts, extract parent molecule |
| **Neutralization** | Neutralize charged species |
| **Tautomer Canonicalization** | Canonicalize tautomeric forms |
| **Duplicate Removal** | Remove duplicates by InChIKey |

### How to Use

**Web Interface:**
1. Navigate to **QSAR-Ready** under the Data Preparation dropdown
2. Paste SMILES (one per line) or upload a CSV/SDF file
3. Configure pipeline options
4. Click **Process**
5. Monitor progress via real-time status updates
6. Review results with per-molecule step details
7. Download curated dataset in CSV, SDF, or JSON format

**API — Single molecule:**
```bash
curl -X POST http://localhost:8001/api/v1/qsar-ready/single \
  -H "Content-Type: application/json" \
  -d '{
    "smiles": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",
    "config": {}
  }'
```

**API — Batch upload:**
```bash
curl -X POST http://localhost:8001/api/v1/qsar-ready/batch/upload \
  -F "file=@molecules.csv" \
  -F 'config={}'
```

**API — Check status and download:**
```bash
# Check status
curl http://localhost:8001/api/v1/qsar-ready/batch/{job_id}/status

# Get paginated results
curl "http://localhost:8001/api/v1/qsar-ready/batch/{job_id}/results?page=1&per_page=50"

# Download as CSV
curl http://localhost:8001/api/v1/qsar-ready/batch/{job_id}/download/csv -o curated.csv
```

### Result Status

Each molecule receives one of these statuses:

| Status | Meaning |
|--------|---------|
| **ok** | Successfully curated |
| **rejected** | Failed a pipeline step (see `rejection_reason`) |
| **duplicate** | Duplicate of another molecule (by InChIKey) |
| **error** | Processing error |

---

## Structure Filter

The Structure Filter (`/structure-filter`) provides multi-stage funnel filtering for generative chemistry outputs, screening libraries, or any SMILES collection.

### Filter Stages

Molecules pass through a configurable sequence of filter stages. Each stage can be enabled/disabled:

- **Property filters** — MW, LogP, TPSA, HBD, HBA, rotatable bonds, aromatic rings
- **Substructure matching** — Custom SMARTS patterns for inclusion/exclusion
- **Drug-likeness presets** — Pre-configured filter chains

### Presets

| Preset | Description |
|--------|-------------|
| **drug_like** | Lipinski-based drug-likeness criteria |
| **lead_like** | Lead-likeness criteria (MW 200–450, LogP −1 to 4) |
| **fragment_like** | Rule of Three fragment criteria |
| **permissive** | Minimal filtering (basic validity only) |

### How to Use

**Web Interface:**
1. Navigate to **Structure Filter** under the Library dropdown
2. Paste SMILES or upload a file
3. Select a preset or configure custom filters
4. Click **Filter**
5. View the funnel visualization showing molecules passing/failing each stage
6. Download passing molecules

**API — Synchronous (≤1,000 molecules):**
```bash
curl -X POST http://localhost:8001/api/v1/structure-filter/filter \
  -H "Content-Type: application/json" \
  -d '{
    "smiles_list": ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"],
    "preset": "drug_like"
  }'
```

**API — Batch upload (>1,000 molecules):**
```bash
curl -X POST http://localhost:8001/api/v1/structure-filter/batch/upload \
  -F "file=@molecules.csv" \
  -F "preset=drug_like"
```

### Scoring Mode

The Structure Filter also provides a 0–1 scoring mode for integration with generative models:

```bash
curl -X POST http://localhost:8001/api/v1/structure-filter/score \
  -H "Content-Type: application/json" \
  -d '{
    "smiles_list": ["CCO", "c1ccccc1"],
    "preset": "drug_like"
  }'
```

A REINVENT-compatible scoring endpoint is also available at `POST /structure-filter/reinvent-score`.

---

## Dataset Audit

The Dataset Audit (`/dataset-audit`) provides comprehensive dataset health auditing for chemical datasets.

### Features

| Feature | Description |
|---------|-------------|
| **Health Score** | Overall dataset quality score with sub-scores for validity, diversity, property distributions, and standardization |
| **Contradictory Labels** | Detects molecules with contradictory activity labels (same structure, different labels) |
| **Dataset Diff** | Compare two dataset versions to find added, removed, and modified molecules |
| **Curation Report** | Detailed curation recommendations and downloadable curated CSV |
| **Treemap Drill-down** | Interactive treemap visualization of dataset issues |

### How to Use

**Web Interface:**
1. Navigate to **Dataset Audit** under the Data Preparation dropdown
2. Upload a CSV or SDF file (with optional activity column)
3. Monitor audit progress (health audit, contradiction detection, curation)
4. Review health score and sub-scores
5. Explore contradictory labels, property distributions, and issues
6. Optionally upload a second file for dataset diff comparison
7. Download the curation report or curated CSV

**API:**
```bash
# Upload dataset for auditing
curl -X POST http://localhost:8001/api/v1/dataset/upload \
  -F "file=@dataset.csv" \
  -F "smiles_column=SMILES" \
  -F "activity_column=Activity"

# Check audit status
curl http://localhost:8001/api/v1/dataset/{job_id}/status

# Get audit results
curl http://localhost:8001/api/v1/dataset/{job_id}/results

# Compare with another dataset version
curl -X POST http://localhost:8001/api/v1/dataset/{job_id}/diff \
  -F "file=@dataset_v2.csv"

# Download curation report
curl http://localhost:8001/api/v1/dataset/{job_id}/download/report -o report.json

# Download curated CSV
curl http://localhost:8001/api/v1/dataset/{job_id}/download/csv -o curated.csv
```

---

## Diagnostics

The Diagnostics page (`/diagnostics`, under the Data Preparation dropdown) provides low-level chemical structure analysis tools:

| Tool | Description |
|------|-------------|
| **SMILES Diagnostics** | Parse and analyze SMILES string validity, atom-by-atom |
| **InChI Layer Diff** | Compare InChI layers between two molecules to identify exact structural differences |
| **Round-Trip Validation** | Validate SMILES→MOL→SMILES round-trip fidelity to detect representation loss |
| **File Pre-Validation** | Check a file for parseable SMILES before uploading to batch processing |
| **Coordinate Dimension Analysis** | Detect 2D, 3D, or missing coordinate data in MOL blocks/SDF files |

These tools are useful for debugging parsing errors, investigating structural discrepancies, and pre-screening files before batch processing.

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

ChemAudit provides comprehensive molecular scoring across 10+ dimensions. For complete calculation formulas, thresholds, and academic references, see the [Scoring Methodology](SCORING_METHODOLOGY.md).

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

### Ligand Efficiency

Measures how efficiently a molecule uses its heavy atoms to achieve binding potency.

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| **LE (Ligand Efficiency)** | pActivity / heavy atom count | ≥ 0.3 is generally considered efficient |
| **LLE (Lipophilic Ligand Efficiency)** | pActivity − LogP | Higher is better; measures potency relative to lipophilicity |

When no experimental activity data is available, ChemAudit uses QED as a proxy pIC50 value.

### Bioavailability Radar

A 6-axis radar chart showing whether a molecule falls within optimal ranges for oral bioavailability:

| Axis | Property | Optimal Range |
|------|----------|---------------|
| **LIPO** | Wildman-Crippen LogP | −0.7 to 5.0 |
| **SIZE** | Molecular Weight | 150–500 Da |
| **POLAR** | TPSA | 20–130 A² |
| **INSOLU** | LogS (ESOL) | −6 to 0 |
| **INSATU** | Fsp3 | ≥ 0.25 |
| **FLEX** | Rotatable Bonds | ≤ 9 |

All axes are normalized to 0–1. A molecule with all 6 axes in range has excellent predicted oral bioavailability.

### Property Breakdown

Per-atom contribution analysis for key molecular properties:

- **TPSA Breakdown** — Per-atom topological polar surface area contributions with functional group summaries
- **LogP Breakdown** — Per-atom Wildman-Crippen LogP contributions with functional group summaries

These breakdowns help identify which atoms or functional groups contribute most to a property, useful for lead optimization.

### Salt Inventory

Detects and inventories salt forms in a molecule:

- **has_salts** — Whether the molecule contains salt fragments
- **parent_smiles** — The parent (desalted) structure
- **fragments** — List of all fragments with classification
- **total_fragments** — Fragment count

**How to Get Scores (API):**
```bash
curl -X POST http://localhost:8001/api/v1/score \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "CCO",
    "include": ["ml_readiness", "druglikeness", "safety_filters", "admet", "np_likeness", "scaffold", "aggregator", "ligand_efficiency", "bioavailability_radar", "salt_inventory"]
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
| **PubChem** | Properties, synonyms, IUPAC name, CID | 30 req/min |
| **ChEMBL** | Bioactivity, targets, clinical phase, ChEMBL ID | 30 req/min |
| **COCONUT** | Natural product data, organism source | 30 req/min |
| **Wikidata** | SMILES (isomeric), InChI, InChIKey, CAS, formula, mass | 30 req/min |

### How to Search

**Web Interface:**
1. Enter your molecule on the Single Validation page
2. Navigate to **Database Lookup** tab
3. View cross-references from PubChem, ChEMBL, COCONUT, and Wikidata
4. The **Cross-Database Comparison** panel runs automatically, comparing structural identifiers across all databases that returned results

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

# Search Wikidata
curl -X POST http://localhost:8001/api/v1/integrations/wikidata/lookup \
  -H "Content-Type: application/json" \
  -d '{"molecule": "CCO", "format": "smiles"}'
```

### Identifier Resolution

ChemAudit can resolve a wide range of chemical identifiers to canonical SMILES with cross-references. Supported identifier types:

| Type | Example |
|------|---------|
| SMILES | `CCO` |
| InChI | `InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3` |
| InChIKey | `LFQSCWFLJHTTHZ-UHFFFAOYSA-N` |
| PubChem CID | `702` |
| ChEMBL ID | `CHEMBL545` |
| CAS Number | `64-17-5` |
| DrugBank ID | `DB00898` |
| ChEBI ID | `CHEBI:16236` |
| UNII | `3K9958V90M` |
| Wikipedia URL | `https://en.wikipedia.org/wiki/Ethanol` |
| Compound name | `aspirin`, `caffeine` |

The identifier type is auto-detected, or you can specify it explicitly.

**API:**
```bash
curl -X POST http://localhost:8001/api/v1/integrations/resolve \
  -H "Content-Type: application/json" \
  -d '{"identifier": "CHEMBL25", "identifier_type": "auto"}'
```

The response includes the resolved SMILES, InChI, InChIKey, molecular formula, weight, IUPAC name, resolution source, resolution chain, and cross-references to PubChem, ChEMBL, DrugBank, ChEBI, and KEGG.

### Database Comparison

Compare how different databases represent the same molecule. ChemAudit fetches the structure from PubChem, ChEMBL, COCONUT, and Wikidata, then compares:

| Identifier | Comparison Method |
|------------|-------------------|
| **SMILES** | RDKit canonicalization, then tautomer-invariant comparison |
| **InChIKey** | Layer-by-layer analysis (connectivity, stereo, version) |
| **InChI** | Layer analysis (formula, connectivity, hydrogen, charge, stereo) |

**Verdict values:** `consistent`, `minor_differences`, `major_discrepancies`, `no_data`

**API:**
```bash
curl -X POST http://localhost:8001/api/v1/integrations/compare \
  -H "Content-Type: application/json" \
  -d '{"smiles": "C[C@H](N)C(=O)O"}'
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
