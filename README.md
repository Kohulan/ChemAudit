<div align="center">

<img src="docs/assets/Banner_logo.png" alt="ChemAudit" width="600" />

### Chemical Structure Validation & Quality Assessment Platform

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18492676.svg)](https://doi.org/10.5281/zenodo.18492676)
[![Release](https://img.shields.io/github/v/release/Kohulan/ChemAudit)](https://github.com/Kohulan/ChemAudit/releases)
[![Tests](https://github.com/Kohulan/ChemAudit/actions/workflows/test.yml/badge.svg)](https://github.com/Kohulan/ChemAudit/actions/workflows/test.yml)
[![Documentation](https://img.shields.io/badge/docs-available-brightgreen.svg)](https://kohulan.github.io/ChemAudit/)
[![Contributors](https://img.shields.io/github/contributors/Kohulan/ChemAudit)](https://github.com/Kohulan/ChemAudit/graphs/contributors)
[![Issues](https://img.shields.io/github/issues/Kohulan/ChemAudit)](https://github.com/Kohulan/ChemAudit/issues)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![React 18](https://img.shields.io/badge/react-18-61DAFB.svg?logo=react)](https://reactjs.org/)
[![FastAPI](https://img.shields.io/badge/FastAPI-009688.svg?logo=fastapi&logoColor=white)](https://fastapi.tiangolo.com/)
[![RDKit](https://img.shields.io/badge/RDKit-3776AB.svg)](https://www.rdkit.org/)
[![Docker](https://img.shields.io/badge/docker-ready-2496ED.svg?logo=docker&logoColor=white)](https://www.docker.com/)

<br />

**Validate &bull; Standardize &bull; Score &bull; Profile &bull; Curate &bull; Analyze**

*A comprehensive web platform for cheminformatics workflows, drug discovery, ML dataset curation, and generative chemistry evaluation*

[Features](#-features) &bull;
[Quick Start](#-quick-start) &bull;
[Documentation](https://www.kohulanr.com/ChemAudit/) &bull;
[API](#-api-reference) &bull;
[CLI](#-command-line-interface) &bull;
[Contributing](#-contributing)

<br />

<img src="docs/assets/screenshot-dashboard.png" alt="ChemAudit Dashboard" width="800" />

</div>

---

## ✨ Features

<table>
<tr>
<td width="50%">

### Structure Validation
Comprehensive chemical structure analysis with 15+ validation checks

- Valence, connectivity & aromaticity errors
- Stereochemistry & ring system analysis
- Atom & bond type verification
- Configurable validation profiles with 8 presets

</td>
<td width="50%">

### Structural Alerts
Screen compounds against known problematic substructures

- **PAINS** - Pan-Assay Interference detection
- **BRENK** - Unwanted group filters
- **Kazius & NIBR** - Additional alert catalogs
- Concern-group deduplication across catalogs

</td>
</tr>
<tr>
<td width="50%">

### ML-Readiness Scoring
Evaluate compound suitability for machine learning

- 4-dimension scientific assessment
- Descriptor calculability & fingerprint validation
- Molecular complexity scoring
- Profile-based scoring with custom presets

</td>
<td width="50%">

### Standardization Pipeline
ChEMBL-compatible molecular standardization

- Salt stripping & neutralization
- Tautomer canonicalization
- Stereochemistry normalization
- Cross-pipeline comparison (RDKit vs ChEMBL)

</td>
</tr>
</table>

### Data Preparation Suite

<table>
<tr>
<td width="50%">

#### QSAR-Ready Pipeline
A 10-step curation pipeline producing ML-ready structures

- Desalting, neutralisation, tautomer canonicalisation
- InChIKey deduplication with change tracking per stage
- Configurable steps with 3 presets (QSAR-2D, QSAR-3D, Minimal)
- Batch processing via Celery with WebSocket progress

</td>
<td width="50%">

#### Structure Filter
Validation funnel for generative model output (REINVENT, etc.)

- 6-stage funnel: parse, valence, alerts, property, SA, dedup
- Optional novelty check via ChEMBL Tanimoto similarity
- REINVENT-compatible scoring endpoint
- Interactive funnel visualisation showing drop-out points

</td>
</tr>
<tr>
<td width="50%">

#### Dataset Audit
Upload a dataset and get a comprehensive health score

- Structural issues, standardisation inconsistencies, property distributions
- Contradictory label detection (same InChIKey, opposite activity)
- Full curation report generation
- Dataset diff to compare before/after curation

</td>
<td width="50%">

#### Diagnostics
Debugging toolkit for chemical data issues

- SMILES error annotation with exact break position
- InChI layer-by-layer diff tables
- SMILES-InChI-SMILES round-trip loss detection
- Cross-pipeline comparison (RDKit vs ChEMBL vs Minimal)
- File pre-validator: encoding, duplicate headers, malformed rows

</td>
</tr>
</table>

### Compound Profiler

<table>
<tr>
<td width="50%">

- **PFI** - Property Forecast Index (Young et al. 2011)
- **Ligand Efficiency** - LE, LLE, LELP, BEI, SEI metrics
- **SA Comparison** - Three independent scorers (SA Score, SCScore, SYBA)
- **3D Shape Analysis** - PMI ternary plots via ETKDGv3/MMFF94

</td>
<td width="50%">

- **Skin Permeation** - Potts-Guy 1992 model
- **CNS MPO** - Wager 2010 multi-parameter optimisation
- **Custom MPO Editor** - Define your own desirability functions with sigmoid/ramp/step shapes per property

</td>
</tr>
</table>

### Safety Assessment

- **CYP Soft Spot Prediction** - 11 SMARTS-based metabolic liability patterns with atom-level highlighting
- **hERG Liability Scoring** - 4-factor amphiphile risk model
- **Beyond Rule of 5 (bRo5)** - Doak et al. 2014 thresholds for extended chemical space
- **REOS Filtering** - Walters 1999 rapid elimination of swill
- **Concern-Group Deduplication** - Same alert not flagged by multiple catalogs

### Database Integrations & Identifier Resolution

- **Universal Identifier Resolution** - Paste SMILES, InChI, InChIKey, CAS, common name, ChEMBL ID, PubChem CID, or more
- **Cross-Reference** across PubChem, ChEMBL, COCONUT, Wikidata, ChEBI, and UniChem
- **Side-by-side comparison** with stereochemistry-aware diffing
- **SureChEMBL** patent presence lookups via UniChem cross-reference

### Batch Analytics

- **Butina Clustering** with configurable distance cutoff
- **Murcko Scaffold Analysis** with interactive treemap visualisation
- **Chemotype Taxonomy** - SMARTS-based classification with 50+ rules
- **Registration Hash Collision Detection** - RDKit RegistrationHash v2 with tautomer support
- **MCS Comparison** - Select any two molecules for side-by-side analysis with Tanimoto similarity, common substructure, and property delta table
- **Click-to-filter drill-down** across all analytics charts
- **Shareable Permalink Reports** with auto-snapshot persistence

### Export System

Full 78-column audit trail across 6 sections (Validation, Deep Validation, Scoring, Safety, Compound Profile, Standardization). One declarative registry drives all five exporters:

| Format | Features |
|--------|----------|
| **CSV** | All 78 audit columns, always included |
| **Excel** | Single-sheet or multi-sheet layout, optional 2D structure depictions, conditional score formatting |
| **JSON** | Nested audit structure |
| **SDF** | Optional full audit data toggle |
| **PDF** | Batch report with optional full audit data toggle |

---

## 🚀 Quick Start

### Using Docker (Recommended)

```bash
# Clone the repository
git clone https://github.com/Kohulan/ChemAudit.git
cd chemaudit

# Start all services (development)
docker-compose up -d

# View logs
docker-compose logs -f
```

**Access Points:**
| Service | URL |
|---------|-----|
| Web UI | http://localhost:3002 |
| API Docs | http://localhost:8001/docs |
| MCP Server | http://localhost:8001/mcp |
| Metrics | http://localhost:9090 |

### Production Deployment

Use the interactive deploy script to select a deployment profile:

```bash
# Interactive mode - shows profile menu
./deploy.sh

# Or specify profile directly
./deploy.sh medium
```

**Available Profiles:**

| Profile | Max Molecules | Max File Size | Workers | Use Case |
|---------|---------------|---------------|---------|----------|
| `small` | 1,000 | 100 MB | 2 | Development |
| `medium` | 10,000 | 500 MB | 4 | Standard production |
| `large` | 50,000 | 500 MB | 8 | High-throughput |
| `xl` | 100,000 | 1 GB | 12 | Enterprise |
| `coconut` | 1,000,000 | 1 GB | 16 | Full COCONUT DB |

See [Deployment Guide](docs/DEPLOYMENT.md) for detailed configuration.

### Manual Installation

<details>
<summary><b>Backend Setup</b></summary>

```bash
cd backend
pip install -e ".[dev]"
uvicorn app.main:app --reload --port 8000
```

</details>

<details>
<summary><b>Frontend Setup</b></summary>

```bash
cd frontend
npm install
npm run dev
```

</details>

---

## 💻 Command-Line Interface

ChemAudit includes a CLI tool with 4 subcommands:

```bash
# Validate a molecule
chemaudit validate --smiles "CCO" --format table

# Score a molecule
chemaudit score --smiles "CC(=O)Oc1ccccc1C(=O)O"

# Standardize a molecule
chemaudit standardize --smiles "CC(=O)[O-].[Na+]"

# Profile a molecule
chemaudit profile --smiles "c1ccccc1"

# Process a file
chemaudit validate --file molecules.csv --format json

# Pipe from stdin
echo "CCO" | chemaudit validate --format table
```

All commands support `--local` (offline mode), `--server` (remote API), and `--format json|table` output.

---

## 📦 Batch Processing

Process large datasets with real-time progress tracking:

| Feature | Specification |
|---------|---------------|
| **Max File Size** | Up to 1 GB (profile-dependent) |
| **Max Molecules** | Up to 1,000,000 per batch (profile-dependent) |
| **Supported Formats** | SDF, CSV |
| **Progress Tracking** | Real-time WebSocket updates |
| **Export Formats** | CSV, JSON, Excel, SDF, PDF Report |

```python
# Python client example
from chemaudit import ChemAuditClient

client = ChemAuditClient("http://localhost:8000")

# Upload and process
job = client.upload_batch("molecules.sdf")

# Monitor progress
for update in client.stream_progress(job.job_id):
    print(f"Progress: {update.progress}%")

# Get results
results = client.get_results(job.job_id)
```

---

## 📸 Screenshots

<div align="center">
<table>
<tr>
<td align="center"><b>Single Molecule Validation</b></td>
<td align="center"><b>Batch Processing</b></td>
</tr>
<tr>
<td><img src="docs/assets/screenshot-single.png" width="400" /></td>
<td><img src="docs/assets/screenshot-batch.png" width="400" /></td>
</tr>
<tr>
<td align="center"><b>Scoring Dashboard</b></td>
<td align="center"><b>Database Lookup</b></td>
</tr>
<tr>
<td><img src="docs/assets/screenshot-scoring.png" width="400" /></td>
<td><img src="docs/assets/screenshot-integrations.png" width="400" /></td>
</tr>
</table>
</div>

---

## 🛠 Tech Stack

<div align="center">

| Layer | Technologies |
|-------|-------------|
| **Frontend** | ![React](https://img.shields.io/badge/React-61DAFB?logo=react&logoColor=black) ![TypeScript](https://img.shields.io/badge/TypeScript-3178C6?logo=typescript&logoColor=white) ![Vite](https://img.shields.io/badge/Vite-646CFF?logo=vite&logoColor=white) ![Tailwind](https://img.shields.io/badge/Tailwind-06B6D4?logo=tailwindcss&logoColor=white) |
| **Backend** | ![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white) ![FastAPI](https://img.shields.io/badge/FastAPI-009688?logo=fastapi&logoColor=white) ![RDKit](https://img.shields.io/badge/RDKit-3776AB) ![Celery](https://img.shields.io/badge/Celery-37814A?logo=celery&logoColor=white) |
| **Database** | ![PostgreSQL](https://img.shields.io/badge/PostgreSQL-4169E1?logo=postgresql&logoColor=white) ![Redis](https://img.shields.io/badge/Redis-DC382D?logo=redis&logoColor=white) |
| **Infrastructure** | ![Docker](https://img.shields.io/badge/Docker-2496ED?logo=docker&logoColor=white) ![Nginx](https://img.shields.io/badge/Nginx-009639?logo=nginx&logoColor=white) |
| **Monitoring** | ![Prometheus](https://img.shields.io/badge/Prometheus-E6522C?logo=prometheus&logoColor=white) ![Grafana](https://img.shields.io/badge/Grafana-F46800?logo=grafana&logoColor=white) |
| **AI Integration** | MCP Server (Model Context Protocol) for AI-assisted workflows |

</div>

---

## 📚 Documentation

| Document | Description |
|----------|-------------|
| [Getting Started](docs/GETTING_STARTED.md) | Installation and first steps |
| [User Guide](docs/USER_GUIDE.md) | Complete usage instructions |
| [API Reference](docs/API_REFERENCE.md) | Full REST API documentation |
| [Deployment](docs/DEPLOYMENT.md) | Production deployment guide |
| [Troubleshooting](docs/TROUBLESHOOTING.md) | Common issues & solutions |

**Interactive API Docs:** http://localhost:8000/docs

---

## 📡 API Reference

### Validate a Molecule

```bash
curl -X POST http://localhost:8000/api/v1/validate \
  -H "Content-Type: application/json" \
  -d '{"molecule": "CC(=O)Oc1ccccc1C(=O)O", "format": "smiles"}'
```

### Resolve an Identifier

```bash
curl -X POST http://localhost:8000/api/v1/integrations/resolve \
  -H "Content-Type: application/json" \
  -d '{"identifier": "aspirin"}'
```

### Run Safety Assessment

```bash
curl -X POST http://localhost:8000/api/v1/safety/assess \
  -H "Content-Type: application/json" \
  -d '{"smiles": "c1ccc2c(c1)nc(n2)Sc3nnnn3C"}'
```

### Profile a Compound

```bash
curl -X POST http://localhost:8000/api/v1/profiler/full \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CC(=O)Oc1ccccc1C(=O)O"}'
```

---

## 🏗 Project Structure

```
chemaudit/
├── backend/
│   ├── app/
│   │   ├── api/routes/          # 22 API route modules
│   │   ├── core/                # Config, security, rate limiting, sessions
│   │   ├── schemas/             # Pydantic v2 models
│   │   ├── services/            # Business logic
│   │   │   ├── validation/      # 15+ validation checks
│   │   │   ├── scoring/         # ML-readiness, NP-likeness
│   │   │   ├── alerts/          # PAINS, BRENK, Kazius, NIBR
│   │   │   ├── profiler/        # PFI, SA, MPO, ligand efficiency, 3D shape
│   │   │   ├── safety/          # CYP, hERG, bRo5, REOS, unified screen
│   │   │   ├── diagnostics/     # SMILES, InChI, round-trip, cross-pipeline
│   │   │   ├── qsar_ready/      # 10-step curation pipeline
│   │   │   ├── structure_filter/ # Generative model validation funnel
│   │   │   ├── dataset_intelligence/ # Health audit, contradictory labels, diff
│   │   │   ├── analytics/       # Clustering, scaffold, taxonomy, MCS, RegHash
│   │   │   ├── integrations/    # PubChem, ChEMBL, COCONUT, Wikidata, SureChEMBL
│   │   │   └── export/          # CSV, Excel, SDF, JSON, PDF (78-column audit)
│   │   ├── cli.py               # CLI tool (4 subcommands)
│   │   └── templates/           # PDF report templates
│   └── tests/
├── frontend/
│   ├── src/
│   │   ├── components/          # React components
│   │   ├── pages/               # 15 route pages
│   │   ├── services/            # API client
│   │   ├── hooks/               # Custom hooks
│   │   └── types/               # TypeScript definitions
│   └── tests/
├── client/                      # Python client library
├── docs-site/                   # Docusaurus documentation site
├── nginx/                       # Reverse proxy config
├── .github/workflows/           # CI/CD (tests, lint, release, secret scanning)
├── SECURITY.md                  # Security policy
└── docker-compose.yml
```

---

## 🔒 Security

ChemAudit implements defense-in-depth security:

- **API Key Authentication** with Redis-backed key management
- **Rate Limiting** per-IP and per-API-key tiers with progressive IP banning
- **Session Isolation** via HttpOnly cookies and PostgreSQL row-level security
- **WebSocket Ownership Checks** preventing cross-session data access
- **CSRF Protection** and Content Security Policy headers
- **Secret Scanning** via Gitleaks in CI

See [SECURITY.md](SECURITY.md) for vulnerability reporting.

---

## 🧪 Running Tests

```bash
# Backend tests
cd backend
pytest tests/ -v

# Frontend tests
cd frontend
npm test

# With coverage
pytest tests/ --cov=app --cov-report=xml
```

---

## 🤝 Contributing

Contributions are welcome! Please read our contributing guidelines before submitting a PR.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

<div align="center">

*Built on the shoulders of these outstanding open-source projects and databases*

<br />

<table>
<tr>
<td align="center" width="140">
<br />
<a href="https://www.rdkit.org/">
<img src="https://img.shields.io/badge/RDKit-3776AB?style=for-the-badge&logoColor=white" alt="RDKit" />
</a>
<br /><br />
<sub>Cheminformatics<br />Toolkit</sub>
<br /><br />
</td>
<td align="center" width="140">
<br />
<a href="https://www.ebi.ac.uk/chembl/">
<img src="frontend/public/assets/logos/chembl.png" alt="ChEMBL" height="50" />
</a>
<br />
<sub><b>ChEMBL</b></sub>
<br />
<sub>Bioactivity<br />Database</sub>
<br /><br />
</td>
<td align="center" width="140">
<br />
<a href="https://pubchem.ncbi.nlm.nih.gov/">
<img src="frontend/public/assets/logos/pubchem.png" alt="PubChem" height="50" />
</a>
<br />
<sub><b>PubChem</b></sub>
<br />
<sub>Chemical<br />Database</sub>
<br /><br />
</td>
<td align="center" width="140">
<br />
<a href="https://coconut.naturalproducts.net/">
<img src="frontend/public/assets/logos/coconut.png" alt="COCONUT" height="50" />
</a>
<br />
<sub><b>COCONUT</b></sub>
<br />
<sub>Natural Products<br />Database</sub>
<br /><br />
</td>
</tr>
<tr>
<td align="center" width="140">
<br />
<a href="https://www.ebi.ac.uk/chebi/">
<img src="https://img.shields.io/badge/ChEBI-5C2D91?style=for-the-badge&logoColor=white" alt="ChEBI" />
</a>
<br /><br />
<sub>Chemical Entities of<br />Biological Interest</sub>
<br /><br />
</td>
<td align="center" width="140">
<br />
<a href="https://www.ebi.ac.uk/unichem/">
<img src="https://img.shields.io/badge/UniChem-2E7D32?style=for-the-badge&logoColor=white" alt="UniChem" />
</a>
<br /><br />
<sub>Cross-Reference<br />Mapping Service</sub>
<br /><br />
</td>
<td align="center" width="140">
<br />
<a href="https://www.wikidata.org/">
<img src="frontend/public/assets/logos/wikidata.svg" alt="Wikidata" height="50" />
</a>
<br />
<sub><b>Wikidata</b></sub>
<br />
<sub>Open<br />Knowledge Base</sub>
<br /><br />
</td>
<td align="center" width="140">
<br />
<a href="https://www.surechembl.org/">
<img src="frontend/public/assets/logos/surechembl.png" alt="SureChEMBL" height="50" />
</a>
<br />
<sub><b>SureChEMBL</b></sub>
<br />
<sub>Patent Chemistry<br />Database</sub>
<br /><br />
</td>
</tr>
</table>

</div>

---

<div align="center">

<a href="https://www.uni-jena.de/">
<img src="frontend/public/cheminf-logo.png" alt="Natural Products Cheminformatics - Friedrich Schiller University Jena" height="80" />
</a>

<br /><br />

Made with ☕ by [Kohulan Rajan](https://github.com/Kohulan) at [Friedrich Schiller University Jena](https://www.uni-jena.de/)

<br />

**[Back to Top](#)**

</div>
