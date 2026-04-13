<div align="center">

<img src="assets/logo.png" alt="ChemAudit" width="100" />

# Documentation

### ChemAudit Documentation Hub

</div>

---

## 🗂️ Quick Navigation

| Document | Description |
|----------|-------------|
| [**Getting Started**](GETTING_STARTED.md) | Installation and first steps |
| [**User Guide**](USER_GUIDE.md) | Complete feature walkthrough |
| [**Scoring Methodology**](SCORING_METHODOLOGY.md) | All formulas, thresholds, and citations |
| [**API Reference**](API_REFERENCE.md) | REST API documentation |
| [**Deployment**](DEPLOYMENT.md) | Production deployment guide |
| [**Troubleshooting**](TROUBLESHOOTING.md) | Common issues and solutions |

---

## 🎯 Quick Start

```bash
# Clone and start
git clone https://github.com/Kohulan/ChemAudit.git
cd chemaudit

# Create .env with required secrets
cp .env.example .env
# Edit .env to set POSTGRES_PASSWORD, SECRET_KEY, API_KEY_ADMIN_SECRET, CSRF_SECRET_KEY

# Development mode
docker-compose up -d
```

**Access (Development):**
- 🌐 Web UI: http://localhost:3002
- 📖 API Docs: http://localhost:8001/api/v1/docs
- 📖 API ReDoc: http://localhost:8001/api/v1/redoc

**Access (Production):**
- 🌐 Web UI + API: http://localhost (behind Nginx)

---

## ✨ Key Features

| Feature | Description |
|---------|-------------|
| **Single Validation** | Validate SMILES, InChI, MOL blocks, or IUPAC names with 27 structural checks |
| **Batch Processing** | Process up to 1M molecules with real-time WebSocket progress |
| **Batch Analytics** | Butina clustering, chemical taxonomy (~50 SMARTS rules), t-SNE/PCA chemical space, scaffold analysis, registration hashing |
| **Structural Alerts** | Screen against PAINS, BRENK, NIH, ZINC, and ChEMBL catalogs (~1,500+ patterns) |
| **Scoring** | 10+ scoring modules: ML-readiness, drug-likeness, ADMET, NP-likeness, scaffold analysis, aggregator likelihood, ligand efficiency, bioavailability radar, property breakdown, salt inventory |
| **QSAR-Ready Pipeline** | Multi-step ML-readiness curation: standardization, salt stripping, neutralization, tautomer canonicalization, duplicate removal |
| **Structure Filter** | Multi-stage funnel filtering with property filters, SMARTS substructure matching, presets (drug-like, lead-like, fragment-like), REINVENT scoring |
| **Dataset Audit** | Dataset health scoring, contradictory label detection, dataset diff/comparison, curation reports |
| **Identifier Resolution** | Resolve 10+ identifier types (SMILES, InChI, InChIKey, CID, ChEMBL, CAS, DrugBank, ChEBI, UNII, Wikipedia, names) with cross-database linking |
| **Database Comparison** | Cross-database structure comparison against PubChem, ChEMBL, COCONUT, Wikidata with consistency verdicts |
| **Diagnostics** | SMILES diagnostics, InChI layer diff, round-trip validation, file pre-validation, coordinate analysis |
| **Standardization** | ChEMBL-compatible pipeline: sanitize, get parent, remove salts, optional tautomer canonicalization |
| **Database Lookup** | Cross-reference PubChem, ChEMBL, COCONUT, and Wikidata |
| **Export** | CSV, Excel, SDF, JSON, PDF, fingerprint matrix, deduplicated set, scaffold-grouped, property matrix |

---

## 📖 Documentation Overview

### For Users

| Guide | What You'll Learn |
|-------|-------------------|
| [Getting Started](GETTING_STARTED.md) | Install, configure, validate your first molecule |
| [User Guide](USER_GUIDE.md) | All features: batch processing, alerts, scoring, standardization, database lookup, export |
| [Scoring Methodology](SCORING_METHODOLOGY.md) | Detailed formulas, thresholds, and academic references for all scoring algorithms |
| [Troubleshooting](TROUBLESHOOTING.md) | Fix common issues |

### For Developers

| Guide | What You'll Learn |
|-------|-------------------|
| [API Reference](API_REFERENCE.md) | All endpoints, parameters, request/response schemas, rate limits, WebSocket, API key auth |
| [Deployment](DEPLOYMENT.md) | Production setup with Docker, Nginx, deployment profiles, SSL, monitoring |

---

## 🔗 External Resources

- [RDKit Documentation](https://www.rdkit.org/docs/)
- [FastAPI Documentation](https://fastapi.tiangolo.com/)
- [React Documentation](https://react.dev/)

---

<div align="center">

**Questions?** [Open an Issue](https://github.com/Kohulan/ChemAudit/issues)

</div>
