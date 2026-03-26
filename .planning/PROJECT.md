# ChemAudit — Advanced Chemical Structure Validation Suite

## What This Is

A comprehensive web-based chemical structure validation, standardization, and assessment platform. React frontend + FastAPI backend for validating, standardizing, scoring, and analyzing chemical structures — both individually and in batch. Used by medicinal chemists, cheminformatics scientists, and ML practitioners to curate compound datasets.

## Core Value

Every chemical structure submitted gets a thorough, transparent, and reproducible quality assessment — from basic validity through ML-readiness — so scientists can trust their molecular data.

## Requirements

### Validated

- ✓ Single molecule validation with 28 checks (parsability, sanitization, valence, aromaticity, connectivity, stereo, representation, deep validation) — v1.0 + v2.0
- ✓ Batch validation up to 10K molecules with Celery async processing — v1.0
- ✓ 4-stage ChEMBL standardization pipeline with full provenance tracking — v1.0 + v2.0
- ✓ ML-readiness scoring (451 descriptors + 7 fingerprint types, 100-point scale) — v1.0
- ✓ Drug-likeness scoring (Lipinski, Veber, QED, Ghose, Egan, Muegge, consensus, lead-likeness) — v1.0 + v2.0
- ✓ NP-likeness scoring with fragment breakdown — v1.0 + v2.0
- ✓ ADMET property scoring (SA score, ESOL, Fsp3, CNS MPO, Pfizer 3/75, GSK 4/400, Golden Triangle) — v1.0 + v2.0
- ✓ Safety filter screening (PAINS, Brenk, NIH, ZINC, 7 ChEMBL catalogs) — v1.0
- ✓ Deep validation: stereo, tautomer, composition, complexity (17 checks) — v2.0
- ✓ Standardization provenance with atom-level diffs — v2.0
- ✓ Batch analytics: dedup, scaffold, chemical space, MMP, statistics — v2.0
- ✓ Scoring expansion: consensus, property breakdowns, bioavailability radar, BOILED-Egg — v2.0
- ✓ 9 interactive visualizations (histograms, scatter, treemaps, radar, timeline) — v2.0
- ✓ Advanced exports (fingerprint, dedup, scaffold, property matrix, enhanced PDF) — v2.0
- ✓ Custom scoring profiles with 8 presets + bookmarks + audit trail — v2.0
- ✓ Webhooks, email notifications, permalinks — v2.0
- ✓ IUPAC name input via OPSIN — v2.0
- ✓ Universal identifier resolver (SMILES, InChI, CAS, ChEMBL, PubChem, DrugBank, etc.) — v2.1
- ✓ Cross-database comparison (PubChem vs ChEMBL vs COCONUT vs Wikidata) — v2.1
- ✓ Multi-format export (CSV, Excel, SDF, JSON, PDF) — v1.0
- ✓ WebSocket real-time batch progress — v1.0
- ✓ API key authentication, CSRF protection, rate limiting — v1.0
- ✓ Docusaurus documentation site — v1.0
- ✓ Docker deployment with configurable profiles — v1.0

### Active

## Current Milestone: v3.0 Advanced Profiling, Safety Intelligence & Dataset Curation

**Goal:** Transform ChemAudit from a validation/scoring suite into a comprehensive compound profiling and dataset curation platform.

**Target features:**

**Phase 7 — Compound Profiling Engine (9 features)**
- [ ] PFI (Property Forecast Index) with risk classification
- [ ] #stars (QikProp-equivalent outlier count)
- [ ] Abbott Bioavailability Score (4-class probability)
- [ ] Consensus LogP (multi-method average)
- [ ] Skin Permeation (Potts-Guy log Kp)
- [ ] 3D Shape Descriptors (PMI/NPR/PBF with shape classification)
- [ ] Custom MPO Framework (desirability functions + CNS MPO preset)
- [ ] Extended Ligand Efficiency (LLE, LELP, SEI on top of existing LE/BEI)
- [ ] SA Comparison (SA_Score + SCScore + SYBA side-by-side)

**Phase 8 — Enhanced Structural Alerts & Safety (8 features)**
- [ ] Custom SMARTS alerts (reactive, chelator, fluorescent, redox, phospholipidosis — 21 patterns)
- [ ] Kazius mutagenicity toxicophores (29 patterns)
- [ ] NIBR Novartis screening deck filters (via medchem library)
- [ ] Complexity percentile filter (vs commercial compound distributions)
- [ ] CYP metabolism soft-spot prediction (SMARTS-based, atom-level)
- [ ] hERG liability flag (rule-based amphiphile assessment)
- [ ] bRo5 beyond-Rule-of-5 (for macrocycles/PROTACs)
- [ ] REOS Rapid Elimination of Swill filter

**Phase 9 — Structure Quality Diagnostics (5 features)**
- [ ] SMILES error diagnostics (position-specific errors with fix suggestions)
- [ ] InChI layer-by-layer diff comparison
- [ ] Format round-trip lossiness checker (SMILES↔InChI↔MOL)
- [ ] Cross-pipeline standardization comparison (RDKit vs ChEMBL vs minimal)
- [ ] File format pre-validator (SDF block integrity, CSV structure)

**Phase 10 — QSAR-Ready Pipeline (3 features)**
- [ ] 10-step configurable curation pipeline (parse→metals→desalt→normalize→neutralize→tautomer→stereo→isotope→filter→canonical)
- [ ] Batch file upload with QSAR-ready processing + InChIKey deduplication
- [ ] Preset configurations (QSAR-2D, QSAR-3D, custom)

**Phase 11 — Generative Chemistry Filter (3 features)**
- [ ] Multi-stage funnel pipeline (parse→valence→alerts→rules→SA→dedup)
- [ ] REINVENT-compatible REST scoring API (SMILES list → 0-1 scores)
- [ ] 4 preset configurations (drug-like, lead-like, fragment-like, permissive)

**Phase 12 — Dataset Intelligence (4 features)**
- [ ] Dataset health audit (composite 0-100 from parsability, stereo, uniqueness, alerts, standardization)
- [ ] Contradictory label detection (same InChIKey with divergent activity)
- [ ] Dataset diff tool (compare two uploads, show added/removed/modified)
- [ ] Reproducible curation report (serializable JSON audit trail)

**Phase 13 — Batch Analytics Extensions (4 features)**
- [ ] Butina compound clustering (Tanimoto distance, configurable cutoff)
- [ ] MCS molecule comparator (maximum common substructure + property deltas)
- [ ] ClassyFire-style chemical taxonomy (SMARTS-based chemotype classification)
- [ ] Registration hash (RDKit RegistrationHash for compound uniqueness)

**Phase 14 — Ecosystem & Workflow Integration (3 features)**
- [ ] MCP Server (via fastapi-mcp, expose validation/scoring/standardization as LLM tools)
- [ ] CLI tool (chemaudit-cli wrapping Python client)
- [ ] SureChEMBL patent presence check (InChIKey-based lookup)

### Out of Scope

- External ML model integration — all features use RDKit/MolVS/ChEMBL pipeline + rule-based methods only
- 3D conformation generation as deliverable — conformer generation used internally for PMI/NPR/PBF computation only (transient, not exported)
- Reaction validation — molecules only, not reactions
- User accounts/multi-tenancy — API key auth is sufficient
- Real-time collaboration — single-user workflow
- RAscore — requires Python 3.7 + TensorFlow 2.5, production-incompatible
- Lilly Medchem Rules (275-rule demerit system) — requires C++ compilation + Ruby driver, no Python package; comparable coverage via RDKit FilterCatalog + medchem NIBR + Kazius + custom SMARTS
- UMAP dimensionality reduction — heavy C++ dependency, outside "zero external ML models" constraint

## Context

- **Existing codebase**: Fully functional application with v1.0 (5 phases) + v2.0 (6 phases) + identifier resolution/cross-db comparison delivered
- **Feature specifications**: Planning docs in `.planning/New_featues/` (4 guides: complete, backend, implementation, UI/UX) with verified API references and test SMILES
- **Architecture**: Plugin-based validation checks via `@CheckRegistry.register()` decorator pattern — 28 checks registered
- **Scoring**: 12 scoring modules operational (druglikeness, admet, ml_readiness, np_likeness, aggregator, bioavailability_radar, property_breakdown, salt_inventory, safety_filters, scaffold, profile_scoring)
- **Alerts**: AlertManager + FilterCatalog with PAINS/Brenk/NIH/ZINC/ChEMBL catalogs
- **Batch processing**: Celery chord pattern with Redis, analytics pipeline, 6 analytics modules
- **Integrations**: Universal resolver, cross-database comparison, PubChem/ChEMBL/COCONUT/Wikidata/UniChem
- **Frontend**: React 18 + TypeScript + Tailwind + Framer Motion + RDKit.js + Recharts, 9 pages, 78 component files
- **New dependencies for v3.0**: medchem (NIBR filters, pip install), vendored SCScore (numpy standalone), SYBA (from GitHub source)

## Constraints

- **Tech stack**: RDKit, MolVS, ChEMBL Pipeline, NumPy/SciPy, medchem — zero external ML model APIs
- **Backward compatibility**: All existing API endpoints, schemas, and behaviors must be preserved. New alerts layer on top of existing safety_filters.py
- **Quality bar**: Production code only — no pseudo code, no bandaids, no quick fixes. Tests per feature.
- **Batch analytics**: Heavy compute features must run as Celery async tasks
- **Full stack**: Each phase delivers backend + frontend together, fully integrated
- **Accuracy priority**: This is a production system used by scientists — correctness of chemical calculations takes priority over speed of delivery

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Follow recommended build order from features roadmap | Respects dependency chains, maximizes impact early | ✓ Good (v2.0) |
| All 79 features in single milestone cycle | User wants comprehensive coverage | ✓ Good (v2.0) |
| Adapt implementations from milestone docs | Use as guidance, improve where better approaches exist | ✓ Good (v2.0) |
| Celery for batch analytics | Consistent with existing async architecture | ✓ Good (v2.0) |
| Full stack per phase | Each phase ships with backend + frontend integrated | ✓ Good (v2.0) |
| Reuse existing scoring modules | PFI/3/75/GSK/etc. already exist — build new profiler on top, don't duplicate | — Pending |
| Alert system as additional layer | New custom/Kazius/NIBR alerts coexist with existing safety_filters.py | — Pending |
| medchem for NIBR filters | Production-ready, pip-installable, Apache-2.0, actively maintained | — Pending |
| SCScore vendored (numpy standalone) | No TF dependency, single file + weights, reliable for production | — Pending |
| Skip RAscore | Requires Python 3.7 + TF 2.5, PyPI name collision, production-incompatible | — Pending |
| Skip Lilly Medchem Rules | C++/Ruby, custom query format, no Python package; ~1,700 patterns from other sources provide comparable coverage | — Pending |
| Include 3D shape descriptors | Conformer generation is transient computation for PMI/NPR/PBF — not a standalone feature | — Pending |
| MCP Server via fastapi-mcp | Low effort, first validation/scoring MCP server in cheminformatics, high strategic value | — Pending |

## Evolution

This document evolves at phase transitions and milestone boundaries.

**After each phase transition** (via `/gsd:transition`):
1. Requirements invalidated? → Move to Out of Scope with reason
2. Requirements validated? → Move to Validated with phase reference
3. New requirements emerged? → Add to Active
4. Decisions to log? → Add to Key Decisions
5. "What This Is" still accurate? → Update if drifted

**After each milestone** (via `/gsd:complete-milestone`):
1. Full review of all sections
2. Core Value check — still the right priority?
3. Audit Out of Scope — reasons still valid?
4. Update Context with current state

---
*Last updated: 2026-03-26 after v3.0 milestone initialization*
