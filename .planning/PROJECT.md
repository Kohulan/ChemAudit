# ChemStructVal

## What This Is

A comprehensive web-based chemical structure validation platform that consolidates disparate validation tools into a unified, visual interface. It serves general chemists who currently validate structures manually with fragmented tools, providing accessible, professional-grade quality assessment without requiring programming expertise.

## Core Value

Any chemist can validate chemical structures comprehensively — from basic parsing to ML-readiness assessment — through an intuitive web interface, with batch processing for large datasets.

## Requirements

### Validated

(None yet — ship to validate)

### Active

**Validation Engine:**
- [ ] Single molecule validation with 8+ comprehensive checks
- [ ] SMILES, InChI, MOL/SDF format support with auto-detection
- [ ] Parsability, sanitization, valence, aromaticity checks
- [ ] Stereochemistry validation (undefined centers, conflicting stereo)
- [ ] Representation consistency (SMILES/InChI roundtrip)
- [ ] Overall quality score calculation (0-100)

**Structural Alerts:**
- [ ] PAINS A/B/C pattern screening (480+ patterns)
- [ ] BRENK alerts (105 patterns)
- [ ] NIH, ZINC, ChEMBL alert sets (BMS, Dundee, Glaxo, etc.)
- [ ] Matched atom highlighting in structure viewer

**Scoring & Assessment:**
- [ ] ML-readiness scoring (descriptor calculability, fingerprints, size)
- [ ] NP-likeness scoring (natural product likelihood)
- [ ] Scaffold analysis (Murcko extraction)

**Standardization:**
- [ ] ChEMBL-compatible standardization pipeline
- [ ] Salt stripping, solvent removal, normalization
- [ ] Tautomer canonicalization
- [ ] Parent structure extraction
- [ ] Before/after comparison view

**Batch Processing:**
- [ ] File upload (SDF, CSV) up to 10K molecules
- [ ] Progress tracking with WebSocket updates
- [ ] Job queue with Celery workers
- [ ] Partial failure handling
- [ ] Result aggregation and statistics

**External Integrations:**
- [ ] DECIMER OCSR output validation
- [ ] COCONUT natural products database lookup
- [ ] PubChem compound cross-reference
- [ ] ChEMBL bioactivity data lookup

**Frontend:**
- [ ] RDKit.js molecule rendering with atom highlighting
- [ ] Interactive validation results dashboard
- [ ] Batch upload with drag-and-drop
- [ ] Export (CSV, Excel, SDF, JSON)
- [ ] PDF report generation

**API:**
- [ ] RESTful API with OpenAPI documentation
- [ ] Rate limiting and optional API key auth
- [ ] Python client library

**Deployment:**
- [ ] Docker Compose for development
- [ ] Production-ready Docker images
- [ ] Self-hosting documentation

### Out of Scope

- Custom ML model training — beyond v1 scope, focus on validation
- Reaction validation — different problem domain
- Protein-ligand analysis — separate tool category
- Commercial database integration — licensing complexity
- Real-time collaborative editing — unnecessary for validation workflow
- Mobile app — web-first, responsive design sufficient

## Context

**Existing Documentation:**
- `docs/PRD.md` — Detailed product requirements
- `docs/PHASES.md` — 14-week development breakdown with checklists
- `docs/ARCHITECTURE.md` — System design, code structure, patterns
- `docs/CLAUDE_CODE_GUIDE.md` — Implementation snippets and patterns
- `docs/SKILLS.md` — Technical dependencies and requirements

**Codebase Map:**
- `.planning/codebase/` — Technology stack, architecture, conventions

**Target Users:**
General chemists (database curators, ML researchers, medicinal chemists, natural product scientists) who need accessible structure validation without CLI expertise.

**Problem Being Solved:**
Chemists currently validate structures manually using fragmented tools. ChemStructVal consolidates validation, alerts, standardization, and scoring into one visual interface.

## Constraints

- **Tech Stack**: React 18 + FastAPI + PostgreSQL + Redis + RDKit (documented, non-negotiable)
- **Performance**: <10 seconds for single validation, 10K molecules per batch
- **Accuracy**: 99%+ agreement with reference implementations (RDKit, ChEMBL pipeline)
- **Browser Support**: Modern browsers (Chrome, Firefox, Safari, Edge)
- **Deployment**: Must support both hosted web service and Docker self-hosting

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| RDKit.js for frontend rendering | Consistent with backend RDKit, no server roundtrip for display | — Pending |
| Celery + Redis for batch processing | Proven stack for async job processing at scale | — Pending |
| ChEMBL structure pipeline for standardization | Industry standard, well-tested | — Pending |
| PostgreSQL for job/result storage | Reliable, good JSON support, familiar | — Pending |

---
*Last updated: 2026-01-20 after initialization*
