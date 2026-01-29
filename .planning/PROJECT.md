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

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| Single molecule validation with 8+ checks | Phase 1 | Pending |
| SMILES, InChI, MOL/SDF format support | Phase 1 | Pending |
| Parsability, sanitization, valence, aromaticity checks | Phase 1 | Pending |
| Stereochemistry validation | Phase 1 | Pending |
| Representation consistency | Phase 1 | Pending |
| Overall quality score (0-100) | Phase 1 | Pending |
| PAINS A/B/C pattern screening | Phase 2 | Pending |
| BRENK alerts | Phase 2 | Pending |
| NIH, ZINC, ChEMBL alert sets | Phase 2 | Pending |
| Matched atom highlighting | Phase 2 | Pending |
| ML-readiness scoring | Phase 2 | Pending |
| NP-likeness scoring | Phase 2 | Pending |
| Scaffold analysis | Phase 2 | Pending |
| ChEMBL-compatible standardization | Phase 2 | Pending |
| Salt stripping, solvent removal, normalization | Phase 2 | Pending |
| Tautomer canonicalization | Phase 2 | Pending |
| Parent structure extraction | Phase 2 | Pending |
| Before/after comparison view | Phase 2 | Pending |
| File upload (SDF, CSV) up to 10K | Phase 2 | Pending |
| Progress tracking with WebSocket | Phase 2 | Pending |
| Job queue with Celery workers | Phase 2 | Pending |
| Partial failure handling | Phase 2 | Pending |
| Result aggregation and statistics | Phase 2 | Pending |
| RDKit.js molecule rendering | Phase 1 | Pending |
| Interactive validation results dashboard | Phase 1 | Pending |
| Batch upload with drag-and-drop | Phase 3 | Pending |
| Export (CSV, Excel, SDF, JSON) | Phase 3 | Pending |
| PDF report generation | Phase 3 | Pending |
| RESTful API with OpenAPI docs | Phase 1 | Pending |
| Rate limiting and API key auth | Phase 3 | Pending |
| Python client library | Phase 3 | Pending |
| DECIMER OCSR output validation | Phase 3 | Pending |
| COCONUT database lookup | Phase 3 | Pending |
| PubChem cross-reference | Phase 3 | Pending |
| ChEMBL bioactivity lookup | Phase 3 | Pending |
| Docker Compose for development | Phase 1 | Pending |
| Production-ready Docker images | Phase 1 | Pending |
| Self-hosting documentation | Phase 1 | Pending |

## Context

**Existing Documentation:**
- `docs/PRD.md` — Detailed product requirements
- `docs/PHASES.md` — 14-week development breakdown with checklists
- `docs/ARCHITECTURE.md` — System design, code structure, patterns
- `docs/CLAUDE_CODE_GUIDE.md` — Implementation snippets and patterns
- `docs/SKILLS.md` — Technical dependencies and requirements

**Planning Files:**
- `.planning/ROADMAP.md` — Phase roadmap with success criteria
- `.planning/STATE.md` — Project state and session continuity
- `.planning/research/` — Research findings (stack, features, architecture, pitfalls)

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
*Last updated: 2026-01-20 after roadmap creation*
