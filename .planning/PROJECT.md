# ChemVault

## What This Is

A comprehensive web-based chemical structure validation platform that consolidates disparate validation tools into a unified, visual interface. It serves general chemists who currently validate structures manually with fragmented tools, providing accessible, professional-grade quality assessment without requiring programming expertise.

**Current State:** v1.0 shipped (2026-01-29) — Production-ready with all planned features

## Core Value

Any chemist can validate chemical structures comprehensively — from basic parsing to ML-readiness assessment — through an intuitive web interface, with batch processing for large datasets.

## Requirements

### Validated

**v1.0 (Shipped 2026-01-29):**

- ✓ Single molecule validation with 11 comprehensive checks — v1.0
- ✓ SMILES, InChI, MOL/SDF format support with auto-detection — v1.0
- ✓ Parsability, sanitization, valence, aromaticity checks — v1.0
- ✓ Stereochemistry validation (undefined centers, conflicting stereo) — v1.0
- ✓ Representation consistency (SMILES/InChI roundtrip) — v1.0
- ✓ Overall quality score calculation (0-100) — v1.0
- ✓ PAINS A/B/C pattern screening (480+ patterns) — v1.0
- ✓ BRENK alerts (105 patterns) — v1.0
- ✓ NIH, ZINC, ChEMBL alert sets (2600+ patterns total) — v1.0
- ✓ Matched atom highlighting in structure viewer — v1.0
- ✓ ML-readiness scoring (descriptor calculability, fingerprints, size) — v1.0
- ✓ NP-likeness scoring (natural product likelihood) — v1.0
- ✓ Scaffold analysis (Murcko extraction) — v1.0
- ✓ ChEMBL-compatible standardization pipeline — v1.0
- ✓ Salt stripping, solvent removal, normalization — v1.0
- ✓ Tautomer canonicalization — v1.0
- ✓ Parent structure extraction — v1.0
- ✓ Before/after comparison view — v1.0
- ✓ File upload (SDF, CSV) up to 10K molecules — v1.0
- ✓ Progress tracking with WebSocket updates — v1.0
- ✓ Job queue with Celery workers — v1.0
- ✓ Partial failure handling — v1.0
- ✓ Result aggregation and statistics — v1.0
- ✓ RDKit.js molecule rendering with atom highlighting — v1.0
- ✓ Interactive validation results dashboard — v1.0
- ✓ Batch upload with drag-and-drop — v1.0
- ✓ Export (CSV, Excel, SDF, JSON) — v1.0
- ✓ PDF report generation — v1.0
- ✓ RESTful API with OpenAPI documentation — v1.0
- ✓ Rate limiting and optional API key auth — v1.0
- ✓ Python client library (chemvault-client) — v1.0
- ✓ DECIMER OCSR output validation — v1.0
- ✓ COCONUT natural products database lookup — v1.0
- ✓ PubChem compound cross-reference — v1.0
- ✓ ChEMBL bioactivity data lookup — v1.0
- ✓ Docker Compose for development — v1.0
- ✓ Production-ready Docker images — v1.0
- ✓ Self-hosting documentation — v1.0
- ✓ Performance <3s single validation (achieved) — v1.0
- ✓ 10K molecules batch processing (achieved) — v1.0
- ✓ 99%+ reference accuracy (achieved) — v1.0
- ✓ 404 Not Found page — v1.0
- ✓ React Error Boundary — v1.0
- ✓ Keyboard shortcuts (Ctrl+Enter) — v1.0
- ✓ Recent molecules history — v1.0
- ✓ URL sharing with encoded SMILES — v1.0
- ✓ Loading skeletons — v1.0
- ✓ Environment-based API URL configuration — v1.0
- ✓ Prometheus monitoring + Grafana dashboards — v1.0

### Active

(None — run `/gsd:new-milestone` to define v1.1 requirements)

### Out of Scope

- Custom ML model training — beyond v1 scope, focus on validation
- Reaction validation — different problem domain
- Protein-ligand analysis — separate tool category
- Commercial database integration — licensing complexity
- Real-time collaborative editing — unnecessary for validation workflow
- Mobile app — web-first, responsive design sufficient
- Frontend UI for API key management — backend-only by design (API keys managed via CLI/API)

## Context

**Codebase:**
- ~27K lines of code (14K Python backend, 13K TypeScript frontend)
- Tech stack: React 18, TypeScript, FastAPI, PostgreSQL, Redis, Celery, RDKit
- Fully Dockerized with production deployment configuration

**Documentation:**
- `docs/` — User guide, deployment guide, API reference, troubleshooting
- `.planning/` — Project planning, phase history, milestone archives

**Target Users:**
General chemists (database curators, ML researchers, medicinal chemists, natural product scientists) who need accessible structure validation without CLI expertise.

## Constraints

- **Tech Stack**: React 18 + FastAPI + PostgreSQL + Redis + RDKit (documented, non-negotiable)
- **Performance**: <3 seconds for single validation, 10K molecules per batch in <5 min
- **Accuracy**: 99%+ agreement with reference implementations (RDKit, ChEMBL pipeline)
- **Browser Support**: Modern browsers (Chrome, Firefox, Safari, Edge)
- **Deployment**: Supports both hosted web service and Docker self-hosting

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| RDKit.js for frontend rendering | Consistent with backend RDKit, no server roundtrip for display | ✓ Good — fast, consistent |
| Celery + Redis for batch processing | Proven stack for async job processing at scale | ✓ Good — handles 10K molecules |
| ChEMBL structure pipeline for standardization | Industry standard, well-tested | ✓ Good — reliable |
| PostgreSQL for job/result storage | Reliable, good JSON support, familiar | ✓ Good |
| Tautomer canonicalization OFF by default | Preserves E/Z stereochemistry | ✓ Good — prevents data loss |
| Stereochemistry tracking in standardization | Users need to know if stereo changed | ✓ Good — explicit warnings |
| 100-molecule chunk size for batch | Balance between progress granularity and overhead | ✓ Good |
| InChIKey as cache key | Canonical identifier, handles SMILES variations | ✓ Good |
| SHA256 for cache hashing | Security improvement over MD5 | ✓ Good |

---
*Last updated: 2026-01-29 after v1.0 milestone*
