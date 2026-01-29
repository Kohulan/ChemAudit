# External Integrations

**Analysis Date:** 2026-01-20

## APIs & External Services

**Cheminformatics:**
- RDKit - Core chemistry operations (local library, no external API)
  - SDK/Client: `rdkit` Python package, `@rdkit/rdkit` WebAssembly
  - Auth: None (local)
- MolVS - Molecule standardization (local library)
  - SDK/Client: `molvs` Python package
  - Auth: None (local)
- ChEMBL Structure Pipeline - Standardization (local library)
  - SDK/Client: `chembl-structure-pipeline` Python package
  - Auth: None (local)

**Planned Integrations (Phase 3):**
- PubChem - Compound lookup by ID
  - SDK/Client: REST API
  - Auth: None (public API, rate limited)
- ChEMBL - Compound lookup and validation comparison
  - SDK/Client: REST API
  - Auth: None (public API, rate limited)
- DECIMER - OCSR output validation
  - SDK/Client: Custom integration
  - Auth: TBD
- COCONUT - Natural products database validation
  - SDK/Client: Custom integration
  - Auth: TBD

## Data Storage

**Databases:**
- PostgreSQL 15+
  - Connection: `DATABASE_URL` env var
  - Client: SQLAlchemy 2.0 + asyncpg
  - Schema: Jobs, results, alerts (see `docs/ARCHITECTURE.md`)
  - ORM: SQLAlchemy with Alembic migrations

**File Storage:**
- Local filesystem for development
- Uploaded files: Temporary storage (7 day retention)
- Generated reports: Temporary storage (30 day retention)
- Production: Consider object storage (S3-compatible)

**Caching:**
- Redis 7+
  - Connection: `REDIS_URL` env var
  - Client: `redis` Python package
  - Usage: Validation result caching, job progress, rate limiting

## Authentication & Identity

**Auth Provider:**
- Phase 1: None (anonymous access)
- Phase 3: Optional API key authentication
  - Implementation: Custom API key generation and validation
  - Storage: PostgreSQL
  - Rate limit tiers: Anonymous (10/min), Registered (60/min), API Key (300/min)

**Planned (Future):**
- OAuth integration (GitHub, Google) for user accounts
- JWT tokens for session management

## Monitoring & Observability

**Error Tracking:**
- Planned: Sentry integration
- Implementation: FastAPI middleware, React error boundaries

**Logs:**
- Structured JSON logging (Python `logging` module)
- Log levels: DEBUG, INFO, WARNING, ERROR
- Correlation IDs for request tracing
- Configuration: `backend/app/core/logging.py`

**Metrics (Planned):**
- Prometheus metrics endpoint
- Grafana dashboards
- Key metrics: Request latency, validation throughput, job queue depth

## CI/CD & Deployment

**Hosting:**
- Docker containers
- Production: Kubernetes or Docker Swarm (flexible)
- CDN: Static frontend assets

**CI Pipeline:**
- GitHub Actions
- Workflows: `ci.yml`, `deploy.yml`, `release.yml`
- Jobs: Backend tests, frontend tests, linting, build, deploy

**Container Registry:**
- Docker Hub or GitHub Container Registry (TBD)

## Environment Configuration

**Required env vars (Backend):**
- `DATABASE_URL` - PostgreSQL connection (required)
- `REDIS_URL` - Redis connection (required)
- `CELERY_BROKER_URL` - Celery broker URL (required for batch)
- `CORS_ORIGINS` - Allowed frontend origins (required)
- `APP_NAME` - Application name (default: ChemVault)
- `DEBUG` - Debug mode (default: false)
- `MAX_BATCH_SIZE` - Max molecules per batch (default: 100000)
- `MAX_FILE_SIZE_MB` - Max upload size (default: 100)

**Required env vars (Frontend):**
- `VITE_API_URL` - Backend API base URL

**Secrets location:**
- Development: `.env` files (git-ignored)
- Production: Environment variables or secrets manager

## Webhooks & Callbacks

**Incoming:**
- None currently planned

**Outgoing:**
- Planned (Phase 2): Email notifications on batch job completion
  - Provider: TBD (SendGrid, SES, or SMTP)

## WebSocket Connections

**Real-time Updates:**
- Endpoint: `ws://api.chemvault.com/ws/jobs/{job_id}`
- Purpose: Batch job progress updates
- Protocol: JSON messages (progress, completed, error)
- Implementation: FastAPI WebSocket support
- Client: React hook `useWebSocket.ts`

## File Format Support

**Input Formats:**
- SMILES (canonical, isomeric, extended/CXSMILES)
- MOL/SDF (V2000, V3000)
- InChI / InChIKey
- CSV with SMILES column
- Excel (.xlsx) with SMILES column (Phase 1.1)
- Compressed files (.zip, .gz) (Phase 2)

**Output Formats:**
- JSON (API responses, validation results)
- CSV (batch results export)
- Excel (.xlsx) with formatting (Phase 3)
- SDF (standardized structures)
- PDF (summary reports, Phase 3)
- PNG/SVG (structure images)

## Rate Limiting

**Implementation:**
- slowapi library with Redis backend
- Per-IP rate limiting
- Configurable limits by tier

**Default Limits:**
| Tier | Single validations/min | Batch jobs/hour |
|------|------------------------|-----------------|
| Anonymous | 10 | 1 |
| Registered | 60 | 10 |
| API Key | 300 | 50 |

---

*Integration audit: 2026-01-20*
