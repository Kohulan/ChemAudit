# Directory Structure

**Analysis Date:** 2026-01-20

## Project Layout

```
chemvault/
├── backend/                    # FastAPI Python backend
│   ├── app/
│   │   ├── __init__.py
│   │   ├── main.py            # FastAPI app entry
│   │   ├── api/
│   │   │   ├── routes/        # HTTP endpoints
│   │   │   │   ├── validation.py
│   │   │   │   ├── batch.py
│   │   │   │   ├── jobs.py
│   │   │   │   ├── alerts.py
│   │   │   │   └── health.py
│   │   │   └── dependencies.py
│   │   ├── core/
│   │   │   ├── config.py      # pydantic-settings
│   │   │   ├── logging.py
│   │   │   └── exceptions.py
│   │   ├── services/
│   │   │   ├── validation/
│   │   │   │   ├── engine.py  # ValidationEngine
│   │   │   │   ├── checks/    # Check implementations
│   │   │   │   └── scorers/   # NP/ML scorers
│   │   │   ├── standardization/
│   │   │   │   ├── pipeline.py
│   │   │   │   └── steps/
│   │   │   ├── parser/
│   │   │   │   └── molecule_parser.py
│   │   │   └── export/
│   │   ├── models/            # SQLAlchemy models
│   │   ├── schemas/           # Pydantic schemas
│   │   ├── tasks/             # Celery tasks
│   │   └── db/                # Database layer
│   ├── data/
│   │   └── alerts/            # Alert pattern JSON files
│   ├── tests/
│   │   ├── conftest.py
│   │   ├── test_validation/
│   │   ├── test_api/
│   │   └── fixtures/
│   ├── alembic/               # DB migrations
│   ├── pyproject.toml
│   └── Dockerfile
│
├── frontend/                   # React TypeScript frontend
│   ├── src/
│   │   ├── main.tsx
│   │   ├── App.tsx
│   │   ├── components/
│   │   │   ├── ui/            # shadcn/ui components
│   │   │   ├── layout/        # Header, Footer, Layout
│   │   │   ├── molecules/     # MoleculeViewer, StructureInput
│   │   │   └── validation/    # ValidationResults, IssueCard
│   │   ├── pages/             # Route components
│   │   ├── hooks/             # React Query hooks
│   │   ├── services/          # API client
│   │   ├── types/             # TypeScript types
│   │   ├── utils/
│   │   └── styles/
│   ├── public/
│   ├── package.json
│   ├── tsconfig.json
│   ├── vite.config.ts
│   ├── tailwind.config.js
│   └── Dockerfile
│
├── docs/                       # Project documentation
│   ├── ARCHITECTURE.md
│   ├── PHASES.md
│   ├── PRD.md
│   ├── SKILLS.md
│   ├── PROJECT.md
│   └── CLAUDE_CODE_GUIDE.md
│
├── .github/
│   └── workflows/
│       ├── ci.yml
│       ├── deploy.yml
│       └── release.yml
│
├── docker-compose.yml          # Development
├── docker-compose.prod.yml     # Production
├── CLAUDE.md                   # Claude Code instructions
└── README.md
```

## Key Locations

**Configuration:**
- `backend/app/core/config.py` - Backend settings (pydantic-settings)
- `frontend/vite.config.ts` - Vite build config
- `frontend/tailwind.config.js` - Tailwind CSS config
- `docker-compose.yml` - Development services

**Entry Points:**
- `backend/app/main.py` - FastAPI application
- `frontend/src/main.tsx` - React application
- `backend/app/tasks/celery_app.py` - Celery worker

**Business Logic:**
- `backend/app/services/validation/engine.py` - Core validation
- `backend/app/services/validation/checks/` - Individual checks
- `backend/app/services/standardization/pipeline.py` - Standardization

**API Layer:**
- `backend/app/api/routes/` - HTTP endpoint handlers
- `backend/app/schemas/` - Request/response models
- `frontend/src/services/api.ts` - Frontend API client

**Database:**
- `backend/app/models/` - SQLAlchemy ORM models
- `backend/app/db/repositories/` - Data access layer
- `backend/alembic/` - Database migrations

**Tests:**
- `backend/tests/` - pytest tests
- `frontend/` - Vitest tests (co-located or in `__tests__/`)

## Naming Conventions

**Files:**
- Python: `snake_case.py`
- TypeScript: `PascalCase.tsx` (components), `camelCase.ts` (utils)
- Tests: `test_*.py` (pytest), `*.test.ts` (Vitest)

**Classes:**
- Python: `PascalCase` (e.g., `ValidationEngine`, `BaseCheck`)
- TypeScript: `PascalCase` (e.g., `MoleculeViewer`)

**Functions:**
- Python: `snake_case` (e.g., `validate_molecule`)
- TypeScript: `camelCase` (e.g., `useValidation`)

**Variables:**
- Python: `snake_case`, `UPPER_CASE` for constants
- TypeScript: `camelCase`, `UPPER_CASE` for constants

**API Endpoints:**
- RESTful: `/api/v1/{resource}` (e.g., `/api/v1/validate`)
- Kebab-case for multi-word resources

---

*Structure analysis: 2026-01-20*
