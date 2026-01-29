# ChemStructVal Project

## Overview
ChemStructVal is a comprehensive web-based chemical structure validation suite.
React frontend + FastAPI backend for validating, standardizing, and assessing ML-readiness of chemical structures.

## Documentation Location
All project documentation is in `docs/`:
- `docs/PHASES.md` - Development phases and tasks (START HERE for implementation)
- `docs/ARCHITECTURE.md` - System design and code structure
- `docs/PRD.md` - Product requirements and specifications
- `docs/CLAUDE_CODE_GUIDE.md` - Code snippets and implementation patterns
- `docs/SKILLS.md` - Technical requirements and dependencies

## Tech Stack
- **Frontend**: React 18, TypeScript, Vite, Tailwind CSS, shadcn/ui, RDKit.js
- **Backend**: Python 3.11+, FastAPI, RDKit, MolVS, chembl-structure-pipeline
- **Database**: PostgreSQL with asyncpg
- **Queue**: Redis + Celery for batch processing
- **Testing**: pytest (backend), Vitest (frontend)

## Commands
```bash
# Backend
cd backend
poetry install              # Install dependencies
poetry run uvicorn app.main:app --reload  # Start dev server
poetry run pytest           # Run tests

# Frontend
cd frontend
npm install                 # Install dependencies
npm run dev                 # Start dev server
npm test                    # Run tests

# Docker
docker-compose up -d        # Start all services
docker-compose logs -f      # View logs
```

## Code Style
- Python: Black formatter, isort for imports, type hints required
- TypeScript: ESLint + Prettier, strict mode enabled
- Always include docstrings for public functions
- Write tests for new features

## Current Phase
Phase 1 - Foundation (Weeks 1-4)
Currently working on: [UPDATE AS YOU PROGRESS]

## Important Patterns
- See `docs/ARCHITECTURE.md` for ValidationEngine class structure
- All validation checks inherit from BaseCheck abstract class
- Use Pydantic models for request/response schemas
- Follow the file structure in docs/ARCHITECTURE.md