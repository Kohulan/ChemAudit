# ChemStructVal: Claude Code Implementation Guide

This guide provides specific instructions for implementing ChemStructVal using Claude Code, optimized for efficient agentic development.

---

## Project Initialization Commands

### 1. Create Project Structure

```bash
# Create root directory
mkdir -p chemstructval
cd chemstructval

# Initialize git
git init
echo "# ChemStructVal" > README.md

# Create directory structure
mkdir -p backend/app/{api/routes,core,services/{validation/checks,standardization/steps,parser,export},models,schemas,tasks,db/repositories}
mkdir -p backend/{tests/{test_validation,test_api,fixtures},data/alerts,scripts,alembic/versions}
mkdir -p frontend/src/{components/{ui,layout,molecules,validation},pages,hooks,services,types,utils,styles}
mkdir -p docs/{api,user-guide,deployment}

# Create initial files
touch backend/app/__init__.py
touch backend/app/main.py
touch backend/app/api/__init__.py
touch backend/app/core/__init__.py
touch backend/app/services/__init__.py
```

### 2. Backend Setup

```bash
cd backend

# Create pyproject.toml
cat > pyproject.toml << 'EOF'
[project]
name = "chemstructval"
version = "0.1.0"
description = "Chemical Structure Validation Suite"
requires-python = ">=3.10"
dependencies = [
    "fastapi>=0.100.0",
    "uvicorn[standard]>=0.22.0",
    "python-multipart>=0.0.6",
    "pydantic>=2.0.0",
    "pydantic-settings>=2.0.0",
    "rdkit>=2023.9.0",
    "molvs>=0.1.1",
    "chembl-structure-pipeline>=1.2.0",
    "asyncpg>=0.28.0",
    "redis>=4.6.0",
    "sqlalchemy>=2.0.0",
    "alembic>=1.11.0",
    "celery>=5.3.0",
    "pandas>=2.0.0",
    "aiofiles>=23.0.0",
    "python-jose>=3.3.0",
    "slowapi>=0.1.8",
]

[project.optional-dependencies]
dev = [
    "pytest>=7.4.0",
    "pytest-asyncio>=0.21.0",
    "pytest-cov>=4.1.0",
    "httpx>=0.24.0",
    "black>=23.0.0",
    "isort>=5.12.0",
    "mypy>=1.4.0",
]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"
EOF

# Install dependencies (use pip or uv)
pip install -e ".[dev]"
```

### 3. Frontend Setup

```bash
cd ../frontend

# Create package.json
cat > package.json << 'EOF'
{
  "name": "chemstructval-frontend",
  "version": "0.1.0",
  "type": "module",
  "scripts": {
    "dev": "vite",
    "build": "tsc && vite build",
    "preview": "vite preview",
    "lint": "eslint . --ext ts,tsx",
    "test": "vitest"
  },
  "dependencies": {
    "react": "^18.2.0",
    "react-dom": "^18.2.0",
    "react-router-dom": "^6.14.0",
    "@tanstack/react-query": "^4.32.0",
    "axios": "^1.4.0",
    "@rdkit/rdkit": "^2023.9.1",
    "lucide-react": "^0.263.0",
    "recharts": "^2.7.0",
    "react-dropzone": "^14.2.0",
    "clsx": "^2.0.0",
    "tailwind-merge": "^1.14.0"
  },
  "devDependencies": {
    "@types/react": "^18.2.0",
    "@types/react-dom": "^18.2.0",
    "@vitejs/plugin-react": "^4.0.0",
    "typescript": "^5.1.0",
    "vite": "^4.4.0",
    "tailwindcss": "^3.3.0",
    "postcss": "^8.4.0",
    "autoprefixer": "^10.4.0",
    "vitest": "^0.34.0",
    "@testing-library/react": "^14.0.0"
  }
}
EOF

npm install
```

---

## Core Implementation Files

### Backend: Main Application

```python
# backend/app/main.py
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager

from app.api.routes import validation, batch, jobs, alerts, health
from app.core.config import settings
from app.core.logging import setup_logging

@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup
    setup_logging()
    yield
    # Shutdown

app = FastAPI(
    title="ChemStructVal API",
    description="Chemical Structure Validation Suite",
    version="0.1.0",
    lifespan=lifespan,
)

# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Routes
app.include_router(health.router, prefix="/api/v1", tags=["health"])
app.include_router(validation.router, prefix="/api/v1", tags=["validation"])
app.include_router(batch.router, prefix="/api/v1", tags=["batch"])
app.include_router(jobs.router, prefix="/api/v1", tags=["jobs"])
app.include_router(alerts.router, prefix="/api/v1", tags=["alerts"])
```

### Backend: Configuration

```python
# backend/app/core/config.py
from pydantic_settings import BaseSettings
from typing import List

class Settings(BaseSettings):
    # App
    APP_NAME: str = "ChemStructVal"
    DEBUG: bool = False
    
    # Database
    DATABASE_URL: str = "postgresql+asyncpg://user:pass@localhost:5432/chemstructval"
    
    # Redis
    REDIS_URL: str = "redis://localhost:6379"
    
    # Celery
    CELERY_BROKER_URL: str = "redis://localhost:6379/0"
    
    # CORS
    CORS_ORIGINS: List[str] = ["http://localhost:3000"]
    
    # Limits
    MAX_BATCH_SIZE: int = 100000
    MAX_FILE_SIZE_MB: int = 100
    
    class Config:
        env_file = ".env"

settings = Settings()
```

### Backend: Validation Schemas

```python
# backend/app/schemas/validation.py
from pydantic import BaseModel, Field, validator
from typing import List, Optional, Dict, Any
from enum import Enum

class Severity(str, Enum):
    CRITICAL = "critical"
    ERROR = "error"
    WARNING = "warning"
    INFO = "info"
    PASS = "pass"

class CheckResult(BaseModel):
    check_name: str
    passed: bool
    severity: Severity
    message: str
    affected_atoms: List[int] = []
    details: Dict[str, Any] = {}

class MoleculeInfo(BaseModel):
    input_smiles: str
    canonical_smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    num_atoms: Optional[int] = None
    num_bonds: Optional[int] = None

class AlertMatch(BaseModel):
    alert_name: str
    alert_set: str
    matched_atoms: List[int]
    severity: Severity = Severity.WARNING

class ValidationRequest(BaseModel):
    molecule: str = Field(..., min_length=1, max_length=10000)
    format: str = Field(default="auto", pattern="^(auto|smiles|inchi|mol)$")
    checks: List[str] = Field(default=["all"])
    options: Optional[Dict[str, Any]] = None
    
    @validator('molecule')
    def sanitize_molecule(cls, v):
        # Basic sanitization
        dangerous = ['<', '>', '&', ';', '|', '$', '`']
        if any(c in v for c in dangerous):
            raise ValueError("Invalid characters in molecule string")
        return v

class ValidationResponse(BaseModel):
    status: str = "completed"
    molecule_info: MoleculeInfo
    overall_score: int = Field(ge=0, le=100)
    ml_readiness_score: int = Field(ge=0, le=100)
    np_likeness_score: Optional[float] = None
    issues: List[CheckResult]
    alerts: List[AlertMatch] = []
    standardized_smiles: Optional[str] = None
    execution_time_ms: int
```

### Backend: Validation Engine (Core)

```python
# backend/app/services/validation/engine.py
from abc import ABC, abstractmethod
from typing import List, Optional, Dict, Type
from rdkit import Chem
import time

from app.schemas.validation import (
    CheckResult, ValidationResponse, MoleculeInfo, Severity
)

class BaseCheck(ABC):
    """Abstract base class for validation checks."""
    
    name: str
    description: str
    category: str
    
    @abstractmethod
    def run(self, mol: Chem.Mol) -> CheckResult:
        pass

class ValidationEngine:
    """Main validation orchestrator."""
    
    def __init__(self):
        self._checks: Dict[str, BaseCheck] = {}
        self._categories: Dict[str, List[str]] = {}
    
    def register_check(self, check: BaseCheck) -> None:
        self._checks[check.name] = check
        cat = check.category
        if cat not in self._categories:
            self._categories[cat] = []
        self._categories[cat].append(check.name)
    
    def get_check(self, name: str) -> Optional[BaseCheck]:
        return self._checks.get(name)
    
    def list_checks(self) -> Dict[str, List[str]]:
        return self._categories.copy()
    
    def validate(
        self,
        mol: Chem.Mol,
        checks: Optional[List[str]] = None,
        options: Optional[Dict] = None
    ) -> ValidationResponse:
        """Run validation on a molecule."""
        start_time = time.time()
        
        # Determine which checks to run
        if checks is None or "all" in checks:
            checks_to_run = list(self._checks.keys())
        else:
            checks_to_run = [c for c in checks if c in self._checks]
        
        # Run checks
        results: List[CheckResult] = []
        for check_name in checks_to_run:
            check = self._checks[check_name]
            try:
                result = check.run(mol)
                results.append(result)
            except Exception as e:
                results.append(CheckResult(
                    check_name=check_name,
                    passed=False,
                    severity=Severity.ERROR,
                    message=f"Check failed: {str(e)}",
                    affected_atoms=[]
                ))
        
        # Calculate score
        overall_score = self._calculate_score(results)
        
        # Extract molecule info
        mol_info = self._extract_mol_info(mol)
        
        execution_time = int((time.time() - start_time) * 1000)
        
        return ValidationResponse(
            molecule_info=mol_info,
            overall_score=overall_score,
            ml_readiness_score=self._calculate_ml_score(results),
            issues=[r for r in results if not r.passed],
            execution_time_ms=execution_time
        )
    
    def _calculate_score(self, results: List[CheckResult]) -> int:
        score = 100
        severity_impacts = {
            Severity.CRITICAL: -50,
            Severity.ERROR: -20,
            Severity.WARNING: -5,
            Severity.INFO: 0,
            Severity.PASS: 0,
        }
        for r in results:
            if not r.passed:
                score += severity_impacts.get(r.severity, 0)
        return max(0, min(100, score))
    
    def _calculate_ml_score(self, results: List[CheckResult]) -> int:
        # Simple ML score based on basic checks
        ml_relevant = ["parsability", "sanitization", "valence"]
        ml_results = [r for r in results if r.check_name in ml_relevant]
        if not ml_results:
            return 100
        passed = sum(1 for r in ml_results if r.passed)
        return int(100 * passed / len(ml_results))
    
    def _extract_mol_info(self, mol: Chem.Mol) -> MoleculeInfo:
        from rdkit.Chem import Descriptors, inchi
        
        try:
            canonical = Chem.MolToSmiles(mol)
            mol_inchi = inchi.MolToInchi(mol)
            mol_inchikey = inchi.MolToInchiKey(mol)
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            mw = Descriptors.MolWt(mol)
        except:
            canonical = None
            mol_inchi = None
            mol_inchikey = None
            formula = None
            mw = None
        
        return MoleculeInfo(
            input_smiles=Chem.MolToSmiles(mol) if mol else "",
            canonical_smiles=canonical,
            inchi=mol_inchi,
            inchikey=mol_inchikey,
            molecular_formula=formula,
            molecular_weight=mw,
            num_atoms=mol.GetNumAtoms() if mol else 0,
            num_bonds=mol.GetNumBonds() if mol else 0
        )


# Create global engine instance
validation_engine = ValidationEngine()
```

### Backend: Basic Checks Implementation

```python
# backend/app/services/validation/checks/basic.py
from rdkit import Chem
from app.services.validation.engine import BaseCheck
from app.schemas.validation import CheckResult, Severity

class ParsabilityCheck(BaseCheck):
    name = "parsability"
    description = "Check if molecule can be parsed"
    category = "basic"
    
    def run(self, mol: Chem.Mol) -> CheckResult:
        # If we got here with a mol, it's parsable
        return CheckResult(
            check_name=self.name,
            passed=True,
            severity=Severity.PASS,
            message="Molecule parsed successfully"
        )

class SanitizationCheck(BaseCheck):
    name = "sanitization"
    description = "Check if molecule passes RDKit sanitization"
    category = "basic"
    
    def run(self, mol: Chem.Mol) -> CheckResult:
        try:
            # Try to sanitize a copy
            test_mol = Chem.RWMol(mol)
            Chem.SanitizeMol(test_mol)
            return CheckResult(
                check_name=self.name,
                passed=True,
                severity=Severity.PASS,
                message="Molecule sanitizes correctly"
            )
        except Exception as e:
            return CheckResult(
                check_name=self.name,
                passed=False,
                severity=Severity.ERROR,
                message=f"Sanitization failed: {str(e)}"
            )

class ValenceCheck(BaseCheck):
    name = "valence"
    description = "Check for valence errors"
    category = "basic"
    
    def run(self, mol: Chem.Mol) -> CheckResult:
        problems = Chem.DetectChemistryProblems(mol)
        valence_problems = [p for p in problems if "valence" in p.GetType().lower()]
        
        if not valence_problems:
            return CheckResult(
                check_name=self.name,
                passed=True,
                severity=Severity.PASS,
                message="No valence errors"
            )
        
        affected = []
        for p in valence_problems:
            if hasattr(p, 'GetAtomIdx'):
                affected.append(p.GetAtomIdx())
        
        return CheckResult(
            check_name=self.name,
            passed=False,
            severity=Severity.CRITICAL,
            message=f"Found {len(valence_problems)} valence error(s)",
            affected_atoms=affected
        )

class ConnectivityCheck(BaseCheck):
    name = "connectivity"
    description = "Check for disconnected fragments"
    category = "basic"
    
    def run(self, mol: Chem.Mol) -> CheckResult:
        frags = Chem.GetMolFrags(mol)
        
        if len(frags) == 1:
            return CheckResult(
                check_name=self.name,
                passed=True,
                severity=Severity.PASS,
                message="Single connected component"
            )
        
        return CheckResult(
            check_name=self.name,
            passed=False,
            severity=Severity.WARNING,
            message=f"Found {len(frags)} disconnected fragments",
            details={"fragment_count": len(frags)}
        )
```

### Backend: Validation Route

```python
# backend/app/api/routes/validation.py
from fastapi import APIRouter, HTTPException
from rdkit import Chem

from app.schemas.validation import ValidationRequest, ValidationResponse
from app.services.validation.engine import validation_engine
from app.services.parser.molecule_parser import parse_molecule

router = APIRouter()

@router.post("/validate", response_model=ValidationResponse)
async def validate_molecule(request: ValidationRequest):
    """Validate a single molecule."""
    
    # Parse molecule
    try:
        mol = parse_molecule(request.molecule, request.format)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    
    if mol is None:
        raise HTTPException(
            status_code=400, 
            detail="Could not parse molecule"
        )
    
    # Run validation
    result = validation_engine.validate(
        mol,
        checks=request.checks,
        options=request.options
    )
    
    return result

@router.get("/checks")
async def list_checks():
    """List all available validation checks."""
    return validation_engine.list_checks()
```

### Frontend: API Service

```typescript
// frontend/src/services/api.ts
import axios from 'axios';

const API_BASE_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000/api/v1';

export const api = axios.create({
  baseURL: API_BASE_URL,
  timeout: 30000,
  headers: {
    'Content-Type': 'application/json',
  },
});

// Types
export interface ValidationRequest {
  molecule: string;
  format?: 'auto' | 'smiles' | 'inchi' | 'mol';
  checks?: string[];
  options?: Record<string, unknown>;
}

export interface CheckResult {
  check_name: string;
  passed: boolean;
  severity: 'critical' | 'error' | 'warning' | 'info' | 'pass';
  message: string;
  affected_atoms: number[];
  details?: Record<string, unknown>;
}

export interface MoleculeInfo {
  input_smiles: string;
  canonical_smiles?: string;
  inchi?: string;
  inchikey?: string;
  molecular_formula?: string;
  molecular_weight?: number;
}

export interface ValidationResponse {
  status: string;
  molecule_info: MoleculeInfo;
  overall_score: number;
  ml_readiness_score: number;
  np_likeness_score?: number;
  issues: CheckResult[];
  alerts: Array<{
    alert_name: string;
    alert_set: string;
    matched_atoms: number[];
  }>;
  standardized_smiles?: string;
  execution_time_ms: number;
}

// API functions
export const validationApi = {
  validate: async (request: ValidationRequest): Promise<ValidationResponse> => {
    const response = await api.post<ValidationResponse>('/validate', request);
    return response.data;
  },
  
  getChecks: async (): Promise<Record<string, string[]>> => {
    const response = await api.get<Record<string, string[]>>('/checks');
    return response.data;
  },
};
```

### Frontend: Validation Hook

```typescript
// frontend/src/hooks/useValidation.ts
import { useMutation, useQuery } from '@tanstack/react-query';
import { validationApi, ValidationRequest, ValidationResponse } from '../services/api';

export function useValidation() {
  const mutation = useMutation({
    mutationFn: (request: ValidationRequest) => validationApi.validate(request),
  });

  return {
    validate: mutation.mutate,
    validateAsync: mutation.mutateAsync,
    isLoading: mutation.isPending,
    error: mutation.error,
    result: mutation.data,
    reset: mutation.reset,
  };
}

export function useAvailableChecks() {
  return useQuery({
    queryKey: ['checks'],
    queryFn: validationApi.getChecks,
    staleTime: Infinity, // Checks don't change
  });
}
```

### Frontend: Main Validation Component

```tsx
// frontend/src/components/validation/ValidationResults.tsx
import React from 'react';
import { ValidationResponse } from '../../services/api';

interface Props {
  result: ValidationResponse;
}

export function ValidationResults({ result }: Props) {
  const scoreColor = result.overall_score >= 80 
    ? 'text-green-600' 
    : result.overall_score >= 50 
    ? 'text-yellow-600' 
    : 'text-red-600';

  return (
    <div className="space-y-6">
      {/* Score */}
      <div className="flex items-center justify-between p-4 bg-white rounded-lg shadow">
        <div>
          <h3 className="text-lg font-medium">Overall Score</h3>
          <p className="text-sm text-gray-500">
            Based on {result.issues.length + (result.overall_score === 100 ? 0 : 1)} checks
          </p>
        </div>
        <div className={`text-4xl font-bold ${scoreColor}`}>
          {result.overall_score}
        </div>
      </div>

      {/* Molecule Info */}
      <div className="p-4 bg-white rounded-lg shadow">
        <h3 className="text-lg font-medium mb-3">Molecule Info</h3>
        <dl className="grid grid-cols-2 gap-2 text-sm">
          <dt className="text-gray-500">Canonical SMILES</dt>
          <dd className="font-mono">{result.molecule_info.canonical_smiles}</dd>
          <dt className="text-gray-500">InChIKey</dt>
          <dd className="font-mono text-xs">{result.molecule_info.inchikey}</dd>
          <dt className="text-gray-500">Formula</dt>
          <dd>{result.molecule_info.molecular_formula}</dd>
          <dt className="text-gray-500">Molecular Weight</dt>
          <dd>{result.molecule_info.molecular_weight?.toFixed(2)}</dd>
        </dl>
      </div>

      {/* Issues */}
      {result.issues.length > 0 && (
        <div className="p-4 bg-white rounded-lg shadow">
          <h3 className="text-lg font-medium mb-3">Issues Found</h3>
          <ul className="space-y-2">
            {result.issues.map((issue, idx) => (
              <li 
                key={idx}
                className={`p-3 rounded ${
                  issue.severity === 'critical' ? 'bg-red-50 border-l-4 border-red-500' :
                  issue.severity === 'error' ? 'bg-orange-50 border-l-4 border-orange-500' :
                  issue.severity === 'warning' ? 'bg-yellow-50 border-l-4 border-yellow-500' :
                  'bg-blue-50 border-l-4 border-blue-500'
                }`}
              >
                <div className="flex justify-between">
                  <span className="font-medium">{issue.check_name}</span>
                  <span className="text-xs uppercase">{issue.severity}</span>
                </div>
                <p className="text-sm text-gray-600 mt-1">{issue.message}</p>
              </li>
            ))}
          </ul>
        </div>
      )}

      {/* Execution Time */}
      <p className="text-xs text-gray-400 text-right">
        Completed in {result.execution_time_ms}ms
      </p>
    </div>
  );
}
```

---

## Docker Configuration

```yaml
# docker-compose.yml
version: '3.8'

services:
  backend:
    build:
      context: ./backend
      dockerfile: Dockerfile
    ports:
      - "8000:8000"
    environment:
      - DATABASE_URL=postgresql+asyncpg://chemstructval:chemstructval@postgres:5432/chemstructval
      - REDIS_URL=redis://redis:6379
      - CELERY_BROKER_URL=redis://redis:6379/0
    volumes:
      - ./backend:/app
    depends_on:
      - postgres
      - redis
    command: uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload

  frontend:
    build:
      context: ./frontend
      dockerfile: Dockerfile.dev
    ports:
      - "3000:3000"
    environment:
      - VITE_API_URL=http://localhost:8000/api/v1
    volumes:
      - ./frontend:/app
      - /app/node_modules
    command: npm run dev -- --host

  postgres:
    image: postgres:15-alpine
    environment:
      - POSTGRES_USER=chemstructval
      - POSTGRES_PASSWORD=chemstructval
      - POSTGRES_DB=chemstructval
    volumes:
      - postgres_data:/var/lib/postgresql/data
    ports:
      - "5432:5432"

  redis:
    image: redis:7-alpine
    ports:
      - "6379:6379"

volumes:
  postgres_data:
```

```dockerfile
# backend/Dockerfile
FROM python:3.11-slim

WORKDIR /app

# Install system dependencies for RDKit
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY pyproject.toml .
RUN pip install --no-cache-dir -e .

COPY . .

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
```

---

## Testing Commands

```bash
# Run backend tests
cd backend
pytest -v --cov=app

# Run specific test file
pytest tests/test_validation/test_basic_checks.py -v

# Run frontend tests
cd frontend
npm test

# Run all with coverage
pytest --cov=app --cov-report=html
```

---

## Development Workflow

1. **Start services**: `docker-compose up -d`
2. **Backend dev**: Code changes auto-reload with uvicorn
3. **Frontend dev**: Code changes auto-reload with Vite
4. **API docs**: Visit `http://localhost:8000/docs`
5. **Run tests**: `pytest` (backend) / `npm test` (frontend)

---

## Common Claude Code Commands

```
# Add new validation check
Create a new file backend/app/services/validation/checks/stereo.py with UndefinedStereoCentersCheck class

# Add new API endpoint
Add a new route to backend/app/api/routes/validation.py for getting NP-likeness score

# Create React component
Create frontend/src/components/molecules/MoleculeViewer.tsx that renders RDKit.js molecules

# Fix import error
The import 'xyz' is failing in file abc.py, fix it

# Add tests
Write pytest tests for the ValenceCheck class in backend/tests/test_validation/test_basic.py
```

---

*This guide is optimized for Claude Code agentic development.*
