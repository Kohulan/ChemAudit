# ChemStructVal: Technical Architecture

## System Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              ChemStructVal                                   │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌─────────────────────────────┐    ┌──────────────────────────────────┐   │
│  │     React Frontend          │    │         FastAPI Backend          │   │
│  │  ┌─────────────────────┐   │    │  ┌────────────────────────────┐  │   │
│  │  │  Structure Input    │   │    │  │    Validation Engine       │  │   │
│  │  ├─────────────────────┤   │    │  │  ┌──────────────────────┐  │  │   │
│  │  │  RDKit.js Viewer    │   │◄───┼──┤  │  Molecule Parser     │  │  │   │
│  │  ├─────────────────────┤   │    │  │  ├──────────────────────┤  │  │   │
│  │  │  Results Dashboard  │   │    │  │  │  Basic Checks        │  │  │   │
│  │  ├─────────────────────┤   │    │  │  ├──────────────────────┤  │  │   │
│  │  │  Batch Manager      │   │    │  │  │  Stereo Checks       │  │  │   │
│  │  └─────────────────────┘   │    │  │  ├──────────────────────┤  │  │   │
│  └─────────────────────────────┘    │  │  │  Alert Engine        │  │  │   │
│                                      │  │  ├──────────────────────┤  │  │   │
│  ┌─────────────────────────────┐    │  │  │  NP/ML Scorers       │  │  │   │
│  │        Job Queue            │    │  │  ├──────────────────────┤  │  │   │
│  │  ┌─────────────────────┐   │    │  │  │  Standardizer        │  │  │   │
│  │  │    Celery Worker    │◄──┼────┤  │  └──────────────────────┘  │  │   │
│  │  └─────────────────────┘   │    │  └────────────────────────────┘  │   │
│  │  ┌─────────────────────┐   │    └──────────────────────────────────┘   │
│  │  │       Redis         │   │                                           │
│  │  └─────────────────────┘   │    ┌──────────────────────────────────┐   │
│  └─────────────────────────────┘    │        PostgreSQL               │   │
│                                      │  ┌────────────────────────────┐  │   │
│                                      │  │  Jobs, Results, Alerts     │  │   │
│                                      │  └────────────────────────────┘  │   │
│                                      └──────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Directory Structure

```
chemstructval/
├── backend/
│   ├── app/
│   │   ├── __init__.py
│   │   ├── main.py                 # FastAPI application entry
│   │   ├── api/
│   │   │   ├── __init__.py
│   │   │   ├── routes/
│   │   │   │   ├── __init__.py
│   │   │   │   ├── validation.py   # /api/v1/validate
│   │   │   │   ├── batch.py        # /api/v1/batch
│   │   │   │   ├── jobs.py         # /api/v1/jobs
│   │   │   │   ├── alerts.py       # /api/v1/alerts
│   │   │   │   └── health.py       # /api/v1/health
│   │   │   └── dependencies.py     # Shared dependencies
│   │   ├── core/
│   │   │   ├── __init__.py
│   │   │   ├── config.py           # Settings management
│   │   │   ├── logging.py          # Logging configuration
│   │   │   └── exceptions.py       # Custom exceptions
│   │   ├── services/
│   │   │   ├── __init__.py
│   │   │   ├── validation/
│   │   │   │   ├── __init__.py
│   │   │   │   ├── engine.py       # ValidationEngine
│   │   │   │   ├── checks/
│   │   │   │   │   ├── __init__.py
│   │   │   │   │   ├── base.py     # BaseCheck abstract class
│   │   │   │   │   ├── basic.py    # Valence, parsability, etc.
│   │   │   │   │   ├── stereo.py   # Stereochemistry checks
│   │   │   │   │   ├── alerts.py   # Structural alerts
│   │   │   │   │   ├── representation.py  # Round-trip checks
│   │   │   │   │   └── ml_readiness.py    # ML-readiness
│   │   │   │   └── scorers/
│   │   │   │       ├── __init__.py
│   │   │   │       ├── np_likeness.py
│   │   │   │       └── ml_readiness.py
│   │   │   ├── standardization/
│   │   │   │   ├── __init__.py
│   │   │   │   ├── pipeline.py     # StandardizationPipeline
│   │   │   │   └── steps/
│   │   │   │       ├── __init__.py
│   │   │   │       ├── salt_stripper.py
│   │   │   │       ├── normalizer.py
│   │   │   │       └── tautomer.py
│   │   │   ├── parser/
│   │   │   │   ├── __init__.py
│   │   │   │   └── molecule_parser.py
│   │   │   └── export/
│   │   │       ├── __init__.py
│   │   │       ├── csv_exporter.py
│   │   │       ├── excel_exporter.py
│   │   │       ├── sdf_exporter.py
│   │   │       └── pdf_report.py
│   │   ├── models/
│   │   │   ├── __init__.py
│   │   │   ├── job.py              # SQLAlchemy models
│   │   │   ├── alert.py
│   │   │   └── result.py
│   │   ├── schemas/
│   │   │   ├── __init__.py
│   │   │   ├── validation.py       # Pydantic schemas
│   │   │   ├── batch.py
│   │   │   ├── job.py
│   │   │   └── common.py
│   │   ├── tasks/
│   │   │   ├── __init__.py
│   │   │   ├── celery_app.py
│   │   │   └── batch_tasks.py
│   │   └── db/
│   │       ├── __init__.py
│   │       ├── database.py         # DB connection
│   │       └── repositories/
│   │           ├── __init__.py
│   │           ├── job_repo.py
│   │           └── result_repo.py
│   ├── data/
│   │   └── alerts/
│   │       ├── pains_a.json
│   │       ├── pains_b.json
│   │       ├── pains_c.json
│   │       ├── brenk.json
│   │       ├── nih.json
│   │       └── chembl/
│   │           └── *.json
│   ├── tests/
│   │   ├── __init__.py
│   │   ├── conftest.py
│   │   ├── test_validation/
│   │   ├── test_api/
│   │   └── fixtures/
│   │       └── molecules.py
│   ├── scripts/
│   │   └── load_alerts.py
│   ├── pyproject.toml
│   ├── Dockerfile
│   └── alembic/                    # DB migrations
│       ├── versions/
│       └── env.py
│
├── frontend/
│   ├── src/
│   │   ├── main.tsx
│   │   ├── App.tsx
│   │   ├── components/
│   │   │   ├── ui/                 # shadcn components
│   │   │   │   ├── button.tsx
│   │   │   │   ├── card.tsx
│   │   │   │   ├── dialog.tsx
│   │   │   │   └── ...
│   │   │   ├── layout/
│   │   │   │   ├── Header.tsx
│   │   │   │   ├── Footer.tsx
│   │   │   │   └── Layout.tsx
│   │   │   ├── molecules/
│   │   │   │   ├── MoleculeViewer.tsx
│   │   │   │   ├── StructureInput.tsx
│   │   │   │   └── MoleculeComparison.tsx
│   │   │   └── validation/
│   │   │       ├── ValidationResults.tsx
│   │   │       ├── IssueCard.tsx
│   │   │       ├── ScoreGauge.tsx
│   │   │       ├── AlertResults.tsx
│   │   │       └── BatchResults.tsx
│   │   ├── pages/
│   │   │   ├── Home.tsx
│   │   │   ├── SingleValidation.tsx
│   │   │   ├── BatchValidation.tsx
│   │   │   ├── JobHistory.tsx
│   │   │   ├── Documentation.tsx
│   │   │   └── ApiDocs.tsx
│   │   ├── hooks/
│   │   │   ├── useValidation.ts
│   │   │   ├── useBatchJob.ts
│   │   │   ├── useMolecule.ts
│   │   │   └── useWebSocket.ts
│   │   ├── services/
│   │   │   ├── api.ts              # Axios instance
│   │   │   ├── validation.ts       # Validation API calls
│   │   │   └── batch.ts            # Batch API calls
│   │   ├── types/
│   │   │   ├── validation.ts
│   │   │   ├── molecule.ts
│   │   │   └── api.ts
│   │   ├── utils/
│   │   │   ├── formatting.ts
│   │   │   └── chemistry.ts
│   │   └── styles/
│   │       └── globals.css
│   ├── public/
│   ├── package.json
│   ├── tsconfig.json
│   ├── vite.config.ts
│   ├── tailwind.config.js
│   └── Dockerfile
│
├── docker-compose.yml
├── docker-compose.prod.yml
├── .github/
│   └── workflows/
│       ├── ci.yml
│       ├── deploy.yml
│       └── release.yml
├── docs/
│   ├── api/
│   ├── user-guide/
│   └── deployment/
└── README.md
```

---

## Core Components

### 1. Validation Engine

The central component that orchestrates all validation checks.

```python
# backend/app/services/validation/engine.py

from abc import ABC, abstractmethod
from typing import List, Optional
from rdkit import Chem

class BaseCheck(ABC):
    """Abstract base class for all validation checks."""
    
    name: str
    description: str
    severity_default: Severity
    
    @abstractmethod
    def run(self, mol: Chem.Mol) -> CheckResult:
        """Execute the check on a molecule."""
        pass
    
    @abstractmethod
    def get_affected_atoms(self, mol: Chem.Mol) -> List[int]:
        """Return indices of atoms affected by issues."""
        pass


class ValidationEngine:
    """Main validation orchestrator."""
    
    def __init__(self):
        self._checks: Dict[str, BaseCheck] = {}
        self._load_default_checks()
    
    def register_check(self, check: BaseCheck):
        """Register a new check."""
        self._checks[check.name] = check
    
    def validate(
        self, 
        mol: Chem.Mol, 
        checks: Optional[List[str]] = None,
        options: Optional[ValidationOptions] = None
    ) -> ValidationResult:
        """
        Run validation checks on a molecule.
        
        Args:
            mol: RDKit molecule object
            checks: List of check names to run (None = all)
            options: Validation options
            
        Returns:
            ValidationResult with all check outcomes
        """
        results = []
        checks_to_run = checks or list(self._checks.keys())
        
        for check_name in checks_to_run:
            check = self._checks.get(check_name)
            if check:
                result = check.run(mol)
                results.append(result)
        
        return ValidationResult(
            checks=results,
            overall_score=self._calculate_score(results),
            molecule_info=self._extract_mol_info(mol)
        )
    
    def _calculate_score(self, results: List[CheckResult]) -> int:
        """Calculate overall score from individual check results."""
        score = 100
        for result in results:
            if not result.passed:
                score += result.severity.score_impact
        return max(0, score)
```

### 2. Check Implementation Pattern

```python
# backend/app/services/validation/checks/stereo.py

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

class UndefinedStereoCentersCheck(BaseCheck):
    """Check for undefined stereocenters."""
    
    name = "undefined_stereocenters"
    description = "Identifies chiral centers without defined stereochemistry"
    severity_default = Severity.WARNING
    
    def __init__(self, max_undefined: int = 0):
        self.max_undefined = max_undefined
    
    def run(self, mol: Chem.Mol) -> CheckResult:
        # Find all chiral centers
        chiral_centers = Chem.FindMolChiralCenters(
            mol, 
            includeUnassigned=True,
            useLegacyImplementation=False
        )
        
        # Count undefined (marked as '?')
        undefined = [c for c in chiral_centers if c[1] == '?']
        count = len(undefined)
        
        passed = count <= self.max_undefined
        
        return CheckResult(
            check_name=self.name,
            passed=passed,
            severity=Severity.WARNING if not passed else Severity.PASS,
            message=f"Found {count} undefined stereocenter(s)" if not passed else "All stereocenters defined",
            affected_atoms=[c[0] for c in undefined],
            details={
                "total_stereocenters": len(chiral_centers),
                "undefined_count": count,
                "undefined_atoms": [c[0] for c in undefined]
            }
        )
    
    def get_affected_atoms(self, mol: Chem.Mol) -> List[int]:
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        return [c[0] for c in chiral_centers if c[1] == '?']
```

### 3. Alert Engine

```python
# backend/app/services/validation/checks/alerts.py

import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import FilterCatalog

class AlertEngine:
    """Engine for structural alert screening."""
    
    BUILTIN_SETS = {
        "PAINS_A": FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A,
        "PAINS_B": FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B,
        "PAINS_C": FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C,
        "BRENK": FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK,
        "NIH": FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH,
        "ZINC": FilterCatalog.FilterCatalogParams.FilterCatalogs.ZINC,
    }
    
    def __init__(self, alert_sets: List[str] = None):
        self.alert_sets = alert_sets or ["PAINS_A", "PAINS_B", "PAINS_C"]
        self._catalog = self._build_catalog()
    
    def _build_catalog(self) -> FilterCatalog.FilterCatalog:
        """Build filter catalog from selected alert sets."""
        params = FilterCatalog.FilterCatalogParams()
        
        for alert_set in self.alert_sets:
            if alert_set in self.BUILTIN_SETS:
                params.AddCatalog(self.BUILTIN_SETS[alert_set])
        
        return FilterCatalog.FilterCatalog(params)
    
    def screen(self, mol: Chem.Mol) -> List[AlertMatch]:
        """Screen molecule against alert patterns."""
        matches = []
        
        for entry in self._catalog.GetMatches(mol):
            match = AlertMatch(
                alert_name=entry.GetDescription(),
                alert_set=entry.GetProp("Scope") if entry.HasProp("Scope") else "Unknown",
                matched_atoms=self._get_matched_atoms(mol, entry),
                severity=Severity.WARNING
            )
            matches.append(match)
        
        return matches
    
    def _get_matched_atoms(self, mol: Chem.Mol, entry) -> List[int]:
        """Get atom indices that matched the alert pattern."""
        matches = entry.GetFilterMatches(mol)
        atoms = set()
        for match in matches:
            for atom_pair in match.atomPairs:
                atoms.add(atom_pair[1])
        return list(atoms)
```

### 4. Standardization Pipeline

```python
# backend/app/services/standardization/pipeline.py

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from chembl_structure_pipeline import standardizer as chembl_std

class StandardizationPipeline:
    """Configurable molecular standardization pipeline."""
    
    def __init__(self, config: StandardizationConfig = None):
        self.config = config or StandardizationConfig()
        self._steps = self._build_steps()
    
    def _build_steps(self) -> List[StandardizationStep]:
        """Build pipeline steps based on configuration."""
        steps = []
        
        if self.config.remove_salts:
            steps.append(SaltStripper())
        
        if self.config.remove_solvents:
            steps.append(SolventRemover())
        
        if self.config.normalize:
            steps.append(Normalizer())
        
        if self.config.reionize:
            steps.append(Reionizer())
        
        if self.config.canonicalize_tautomer:
            steps.append(TautomerCanonicalizer())
        
        if self.config.remove_isotopes:
            steps.append(IsotopeRemover())
        
        return steps
    
    def standardize(self, mol: Chem.Mol) -> StandardizationResult:
        """
        Apply standardization pipeline to molecule.
        
        Returns:
            StandardizationResult with original and standardized forms
        """
        original_smiles = Chem.MolToSmiles(mol)
        current_mol = Chem.RWMol(mol)
        operations = []
        
        for step in self._steps:
            try:
                result = step.apply(current_mol)
                if result.modified:
                    operations.append(result.operation_name)
                    current_mol = result.mol
            except Exception as e:
                operations.append(f"{step.name}: ERROR - {str(e)}")
        
        standardized_smiles = Chem.MolToSmiles(current_mol)
        
        return StandardizationResult(
            original_smiles=original_smiles,
            standardized_smiles=standardized_smiles,
            standardized_mol=current_mol.GetMol(),
            operations_applied=operations,
            changed=(original_smiles != standardized_smiles)
        )


class SaltStripper(StandardizationStep):
    """Remove salt/counterion fragments."""
    
    name = "salt_stripper"
    
    def apply(self, mol: Chem.RWMol) -> StepResult:
        # Use ChEMBL salt stripping
        parent, excluded = chembl_std.get_parent_molblock(
            Chem.MolToMolBlock(mol)
        )
        
        if excluded:
            return StepResult(
                mol=Chem.MolFromMolBlock(parent),
                modified=True,
                operation_name="Removed salt components"
            )
        
        return StepResult(mol=mol, modified=False)
```

---

## Database Schema

```sql
-- Jobs table
CREATE TABLE validation_jobs (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    status VARCHAR(20) NOT NULL DEFAULT 'pending',
    input_format VARCHAR(20) NOT NULL,
    molecule_count INTEGER NOT NULL,
    options JSONB,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    started_at TIMESTAMP WITH TIME ZONE,
    completed_at TIMESTAMP WITH TIME ZONE,
    error_message TEXT,
    progress_current INTEGER DEFAULT 0,
    progress_total INTEGER DEFAULT 0
);

-- Results table (for batch jobs)
CREATE TABLE validation_results (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    job_id UUID REFERENCES validation_jobs(id) ON DELETE CASCADE,
    molecule_index INTEGER NOT NULL,
    input_smiles TEXT NOT NULL,
    canonical_smiles TEXT,
    inchi TEXT,
    inchikey VARCHAR(27),
    overall_score INTEGER,
    ml_readiness_score INTEGER,
    np_likeness_score FLOAT,
    issues JSONB,
    standardized_smiles TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    
    UNIQUE(job_id, molecule_index)
);

-- Alerts configuration
CREATE TABLE alert_sets (
    id SERIAL PRIMARY KEY,
    name VARCHAR(50) UNIQUE NOT NULL,
    description TEXT,
    pattern_count INTEGER,
    enabled BOOLEAN DEFAULT true,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);

CREATE TABLE alert_patterns (
    id SERIAL PRIMARY KEY,
    alert_set_id INTEGER REFERENCES alert_sets(id),
    pattern_name VARCHAR(100) NOT NULL,
    smarts TEXT NOT NULL,
    description TEXT,
    severity VARCHAR(20) DEFAULT 'warning'
);

-- Indexes
CREATE INDEX idx_results_job_id ON validation_results(job_id);
CREATE INDEX idx_results_inchikey ON validation_results(inchikey);
CREATE INDEX idx_jobs_status ON validation_jobs(status);
CREATE INDEX idx_jobs_created ON validation_jobs(created_at);
```

---

## API Endpoints Detail

### Single Molecule Validation

```yaml
POST /api/v1/validate
  description: Validate a single molecule
  request:
    body:
      molecule: string (SMILES, InChI, or MOL block)
      format: string (auto|smiles|inchi|mol)
      checks: array[string] (optional, default: all)
      options:
        standardize: boolean (default: false)
        alert_sets: array[string] (default: ["PAINS"])
        include_3d: boolean (default: false)
  response:
    200:
      status: "completed"
      molecule_info:
        input_smiles: string
        canonical_smiles: string
        inchi: string
        inchikey: string
        molecular_formula: string
        molecular_weight: float
      overall_score: integer (0-100)
      ml_readiness_score: integer (0-100)
      np_likeness_score: float (-5 to +5)
      issues: array[Issue]
      alerts: array[Alert]
      standardized_smiles: string (if standardize=true)
      execution_time_ms: integer
    400:
      error: "Invalid molecule"
      details: string
    422:
      error: "Validation error"
      details: array[ValidationError]
```

### Batch Processing

```yaml
POST /api/v1/validate/batch
  description: Start batch validation job
  request:
    body:
      file: File (multipart/form-data)
      format: string (sdf|csv)
      smiles_column: string (for CSV)
      options:
        checks: array[string]
        standardize: boolean
        alert_sets: array[string]
  response:
    202:
      job_id: string (UUID)
      status: "pending"
      molecule_count: integer
      estimated_time_seconds: integer

GET /api/v1/jobs/{job_id}
  description: Get job status and progress
  response:
    200:
      job_id: string
      status: string (pending|running|completed|failed|cancelled)
      progress:
        current: integer
        total: integer
        percentage: float
      started_at: datetime
      estimated_completion: datetime

GET /api/v1/jobs/{job_id}/results
  description: Get job results
  query:
    page: integer (default: 1)
    per_page: integer (default: 100, max: 1000)
    filter_severity: string (critical|error|warning)
    sort_by: string (score|molecule_index)
  response:
    200:
      job_id: string
      total_molecules: integer
      summary:
        passed: integer
        critical: integer
        errors: integer
        warnings: integer
      results: array[ValidationResult]
      pagination:
        page: integer
        per_page: integer
        total_pages: integer

DELETE /api/v1/jobs/{job_id}
  description: Cancel or delete job
  response:
    200:
      message: "Job cancelled"
```

---

## WebSocket Protocol

For real-time progress updates during batch processing:

```typescript
// Connect
ws://api.chemstructval.com/ws/jobs/{job_id}

// Server -> Client messages
interface ProgressMessage {
  type: "progress";
  job_id: string;
  current: number;
  total: number;
  current_molecule?: string;
  estimated_remaining_seconds: number;
}

interface CompletedMessage {
  type: "completed";
  job_id: string;
  summary: {
    passed: number;
    failed: number;
    total_time_seconds: number;
  };
}

interface ErrorMessage {
  type: "error";
  job_id: string;
  error: string;
}

// Client -> Server messages
interface CancelMessage {
  type: "cancel";
}
```

---

## Frontend State Management

```typescript
// types/validation.ts
interface ValidationState {
  // Single validation
  currentMolecule: string | null;
  validationResult: ValidationResult | null;
  isValidating: boolean;
  error: string | null;
  
  // Batch processing
  activeJob: BatchJob | null;
  jobProgress: Progress | null;
  batchResults: BatchResults | null;
  
  // History
  recentValidations: ValidationResult[];
  savedJobs: BatchJob[];
}

// React Query for API calls
function useValidation(molecule: string) {
  return useQuery({
    queryKey: ['validation', molecule],
    queryFn: () => validationApi.validate(molecule),
    enabled: !!molecule,
    staleTime: 5 * 60 * 1000, // Cache for 5 minutes
  });
}

function useBatchJob(jobId: string) {
  const queryClient = useQueryClient();
  
  // WebSocket for real-time updates
  useEffect(() => {
    const ws = new WebSocket(`${WS_URL}/jobs/${jobId}`);
    
    ws.onmessage = (event) => {
      const message = JSON.parse(event.data);
      
      if (message.type === 'progress') {
        queryClient.setQueryData(['job', jobId], (old) => ({
          ...old,
          progress: message,
        }));
      }
    };
    
    return () => ws.close();
  }, [jobId]);
  
  return useQuery({
    queryKey: ['job', jobId],
    queryFn: () => batchApi.getJob(jobId),
    refetchInterval: (data) => 
      data?.status === 'completed' ? false : 5000,
  });
}
```

---

## Deployment Architecture

### Docker Compose (Development)

```yaml
version: '3.8'

services:
  backend:
    build:
      context: ./backend
      dockerfile: Dockerfile
    ports:
      - "8000:8000"
    environment:
      - DATABASE_URL=postgresql://user:pass@postgres:5432/chemstructval
      - REDIS_URL=redis://redis:6379
      - CELERY_BROKER_URL=redis://redis:6379/0
    volumes:
      - ./backend:/app
    depends_on:
      - postgres
      - redis
    command: uvicorn app.main:app --host 0.0.0.0 --reload

  celery_worker:
    build:
      context: ./backend
      dockerfile: Dockerfile
    environment:
      - DATABASE_URL=postgresql://user:pass@postgres:5432/chemstructval
      - REDIS_URL=redis://redis:6379
      - CELERY_BROKER_URL=redis://redis:6379/0
    volumes:
      - ./backend:/app
    depends_on:
      - postgres
      - redis
    command: celery -A app.tasks.celery_app worker --loglevel=info

  frontend:
    build:
      context: ./frontend
      dockerfile: Dockerfile.dev
    ports:
      - "3000:3000"
    volumes:
      - ./frontend:/app
      - /app/node_modules
    command: npm run dev

  postgres:
    image: postgres:15
    environment:
      - POSTGRES_USER=user
      - POSTGRES_PASSWORD=pass
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

### Production Deployment (Kubernetes/Docker Swarm)

```
                    ┌─────────────────┐
                    │  Load Balancer  │
                    │   (nginx/HAProxy)│
                    └────────┬────────┘
                             │
           ┌─────────────────┼─────────────────┐
           │                 │                 │
    ┌──────▼──────┐   ┌──────▼──────┐   ┌──────▼──────┐
    │  Frontend   │   │  Frontend   │   │   Backend   │
    │  (Static)   │   │  (Static)   │   │   API x3    │
    │   via CDN   │   │   via CDN   │   │             │
    └─────────────┘   └─────────────┘   └──────┬──────┘
                                               │
                    ┌──────────────────────────┼──────┐
                    │                          │      │
             ┌──────▼──────┐            ┌──────▼──────┼──────┐
             │   Redis     │            │  PostgreSQL │      │
             │  (Cache +   │            │  (Primary)  │Replica│
             │   Queue)    │            │             │      │
             └─────────────┘            └─────────────┴──────┘
                    │
          ┌─────────┼─────────┐
          │         │         │
   ┌──────▼──┐ ┌────▼───┐ ┌───▼────┐
   │ Worker  │ │ Worker │ │ Worker │
   │   1     │ │   2    │ │   3    │
   └─────────┘ └────────┘ └────────┘
```

---

## Performance Considerations

### Caching Strategy

```python
# Redis caching for validation results
CACHE_CONFIG = {
    "validation_result": {
        "ttl": 3600,  # 1 hour
        "key_format": "val:{inchikey}:{checks_hash}"
    },
    "alert_patterns": {
        "ttl": 86400,  # 24 hours
        "key_format": "alerts:{set_name}"
    },
    "job_progress": {
        "ttl": 3600,
        "key_format": "job:{job_id}:progress"
    }
}

async def get_cached_validation(inchikey: str, checks: List[str]) -> Optional[ValidationResult]:
    """Get cached validation result if available."""
    checks_hash = hash(tuple(sorted(checks)))
    key = f"val:{inchikey}:{checks_hash}"
    
    cached = await redis.get(key)
    if cached:
        return ValidationResult.parse_raw(cached)
    return None
```

### Batch Processing Optimization

```python
# Parallel processing with controlled concurrency
async def process_batch(molecules: List[str], options: dict, job_id: str):
    """Process batch with parallel workers."""
    
    # Chunk molecules for processing
    chunk_size = 100
    chunks = [molecules[i:i + chunk_size] 
              for i in range(0, len(molecules), chunk_size)]
    
    # Process chunks in parallel with semaphore
    semaphore = asyncio.Semaphore(4)  # Max 4 concurrent chunks
    
    async def process_chunk(chunk, start_idx):
        async with semaphore:
            results = []
            for i, mol in enumerate(chunk):
                result = await validate_single(mol, options)
                results.append(result)
                
                # Update progress
                await update_progress(job_id, start_idx + i + 1, len(molecules))
            
            return results
    
    tasks = [process_chunk(chunk, i * chunk_size) 
             for i, chunk in enumerate(chunks)]
    
    all_results = await asyncio.gather(*tasks)
    return [r for chunk_results in all_results for r in chunk_results]
```

---

## Security Measures

```python
# Input validation
from pydantic import BaseModel, validator
from rdkit import Chem

class MoleculeInput(BaseModel):
    molecule: str
    format: str = "auto"
    
    @validator('molecule')
    def validate_molecule(cls, v):
        # Size limit
        if len(v) > 10000:
            raise ValueError("Molecule string too long")
        
        # Prevent code injection
        if any(c in v for c in ['<', '>', '&', ';', '|', '$', '`']):
            raise ValueError("Invalid characters in molecule")
        
        return v

# Rate limiting
from slowapi import Limiter
from slowapi.util import get_remote_address

limiter = Limiter(key_func=get_remote_address)

@app.post("/api/v1/validate")
@limiter.limit("60/minute")
async def validate(request: Request, input: MoleculeInput):
    ...

# File upload validation
MAX_FILE_SIZE = 100 * 1024 * 1024  # 100MB
ALLOWED_EXTENSIONS = {'.sdf', '.mol', '.csv', '.xlsx'}

async def validate_upload(file: UploadFile):
    # Check extension
    ext = Path(file.filename).suffix.lower()
    if ext not in ALLOWED_EXTENSIONS:
        raise HTTPException(400, f"File type {ext} not allowed")
    
    # Check size
    contents = await file.read()
    if len(contents) > MAX_FILE_SIZE:
        raise HTTPException(400, "File too large")
    
    await file.seek(0)  # Reset for reading
    return contents
```

---

*This architecture document should be updated as the system evolves.*
