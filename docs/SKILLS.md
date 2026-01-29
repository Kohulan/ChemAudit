# ChemStructVal: Technical Skills & Requirements

## Overview

This document outlines the technical skills, libraries, and knowledge required to build ChemStructVal. It serves as a reference for developers and a checklist for project readiness.

---

## Core Technology Skills

### Frontend Development

#### React & TypeScript
```
Required Knowledge:
├── React 18+ fundamentals
│   ├── Functional components
│   ├── Hooks (useState, useEffect, useCallback, useMemo, useReducer)
│   ├── Context API for state management
│   └── React Router for navigation
├── TypeScript
│   ├── Type definitions and interfaces
│   ├── Generics
│   └── Type guards
└── Build tools
    ├── Vite or Create React App
    └── npm/yarn package management
```

#### UI Framework
```
Tailwind CSS:
├── Utility-first CSS
├── Custom configuration
├── Responsive design
└── Dark mode support

Component Libraries:
├── shadcn/ui (primary)
├── Headless UI
└── Radix UI primitives
```

#### Chemical Structure Visualization
```
Required Libraries:
├── RDKit.js
│   ├── @rdkit/rdkit (WebAssembly build)
│   ├── Molecule rendering
│   ├── Substructure highlighting
│   └── SVG generation
├── Ketcher (optional, for structure editing)
│   └── ketcher-react
└── 3Dmol.js (optional, for 3D viewing)
    └── 3dmol/3dmol.js
```

---

### Backend Development

#### Python Framework
```
FastAPI:
├── Async request handling
├── Pydantic models for validation
├── Dependency injection
├── Background tasks
├── WebSocket support (for progress updates)
└── OpenAPI/Swagger documentation
```

#### Required Python Skills
```python
# Key concepts needed:
├── Async/await patterns
├── Type hints
├── Context managers
├── Decorators
├── Generators (for batch processing)
└── Exception handling patterns
```

---

### Cheminformatics Libraries

#### RDKit (Primary)
```python
# Core modules needed:
from rdkit import Chem
from rdkit.Chem import (
    AllChem,              # Coordinates, conformers
    Descriptors,          # Molecular descriptors
    rdMolDescriptors,     # Advanced descriptors
    FilterCatalog,        # PAINS, BRENK, etc.
    MolStandardize,       # Standardization
    Draw,                 # 2D depiction
    inchi,                # InChI generation
    rdDepictor,           # 2D coordinate generation
    rdMolEnumerator,      # Tautomer enumeration
)
from rdkit.Chem.MolStandardize import rdMolStandardize

# Key operations:
├── Mol parsing (SMILES, MOL, SDF)
├── Sanitization and validation
├── Standardization pipeline
├── Property calculation
├── Substructure matching (SMARTS)
├── Fingerprint generation
└── 2D coordinate generation
```

#### CDK (Chemistry Development Kit)
```java
// Via Python wrapper (jpype) or REST service
org.openscience.cdk:
├── AtomContainer operations
├── Aromaticity detection
├── Stereochemistry handling
├── Descriptor calculation
└── InChI generation (for comparison)
```

#### Additional Chemistry Libraries
```python
# MolVS (built on RDKit)
from molvs import Standardizer, Validator

# ChEMBL Structure Pipeline
from chembl_structure_pipeline import checker, standardizer

# Open Babel (optional, for format conversion)
from openbabel import openbabel
```

---

### Validation Rule Knowledge

#### Structural Validation Rules

| Category | Knowledge Required |
|----------|-------------------|
| Valence Rules | Electron counting, hypervalent atoms, radicals |
| Stereochemistry | CIP rules (R/S, E/Z), pseudochiral centers |
| Aromaticity | Hückel's rule, heterocycles, fused systems |
| Ring Strain | Small rings, bridged systems, spiro centers |
| Charge States | Zwitterions, betaines, resonance forms |
| Tautomerism | Keto-enol, amide-imidol, ring-chain |

#### Structural Alert Patterns (SMARTS)

```
Required SMARTS knowledge:
├── Basic atom matching [C], [N], [O]
├── Bond type specification -, =, #, :
├── Logical operators [!], [,], [;]
├── Ring membership specification
├── Recursive SMARTS
└── Reaction SMARTS (optional)
```

**Key Alert Categories:**
```
PAINS (480 patterns):
├── Quinones and related
├── Catechols
├── Rhodanines
├── Enones
└── ...

BRENK (105 patterns):
├── Alkyl halides
├── Michael acceptors
├── Epoxides
└── ...
```

---

### Database Skills

#### PostgreSQL
```sql
-- Required knowledge:
├── Schema design for chemical data
├── JSONB for flexible metadata
├── Full-text search
├── Indexing strategies
└── Connection pooling (asyncpg)

-- Example schema pattern:
CREATE TABLE validation_jobs (
    id UUID PRIMARY KEY,
    status VARCHAR(20),
    input_format VARCHAR(10),
    results JSONB,
    created_at TIMESTAMP,
    completed_at TIMESTAMP
);
```

#### Redis (Caching)
```python
# Caching strategies:
├── Job status caching
├── Validation result caching (by InChIKey)
├── Rate limiting
└── Session management
```

---

### DevOps & Infrastructure

#### Docker
```dockerfile
# Skills needed:
├── Multi-stage builds
├── Layer optimization
├── docker-compose for local dev
├── Health checks
└── Volume management
```

#### CI/CD (GitHub Actions)
```yaml
# Workflow knowledge:
├── Testing workflows
├── Build and publish
├── Deployment automation
└── Dependency caching
```

---

## Library Reference

### Python Dependencies

```toml
# pyproject.toml or requirements.txt
[dependencies]
# Core
fastapi = "^0.100.0"
uvicorn = "^0.22.0"
python-multipart = "^0.0.6"
pydantic = "^2.0.0"

# Chemistry
rdkit = "^2023.9.0"
molvs = "^0.1.1"
chembl-structure-pipeline = "^1.2.0"

# Database
asyncpg = "^0.28.0"
redis = "^4.6.0"
sqlalchemy = "^2.0.0"

# Utilities
pandas = "^2.0.0"
numpy = "^1.24.0"
aiofiles = "^23.0.0"

# Testing
pytest = "^7.4.0"
pytest-asyncio = "^0.21.0"
httpx = "^0.24.0"
```

### JavaScript Dependencies

```json
{
  "dependencies": {
    "react": "^18.2.0",
    "react-dom": "^18.2.0",
    "react-router-dom": "^6.14.0",
    "@tanstack/react-query": "^4.32.0",
    "axios": "^1.4.0",
    "@rdkit/rdkit": "^2023.9.1",
    "tailwindcss": "^3.3.0",
    "lucide-react": "^0.263.0",
    "recharts": "^2.7.0",
    "react-dropzone": "^14.2.0"
  },
  "devDependencies": {
    "typescript": "^5.1.0",
    "vite": "^4.4.0",
    "@types/react": "^18.2.0",
    "vitest": "^0.34.0",
    "@testing-library/react": "^14.0.0"
  }
}
```

---

## API Design Patterns

### RESTful Endpoints
```
Required patterns:
├── Resource-based URLs
├── HTTP method semantics (GET, POST, PUT, DELETE)
├── Status code usage
├── Pagination
├── Filtering
└── Error response format
```

### WebSocket (Progress Updates)
```python
# For batch processing progress
├── Connection management
├── Message protocols (JSON)
├── Heartbeat/keepalive
└── Error recovery
```

---

## Testing Knowledge

### Unit Testing
```python
# pytest patterns:
├── Fixtures for test molecules
├── Parametrized tests for validation rules
├── Mocking external services
└── Async test handling
```

### Integration Testing
```python
# API testing:
├── TestClient usage
├── Database fixtures
├── File upload testing
└── WebSocket testing
```

### Frontend Testing
```typescript
// React Testing Library:
├── Component rendering
├── User interaction simulation
├── Async state testing
└── Accessibility testing
```

---

## File Format Knowledge

### Input Formats
```
Must handle:
├── SMILES (various flavors)
│   ├── Canonical SMILES
│   ├── Isomeric SMILES
│   └── Extended SMILES (CXSMILES)
├── SDF/MOL (V2000, V3000)
├── InChI / InChIKey
├── CSV with SMILES column
└── Excel (.xlsx)
```

### Output Formats
```
Must generate:
├── JSON (validation results)
├── CSV (batch results)
├── Excel (formatted reports)
├── PDF (summary reports)
├── SDF (standardized structures)
└── PNG/SVG (structure images)
```

---

## Security Considerations

```
Knowledge required:
├── Input sanitization (file uploads)
├── Rate limiting
├── CORS configuration
├── File size limits
├── Temporary file cleanup
└── No arbitrary code execution from SMILES
```

---

## Performance Optimization

```
Techniques needed:
├── Batch processing with workers
├── Async I/O for file operations
├── Caching validation results
├── Lazy loading for frontend
├── WebAssembly for client-side chemistry
└── Database query optimization
```

---

## Learning Resources

### Cheminformatics
- RDKit Documentation: https://www.rdkit.org/docs/
- RDKit Cookbook: https://www.rdkit.org/docs/Cookbook.html
- TeachOpenCADD: https://github.com/volkamerlab/teachopencadd

### Web Development
- FastAPI Tutorial: https://fastapi.tiangolo.com/tutorial/
- React Documentation: https://react.dev/
- TypeScript Handbook: https://www.typescriptlang.org/docs/

### Chemical Validation
- ChEMBL Curation Paper: https://doi.org/10.1186/s13321-020-00456-1
- MolVS Documentation: https://molvs.readthedocs.io/
- PAINS Original Paper: https://doi.org/10.1021/jm901137j

---

## Skill Assessment Checklist

Use this checklist to assess team readiness:

### Essential (Must Have)
- [ ] Python (intermediate+)
- [ ] FastAPI basics
- [ ] RDKit molecule handling
- [ ] React functional components
- [ ] TypeScript basics
- [ ] REST API design
- [ ] Git version control

### Important (Should Have)
- [ ] RDKit advanced (standardization, fingerprints)
- [ ] PostgreSQL
- [ ] Docker
- [ ] Tailwind CSS
- [ ] Testing (pytest, Jest)

### Nice to Have
- [ ] CDK / Open Babel
- [ ] Redis
- [ ] WebSocket implementation
- [ ] GitHub Actions
- [ ] RDKit.js

---

*This document should be updated as the project evolves and new requirements emerge.*
