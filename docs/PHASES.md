# ChemStructVal: Development Phases & To-Do Lists

## Overview

This document breaks down the ChemStructVal development into 4 phases over approximately 14 weeks. Each phase has specific goals, deliverables, and detailed task lists.

---

## Phase Summary

| Phase | Duration | Focus | Deliverable |
|-------|----------|-------|-------------|
| Phase 1 | 4 weeks | Foundation | MVP with basic validation |
| Phase 2 | 4 weeks | Enhancement | Full validation suite + batch |
| Phase 3 | 3 weeks | Integration | API + external integrations |
| Phase 4 | 3 weeks | Polish | Reports, docs, deployment |

---

# Phase 1: Foundation (Weeks 1-4)

## Goals
- Establish project structure and development environment
- Implement core validation engine
- Create basic React UI for single-molecule validation
- Set up CI/CD pipeline

## Week 1: Project Setup & Architecture

### Backend Setup
- [ ] Initialize Python project with Poetry/uv
  - [ ] Create `pyproject.toml` with dependencies
  - [ ] Set up virtual environment
  - [ ] Configure pre-commit hooks (black, isort, mypy)
- [ ] Create FastAPI application structure
  ```
  backend/
  ├── app/
  │   ├── __init__.py
  │   ├── main.py
  │   ├── api/
  │   │   ├── __init__.py
  │   │   ├── routes/
  │   │   └── dependencies.py
  │   ├── core/
  │   │   ├── config.py
  │   │   └── logging.py
  │   ├── services/
  │   │   └── validation/
  │   ├── models/
  │   └── schemas/
  ├── tests/
  └── Dockerfile
  ```
- [ ] Set up configuration management (pydantic-settings)
- [ ] Configure logging (structured JSON logs)
- [ ] Create health check endpoint
- [ ] Write initial tests for app startup

### Frontend Setup
- [ ] Initialize React project with Vite + TypeScript
  ```bash
  npm create vite@latest frontend -- --template react-ts
  ```
- [ ] Install and configure Tailwind CSS
- [ ] Set up shadcn/ui component library
- [ ] Create folder structure
  ```
  frontend/
  ├── src/
  │   ├── components/
  │   │   ├── ui/           # shadcn components
  │   │   ├── molecules/    # structure viewers
  │   │   └── validation/   # validation result displays
  │   ├── pages/
  │   ├── hooks/
  │   ├── services/
  │   ├── types/
  │   └── utils/
  ├── public/
  └── package.json
  ```
- [ ] Configure React Router
- [ ] Set up Axios for API calls
- [ ] Create basic layout component (header, footer, main)

### DevOps Setup
- [ ] Create Docker Compose for local development
  ```yaml
  services:
    backend:
      build: ./backend
      ports: ["8000:8000"]
    frontend:
      build: ./frontend
      ports: ["3000:3000"]
    postgres:
      image: postgres:15
    redis:
      image: redis:7
  ```
- [ ] Create GitHub repository with branch protection
- [ ] Set up GitHub Actions for CI
  - [ ] Backend tests workflow
  - [ ] Frontend tests workflow
  - [ ] Linting workflow
- [ ] Create development environment documentation

### Deliverable: Running dev environment with "Hello World" API and React app

---

## Week 2: Core Validation Engine

### Molecule Parsing Service
- [ ] Create `MoleculeParser` class
  ```python
  class MoleculeParser:
      def parse_smiles(self, smiles: str) -> ParseResult
      def parse_molblock(self, molblock: str) -> ParseResult
      def parse_sdf(self, sdf_content: str) -> List[ParseResult]
      def detect_format(self, content: str) -> MoleculeFormat
  ```
- [ ] Implement SMILES parsing with error handling
- [ ] Implement MOL/SDF parsing
- [ ] Create format auto-detection
- [ ] Add support for multiple SMILES flavors
- [ ] Write comprehensive parser tests (50+ test molecules)

### Basic Validation Checks
- [ ] Create `ValidationEngine` class
  ```python
  class ValidationEngine:
      def validate(self, mol: Mol, checks: List[str]) -> ValidationResult
      def get_available_checks(self) -> List[CheckInfo]
  ```
- [ ] Implement individual check classes:
  - [ ] `ParsabilityCheck` - Can molecule be parsed
  - [ ] `SanitizationCheck` - RDKit sanitization passes
  - [ ] `ValenceCheck` - All atoms have valid valence
  - [ ] `AromaticityCheck` - Kekulization succeeds
  - [ ] `ConnectivityCheck` - Single connected component
- [ ] Create check result schema
  ```python
  class CheckResult(BaseModel):
      check_name: str
      passed: bool
      severity: Severity
      message: str
      affected_atoms: List[int] = []
  ```
- [ ] Implement overall score calculation
- [ ] Write unit tests for each check (10+ test cases each)

### API Endpoint
- [ ] Create `/api/v1/validate` POST endpoint
- [ ] Define request/response schemas
- [ ] Add input validation (Pydantic)
- [ ] Implement basic error handling
- [ ] Add request logging
- [ ] Write API integration tests

### Deliverable: Working validation API for single molecules with 5 basic checks

---

## Week 3: React UI - Single Molecule

### Structure Input Component
- [ ] Create `StructureInput` component
  - [ ] SMILES text input with debounced validation
  - [ ] "Validate" button
  - [ ] Format indicator
  - [ ] Loading state
- [ ] Add input validation feedback
- [ ] Create example molecules dropdown
- [ ] Style with Tailwind (clean, professional)

### Structure Visualization
- [ ] Install and configure RDKit.js
  ```typescript
  import { RDKitModule } from '@rdkit/rdkit';
  ```
- [ ] Create `MoleculeViewer` component
  - [ ] SVG rendering from SMILES
  - [ ] Atom highlighting capability
  - [ ] Zoom controls (optional)
- [ ] Handle rendering errors gracefully
- [ ] Add loading skeleton

### Validation Results Display
- [ ] Create `ValidationResults` component
  - [ ] Overall score with circular progress
  - [ ] Color-coded status (green/yellow/red)
  - [ ] Issue list with severity icons
  - [ ] Expandable issue details
- [ ] Create `IssueCard` sub-component
  - [ ] Severity badge
  - [ ] Description
  - [ ] Affected atoms display
- [ ] Add "Copy Results" functionality

### Page Assembly
- [ ] Create `SingleValidationPage` component
- [ ] Wire up API calls
- [ ] Add error handling UI
- [ ] Implement loading states
- [ ] Add success animations

### Deliverable: Functional single-molecule validation page

---

## Week 4: Stereochemistry & Representation Checks

### Stereochemistry Validation
- [ ] Implement stereochemistry check classes:
  - [ ] `UndefinedStereoCentersCheck`
    - Count chiral atoms without defined stereo
    - Report specific atom indices
  - [ ] `UndefinedDoubleBondStereoCheck`
    - Identify E/Z undefined double bonds
  - [ ] `ExcessiveStereoCentersCheck`
    - Flag molecules with >N stereocenters (configurable)
  - [ ] `ConflictingStereoCheck`
    - Detect impossible combinations
- [ ] Add stereo summary to results
- [ ] Write stereo test suite (known molecules with stereo issues)

### Representation Consistency
- [ ] Implement representation checks:
  - [ ] `SmilesRoundtripCheck`
    ```python
    original -> Mol -> SMILES -> Mol -> SMILES
    # Compare InChIKeys
    ```
  - [ ] `InchiGenerationCheck`
    - Can valid InChI be generated
  - [ ] `InchiRoundtripCheck`
    - InChI -> Mol -> InChI consistency
  - [ ] `CanonicalSmilesCheck`
    - Generate and verify canonical form
- [ ] Create representation diff display
- [ ] Add tests with known problematic molecules

### UI Enhancements
- [ ] Add stereochemistry section to results
- [ ] Create stereo visualization overlay
  - [ ] Show stereo bonds
  - [ ] Highlight undefined centers
- [ ] Add "View Canonical" toggle
- [ ] Create before/after comparison view

### Testing & Bug Fixes
- [ ] Add integration tests for full validation flow
- [ ] Fix bugs found in Week 3
- [ ] Performance profiling
- [ ] Code review and refactoring

### Deliverable: Complete Phase 1 MVP with stereo and representation checks

---

## Phase 1 Checklist Summary

### Must Have ✓
- [ ] Project structure established
- [ ] FastAPI backend running
- [ ] React frontend running
- [ ] Docker Compose development environment
- [ ] GitHub CI/CD pipeline
- [ ] 8 basic validation checks implemented
- [ ] Single molecule validation API
- [ ] Working UI with structure rendering
- [ ] Validation results display
- [ ] Basic error handling

### Nice to Have
- [ ] Dark mode toggle
- [ ] Mobile responsive design
- [ ] Keyboard shortcuts

---

# Phase 2: Enhancement (Weeks 5-8)

## Goals
- Implement all structural alert checks
- Add NP-likeness and ML-readiness scoring
- Build standardization pipeline
- Create batch processing system

## Week 5: Structural Alerts

### Alert System Architecture
- [ ] Create `AlertEngine` class
  ```python
  class AlertEngine:
      def __init__(self, alert_sets: List[str])
      def screen(self, mol: Mol) -> List[AlertMatch]
      def get_alert_sets(self) -> List[AlertSetInfo]
  ```
- [ ] Define alert pattern storage schema
- [ ] Create alert loading from JSON/CSV
- [ ] Implement SMARTS matching engine

### Alert Set Implementation
- [ ] PAINS Alerts (480 patterns)
  - [ ] Load PAINS A patterns
  - [ ] Load PAINS B patterns
  - [ ] Load PAINS C patterns
  - [ ] Map to descriptive names
- [ ] BRENK Alerts (105 patterns)
- [ ] NIH Alerts
- [ ] ZINC Filters
- [ ] ChEMBL Alert Sets
  - [ ] BMS (180)
  - [ ] Dundee (105)
  - [ ] Glaxo (55)
  - [ ] Inpharmatica (91)
  - [ ] LINT (57)
  - [ ] MLSMR (116)
  - [ ] SureChEMBL (166)

### Alert API & UI
- [ ] Add `/api/v1/alerts` endpoint
- [ ] Update `/api/v1/validate` to include alerts
- [ ] Create `AlertResults` component
  - [ ] Alert set selector
  - [ ] Matched patterns display
  - [ ] Structure highlighting for matches
- [ ] Add alert filtering options
- [ ] Write alert system tests

### Deliverable: Full structural alerts screening

---

## Week 6: NP-Likeness & ML-Readiness

### NP-Likeness Implementation
- [ ] Integrate RDKit NP_Score
  ```python
  from rdkit.Contrib.NP_Score import npscorer
  ```
- [ ] Create `NPLikenessScorer` class
  - [ ] Score calculation
  - [ ] Confidence assessment
  - [ ] Classification (NP-like/Synthetic-like)
- [ ] Add scaffold analysis
  - [ ] Murcko scaffold extraction
  - [ ] Core/side chain analysis
- [ ] Implement sugar moiety detection
- [ ] Create NP results section in UI

### ML-Readiness Scoring
- [ ] Create `MLReadinessScorer` class
  ```python
  class MLReadinessScorer:
      def score(self, mol: Mol) -> MLReadinessResult
      
  class MLReadinessResult:
      overall_score: int  # 0-100
      parsability: bool
      sanitizes: bool
      descriptors_ok: bool
      fingerprints_ok: bool
      size_ok: bool
      elements_ok: bool
      issues: List[str]
  ```
- [ ] Implement component checks:
  - [ ] Descriptor calculability (RDKit 200 descriptors)
  - [ ] Fingerprint generation (Morgan, MACCS)
  - [ ] Size check (MW, atom count ranges)
  - [ ] Element check (common ML elements)
- [ ] Create weighted score calculation
- [ ] Add ML-readiness UI section
- [ ] Write ML readiness tests

### API Updates
- [ ] Update validation response schema
- [ ] Add NP and ML scores to results
- [ ] Create dedicated endpoints:
  - [ ] `/api/v1/np-likeness`
  - [ ] `/api/v1/ml-readiness`

### Deliverable: NP-likeness and ML-readiness scoring

---

## Week 7: Standardization Pipeline

### Standardizer Implementation
- [ ] Create `StandardizationPipeline` class
  ```python
  class StandardizationPipeline:
      def __init__(self, config: StandardizationConfig)
      def standardize(self, mol: Mol) -> StandardizationResult
      
  class StandardizationResult:
      standardized_mol: Mol
      operations_applied: List[str]
      original_smiles: str
      standardized_smiles: str
      changes_made: List[Change]
  ```
- [ ] Implement standardization steps:
  - [ ] `SaltStripper` - Remove salts/counterions
  - [ ] `SolventRemover` - Remove solvent molecules
  - [ ] `Normalizer` - Functional group normalization
  - [ ] `Reionizer` - Charge normalization
  - [ ] `TautomerCanonicalizer` - Canonical tautomer
  - [ ] `IsotopeRemover` - Strip isotopes (optional)
  - [ ] `StereoRemover` - Strip stereo (optional)
- [ ] Create configurable pipeline
- [ ] Add before/after tracking

### Parent Extraction
- [ ] Implement `ParentExtractor` class
- [ ] Handle multi-component molecules
- [ ] Keep largest fragment logic
- [ ] Handle salts database

### Standardization UI
- [ ] Create `StandardizationResults` component
  - [ ] Side-by-side comparison
  - [ ] Operations applied list
  - [ ] Difference highlighting
- [ ] Add standardization options panel
- [ ] Create "Apply Standardization" workflow

### API Endpoints
- [ ] Create `/api/v1/standardize` endpoint
- [ ] Add standardization to validation flow (optional)
- [ ] Support batch standardization

### Deliverable: Complete standardization pipeline

---

## Week 8: Batch Processing

### Job Queue System
- [ ] Set up Celery with Redis
  ```python
  # tasks.py
  @celery.task
  def validate_batch(job_id: str, molecules: List[str], options: dict)
  ```
- [ ] Create job status management
- [ ] Implement progress tracking
- [ ] Add job persistence (PostgreSQL)
- [ ] Create job cleanup task

### Batch Validation Logic
- [ ] Create `BatchValidator` class
  ```python
  class BatchValidator:
      def start_batch(self, molecules, options) -> JobId
      def get_progress(self, job_id) -> Progress
      def get_results(self, job_id) -> BatchResults
      def cancel(self, job_id) -> bool
  ```
- [ ] Implement chunked processing
- [ ] Add partial failure handling
- [ ] Create result aggregation
- [ ] Implement statistics calculation

### File Upload Handling
- [ ] Create file upload endpoint
- [ ] Implement SDF parsing (multi-mol)
- [ ] Implement CSV parsing
- [ ] Add file validation (size, format)
- [ ] Create temp file management

### Batch UI
- [ ] Create `BatchUpload` component
  - [ ] Drag-and-drop zone
  - [ ] File format selector
  - [ ] Column mapping (CSV)
  - [ ] Upload progress
- [ ] Create `BatchProgress` component
  - [ ] Progress bar
  - [ ] ETA display
  - [ ] Current status
  - [ ] Cancel button
- [ ] Create `BatchResults` component
  - [ ] Summary statistics
  - [ ] Sortable/filterable table
  - [ ] Row selection
  - [ ] Export options
- [ ] Create `JobHistory` component

### WebSocket Progress Updates
- [ ] Implement WebSocket endpoint
- [ ] Create progress message protocol
- [ ] Add frontend WebSocket hook
- [ ] Handle reconnection

### Deliverable: Working batch processing up to 10,000 molecules

---

## Phase 2 Checklist Summary

### Must Have ✓
- [ ] All structural alert sets (PAINS, BRENK, etc.)
- [ ] NP-likeness scoring
- [ ] ML-readiness scoring
- [ ] Standardization pipeline
- [ ] Batch processing (10K molecules)
- [ ] File upload (SDF, CSV)
- [ ] Job progress tracking
- [ ] Batch results export

### Nice to Have
- [ ] Custom alert patterns
- [ ] Excel file support
- [ ] Email notifications

---

# Phase 3: Integration (Weeks 9-11)

## Goals
- Create production-ready REST API
- Build integrations with external services
- Implement export/report generation
- Add user preferences

## Week 9: Production API

### API Refinement
- [ ] Complete OpenAPI documentation
- [ ] Add request validation
- [ ] Implement rate limiting
  ```python
  from slowapi import Limiter
  limiter = Limiter(key_func=get_remote_address)
  ```
- [ ] Add API versioning
- [ ] Create API error codes
- [ ] Implement response caching

### Authentication (Basic)
- [ ] Add optional API key authentication
- [ ] Create API key management
- [ ] Implement usage tracking
- [ ] Add rate limit tiers

### SDK/Client Library
- [ ] Create Python client package
  ```python
  from chemstructval import Client
  
  client = Client(api_key="...")
  result = client.validate("CCO")
  ```
- [ ] Create JavaScript/TypeScript client
- [ ] Write client documentation
- [ ] Add client examples

### API Documentation Page
- [ ] Create interactive API docs page
- [ ] Add code examples (Python, curl, JavaScript)
- [ ] Create Postman collection
- [ ] Write quickstart guide

### Deliverable: Production-ready API with documentation

---

## Week 10: External Integrations

### DECIMER Integration
- [ ] Create DECIMER validation workflow
  ```python
  class DECIMERValidator:
      def validate_ocsr_output(self, smiles: str, confidence: float) -> ValidationResult
      def flag_low_confidence(self, result: DECIMERResult) -> List[Warning]
  ```
- [ ] Add OCSR-specific checks
- [ ] Create DECIMER-focused UI mode

### COCONUT Integration
- [ ] Create COCONUT validation endpoint
- [ ] Add NP database comparison
- [ ] Implement structure lookup
- [ ] Create COCONUT report format

### Cheminformatics Microservice
- [ ] Create microservice module
- [ ] Define API contract
- [ ] Implement validation endpoint
- [ ] Add documentation

### PubChem/ChEMBL Lookup (Optional)
- [ ] Add compound lookup by ID
- [ ] Fetch and validate external structures
- [ ] Create comparison view

### Deliverable: Working integrations with DECIMER, COCONUT

---

## Week 11: Export & Reports

### Export Formats
- [ ] Implement CSV export
  ```python
  def export_csv(results: List[ValidationResult]) -> bytes
  ```
- [ ] Implement Excel export (openpyxl)
  - [ ] Formatted headers
  - [ ] Conditional formatting
  - [ ] Summary sheet
- [ ] Implement SDF export (standardized structures)
- [ ] Implement JSON export

### PDF Report Generation
- [ ] Create report template
  ```
  Report Sections:
  ├── Summary Statistics
  ├── Issue Distribution Chart
  ├── Critical Issues List
  ├── Molecule Gallery (flagged)
  └── Appendix: Full Results
  ```
- [ ] Implement PDF generation (ReportLab or WeasyPrint)
- [ ] Add charts (matplotlib/plotly)
- [ ] Create structure images for report

### Export UI
- [ ] Create `ExportDialog` component
  - [ ] Format selector
  - [ ] Options (include structures, etc.)
  - [ ] Download button
- [ ] Add export to batch results
- [ ] Create "Generate Report" button

### Deliverable: Multi-format export and PDF reports

---

## Phase 3 Checklist Summary

### Must Have ✓
- [ ] OpenAPI documentation
- [ ] Rate limiting
- [ ] API key authentication
- [ ] Python client library
- [ ] CSV/Excel export
- [ ] SDF export
- [ ] Basic PDF report

### Nice to Have
- [ ] DECIMER integration
- [ ] COCONUT integration
- [ ] PubChem lookup
- [ ] JavaScript client

---

# Phase 4: Polish (Weeks 12-14)

## Goals
- Complete documentation
- Performance optimization
- Production deployment
- Launch preparation

## Week 12: Documentation & Testing

### User Documentation
- [ ] Write user guide
  - [ ] Getting started
  - [ ] Single molecule validation
  - [ ] Batch processing
  - [ ] Understanding results
  - [ ] API usage
- [ ] Create video tutorials (optional)
- [ ] Write FAQ
- [ ] Create troubleshooting guide

### Technical Documentation
- [ ] Complete API reference
- [ ] Write architecture docs
- [ ] Create deployment guide
- [ ] Document validation rules
  - [ ] Each check explained
  - [ ] Severity rationale
  - [ ] Examples

### In-App Help
- [ ] Add tooltips to UI elements
- [ ] Create help modal
- [ ] Add example molecules with explanations

### Testing
- [ ] Achieve 80% code coverage
- [ ] Add end-to-end tests (Playwright/Cypress)
- [ ] Create performance benchmarks
- [ ] Run security scan (OWASP ZAP)

### Deliverable: Complete documentation

---

## Week 13: Performance & Optimization

### Backend Optimization
- [ ] Profile validation engine
- [ ] Optimize SMARTS matching
- [ ] Add result caching (Redis)
  ```python
  @cache(ttl=3600, key="validation:{inchikey}")
  def validate(mol):
      ...
  ```
- [ ] Optimize database queries
- [ ] Add connection pooling

### Frontend Optimization
- [ ] Implement code splitting
- [ ] Optimize bundle size
- [ ] Add lazy loading
- [ ] Implement virtual scrolling for large tables
- [ ] Add service worker (optional)

### Batch Processing Optimization
- [ ] Tune worker configuration
- [ ] Optimize chunk sizes
- [ ] Add parallel processing
- [ ] Profile memory usage

### Infrastructure
- [ ] Set up production database
- [ ] Configure Redis cluster (if needed)
- [ ] Set up CDN for static assets
- [ ] Configure auto-scaling (if cloud)

### Deliverable: Optimized application meeting performance targets

---

## Week 14: Deployment & Launch

### Production Deployment
- [ ] Create production Docker images
- [ ] Set up production environment
  - [ ] Database migrations
  - [ ] Environment variables
  - [ ] SSL certificates
- [ ] Configure monitoring (Prometheus/Grafana)
- [ ] Set up error tracking (Sentry)
- [ ] Configure backups

### Pre-Launch Checklist
- [ ] Security review
- [ ] Load testing
- [ ] Cross-browser testing
- [ ] Mobile testing
- [ ] Accessibility audit
- [ ] GDPR compliance check

### Launch Preparation
- [ ] Create landing page
- [ ] Write launch blog post
- [ ] Prepare social media
- [ ] Notify early users
- [ ] Create feedback mechanism

### Post-Launch
- [ ] Monitor metrics
- [ ] Address critical bugs
- [ ] Gather user feedback
- [ ] Plan v1.1 features

### Deliverable: Deployed, production-ready application

---

## Phase 4 Checklist Summary

### Must Have ✓
- [ ] User documentation
- [ ] API documentation
- [ ] 80% test coverage
- [ ] Performance optimization
- [ ] Production deployment
- [ ] Monitoring setup
- [ ] Error tracking

### Nice to Have
- [ ] Video tutorials
- [ ] Custom domain
- [ ] Analytics dashboard

---

# Risk Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| RDKit.js performance | Medium | Fallback to server-side rendering |
| Large batch processing | High | Implement chunking and progress saving |
| Alert pattern accuracy | Medium | Test against known datasets |
| API abuse | Medium | Rate limiting and authentication |
| Browser compatibility | Low | Progressive enhancement |

---

# Success Criteria

## Phase 1 (Week 4)
- [ ] Can validate any valid SMILES
- [ ] Displays 8 check results
- [ ] Structure renders correctly
- [ ] Response time < 3 seconds

## Phase 2 (Week 8)
- [ ] All alert sets working
- [ ] NP/ML scores accurate
- [ ] Standardization matches ChEMBL pipeline
- [ ] Batch processes 1000 molecules in < 2 minutes

## Phase 3 (Week 11)
- [ ] API documentation complete
- [ ] Export formats working
- [ ] Integrations operational
- [ ] PDF reports generated

## Phase 4 (Week 14)
- [ ] All documentation complete
- [ ] Performance targets met
- [ ] Successfully deployed
- [ ] No critical bugs

---

*This document should be reviewed and updated at the end of each week.*
