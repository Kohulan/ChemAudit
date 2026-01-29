# ChemStructVal: Product Requirements Document (PRD)

**Version**: 1.0  
**Last Updated**: January 2026  
**Status**: Draft  

---

## Executive Summary

ChemStructVal is a web-based chemical structure validation platform that consolidates disparate validation tools into a unified, visual interface. It addresses the critical need for accessible, comprehensive chemical data quality assessment serving database curators, ML researchers, medicinal chemists, and natural product scientists.

---

## Table of Contents

1. [Product Vision](#1-product-vision)
2. [User Personas](#2-user-personas)
3. [Functional Requirements](#3-functional-requirements)
4. [Non-Functional Requirements](#4-non-functional-requirements)
5. [User Interface Requirements](#5-user-interface-requirements)
6. [API Requirements](#6-api-requirements)
7. [Data Requirements](#7-data-requirements)
8. [Security Requirements](#8-security-requirements)
9. [Performance Requirements](#9-performance-requirements)
10. [Integration Requirements](#10-integration-requirements)

---

## 1. Product Vision

### 1.1 Vision Statement

ChemStructVal democratizes chemical structure validation by providing an intuitive web interface that makes professional-grade quality assessment accessible to all chemists, regardless of programming expertise.

### 1.2 Goals

| Goal | Success Criteria |
|------|------------------|
| Accessibility | 90% of features usable without coding |
| Comprehensiveness | Cover all major validation categories |
| Speed | <5s for single molecule, <5min for 10K batch |
| Accuracy | 99%+ agreement with reference implementations |
| Adoption | 500 MAU within 12 months |

### 1.3 Scope

**In Scope:**
- Single molecule validation
- Batch validation (up to 100,000 molecules)
- Visual error reporting
- Standardization pipeline
- Structural alerts screening
- NP-likeness assessment
- ML-readiness scoring
- Report generation
- REST API access

**Out of Scope (v1.0):**
- Custom ML model training
- Reaction validation
- Protein-ligand analysis
- Commercial database integration
- Real-time collaborative editing

---

## 2. User Personas

### 2.1 Primary Persona: Dr. Maria Chen - Database Curator

**Background:**
- PhD in Chemistry
- 5 years curating a public compound database
- Moderate Python skills
- Uses ChEMBL pipeline via CLI

**Goals:**
- Validate 500+ new submissions weekly
- Generate reports for internal review
- Identify compounds requiring manual attention

**Pain Points:**
- CLI tools lack visual feedback
- Combining multiple tools is time-consuming
- No easy way to share results with non-technical colleagues

**Feature Priorities:**
1. Batch upload (CSV, SDF)
2. Downloadable reports (Excel, PDF)
3. Severity-based filtering
4. Visual structure comparison

---

### 2.2 Primary Persona: Alex Rodriguez - ML Researcher

**Background:**
- PhD Student in Computational Chemistry
- Building QSAR models
- Strong Python, moderate chemistry

**Goals:**
- Prepare clean datasets for model training
- Ensure molecules are featurizable
- Identify dataset quality issues early

**Pain Points:**
- Days lost to data cleaning
- Models fail silently on bad molecules
- Hard to explain data issues to PI

**Feature Priorities:**
1. ML-readiness score
2. Automated filtering/cleaning
3. SDF/CSV export of clean data
4. Quality metrics summary

---

### 2.3 Primary Persona: Dr. Sarah Okonkwo - Natural Products Researcher

**Background:**
- Professor specializing in marine natural products
- Limited computational background
- Works with NP databases (COCONUT)

**Goals:**
- Verify NP authenticity in datasets
- Identify synthetic contaminants
- Analyze structural novelty

**Pain Points:**
- No easy tool for NP validation
- Undefined stereochemistry is pervasive
- Can't visualize what's wrong

**Feature Priorities:**
1. NP-likeness scoring
2. Stereochemistry assessment
3. Scaffold analysis
4. Simple, visual interface

---

### 2.4 Secondary Persona: James Park - Medicinal Chemist

**Background:**
- Industry chemist at biotech
- Uses commercial tools (Dotmatics)
- Needs quick checks for personal projects

**Goals:**
- Screen compound libraries for liabilities
- Check PAINS/structural alerts
- Validate vendor compounds

**Feature Priorities:**
1. PAINS/BRENK screening
2. Drug-likeness filters
3. Fast single-molecule checks
4. Mobile-friendly interface

---

## 3. Functional Requirements

### 3.1 Input Handling

| ID | Requirement | Priority |
|----|-------------|----------|
| F-IN-01 | Accept SMILES string input | P0 |
| F-IN-02 | Accept MOL/SDF file upload (single) | P0 |
| F-IN-03 | Accept SDF file upload (multi-molecule) | P0 |
| F-IN-04 | Accept CSV with SMILES column | P0 |
| F-IN-05 | Accept Excel with SMILES column | P1 |
| F-IN-06 | Accept InChI input | P1 |
| F-IN-07 | Accept drawn structure (Ketcher integration) | P2 |
| F-IN-08 | Accept PubChem/ChEMBL ID lookup | P2 |
| F-IN-09 | Accept compressed files (.zip, .gz) | P2 |
| F-IN-10 | Drag-and-drop file upload | P0 |

### 3.2 Validation Checks

#### 3.2.1 Basic Structural Validation

| ID | Check | Description | Priority |
|----|-------|-------------|----------|
| F-VAL-01 | Parsability | Can structure be parsed by RDKit | P0 |
| F-VAL-02 | Sanitization | Does structure pass RDKit sanitization | P0 |
| F-VAL-03 | Valence | Atoms have valid valence states | P0 |
| F-VAL-04 | Charge Balance | Overall charge is chemically reasonable | P1 |
| F-VAL-05 | Radical Check | Identify unexpected radicals | P1 |
| F-VAL-06 | Kekulization | Can aromatic bonds be kekulized | P1 |

#### 3.2.2 Stereochemistry Validation

| ID | Check | Description | Priority |
|----|-------|-------------|----------|
| F-STE-01 | Undefined Stereocenters | Count undefined chiral centers | P0 |
| F-STE-02 | Undefined Double Bonds | Count undefined E/Z bonds | P0 |
| F-STE-03 | Conflicting Stereo | Detect impossible stereo combinations | P1 |
| F-STE-04 | Excessive Stereocenters | Flag molecules with >10 stereocenters | P1 |
| F-STE-05 | Meso Compounds | Identify meso compounds | P2 |
| F-STE-06 | Pseudo-asymmetric | Detect pseudo-asymmetric centers | P2 |

#### 3.2.3 Structural Alerts

| ID | Check | Description | Priority |
|----|-------|-------------|----------|
| F-ALT-01 | PAINS A/B/C | All 480 PAINS patterns | P0 |
| F-ALT-02 | BRENK | Toxicity/unfavorable pharma | P0 |
| F-ALT-03 | NIH | Problematic functional groups | P1 |
| F-ALT-04 | ZINC | Drug-likeness unwanted groups | P1 |
| F-ALT-05 | ChEMBL Alerts | BMS, Dundee, Glaxo, Inpharmatica, LINT, MLSMR, SureChEMBL | P1 |
| F-ALT-06 | Custom Alerts | User-defined SMARTS patterns | P2 |

#### 3.2.4 Representation Validation

| ID | Check | Description | Priority |
|----|-------|-------------|----------|
| F-REP-01 | SMILES Roundtrip | SMILES → Mol → SMILES consistency | P0 |
| F-REP-02 | InChI Roundtrip | InChI → Mol → InChI consistency | P0 |
| F-REP-03 | InChI Generation | Can valid InChI be generated | P0 |
| F-REP-04 | Canonical Form | Generate canonical SMILES | P0 |
| F-REP-05 | Tautomer Enumeration | Enumerate possible tautomers | P1 |
| F-REP-06 | InChIKey Uniqueness | Generate InChIKey for deduplication | P1 |

#### 3.2.5 Natural Products Validation

| ID | Check | Description | Priority |
|----|-------|-------------|----------|
| F-NP-01 | NP-Likeness Score | Ertl score (-5 to +5) | P0 |
| F-NP-02 | NP Classification | Likely NP vs. synthetic | P1 |
| F-NP-03 | Sugar Moieties | Detect glycosylation | P1 |
| F-NP-04 | Scaffold Analysis | Murcko scaffold extraction | P1 |
| F-NP-05 | Biosynthetic Likelihood | Assess biosynthetic plausibility | P2 |

#### 3.2.6 ML-Readiness Assessment

| ID | Check | Description | Priority |
|----|-------|-------------|----------|
| F-ML-01 | Descriptor Calculability | Can standard descriptors be computed | P0 |
| F-ML-02 | Fingerprint Generation | Can fingerprints be generated | P0 |
| F-ML-03 | Size Check | MW, atom count within typical ML ranges | P1 |
| F-ML-04 | Element Check | Contains only common ML elements | P1 |
| F-ML-05 | Overall ML Score | Composite 0-100 readiness score | P0 |
| F-ML-06 | Feature Issues | List specific featurization problems | P1 |

### 3.3 Standardization Pipeline

| ID | Action | Description | Priority |
|----|--------|-------------|----------|
| F-STD-01 | Salt Stripping | Remove salt/counterion fragments | P0 |
| F-STD-02 | Solvent Removal | Remove solvent molecules | P0 |
| F-STD-03 | Charge Normalization | Standardize charge states | P0 |
| F-STD-04 | Tautomer Canonicalization | Convert to canonical tautomer | P1 |
| F-STD-05 | Isotope Removal | Strip isotope labels (optional) | P1 |
| F-STD-06 | Stereo Removal | Strip stereochemistry (optional) | P2 |
| F-STD-07 | Aromatization | Convert to aromatic form | P1 |
| F-STD-08 | Parent Extraction | Get parent structure | P0 |

### 3.4 Output & Reporting

| ID | Requirement | Priority |
|----|-------------|----------|
| F-OUT-01 | JSON validation result | P0 |
| F-OUT-02 | Visual summary dashboard | P0 |
| F-OUT-03 | CSV export | P0 |
| F-OUT-04 | Excel export with formatting | P1 |
| F-OUT-05 | PDF report generation | P1 |
| F-OUT-06 | SDF export (standardized structures) | P0 |
| F-OUT-07 | Structure images (PNG/SVG) | P1 |
| F-OUT-08 | Shareable result URL | P2 |

### 3.5 Batch Processing

| ID | Requirement | Priority |
|----|-------------|----------|
| F-BAT-01 | Process up to 100,000 molecules | P0 |
| F-BAT-02 | Progress indicator | P0 |
| F-BAT-03 | Background job processing | P0 |
| F-BAT-04 | Email notification on completion | P2 |
| F-BAT-05 | Pause/resume capability | P2 |
| F-BAT-06 | Parallel processing | P1 |
| F-BAT-07 | Partial results download | P1 |

---

## 4. Non-Functional Requirements

### 4.1 Performance

| ID | Requirement | Target |
|----|-------------|--------|
| NF-PERF-01 | Single molecule validation | <3 seconds |
| NF-PERF-02 | Batch (1,000 molecules) | <60 seconds |
| NF-PERF-03 | Batch (10,000 molecules) | <10 minutes |
| NF-PERF-04 | File upload (100MB) | <30 seconds |
| NF-PERF-05 | UI responsiveness | <100ms interaction |
| NF-PERF-06 | API response time | <500ms (95th percentile) |

### 4.2 Reliability

| ID | Requirement | Target |
|----|-------------|--------|
| NF-REL-01 | Uptime | 99.5% monthly |
| NF-REL-02 | Data loss | Zero for completed jobs |
| NF-REL-03 | Job completion | 99% of submitted jobs complete |
| NF-REL-04 | Error recovery | Graceful degradation on partial failures |

### 4.3 Scalability

| ID | Requirement | Target |
|----|-------------|--------|
| NF-SCL-01 | Concurrent users | 100 simultaneous |
| NF-SCL-02 | Daily throughput | 1M molecules |
| NF-SCL-03 | Storage | 1 year of job history |
| NF-SCL-04 | Horizontal scaling | Support for multiple workers |

### 4.4 Usability

| ID | Requirement | Target |
|----|-------------|--------|
| NF-USE-01 | Time to first validation | <60 seconds from landing |
| NF-USE-02 | Error clarity | All errors have actionable messages |
| NF-USE-03 | Accessibility | WCAG 2.1 AA compliance |
| NF-USE-04 | Mobile support | Responsive design for tablets |
| NF-USE-05 | Browser support | Chrome, Firefox, Safari, Edge (latest 2 versions) |

---

## 5. User Interface Requirements

### 5.1 Pages

| Page | Description | Priority |
|------|-------------|----------|
| Landing | Hero, quick input, feature overview | P0 |
| Single Validation | Input + results for one molecule | P0 |
| Batch Validation | Upload + job management + results | P0 |
| Job History | Past validation jobs | P1 |
| Documentation | Help, tutorials, API docs | P1 |
| Settings | User preferences | P2 |

### 5.2 Core UI Components

#### 5.2.1 Structure Input
- SMILES text input with validation
- File drag-and-drop zone
- Format auto-detection
- Preview of parsed structure

#### 5.2.2 Validation Results Display
- Overall score (0-100) with color coding
- Category breakdown (collapsible sections)
- Issue severity indicators (Critical/Warning/Info)
- Atom highlighting for specific issues
- Original vs. Standardized comparison

#### 5.2.3 Batch Results Table
- Sortable columns
- Filterable by issue type/severity
- Row selection for bulk actions
- Pagination (client-side for <1000, server-side for more)
- Export selected rows

### 5.3 Visual Design Principles

1. **Clean and Professional**: Scientific tool aesthetic, not flashy
2. **Information Dense**: Show relevant data without scrolling
3. **Color Coding**: Consistent severity colors (Red/Orange/Yellow/Green)
4. **Accessibility**: High contrast, keyboard navigation
5. **Dark Mode**: Optional dark theme

---

## 6. API Requirements

### 6.1 REST Endpoints

```
POST   /api/v1/validate          # Single molecule validation
POST   /api/v1/validate/batch    # Start batch job
GET    /api/v1/jobs/{job_id}     # Job status
GET    /api/v1/jobs/{job_id}/results    # Job results
DELETE /api/v1/jobs/{job_id}     # Cancel/delete job
POST   /api/v1/standardize       # Standardize molecule(s)
GET    /api/v1/alerts            # List available alert sets
GET    /api/v1/health            # Health check
```

### 6.2 Request/Response Formats

```json
// POST /api/v1/validate request
{
  "molecule": "CCO",
  "format": "smiles",
  "checks": ["all"],
  "options": {
    "standardize": true,
    "alert_sets": ["PAINS", "BRENK"]
  }
}

// Response
{
  "status": "completed",
  "molecule_id": "abc123",
  "input_smiles": "CCO",
  "canonical_smiles": "CCO",
  "inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
  "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
  "overall_score": 98,
  "ml_readiness_score": 100,
  "np_likeness_score": -1.2,
  "issues": [],
  "standardized_smiles": "CCO",
  "checks_performed": [...],
  "execution_time_ms": 45
}
```

### 6.3 Rate Limiting

| Tier | Single/min | Batch jobs/hour |
|------|------------|-----------------|
| Anonymous | 10 | 1 |
| Registered | 60 | 10 |
| API Key | 300 | 50 |

---

## 7. Data Requirements

### 7.1 Data Storage

| Data Type | Storage | Retention |
|-----------|---------|-----------|
| Job metadata | PostgreSQL | 1 year |
| Validation results | PostgreSQL (JSONB) | 1 year |
| Uploaded files | Object storage | 7 days |
| Generated reports | Object storage | 30 days |
| Alert definitions | PostgreSQL | Permanent |

### 7.2 Data Export

All user data must be exportable in standard formats:
- Validation results: JSON, CSV, Excel
- Structures: SMILES, SDF
- Reports: PDF

---

## 8. Security Requirements

### 8.1 Data Protection

| ID | Requirement |
|----|-------------|
| SEC-01 | No storage of molecular structures beyond job completion (+7 days) |
| SEC-02 | HTTPS only |
| SEC-03 | Input sanitization for all file uploads |
| SEC-04 | No execution of uploaded code |
| SEC-05 | Rate limiting to prevent abuse |

### 8.2 Authentication (Phase 2)

| ID | Requirement |
|----|-------------|
| SEC-06 | Optional user accounts |
| SEC-07 | API key authentication |
| SEC-08 | OAuth integration (GitHub, Google) |

---

## 9. Performance Requirements

### 9.1 Benchmarks

The following benchmark molecules should validate within specified times:

| Molecule Type | Example | Max Time |
|---------------|---------|----------|
| Simple organic | Ethanol | 500ms |
| Drug-like | Aspirin | 1s |
| Complex NP | Taxol | 3s |
| Large (MW>1000) | Cyclosporin | 5s |

### 9.2 Batch Processing Benchmarks

| Batch Size | Max Processing Time |
|------------|---------------------|
| 100 | 30 seconds |
| 1,000 | 2 minutes |
| 10,000 | 15 minutes |
| 100,000 | 2 hours |

---

## 10. Integration Requirements

### 10.1 External Service Integration

| Service | Integration Type | Priority |
|---------|------------------|----------|
| PubChem | Compound lookup | P2 |
| ChEMBL | Compound lookup | P2 |
| COCONUT | Database validation | P1 |
| DECIMER | OCSR validation | P1 |
| Cheminformatics Microservice | Backend module | P1 |

### 10.2 Export Integrations

| Format | Standard | Priority |
|--------|----------|----------|
| SDF | V2000/V3000 | P0 |
| SMILES | Canonical | P0 |
| InChI | v1.06 | P0 |
| CSV | RFC 4180 | P0 |
| Excel | XLSX | P1 |
| PDF | PDF/A | P1 |

---

## Appendix A: Glossary

| Term | Definition |
|------|------------|
| PAINS | Pan-Assay Interference Compounds - structures that give false positives in assays |
| NP-Likeness | Score indicating how "natural product-like" a structure is |
| InChI | International Chemical Identifier - unique chemical structure identifier |
| SMARTS | SMILES Arbitrary Target Specification - pattern matching language |
| Tautomer | Isomers that readily interconvert via hydrogen shift |
| Stereocenters | Atoms where substituent arrangement creates stereoisomers |

---

## Appendix B: Validation Scoring System

### Severity Levels

| Level | Score Impact | Color | Description |
|-------|--------------|-------|-------------|
| Critical | -50 | Red | Structure is invalid/unusable |
| Error | -20 | Orange | Significant issue, may affect downstream |
| Warning | -5 | Yellow | Minor issue, should review |
| Info | 0 | Blue | FYI, no action needed |
| Pass | 0 | Green | Check passed |

### Overall Score Calculation

```
base_score = 100
for issue in issues:
    base_score += issue.score_impact
overall_score = max(0, base_score)
```

### ML-Readiness Score Components

| Component | Weight | Criteria |
|-----------|--------|----------|
| Parsability | 30% | Successfully parsed |
| Sanitization | 20% | Passes RDKit sanitization |
| Descriptors | 20% | All standard descriptors calculable |
| Fingerprints | 15% | Morgan FP generates without error |
| Size/Elements | 15% | Within typical ML ranges |

---

## Appendix C: Acceptance Criteria Examples

### AC-1: Single Molecule Validation

**Given** a valid SMILES "CCO"  
**When** user submits for validation  
**Then** result shows:
- Overall score ≥ 95
- No critical/error issues
- Valid canonical SMILES
- Valid InChI/InChIKey
- Response time < 3 seconds

### AC-2: PAINS Detection

**Given** a known PAINS compound (rhodanine: "O=C1NC(=S)SC1")  
**When** user submits with PAINS alerts enabled  
**Then** result shows:
- PAINS alert flagged
- Pattern name displayed
- Highlighted atoms in structure
- Severity = Warning

### AC-3: Batch Processing

**Given** an SDF file with 1000 molecules (10% with issues)  
**When** user uploads for batch validation  
**Then**:
- Job starts within 5 seconds
- Progress updates every 10 seconds
- Completes within 2 minutes
- Results show ~100 flagged molecules
- CSV export available

---

*Document maintained by: Project Team*  
*Review schedule: Bi-weekly during development*
