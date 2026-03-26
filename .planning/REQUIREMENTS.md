# Requirements: ChemAudit v3.0

**Defined:** 2026-03-26
**Core Value:** Every chemical structure submitted gets a thorough, transparent, and reproducible quality assessment — from basic validity through ML-readiness — so scientists can trust their molecular data.

## v3.0 Requirements

Requirements for v3.0 release. 39 features across 8 phases. All existing v2.0 functionality preserved.

### Compound Profiling — Core Metrics

- [ ] **PROF-01**: User can view PFI (Property Forecast Index = cLogP + #AromaticRings) with low/moderate/high risk classification for any molecule
- [ ] **PROF-02**: User can view #stars count (properties outside 95th-percentile drug ranges) with per-property detail showing which ranges were exceeded
- [ ] **PROF-03**: User can view Abbott Bioavailability Score (4-class probability: 11%, 17%, 56%, 85%) based on TPSA and Lipinski violations
- [ ] **PROF-04**: User can view Consensus LogP (average of Wildman-Crippen + XLOGP3 approximation) showing each method's individual value
- [ ] **PROF-05**: User can view Skin Permeation (Potts-Guy log Kp) with low/moderate/high classification

### Compound Profiling — Advanced

- [ ] **PROF-06**: User can view 3D Shape Descriptors (PMI/NPR/PBF) with rod/disc/sphere classification, gracefully handling conformer generation failures
- [ ] **PROF-07**: User can configure and compute Custom MPO scores with desirability functions per property, including CNS MPO preset (0-6 scale, Wager 2010)
- [ ] **PROF-08**: User can compute Extended Ligand Efficiency metrics (LLE, LELP, SEI) from supplied activity data, building on existing LE/BEI
- [ ] **PROF-09**: User can view SA Comparison showing SA_Score alongside SCScore and SYBA scores, with graceful fallback for unavailable scorers

### Enhanced Structural Alerts

- [ ] **ALERT-01**: User can screen molecules against 21 custom SMARTS patterns (reactive warheads, metal chelators, redox-active, fluorescent interferents, phospholipidosis risk) with matched atom highlighting
- [ ] **ALERT-02**: User can screen molecules against 29 Kazius mutagenicity toxicophores with per-pattern match details and affected atom indices
- [ ] **ALERT-03**: User can screen molecules against NIBR Novartis screening deck filters (via medchem library) with source attribution
- [ ] **ALERT-04**: User can view complexity percentile assessment flagging molecules outside 5th-95th percentile boundaries vs commercial compound distributions

### Safety Flags

- [ ] **SAFE-01**: User can view CYP metabolism soft-spot predictions (SMARTS-based) with affected atom indices for highlighting and reaction type annotation
- [ ] **SAFE-02**: User can view hERG liability assessment (rule-based amphiphile check) with risk score and contributing factors
- [ ] **SAFE-03**: User can evaluate molecules against bRo5 (beyond-Rule-of-5) oral drug space for macrocycles and PROTACs (MW > 500 only)
- [ ] **SAFE-04**: User can evaluate molecules against REOS (Rapid Elimination of Swill) filter with per-property violation detail

### Structure Quality Diagnostics

- [ ] **DIAG-01**: User receives position-specific SMILES error diagnostics with fix suggestions when parsing fails, instead of generic "invalid SMILES"
- [ ] **DIAG-02**: User can compare two InChI strings layer-by-layer (formula, connections, hydrogens, charge, stereo, isotope) and see exactly which layers differ
- [ ] **DIAG-03**: User can check format round-trip lossiness (SMILES→InChI→SMILES, SMILES→MOL→SMILES) with specific identification of lost stereo, charge, or isotope information
- [ ] **DIAG-04**: User can compare standardization output across 3 pipelines (RDKit MolStandardize, ChEMBL-style, minimal sanitize) and see where they disagree
- [ ] **DIAG-05**: User can pre-validate SDF and CSV files for structural issues (missing M END, malformed counts lines, encoding problems) before molecule parsing

### QSAR-Ready Pipeline

- [ ] **QSAR-01**: User can process molecules through a 10-step configurable curation pipeline (parse→metals→desalt→normalize→neutralize→tautomer→stereo→isotope→filter→canonical) with per-step provenance
- [ ] **QSAR-02**: User can upload CSV/SDF files for batch QSAR-ready processing with InChIKey-based deduplication and summary statistics (ok/rejected/duplicate/error counts)
- [ ] **QSAR-03**: User can select from preset pipeline configurations (QSAR-2D with stereo stripping, QSAR-3D preserving stereo, custom with per-step toggles)

### Generative Chemistry Filter

- [ ] **GCHEM-01**: User can process SMILES lists from generative models through a multi-stage funnel (parse→valence→alerts→property rules→SA threshold→dedup) with per-stage counts and rejection reasons
- [ ] **GCHEM-02**: User can score SMILES via a REINVENT-compatible REST endpoint that returns 0-1 scores (composite of validity, drug-likeness, alert-free, SA)
- [ ] **GCHEM-03**: User can select from 4 preset filter configurations (drug-like, lead-like, fragment-like, permissive) with different thresholds per stage

### Dataset Intelligence

- [ ] **DSET-01**: User can upload a dataset and receive a composite 0-100 health score from 5 sub-scores (parsability, stereo completeness, uniqueness, alert prevalence, standardization consistency)
- [ ] **DSET-02**: User can detect contradictory labels (same InChIKey with >10-fold activity difference) with affected entries and fold-difference detail
- [ ] **DSET-03**: User can diff two CSV/SDF uploads showing added/removed/modified molecules identified by InChIKey
- [ ] **DSET-04**: User can download a reproducible curation report (serializable JSON) recording every standardization, filter, and dedup step applied

### Batch Analytics Extensions

- [ ] **BEXT-01**: User can cluster batch molecules using Butina algorithm (Tanimoto distance, configurable cutoff) with cluster IDs, sizes, and singleton count
- [ ] **BEXT-02**: User can compare any two molecules via MCS (maximum common substructure) with SMARTS, atom/bond counts, property deltas, and Tanimoto similarity
- [ ] **BEXT-03**: User can view ClassyFire-style chemical taxonomy classification for molecules using SMARTS-based chemotype rules covering ~50 drug-relevant categories
- [ ] **BEXT-04**: User can compute Registration Hash (RDKit RegistrationHash) for compound uniqueness tracking with explicit tautomer hash version control

### Ecosystem & Workflow Integration

- [ ] **ECO-01**: ChemAudit exposes validation, scoring, standardization, and integration tools as an MCP server (via fastapi-mcp) with tag-based allowlisting excluding admin endpoints
- [ ] **ECO-02**: User can validate and score molecules from the terminal via a CLI tool (chemaudit-cli) wrapping the Python client library
- [ ] **ECO-03**: User can check patent literature presence for any molecule via SureChEMBL InChIKey lookup

## Future Requirements

Deferred beyond v3.0. Tracked for potential inclusion in later milestones.

### Lilly Medchem Rules
- **FUT-01**: Lilly 275-rule demerit scoring (requires C++ compilation or LillyMol pybind11 build)

### RAscore
- **FUT-02**: Retrosynthetic Accessibility Score (requires Python 3.7 + TensorFlow 2.5)

### Extended SA Comparison
- **FUT-03**: RAscore integration into SA Comparison panel (blocked by FUT-02)

### Matched Molecular Series
- **FUT-04**: Extend existing MMP pairs into SAR chains (deferred — lower priority than other batch extensions)

### Jupyter Notebook Export
- **FUT-05**: Generate .ipynb templates from validation results

## Out of Scope

| Feature | Reason |
|---------|--------|
| Lilly Medchem Rules (v3.0) | Requires C++ compilation + Ruby driver; no Python package; ~1,700 patterns from RDKit + medchem + Kazius provide comparable coverage |
| RAscore | Requires Python 3.7 + TensorFlow 2.5; PyPI name collision; production-incompatible |
| External ML model APIs | All features use RDKit/MolVS/ChEMBL pipeline + rule-based methods only |
| 3D conformation export | Conformer generation is transient for PMI/NPR/PBF — not exported |
| Reaction validation | Molecules only, not reactions |
| User accounts/multi-tenancy | API key auth is sufficient |
| UMAP dimensionality reduction | Heavy C++ dependency |

## Traceability

Validated during roadmap creation (2026-03-26). All 39 requirements mapped to Phases 7-14.

| Requirement | Phase | Status |
|-------------|-------|--------|
| PROF-01 | Phase 7 | Pending |
| PROF-02 | Phase 7 | Pending |
| PROF-03 | Phase 7 | Pending |
| PROF-04 | Phase 7 | Pending |
| PROF-05 | Phase 7 | Pending |
| PROF-06 | Phase 7 | Pending |
| PROF-07 | Phase 7 | Pending |
| PROF-08 | Phase 7 | Pending |
| PROF-09 | Phase 7 | Pending |
| ALERT-01 | Phase 8 | Pending |
| ALERT-02 | Phase 8 | Pending |
| ALERT-03 | Phase 8 | Pending |
| ALERT-04 | Phase 8 | Pending |
| SAFE-01 | Phase 8 | Pending |
| SAFE-02 | Phase 8 | Pending |
| SAFE-03 | Phase 8 | Pending |
| SAFE-04 | Phase 8 | Pending |
| DIAG-01 | Phase 9 | Pending |
| DIAG-02 | Phase 9 | Pending |
| DIAG-03 | Phase 9 | Pending |
| DIAG-04 | Phase 9 | Pending |
| DIAG-05 | Phase 9 | Pending |
| QSAR-01 | Phase 10 | Pending |
| QSAR-02 | Phase 10 | Pending |
| QSAR-03 | Phase 10 | Pending |
| GCHEM-01 | Phase 11 | Pending |
| GCHEM-02 | Phase 11 | Pending |
| GCHEM-03 | Phase 11 | Pending |
| DSET-01 | Phase 12 | Pending |
| DSET-02 | Phase 12 | Pending |
| DSET-03 | Phase 12 | Pending |
| DSET-04 | Phase 12 | Pending |
| BEXT-01 | Phase 13 | Pending |
| BEXT-02 | Phase 13 | Pending |
| BEXT-03 | Phase 13 | Pending |
| BEXT-04 | Phase 13 | Pending |
| ECO-01 | Phase 14 | Pending |
| ECO-02 | Phase 14 | Pending |
| ECO-03 | Phase 14 | Pending |

**Coverage:**
- v3.0 requirements: 39 total
- Mapped to phases: 39
- Unmapped: 0 ✓

---
*Requirements defined: 2026-03-26*
*Last updated: 2026-03-26 after roadmap creation — traceability validated, 39/39 requirements mapped*
