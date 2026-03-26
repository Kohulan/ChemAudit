# Roadmap: ChemAudit

## Milestones

- ✅ **v1.0 Production-Ready Validation Platform** — Phases 1-5 (shipped 2026-01-29)
- ✅ **v2.0 Analytics, Scoring & Workflow** — Phases 1-6 (shipped 2026-02-24)
- ✅ **v2.1 Identifier Resolution & Cross-DB Comparison** — No new phases (shipped 2026-03-26)
- 📋 **v3.0 Advanced Profiling, Safety Intelligence & Dataset Curation** — Phases 7-14 (in progress)

## Phases

<details>
<summary>✅ v2.0 Analytics, Scoring & Workflow (Phases 1-6) — SHIPPED 2026-02-24</summary>

- [x] **Phase 1: Deep Validation** - 17 new validation checks covering stereo, tautomers, composition, and structural complexity (completed 2026-02-23)
- [x] **Phase 2: Standardization Intelligence** - Provenance tracking for all 4 ChEMBL pipeline stages with atom-level change reporting (completed 2026-02-23)
- [x] **Phase 3: Batch Analytics** - Multi-level deduplication, scaffold analysis, chemical space, MMP detection, and statistics (completed 2026-02-23)
- [x] **Phase 4: Scoring Expansion** - Drug-likeness profiles, property breakdowns, bioavailability radar, and BOILED-Egg plot (completed 2026-02-23)
- [x] **Phase 5: Visualizations** - Batch visualization suite and single-molecule deep view (completed 2026-02-24)
- [x] **Phase 6: Export, API & Workflow** - Advanced exports, custom profiles, audit trail, webhooks, and IUPAC input (completed 2026-02-24)

</details>

### 📋 v3.0 Advanced Profiling, Safety Intelligence & Dataset Curation (In Progress)

**Milestone Goal:** Transform ChemAudit from a validation/scoring suite into a comprehensive compound profiling and dataset curation platform — adding compound profiling metrics, enhanced safety intelligence, structure quality diagnostics, a QSAR-ready pipeline, generative chemistry filtering, dataset intelligence, batch analytics extensions, and ecosystem integration.

- [ ] **Phase 7: Compound Profiling Engine** - 9 new profiling metrics (PFI, bioavailability, consensus LogP, skin permeation, 3D shape, extended LE, SA comparison, custom MPO) on a new CompoundProfiler page
- [ ] **Phase 8: Enhanced Structural Alerts & Safety** - 8 new alert and safety features (custom SMARTS, Kazius toxicophores, NIBR filters, CYP soft-spots, hERG, bRo5, REOS, complexity percentile)
- [ ] **Phase 9: Structure Quality Diagnostics** - 5 diagnostic tools (SMILES error position, InChI layer diff, round-trip lossiness, cross-pipeline comparison, file pre-validation)
- [ ] **Phase 10: QSAR-Ready Pipeline** - 10-step configurable curation pipeline with batch file upload and preset configurations
- [ ] **Phase 11: Generative Chemistry Filter** - 6-stage funnel, REINVENT-compatible scoring API, and 4 preset configurations
- [ ] **Phase 12: Dataset Intelligence** - Composite dataset health score, contradictory label detection, dataset diff, and reproducible curation report
- [ ] **Phase 13: Batch Analytics Extensions** - Butina clustering, MCS comparator, chemical taxonomy classification, and registration hash
- [ ] **Phase 14: Ecosystem & Workflow Integration** - MCP server, CLI tool, and SureChEMBL patent check

## Phase Details

### Phase 7: Compound Profiling Engine
**Goal**: Scientists can view a comprehensive compound profile for any molecule — covering PFI risk classification, bioavailability probability, consensus LogP, skin permeation, 3D shape descriptors, extended ligand efficiency metrics, synthesizability comparison (SA/SCScore/SYBA), and user-configurable MPO scoring — all on a new CompoundProfiler page that establishes the shared UI primitives reused throughout v3.0
**Depends on**: Nothing (fully independent; establishes shared frontend components and vendored SCScore/SYBA isolation used by Phases 8 and 11)
**Requirements**: PROF-01, PROF-02, PROF-03, PROF-04, PROF-05, PROF-06, PROF-07, PROF-08, PROF-09
**Success Criteria** (what must be TRUE):
  1. User can enter a SMILES on the CompoundProfiler page and see PFI score with low/moderate/high risk class, #stars outlier count with per-property breakdown, and Abbott bioavailability class (11/17/56/85%) in a single panel
  2. User can view consensus LogP (Wildman-Crippen + XLOGP3 approximation) with both method values, and Skin Permeation (Potts-Guy log Kp) with low/moderate/high classification
  3. User can view 3D shape descriptors (PMI1/2/3, NPR1/2, PBF) with rod/disc/sphere classification on a PMI ternary plot; conformer failures show a graceful fallback message instead of an error
  4. User can configure a Custom MPO with per-property desirability functions (sigmoid/ramp/step), select the CNS MPO preset (Wager 2010, 0-6 scale), and compute the composite score; user can enter an activity value (pIC50/pKd) to compute LLE, LELP, and SEI alongside existing LE/BEI
  5. User can view SA Score, SCScore, and SYBA side-by-side in a synthesizability comparison panel; if SYBA subprocess is unavailable, the panel shows SA Score and SCScore with a SYBA-unavailable notice rather than failing
**Plans**: TBD
**UI hint**: yes

### Phase 8: Enhanced Structural Alerts & Safety
**Goal**: Scientists can screen any molecule against an expanded and categorized alert library — 21 custom SMARTS patterns, 29 Kazius mutagenicity toxicophores, 329 NIBR screening rules, and a complexity percentile filter — alongside new safety flags for CYP soft-spots, hERG liability, bRo5 eligibility, and REOS criteria, all surfaced through a unified `/api/v1/alerts/screen` endpoint with matched atom highlighting
**Depends on**: Nothing (fully independent; medchem installed in Phase 7 is the only shared dependency)
**Requirements**: ALERT-01, ALERT-02, ALERT-03, ALERT-04, SAFE-01, SAFE-02, SAFE-03, SAFE-04
**Success Criteria** (what must be TRUE):
  1. User can screen a molecule against the 21 custom SMARTS alert categories and see which patterns matched, with matched atom indices for highlighting and pattern descriptions
  2. User can see Kazius 29 mutagenicity toxicophore results grouped by concern group (not raw alert count) so a nitro group generates one concern entry, not multiple overlapping alerts from PAINS/Brenk/Kazius simultaneously
  3. User can see NIBR Novartis filter results with source attribution, CYP soft-spot predictions with affected atom indices and reaction type, and hERG liability risk score with contributing factors
  4. User can evaluate any molecule against bRo5 oral drug space (only shown for MW > 500) and REOS filter with per-property violation detail showing the threshold exceeded
  5. User can view a complexity percentile score flagging molecules outside the 5th-95th percentile vs commercial compound distributions
**Plans**: TBD
**UI hint**: yes

### Phase 9: Structure Quality Diagnostics
**Goal**: Scientists get actionable, specific feedback on malformed or ambiguous structures — position-specific SMILES error messages with fix suggestions, InChI layer-by-layer diff for two strings, format round-trip lossiness detection with specific lost-information fields, cross-pipeline standardization comparison across three algorithms, and a file pre-validator that catches SDF/CSV structural issues before parsing begins
**Depends on**: Nothing (fully independent; pure RDKit, no new dependencies)
**Requirements**: DIAG-01, DIAG-02, DIAG-03, DIAG-04, DIAG-05
**Success Criteria** (what must be TRUE):
  1. Submitting an invalid SMILES returns position-specific error information (character position, error type, description) and at least one fix suggestion rather than a generic "invalid SMILES" message
  2. User can compare two InChI strings and see which layers differ (formula, connections, hydrogens, charge, stereo, isotope) with the specific values from each string shown side by side
  3. User can run a round-trip check (SMILES→InChI→SMILES, SMILES→MOL→SMILES) and see exactly which information was lost — stereo centers, charge, or isotope labels — with before/after atom counts
  4. User can compare standardization output from three pipelines (RDKit MolStandardize, ChEMBL-style, minimal sanitize) and identify where they disagree on the canonical form
  5. User can upload an SDF or CSV and receive pre-validation results (missing M END blocks, malformed count lines, encoding issues) before molecule parsing runs, preventing parser hangs on malformed files
**Plans**: TBD
**UI hint**: yes

### Phase 10: QSAR-Ready Pipeline
**Goal**: Scientists can run any molecule or file through a configurable 10-step curation pipeline (parse, desalt, normalize, neutralize, tautomer, stereo, isotope, filter, canonical, dedup) with per-step provenance, InChIKey change tracking, and preset configurations for QSAR-2D, QSAR-3D, and custom workflows
**Depends on**: Phase 9 (file pre-validator must run before SDF parsing in batch upload; soft dependency — Phase 10 ships after Phase 9 to ensure safe file handling)
**Requirements**: QSAR-01, QSAR-02, QSAR-03
**Success Criteria** (what must be TRUE):
  1. User can process a single molecule through the 10-step pipeline and see which steps were applied, what changed at each step, and both the original InChIKey and the standardized InChIKey (with an `inchikey_changed` flag)
  2. User can upload a CSV or SDF file, select a preset (QSAR-2D, QSAR-3D, or custom), and receive a processed output with summary statistics showing counts of ok/rejected/duplicate/error molecules
  3. User can configure a custom pipeline by toggling individual steps on or off (per-step toggles) and the pipeline executes only the enabled steps in the canonical sequence
**Plans**: TBD
**UI hint**: yes

### Phase 11: Generative Chemistry Filter
**Goal**: Generative chemistry researchers can filter SMILES lists from generative models through a multi-stage funnel with per-stage rejection counts, score individual SMILES via a REINVENT-compatible REST endpoint returning 0-1 scores, and select from 4 preset filter configurations matching common generative model target spaces
**Depends on**: Phase 8 (alert stage is most complete with Kazius + NIBR patterns; functions with PAINS/Brenk alone if Phase 8 is delayed), Phase 7 (SA scorer for SA threshold stage)
**Requirements**: GCHEM-01, GCHEM-02, GCHEM-03
**Success Criteria** (what must be TRUE):
  1. User can submit a list of SMILES and see per-stage funnel counts (input → parse → valence → alerts → rules → SA threshold → dedup → output) with rejection reasons for each stage
  2. A POST to the REINVENT-compatible scoring endpoint with a SMILES list returns 0-1 scores where invalid SMILES return `null` (not `0.0`) and valid molecules return a composite score of validity, drug-likeness, alert-free status, and SA
  3. User can select from 4 preset configurations (drug-like, lead-like, fragment-like, permissive) with different property thresholds and SA score cutoffs pre-filled for each
**Plans**: TBD
**UI hint**: yes

### Phase 12: Dataset Intelligence
**Goal**: Scientists curating ML datasets can upload a compound file and receive a composite 0-100 health score across 5 sub-scores, detect contradictory activity labels for the same compound, diff two dataset versions by InChIKey, and download a reproducible JSON curation report recording every standardization, filter, and dedup decision
**Depends on**: Phase 10 (standardization consistency sub-score requires the QSAR pipeline; contradictory label detection and dataset diff are independent of Phase 10)
**Requirements**: DSET-01, DSET-02, DSET-03, DSET-04
**Success Criteria** (what must be TRUE):
  1. User can upload a dataset and see a composite 0-100 health score broken into 5 sub-scores: parsability rate, stereo completeness, uniqueness rate, alert prevalence, and standardization consistency
  2. User can detect contradictory labels by uploading a CSV with InChIKey and activity columns; molecules with the same InChIKey and more than 10-fold activity difference are returned with affected rows and fold-difference values
  3. User can upload two dataset files and receive a diff report showing molecules added, removed, and modified (identified by InChIKey), with the specific property changes for modified molecules
  4. User can download a reproducible JSON curation report that serializes every standardization step, filter applied, and dedup decision made during a session, enabling re-execution of the same workflow
**Plans**: TBD
**UI hint**: yes

### Phase 13: Batch Analytics Extensions
**Goal**: Batch results gain four new analytics capabilities: Butina compound clustering with configurable Tanimoto cutoff, pairwise MCS comparison of any two molecules with property deltas, ClassyFire-style chemical taxonomy classification covering ~50 curated drug-relevant chemotypes, and RDKit RegistrationHash per molecule for compound uniqueness tracking with version control
**Depends on**: Nothing (fully independent of other v3.0 phases; extends existing batch analytics infrastructure)
**Requirements**: BEXT-01, BEXT-02, BEXT-03, BEXT-04
**Success Criteria** (what must be TRUE):
  1. User can trigger Butina clustering on a batch result (capped at 1,000 molecules with an explicit API error above that) and see cluster IDs, cluster sizes, and singleton count with a configurable Tanimoto distance cutoff
  2. User can select any two molecules from a batch and run MCS comparison, receiving the maximum common substructure SMARTS, matched atom/bond counts, Tanimoto similarity, and a property delta table comparing the two molecules
  3. User can view chemical taxonomy classification for batch molecules showing assignment to chemotype categories (covering ~50 drug-relevant categories), with each category backed by a SMARTS rule match
  4. User can compute a RegistrationHash for any molecule — the hash includes the RDKit version used and uses `enable_tautomer_hash_v2=True` explicitly; meso compounds produce the correct hash verified against meso-tartaric acid
**Plans**: TBD
**UI hint**: yes

### Phase 14: Ecosystem & Workflow Integration
**Goal**: ChemAudit becomes accessible to LLM agents via a validated MCP server, to terminal users via a CLI tool, and to patent researchers via a SureChEMBL InChIKey lookup — delivering maximum ecosystem coverage with the full set of v3.0 features available as MCP tools at mount time
**Depends on**: Phases 7-13 (MCP server is most valuable when all v3.0 features are registered; CLI and SureChEMBL lookup are independent)
**Requirements**: ECO-01, ECO-02, ECO-03
**Success Criteria** (what must be TRUE):
  1. The MCP server mounts at `/mcp` and exposes validation, scoring, standardization, profiling, and integration tools as LLM-callable tools; a startup assertion confirms that admin endpoints (`/api/keys`, `/api/admin/config`) are absent from the MCP tool list
  2. User can run `chemaudit-cli validate --smiles "CCO"` from the terminal and receive the same structured validation result as the web UI, using the Python client library as the transport
  3. User can look up any molecule's patent presence via SureChEMBL InChIKey search and receive a binary present/absent result with a direct link to the SureChEMBL entry when found
**Plans**: TBD

## Progress

**Execution Order:**
Phases 7, 8, 9 are fully independent (run in sequence). Phase 10 runs after Phase 9. Phases 11 and 12 run after Phase 10. Phase 13 is independent. Phase 14 runs last.

Recommended sequence: 7 → 8 → 9 → 10 → 11 → 12 → 13 → 14

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 7. Compound Profiling Engine | v3.0 | 0/TBD | Not started | - |
| 8. Enhanced Structural Alerts & Safety | v3.0 | 0/TBD | Not started | - |
| 9. Structure Quality Diagnostics | v3.0 | 0/TBD | Not started | - |
| 10. QSAR-Ready Pipeline | v3.0 | 0/TBD | Not started | - |
| 11. Generative Chemistry Filter | v3.0 | 0/TBD | Not started | - |
| 12. Dataset Intelligence | v3.0 | 0/TBD | Not started | - |
| 13. Batch Analytics Extensions | v3.0 | 0/TBD | Not started | - |
| 14. Ecosystem & Workflow Integration | v3.0 | 0/TBD | Not started | - |
