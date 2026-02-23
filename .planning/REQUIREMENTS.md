# Requirements: ChemAudit v2.0

**Defined:** 2026-02-23
**Core Value:** Every chemical structure submitted gets a thorough, transparent, and reproducible quality assessment — from basic validity through ML-readiness — so scientists can trust their molecular data.

## v1 Requirements

Requirements for v2.0 release. Each maps to roadmap phases. 79 features + 1 prerequisite.

### Deep Validation — Stereo & Tautomer (Milestone 1.1)

- [x] **DVAL-01**: Undefined stereocenter detection with count and enumeration of possible stereoisomers
- [x] **DVAL-02**: Stereoisomer enumeration returning actual SMILES for all possible stereo forms
- [x] **DVAL-03**: Tautomer flagging with canonical form detection and tautomer count
- [x] **DVAL-04**: Aromatic system validation (Kekulization failures, unusual ring sizes, charged aromatics)
- [x] **DVAL-05**: Coordinate dimension check (2D/3D/degenerate detection)

### Deep Validation — Chemical Composition Guards (Milestone 1.2)

- [x] **DVAL-06**: Mixture detection (multi-component vs salt forms with fragment classification)
- [x] **DVAL-07**: Solvent contamination detection (common solvent fragment identification)
- [x] **DVAL-08**: Inorganic/organometallic filter (carbon-free molecule detection)
- [x] **DVAL-09**: Radical electron detection (unpaired electrons flagging)
- [x] **DVAL-10**: Isotope label detection (non-standard isotope identification)
- [x] **DVAL-11**: Empty/trivial molecule check (single atom, no bonds, too small)

### Deep Validation — Structural Complexity Flags (Milestone 1.3)

- [x] **DVAL-12**: Hypervalent atom detection (atoms exceeding normal valence)
- [x] **DVAL-13**: Polymer/repeating unit detection (SGroup markers, MW heuristic)
- [x] **DVAL-14**: Ring strain detection (3/4-membered rings with strain assessment)
- [x] **DVAL-15**: Macrocycle detection (rings > 12 atoms)
- [x] **DVAL-16**: Charged species flag with zwitterion identification
- [x] **DVAL-17**: Explicit hydrogen audit (unusual H count patterns)

### Standardization Intelligence — Provenance (Milestone 2.1)

- [x] **STD-01**: Canonical tautomer generation with provenance (input → canonical mapping)
- [x] **STD-02**: Neutralization report (atom-level charge changes with before/after)
- [x] **STD-03**: Functional group standardization audit (which groups were normalized and how)
- [x] **STD-04**: Parent extraction provenance (removed fragments with names, MW, structure)

### Standardization Intelligence — Normalization (Milestone 2.2)

- [x] **STD-05**: Kekulization/aromaticity normalization report (ring system changes)
- [x] **STD-06**: Stereochemistry normalization tracking (stereo changes during standardization)

### Batch Analytics — Multi-Level Duplicate Detection (Milestone 3.1)

- [ ] **BATCH-01**: Exact duplicate detection (canonical SMILES matching)
- [ ] **BATCH-02**: Tautomeric duplicate detection (canonical tautomer grouping)
- [ ] **BATCH-03**: Stereoisomer-insensitive deduplication (stereo-stripped InChI matching)
- [ ] **BATCH-04**: Salt-form duplicate detection (parent compound matching after desalting)

### Batch Analytics — Scaffold Analysis (Milestone 3.2)

- [ ] **BATCH-05**: Murcko scaffold decomposition and grouping
- [ ] **BATCH-06**: Generic scaffold grouping (element-agnostic framework)
- [ ] **BATCH-07**: Scaffold diversity metrics (unique scaffolds, Shannon entropy, scaffold frequency distribution)
- [ ] **BATCH-08**: R-group decomposition around common core

### Batch Analytics — Chemical Space & Similarity (Milestone 3.3)

- [ ] **BATCH-09**: Chemical space PCA/t-SNE visualization (2D projection of fingerprint space)
- [ ] **BATCH-10**: Similarity search within batch (query molecule → ranked neighbors)
- [ ] **BATCH-11**: Nearest neighbor analysis (isolation score per molecule)
- [ ] **BATCH-12**: Similarity matrix/heatmap data (pairwise Tanimoto)

### Batch Analytics — MMP & Activity Cliffs (Milestone 3.4)

- [ ] **BATCH-13**: Matched molecular pair detection (single-cut BRICS fragmentation)
- [ ] **BATCH-14**: Activity cliff detection (SALI index for structure-activity discontinuities)
- [ ] **BATCH-15**: Lipophilic ligand efficiency (LLE = pIC50 - LogP, requires activity column)

### Batch Analytics — Statistics & Quality (Milestone 3.5)

- [ ] **BATCH-16**: Property distribution statistics (mean, median, std, quartiles per property)
- [ ] **BATCH-17**: Property correlation matrix (pairwise Pearson correlations)
- [ ] **BATCH-18**: Batch quality score (composite metric: validity + diversity + drug-likeness)
- [ ] **BATCH-19**: Outlier detection (IQR-based flagging per property)

### Scoring Expansion — Drug-Likeness Profiles (Milestone 4.1)

- [ ] **SCORE-01**: Salt/counterion inventory report (identified fragments with names and MW)
- [ ] **SCORE-02**: Ligand efficiency (LE = -ΔG / heavy atom count)
- [ ] **SCORE-03**: Drug-likeness consensus score (Lipinski, Veber, Egan, Ghose, Muegge)
- [ ] **SCORE-04**: Lead-likeness assessment (MW 200-350, LogP -1 to 3, RotBonds ≤ 7)
- [ ] **SCORE-05**: Fragment-likeness Rule of 3 (MW < 300, LogP ≤ 3, HBD ≤ 3, HBA ≤ 3)

### Scoring Expansion — Property Breakdowns (Milestone 4.2)

- [ ] **SCORE-06**: Molecular complexity score (Bertz complexity index)
- [ ] **SCORE-07**: TPSA contribution breakdown (per-atom contributions)
- [ ] **SCORE-08**: LogP contribution breakdown (per-atom Crippen contributions)
- [ ] **SCORE-09**: Aggregator likelihood enhancements (extended aggregator detection)
- [ ] **SCORE-10**: Fsp3 visualization (sp3 carbon fraction breakdown)

### Scoring Expansion — NP & Specialized Scoring (Milestone 4.3)

- [ ] **SCORE-11**: NP-likeness breakdown (fragment contribution analysis)
- [ ] **SCORE-12**: Bioavailability radar (6-axis: LIPO, SIZE, POLAR, INSOLU, INSATU, FLEX)
- [ ] **SCORE-13**: BOILED-Egg plot (WLOGP vs TPSA with GI absorption/BBB permeation regions)
- [ ] **SCORE-14**: Property radar comparison (multi-molecule overlay)

### Visualizations — Batch Visualization Suite (Milestone 5.1)

- [ ] **VIZ-01**: Score distribution histogram (configurable bins, multiple score types)
- [ ] **VIZ-02**: Property vs property scatter plots (any two properties, color-coded)
- [ ] **VIZ-03**: Alert frequency bar chart (alert type counts across batch)
- [ ] **VIZ-04**: Validation issue treemap (category → issue → count hierarchy)
- [ ] **VIZ-05**: Scaffold treemap (scaffold frequency visualization)
- [ ] **VIZ-06**: Chemical space scatter plot (PCA/t-SNE with interactive tooltips)

### Visualizations — Single Molecule Deep View (Milestone 5.2)

- [ ] **VIZ-07**: Molecule comparison view (side-by-side property/score comparison)
- [ ] **VIZ-08**: Per-molecule property radar (individual molecule vs dataset average)
- [ ] **VIZ-09**: Batch timeline view (processing timeline with status indicators)

### Export, API & Workflow — Advanced Exports (Milestone 6.1)

- [ ] **WORK-01**: Fingerprint export (Morgan/MACCS/RDKit FP as CSV/numpy matrix)
- [ ] **WORK-02**: Batch deduplication export (deduplicated subset with group representatives)
- [ ] **WORK-03**: Scaffold-grouped export (molecules organized by scaffold family)
- [ ] **WORK-04**: Property matrix export (all computed properties as CSV/Excel)
- [ ] **WORK-05**: Batch summary PDF enhancement (analytics charts + statistics in PDF report)

### Export, API & Workflow — Custom Profiles & Filtering (Milestone 6.2)

- [ ] **WORK-06**: Custom scoring profiles (user-defined property thresholds and weights)
- [ ] **WORK-07**: Preset filter templates (drug-like, lead-like, fragment-like, CNS-penetrant)
- [ ] **WORK-08**: Molecule bookmarking (save/recall individual molecules across sessions)
- [ ] **WORK-09**: Batch subset actions (select → re-validate/re-score/export subsets)
- [ ] **WORK-10**: IUPAC name input support (IUPAC → SMILES conversion via py2opsin)

### Export, API & Workflow — Audit Trail & Notifications (Milestone 6.3)

- [ ] **WORK-11**: Shareable report permalinks (URL-encoded result sharing)
- [ ] **WORK-12**: Webhook on batch complete (configurable HTTP POST callback)
- [ ] **WORK-13**: Email on batch complete (SMTP notification)
- [ ] **WORK-14**: Validation history / audit trail (append-only log of all validations)

### Infrastructure Prerequisite

- [ ] **INFRA-01**: Batch result pagination and Redis TTL policy (prerequisite for Phase 3 analytics)

## v2 Requirements

No v2 requirements at this time. All 79 features are committed to v1.

## Out of Scope

| Feature | Reason |
|---------|--------|
| 3D conformation generation | Adds ETKDG/MMFF overhead, contradicts 2D validation focus |
| External ML model integration | Models go stale, outside RDKit/MolVS/ChEMBL pipeline constraint |
| User accounts / multi-tenancy | API key auth is sufficient; full RBAC is a product pivot |
| Real-time collaboration | Shareable permalinks cover the use case without WebSocket multiplexing |
| Reaction validation | Different domain, different tooling (molecules only) |
| AI/LLM molecule generation | Non-deterministic, incompatible with reproducibility requirement |
| UMAP dimensionality reduction | Heavy C++ dependency, outside "zero external ML models" constraint |

## Traceability

Finalized during roadmap creation — 2026-02-23.

| Requirement | Phase | Plan | Status |
|-------------|-------|------|--------|
| DVAL-01 | Phase 1 (M1.1) | 01-01 | Pending |
| DVAL-02 | Phase 1 (M1.1) | 01-01 | Pending |
| DVAL-03 | Phase 1 (M1.1) | 01-01 | Pending |
| DVAL-04 | Phase 1 (M1.1) | 01-01 | Pending |
| DVAL-05 | Phase 1 (M1.1) | 01-01 | Pending |
| DVAL-06 | Phase 1 (M1.2) | 01-02 | Pending |
| DVAL-07 | Phase 1 (M1.2) | 01-02 | Pending |
| DVAL-08 | Phase 1 (M1.2) | 01-02 | Pending |
| DVAL-09 | Phase 1 (M1.2) | 01-02 | Pending |
| DVAL-10 | Phase 1 (M1.2) | 01-02 | Pending |
| DVAL-11 | Phase 1 (M1.2) | 01-02 | Pending |
| DVAL-12 | Phase 1 (M1.3) | 01-03 | Pending |
| DVAL-13 | Phase 1 (M1.3) | 01-03 | Pending |
| DVAL-14 | Phase 1 (M1.3) | 01-03 | Pending |
| DVAL-15 | Phase 1 (M1.3) | 01-03 | Pending |
| DVAL-16 | Phase 1 (M1.3) | 01-03 | Pending |
| DVAL-17 | Phase 1 (M1.3) | 01-03 | Pending |
| STD-01 | Phase 2 (M2.1) | 02-01 | Complete |
| STD-02 | Phase 2 (M2.1) | 02-01 | Complete |
| STD-03 | Phase 2 (M2.1) | 02-01 | Complete |
| STD-04 | Phase 2 (M2.1) | 02-01 | Complete |
| STD-05 | Phase 2 (M2.2) | 02-02 | Complete |
| STD-06 | Phase 2 (M2.2) | 02-02 | Complete |
| INFRA-01 | Phase 3 (prerequisite) | 03-01 | Pending |
| BATCH-01 | Phase 3 (M3.1) | 03-02 | Pending |
| BATCH-02 | Phase 3 (M3.1) | 03-02 | Pending |
| BATCH-03 | Phase 3 (M3.1) | 03-02 | Pending |
| BATCH-04 | Phase 3 (M3.1) | 03-02 | Pending |
| BATCH-05 | Phase 3 (M3.2) | 03-03 | Pending |
| BATCH-06 | Phase 3 (M3.2) | 03-03 | Pending |
| BATCH-07 | Phase 3 (M3.2) | 03-03 | Pending |
| BATCH-08 | Phase 3 (M3.2) | 03-03 | Pending |
| BATCH-09 | Phase 3 (M3.3) | 03-04 | Pending |
| BATCH-10 | Phase 3 (M3.3) | 03-04 | Pending |
| BATCH-11 | Phase 3 (M3.3) | 03-04 | Pending |
| BATCH-12 | Phase 3 (M3.3) | 03-04 | Pending |
| BATCH-13 | Phase 3 (M3.4) | 03-05 | Pending |
| BATCH-14 | Phase 3 (M3.4) | 03-05 | Pending |
| BATCH-15 | Phase 3 (M3.4) | 03-05 | Pending |
| BATCH-16 | Phase 3 (M3.5) | 03-06 | Pending |
| BATCH-17 | Phase 3 (M3.5) | 03-06 | Pending |
| BATCH-18 | Phase 3 (M3.5) | 03-06 | Pending |
| BATCH-19 | Phase 3 (M3.5) | 03-06 | Pending |
| SCORE-01 | Phase 4 (M4.1) | 04-01 | Pending |
| SCORE-02 | Phase 4 (M4.1) | 04-01 | Pending |
| SCORE-03 | Phase 4 (M4.1) | 04-01 | Pending |
| SCORE-04 | Phase 4 (M4.1) | 04-01 | Pending |
| SCORE-05 | Phase 4 (M4.1) | 04-01 | Pending |
| SCORE-06 | Phase 4 (M4.2) | 04-02 | Pending |
| SCORE-07 | Phase 4 (M4.2) | 04-02 | Pending |
| SCORE-08 | Phase 4 (M4.2) | 04-02 | Pending |
| SCORE-09 | Phase 4 (M4.2) | 04-02 | Pending |
| SCORE-10 | Phase 4 (M4.2) | 04-02 | Pending |
| SCORE-11 | Phase 4 (M4.3) | 04-03 | Pending |
| SCORE-12 | Phase 4 (M4.3) | 04-03 | Pending |
| SCORE-13 | Phase 4 (M4.3) | 04-03 | Pending |
| SCORE-14 | Phase 4 (M4.3) | 04-03 | Pending |
| VIZ-01 | Phase 5 (M5.1) | 05-01 | Pending |
| VIZ-02 | Phase 5 (M5.1) | 05-01 | Pending |
| VIZ-03 | Phase 5 (M5.1) | 05-01 | Pending |
| VIZ-04 | Phase 5 (M5.1) | 05-01 | Pending |
| VIZ-05 | Phase 5 (M5.1) | 05-01 | Pending |
| VIZ-06 | Phase 5 (M5.1) | 05-01 | Pending |
| VIZ-07 | Phase 5 (M5.2) | 05-02 | Pending |
| VIZ-08 | Phase 5 (M5.2) | 05-02 | Pending |
| VIZ-09 | Phase 5 (M5.2) | 05-02 | Pending |
| WORK-01 | Phase 6 (M6.1) | 06-01 | Pending |
| WORK-02 | Phase 6 (M6.1) | 06-01 | Pending |
| WORK-03 | Phase 6 (M6.1) | 06-01 | Pending |
| WORK-04 | Phase 6 (M6.1) | 06-01 | Pending |
| WORK-05 | Phase 6 (M6.1) | 06-01 | Pending |
| WORK-06 | Phase 6 (M6.2) | 06-02 | Pending |
| WORK-07 | Phase 6 (M6.2) | 06-02 | Pending |
| WORK-08 | Phase 6 (M6.2) | 06-02 | Pending |
| WORK-09 | Phase 6 (M6.2) | 06-02 | Pending |
| WORK-10 | Phase 6 (M6.2) | 06-02 | Pending |
| WORK-11 | Phase 6 (M6.3) | 06-03 | Pending |
| WORK-12 | Phase 6 (M6.3) | 06-03 | Pending |
| WORK-13 | Phase 6 (M6.3) | 06-03 | Pending |
| WORK-14 | Phase 6 (M6.3) | 06-03 | Pending |

**Coverage:**
- v1 requirements: 80 total (79 features + 1 infrastructure prerequisite)
- Mapped to phases: 80
- Unmapped: 0

---
*Requirements defined: 2026-02-23*
*Last updated: 2026-02-23 — traceability finalized after roadmap creation*
