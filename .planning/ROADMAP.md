# Roadmap: ChemAudit v2.0

## Overview

ChemAudit v1.0 is production. This roadmap covers the v2.0 enhancement cycle: 80 requirements
(79 features + 1 infrastructure prerequisite) across 6 phases and 18 milestones. The build
extends a fully functional cheminformatics validation platform into a comprehensive analytics
suite — adding deep validation checks, standardization provenance, batch deduplication and
chemical space analysis, expanded scoring, interactive visualizations, and workflow integrations.
Every phase delivers integrated backend + frontend together; no partial feature ships.

## Milestones

- 📋 **v2.0 Deep Validation** — Phase 1 (in progress)
- 📋 **v2.0 Standardization Intelligence** — Phase 2 (planned)
- 📋 **v2.0 Batch Analytics** — Phase 3 (planned)
- 📋 **v2.0 Scoring Expansion** — Phase 4 (planned)
- 📋 **v2.0 Visualizations** — Phase 5 (planned)
- 📋 **v2.0 Export, API & Workflow** — Phase 6 (planned)

## Phases

**Phase Numbering:**
- Integer phases (1–6): Planned milestone work
- Decimal phases (e.g., 2.1): Urgent insertions (marked INSERTED)

- [x] **Phase 1: Deep Validation** - 17 new validation checks covering stereo, tautomers, composition, and structural complexity (completed 2026-02-23)
- [x] **Phase 2: Standardization Intelligence** - Provenance tracking for all 4 ChEMBL pipeline stages with atom-level change reporting (completed 2026-02-23)
- [x] **Phase 3: Batch Analytics** - Multi-level deduplication, scaffold analysis, chemical space, MMP detection, and statistics (preceded by INFRA-01 pagination) (completed 2026-02-23)
- [x] **Phase 4: Scoring Expansion** - Drug-likeness profiles, property breakdowns, bioavailability radar, and BOILED-Egg plot (completed 2026-02-23)
- [x] **Phase 5: Visualizations** - Batch visualization suite and single-molecule deep view built on Phase 3 and 4 data (completed 2026-02-24)
- [x] **Phase 6: Export, API & Workflow** - Advanced exports, custom profiles, audit trail, webhooks, and IUPAC input (completed 2026-02-24)

## Phase Details

### Phase 1: Deep Validation
**Goal**: Users see a richer, more informative validation result — 17 new checks covering stereo completeness, tautomer detection, chemical composition guards, and structural complexity flags, all delivered as plugin checks with zero changes to the existing engine
**Depends on**: Nothing (first phase, zero external dependencies)
**Requirements**: DVAL-01, DVAL-02, DVAL-03, DVAL-04, DVAL-05, DVAL-06, DVAL-07, DVAL-08, DVAL-09, DVAL-10, DVAL-11, DVAL-12, DVAL-13, DVAL-14, DVAL-15, DVAL-16, DVAL-17
**Success Criteria** (what must be TRUE):
  1. A molecule with undefined stereocenters returns a check result with the count of undefined centers and a list of enumerated stereoisomer SMILES
  2. A mixture input (e.g., a SMILES with a dot separator) is flagged as a mixture with each fragment classified as drug component, salt, or solvent
  3. Inorganic, radical, isotope-labeled, and trivially small molecules each return a distinct, named check failure with structured detail fields
  4. Hypervalent atoms, ring strain (3/4-membered rings), macrocycles (>12 atoms), and charged/zwitterionic species each produce a structured warning with affected atom indices
  5. Startup assertion confirms all 17 new check names are registered in CheckRegistry; a CI test cross-checks checks/ file count with registry count
**Plans**: 4 plans (3 backend milestones + 1 frontend integration)

Plans:
- [x] 01-01-PLAN.md — M1.1 Stereo & Tautomer Checks backend (DVAL-01..05) + cache key versioning (COMPLETE: 2 tasks, 31 tests, commits 51bf66e + 58e4b11)
- [x] 01-02-PLAN.md — M1.2 Chemical Composition Guards backend (DVAL-06..11) (COMPLETE: 2 tasks, 133 tests, commits 3b0eab7 + 7af522a)
- [x] 01-03-PLAN.md — M1.3 Structural Complexity Flags backend (DVAL-12..17) + startup assertion (COMPLETE: 2 tasks, 65 tests, commit ab0eb4c)
- [x] 01-04-PLAN.md — Frontend Deep Validation Tab (all DVAL requirements: types, UI, severity config, atom highlighting) (COMPLETE: 3 tasks, 41 tests)

### Phase 2: Standardization Intelligence
**Goal**: Users running standardization see exactly what changed and why — every stage of the ChEMBL pipeline (Checker, Standardizer, GetParent, Tautomer) produces a provenance record with atom-level before/after diffs, removed fragment names, and stereo change tracking
**Depends on**: Phase 1 (stereo and tautomer data feeds into provenance; DVAL-01/DVAL-03 inform STD-01/STD-06)
**Requirements**: STD-01, STD-02, STD-03, STD-04, STD-05, STD-06
**Success Criteria** (what must be TRUE):
  1. A standardized molecule response optionally includes a `provenance` field showing input SMILES, canonical tautomer SMILES, and a `stereo_stripped: true/false` flag
  2. The neutralization report lists each atom that changed charge, its before/after charge state, and the SMARTS pattern that triggered the change
  3. The parent extraction report names each removed fragment (by COUNTERION_NAMES dictionary lookup), its molecular weight, and its SMILES
  4. All provenance fields are additive to the existing StandardizationResponse — existing API consumers receive identical responses if they do not request provenance
**Plans**: 3 plans

Plans:
- [x] 02-01-PLAN.md — M2.1 Provenance Pipeline: ProvenancePipeline wrapper, fragment dictionary, schema extension, route integration, tests (STD-01, STD-02, STD-03, STD-04)
- [x] 02-02-PLAN.md — M2.2 Ring/Stereo Tracking + Frontend: ring aromaticity tracking, per-center stereo detail, ProvenanceTimeline UI, ProvenanceStageCard (STD-05, STD-06)
- [x] 02-03-PLAN.md — Gap closure: fix duplicate fragment dict key, populate DVAL cross-references in stereo provenance (all STD requirements)

### Phase 3: Batch Analytics
**Goal**: A completed batch job exposes a second analytics layer — multi-level deduplication groups, scaffold families, chemical space projections, MMP pairs, and statistical summaries — all computed asynchronously by a post-aggregation Celery chord, accessible via a dedicated analytics endpoint
**Depends on**: Phase 1 (stereo data for dedup), Phase 2 Milestone 2.1 (salt-form dedup reads excluded_fragments from provenance); INFRA-01 pagination must be delivered first (prerequisite)
**Requirements**: INFRA-01, BATCH-01, BATCH-02, BATCH-03, BATCH-04, BATCH-05, BATCH-06, BATCH-07, BATCH-08, BATCH-09, BATCH-10, BATCH-11, BATCH-12, BATCH-13, BATCH-14, BATCH-15, BATCH-16, BATCH-17, BATCH-18, BATCH-19
**Success Criteria** (what must be TRUE):
  1. `GET /api/v1/batch/{job_id}/analytics` returns a status field (pending/complete/failed) and, when complete, deduplication groups across all four levels (exact, tautomeric, stereo-insensitive, salt-form)
  2. Scaffold analysis returns Murcko scaffold SMILES per molecule, generic scaffold SMILES, scaffold frequency distribution, Shannon entropy diversity metric, and R-group decomposition results for a user-specified core
  3. Chemical space endpoint returns PCA 2D coordinates for all molecules and t-SNE coordinates for batches of ≤2000 molecules; similarity search returns ranked neighbors by Tanimoto for a query SMILES
  4. Property distribution statistics return mean, median, std, quartiles, and IQR-based outlier flags per property; batch quality score is a composite 0–100 metric combining validity, diversity, and drug-likeness
  5. All analytics endpoints refuse synchronous computation for large batches; analytics results are cached in Redis with 24-hour TTL; result_storage.get_results() supports page/page_size pagination
**Plans**: 6 plans (INFRA-01 + analytics infrastructure first, then 5 analytics service modules in parallel)

Plans:
- [x] 03-01-PLAN.md — INFRA-01 TTL policy + analytics infrastructure skeleton (storage, schemas, routes, Celery tasks)
- [x] 03-02-PLAN.md — Multi-Level Duplicate Detection: exact, tautomeric, stereo-insensitive, salt-form (BATCH-01..04)
- [x] 03-03-PLAN.md — Scaffold Analysis: Murcko/generic scaffolds, diversity metrics, R-group decomposition (BATCH-05..08)
- [x] 03-04-PLAN.md — Chemical Space & Similarity: PCA, t-SNE, similarity search, nearest neighbors, similarity matrix (BATCH-09..12)
- [x] 03-05-PLAN.md — MMP & Activity Cliffs: BRICS-based MMP, SALI index, LLE (BATCH-13..15)
- [x] 03-06-PLAN.md — Batch Statistics & Quality: property stats, correlations, quality score, outlier detection (BATCH-16..19)

### Phase 4: Scoring Expansion
**Goal**: Single-molecule scoring results include richer drug-likeness profiles (consensus rules, lead-like, fragment-like assessments), per-atom property contribution breakdowns, and visually interpretable radar and BOILED-Egg data — all computed from already-available descriptors, zero new backend dependencies except scoring logic
**Depends on**: Phase 2 Milestone 2.1 (SCORE-01 salt inventory reads excluded_fragments from standardization provenance); Milestones 4.2 and 4.3 are independent
**Requirements**: SCORE-01, SCORE-02, SCORE-03, SCORE-04, SCORE-05, SCORE-06, SCORE-07, SCORE-08, SCORE-09, SCORE-10, SCORE-11, SCORE-12, SCORE-13, SCORE-14
**Success Criteria** (what must be TRUE):
  1. Drug-likeness scoring returns a consensus pass/fail across five rule sets (Lipinski, Veber, Egan, Ghose, Muegge) with per-rule detail and an overall consensus score
  2. Lead-likeness and fragment-likeness (Rule of 3) assessments return structured pass/fail with the specific property values and thresholds that triggered each result
  3. TPSA and LogP contribution breakdowns return per-atom contribution values (atom index, contribution, total), enabling atom-level highlighting in the frontend
  4. The bioavailability radar returns six axis values (LIPO, SIZE, POLAR, INSOLU, INSATU, FLEX) normalized to 0–1 scale; BOILED-Egg returns WLOGP and TPSA values with GI/BBB region membership flags
**Plans**: 3 plans in 2 waves (M4.1 + M4.2 parallel, then M4.3 + frontend)

Plans:
- [x] 04-01-PLAN.md — M4.1 Drug-Likeness Profiles: consensus scoring, lead-likeness, salt inventory, ligand efficiency (SCORE-01..05) [Wave 1] (COMPLETE: 2 tasks, 30 tests)
- [x] 04-02-PLAN.md — M4.2 Property Breakdowns: TPSA/LogP per-atom, Bertz detail, Fsp3 visualization, aggregator extension (SCORE-06..10) [Wave 1] (COMPLETE: 2 tasks, 27 tests)
- [x] 04-03-PLAN.md — M4.3 NP Breakdown + Radars + Frontend: NP fragment analysis, bioavailability radar, BOILED-Egg, property comparison, Scoring Profiles tab (SCORE-11..14) [Wave 2] (COMPLETE: 2 tasks, 44 tests, 6 frontend components)

### Phase 5: Visualizations
**Goal**: Batch results and single-molecule scoring data become visually explorable — 9 new interactive chart components (Recharts-based, zero new npm packages) consuming Phase 3 analytics and Phase 4 scoring endpoints; all are pure presentational components with data fetched by parent containers
**Depends on**: Phase 3 Milestones 3.2, 3.3, 3.5 (scaffold, chemical space, and statistics data for batch viz); Phase 4 Milestones 4.1, 4.2, 4.3 (scoring data for single-molecule deep view)
**Requirements**: VIZ-01, VIZ-02, VIZ-03, VIZ-04, VIZ-05, VIZ-06, VIZ-07, VIZ-08, VIZ-09
**Success Criteria** (what must be TRUE):
  1. The batch results page renders a score distribution histogram, property-vs-property scatter plot, alert frequency bar chart, and validation issue treemap — all from the analytics API, all interactive with tooltips
  2. The scaffold treemap renders scaffold frequency (capped at 50 scaffolds, remainder as "Other") and the chemical space scatter plot renders PCA/t-SNE coordinates (sampled to 2000 visible points for render performance)
  3. The single-molecule detail view shows a side-by-side molecule comparison, a per-molecule property radar overlaid against dataset average, and a batch timeline view with processing status indicators
  4. The similarity matrix heatmap caps rendering at 100x100; canvas rendering is used for larger matrices — no browser crash at scale
**Plans**: 2 plans in 2 waves (M5.1 batch charts in Wave 1, then M5.2 comparison + timeline in Wave 2)

Plans:
- [x] 05-01-PLAN.md — M5.1 Batch Visualization Suite: analytics types, API methods, hooks, tabbed container, 6 chart components (VIZ-01..06) [Wave 1] (COMPLETE: 2 tasks, 12 files, commits a314515 + 3ee40c9)
- [x] 05-02-PLAN.md — M5.2 Single Molecule Deep View: molecule comparison panel, property radar, batch timeline, compare flow integration (VIZ-07..09) [Wave 2] (COMPLETE: 2 tasks, 5 files, commits eaec55f + 4292378)

### Phase 6: Export, API & Workflow
**Goal**: Power users can export richly annotated data in new formats, define custom scoring profiles, bookmark molecules across sessions, receive batch completion notifications, and input structures by IUPAC name — all self-contained features with backward-compatible API additions
**Depends on**: Phase 3 Milestones 3.1 and 3.2 (deduplication and scaffold groups for WORK-02/WORK-03); Phase 5 Milestone 5.1 (analytics charts for batch PDF enhancement in WORK-05); Milestones 6.2 and 6.3 are independent
**Requirements**: WORK-01, WORK-02, WORK-03, WORK-04, WORK-05, WORK-06, WORK-07, WORK-08, WORK-09, WORK-10, WORK-11, WORK-12, WORK-13, WORK-14
**Success Criteria** (what must be TRUE):
  1. Fingerprint export produces a CSV or numpy-compatible matrix of Morgan/MACCS/RDKit fingerprints; scaffold-grouped and deduplicated exports organize molecules by their Phase 3 analytics groupings; property matrix export produces all computed properties in a single CSV/Excel file
  2. The enhanced batch PDF includes analytics charts and statistics from Phase 3 alongside the existing validation summary
  3. Custom scoring profiles allow users to define property thresholds and weights, save them with a name, and apply them to score any molecule; preset filter templates (drug-like, lead-like, fragment-like, CNS-penetrant) are available without configuration
  4. Molecule bookmarks persist across sessions (PostgreSQL-backed); batch subset actions allow selecting molecules from a batch result and re-validating, re-scoring, or exporting the subset
  5. Webhook on batch complete sends an HTTP POST to a configured URL with 3 retries and exponential backoff; email notification uses async SMTP; validation audit trail stores all validation events in an append-only PostgreSQL table with paginated read access; IUPAC name input converts names to SMILES via py2opsin before validation
**Plans**: 8 plans in 3 waves (5 original + 3 gap closure)

Plans:
- [x] 06-01-PLAN.md — Advanced Export Formats: fingerprint, dedup, scaffold, property matrix exporters (WORK-01, WORK-02, WORK-03, WORK-04) [Wave 1] (COMPLETE: 2 tasks, 28 tests, commits 611ed15 + 1c29617)
- [x] 06-02-PLAN.md — ORM Foundation + IUPAC + PDF Enhancement: SQLAlchemy setup, all 4 DB models, Alembic, JPype OPSIN converter, PDF section selection (WORK-05, WORK-10) [Wave 1] (COMPLETE: 2 tasks, 22 tests, commits 0f80f36 + c8a0d57)
- [x] 06-03-PLAN.md — Custom Profiles + Bookmarks + Subset Actions: profiles CRUD with 8 presets, bookmarks with batch-submit, batch subset revalidate/rescore/export (WORK-06, WORK-07, WORK-08, WORK-09) [Wave 2] (COMPLETE: 28 tests, commit 21565f9)
- [x] 06-04-PLAN.md — Notifications + Permalinks + Audit Trail: HMAC webhook, SMTP email, shareable permalinks, append-only audit trail with history API (WORK-11, WORK-12, WORK-13, WORK-14) [Wave 2] (COMPLETE: 24 tests, commit de0c57f)
- [x] 06-05-PLAN.md — Frontend Integration: export dialog, profile builder, preset picker, bookmarks page, history page, IUPAC input badge, subset panel, permalink sharing, tab persistence (all WORK requirements) [Wave 3] (COMPLETE: 17 files, 8 frontend tests, commit 751af96)
- [ ] 06-06-PLAN.md — Gap Closure: Backend Fixes — PubChem REST fallback for IUPAC, batch audit trail wiring, webhook dispatch (WORK-10, WORK-12, WORK-14)
- [ ] 06-07-PLAN.md — Gap Closure: Frontend Wiring — Profiles page + route, SubsetActionPanel integration, Share permalink button (WORK-06, WORK-07, WORK-09, WORK-11)
- [ ] 06-08-PLAN.md — Gap Closure: UI Polish — PDF checkboxes position, input placeholder, bookmark button styling, tab result clearing (WORK-05, WORK-08, WORK-10, WORK-13)

## Progress

**Execution Order:**
Phases execute in dependency order: 1 → 2 → 3 (INFRA-01 first) → 4 (parallel with 3 where possible) → 5 → 6

**Recommended Build Order Within Phases (from Features/ROADMAP.md):**
M1.1 → M1.2 → M1.3 → M2.1 → M3.1 → M4.3 → M3.2 → M3.3 → M3.5 → M4.1 → M4.2 → M5.1 → M5.2 → M3.4 → M2.2 → M6.1 → M6.2 → M6.3

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. Deep Validation | 4/4 | Complete   | 2026-02-23 | - |
| 2. Standardization Intelligence | 3/3 | Complete   | 2026-02-23 | - |
| 3. Batch Analytics | 6/6 | Complete   | 2026-02-23 | - |
| 4. Scoring Expansion | 3/3 | Complete | 2026-02-23 | - |
| 5. Visualizations | 2/2 | Complete | 2026-02-24 | - |
| 6. Export, API & Workflow | 7/8 | In Progress|  | - |

**Total:** 23/26 plans complete (gap closure in progress)
