# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-23)

**Core value:** Every chemical structure submitted gets a thorough, transparent, and reproducible quality assessment — from basic validity through ML-readiness — so scientists can trust their molecular data.
**Current focus:** Phase 5 — Visualizations (plan 05-01 complete, executing 05-02)

## Current Position

Phase: 5 of 6 (Visualizations)
Plan: 1 of 2 in current phase
Status: Plan 05-01 complete — batch visualization suite with 6 interactive charts
Last activity: 2026-02-24 — Completed 05-01: 6 Recharts/Canvas charts, linked brushing, tabbed analytics panel

Progress: [█████████████████] 81%

## Performance Metrics

**Velocity:**
- Total plans completed: 16
- Average duration: 5.1 min
- Total execution time: ~1.35 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 1. Deep Validation | 3/4 | 18 min | 6 min |
| 2. Standardization Intelligence | 3/3 | 17 min | 5.7 min |
| 3. Batch Analytics | 6/6 | 31 min | 5.2 min |
| 4. Scoring Expansion | 3/3 | 15 min | 5.0 min |
| 5. Visualizations | 1/2 | 8 min | 8 min |
| 6. Export, API & Workflow | 0/3 | — | — |

**Recent Trend:**
- Last 7 plans: 03-04 (3.5 min), 03-06 (4 min), 03-05 (8 min), 04-01 (5 min), 04-02 (5 min), 04-03 (5 min), 05-01 (8 min)
- Trend: Stable

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Roadmap: Follow recommended build order from Features/ROADMAP.md (M1.1 → M1.2 → M1.3 → M2.1 → M3.1 → M4.3 → M3.2 → M3.3 → M3.5 → M4.1 → M4.2 → M5.1 → M5.2 → M3.4 → M2.2 → M6.1 → M6.2 → M6.3)
- Phase 3: INFRA-01 (batch pagination + Redis TTL) is plan 03-01 and must ship before any analytics endpoint
- Tautomer canonicalization: Document stereo stripping in all provenance and dedup responses (stereo_stripped flag)
- Validation cache: CHECKS_VERSION=v2 added to cache key in 01-01; key format is now validation:v2:{inchikey}:{checks_hash}
- Batch analytics: Run as post-aggregation Celery chord; frontend polls /batch/{job_id}/analytics separately from batch completion
- 01-01 (M1.1): TautomerDetection uses enumerator.Canonicalize() (GetCanonicalTautomer does not exist); stereoisomer cap=128; multi-fragment checks run on largest fragment only
- 01-02 (M1.2): MolVS REMOVE_FRAGMENTS broad SMARTS ([#7], [#8]) need heavy-atom-count guard to avoid false positives; largest carbon fragment always = drug overriding MW < 50 heuristic
- 01-03 (M1.3): Ring strain heuristic only (size 3/4); macrocycle threshold >12 not >=12; startup check uses logger.warning not hard assert; polymer severity INFO; zwitterion = net 0 with both + and - atoms
- [Phase 02-standardization-intelligence]: ProvenancePipeline wraps StandardizationPipeline externally — no internal modifications; provenance captured by re-calling same ChEMBL functions stage-by-stage
- [Phase 02-standardization-intelligence]: SMILES fragment diffing (set difference) for get_parent provenance — avoids atom-idx comparison pitfall when atom count changes after removal
- [Phase 02-standardization-intelligence]: Separate typed provenance fields per stage (charge_changes, bond_changes, etc.) instead of generic List[dict] to prevent TypeScript any[] on frontend
- [Phase 02-02]: Ring is aromatic iff ALL member atoms report GetIsAromatic()=True; atom count guard prevents ring diff index mismatch on fragment removal
- [Phase 02-02]: StereoProvenance.per_center changed from List[dict] to List[StereoCenterDetailSchema] — fully typed; DVAL cross-refs empty with TODO (optional per research)
- [Phase 02-02]: On-demand MoleculeViewer in ProvenanceStageCard — click "Show structure" to render; avoids RDKit.js cost for all stages by default
- [Phase 02-standardization-intelligence]: Fragment dict O=CO=formic acid: remove duplicate O=C(O)O key (formic acid carbonic), add O=CO for real formic acid; regression test via source inspection
- [Phase 02-standardization-intelligence]: DVAL cross-refs: dval_results kwarg on standardize_with_provenance(); StereoProvenance gains dval_cross_refs field for DVAL-01; tautomer ProvStageRecord.dval_cross_refs for DVAL-03; None default keeps all cross-refs empty (backward compatible)
- [Phase 03-batch-analytics]: BATCH_RESULT_TTL=86400 in config; import-guarded analytics modules in tasks allow incremental plan delivery; analytics_storage key format batch:analytics:{type}:{job_id}
- [Phase 03-02-deduplication]: Per-call TautomerEnumerator (not module-level) — not thread-safe; representative = min(indices) not indices[0] — parallel chunks don't preserve insertion order; stereo InChI stripping uses frozenset(['t','m','s']); total_unique = total_success minus duplicates removed per group
- [Phase 03-04-chemical-space]: PCA uses randomized numpy SVD (no sklearn), seeded np.random.default_rng(42); compute_chemical_space dispatcher added for analytics_tasks.py; pytest slow mark registered in pyproject.toml
- [Phase 03-06-statistics]: compute_all_statistics returns StatisticsResult Pydantic model (not plain dict) — analytics_tasks.py calls .model_dump(); abs < 1e-10 threshold for corrcoef noise (numpy returns -3.4e-16 not NaN for constant arrays); scaffold diversity inline (no scaffold_analysis dep); Lipinski fallback 50% when no data
- [Phase 03-03-scaffold-analysis]: Double-GetScaffoldForMol pattern for generic scaffold: MakeScaffoldGeneric converts exocyclic doubles; second GetScaffoldForMol removes them
- [Phase 03-03-scaffold-analysis]: Acyclic molecules grouped under empty string scaffold key; Shannon entropy edge-case (single scaffold = 0) handled before log2; frequency cap at 50 + Other bucket
- [Phase 03-05-mmp]: BRICS BRICSDecompose for MMP — phenol/aniline have no BRICS bonds, use phenylacetic acid/acetamide for tests; cap applied both inside _detect_mmp_pairs AND compute_mmp_analysis for monkeypatch robustness; compute_mmp returns _MMPResultWrapper with model_dump() for analytics_tasks.py compatibility
- [Phase 01]: Frontend Deep Validation Tab: Framer Motion AnimatePresence for expand/collapse; dynamic verdict from effective severities (not backend score); atom index badges trigger molecule viewer highlighting; fragment table as proper mini-table with colored classification badges
- [Phase 04-01]: ConsensusResult.rule_sets contains ALL properties (not just violations) per rule set; salt parent = fragment with max heavy atom count; LigandEfficiency uses BEI proxy when no external activity value
- [Phase 04-02]: TPSA uses _CalcTPSAContribs directly on mol (NO AddHs); LogP MUST use Chem.AddHs before _GetAtomContribs then fold H contributions back; aggregator confidence = triggered_count / TOTAL_INDICATORS (6)
- [Phase 04-03]: NP breakdown falls back gracefully when npscorer model not found (empty fragments list); ESOL coefficients duplicated inline in bioavailability_radar.py to avoid import cycle; BOILED-Egg: TPSA x-axis, WLOGP y-axis, point-in-ellipse test; radar normalization: in-range=1.0, linearly decreases outside range toward 0
- [Phase 04-03]: Recharts v3.7.0 RadarChart with PolarGrid for bioavailability radar; ScoringProfilesTab auto-fetches via useEffect when smiles prop changes; Framer Motion staggered card animations
- [Phase 05-01]: Canvas 2D for ChemicalSpaceScatter (not Recharts SVG) — required for >800 points; useReducer for brush selection (SET/TOGGLE/ADD_RANGE/CLEAR); SVG serializer for Recharts chart PNG export

### Pending Todos

None yet.

### Blockers/Concerns

- Phase 1: Startup assertion implemented (logger.warning) in 01-03; CI hard assert in test_deep_complexity_checks.py
- Phase 3: Redis memory ceiling under combined analytics load — profile during M3.1 setup
- Phase 6: Java availability in Docker for py2opsin (WORK-10) — validate with spike before committing full implementation

## Session Continuity

Last session: 2026-02-24
Stopped at: Completed 05-01-PLAN.md, executing 05-02
Resume file: .planning/phases/05-visualizations/05-01-SUMMARY.md
