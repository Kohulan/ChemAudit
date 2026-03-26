# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-26)

**Core value:** Every chemical structure submitted gets a thorough, transparent, and reproducible quality assessment — from basic validity through ML-readiness — so scientists can trust their molecular data.
**Current focus:** v3.0 — Phase 7: Compound Profiling Engine

## Current Position

Phase: 7 of 14 (Compound Profiling Engine)
Plan: — (not yet planned)
Status: Ready to plan
Last activity: 2026-03-26 — v3.0 roadmap created (Phases 7-14, 39 requirements)

Progress: [░░░░░░░░░░] 0% (v3.0)

## Performance Metrics

**Velocity (v2.0 baseline):**
- Total plans completed: 27 (v2.0)
- Average duration: ~5 min/plan
- Total execution time: ~2.25 hours (v2.0)

**By Phase (v2.0):**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 1. Deep Validation | 4/4 | 18 min | 4.5 min |
| 2. Standardization Intelligence | 3/3 | 17 min | 5.7 min |
| 3. Batch Analytics | 6/6 | 31 min | 5.2 min |
| 4. Scoring Expansion | 3/3 | 15 min | 5.0 min |
| 5. Visualizations | 2/2 | 16 min | 8.0 min |
| 6. Export, API & Workflow | 9/9 | ~35 min | 3.9 min |

**Recent Trend:** Stable. v3.0 plans not yet started.

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- [v3.0 architecture]: All new services are strictly additive — new directories alongside existing services, no modifications to existing endpoints or schemas
- [v3.0 Phase 7]: SYBA must use subprocess isolation (GPL-3.0 cannot be imported into Apache-2.0 main process); SCScore weights must be converted to .npz at vendor time for NumPy 2.x compatibility
- [v3.0 Phase 7]: ETKDGv3 conformer generation requires 3-attempt fallback cascade (ETKDGv3 with torsion prefs → ETKDGv3 → ETKDG); return `3d_conformer_failed: true` on all failures rather than raising
- [v3.0 Phase 7]: medchem NIBRFilters must be pre-warmed at FastAPI and Celery startup as module-level singleton to avoid 5-10s first-call spike
- [v3.0 Phase 8]: Add `concern_group` field to AlertResult and group by functional group concern before counting — prevents alert overlap inflation from Kazius/BRENK/NIH matching same pattern
- [v3.0 Phase 10]: Pipeline result must carry both `original_inchikey` and `standardized_inchikey` with explicit `inchikey_changed` flag — design in data model before implementation
- [v3.0 Phase 13]: Butina clustering hard-capped at 1,000 molecules; `rdFMCS.FindMCS()` always uses `timeout=10`; RegistrationHash must store RDKit version; `enable_tautomer_hash_v2=True` explicit in every call
- [v3.0 Phase 14]: MCP server must use `include_tags` allowlist in `FastApiMCP()` constructor; startup assertion that `POST /api/keys` is absent from MCP tool list is mandatory

### Pending Todos

None.

### Blockers/Concerns

- [Phase 7 — pre-start]: Validate SCScore .npz weight conversion and SYBA subprocess isolation latency before beginning Phase 7 SA Comparison implementation (research flag from SUMMARY.md)
- [Phase 13 — pre-start]: ClassyFire-style taxonomy needs a curated 50-100 chemotype rule set; evaluate DrugTax package as candidate source during Phase 13 pre-planning (research flag from SUMMARY.md)
- [Phase 8 — pre-start]: Verify medchem ComplexityFilter exact API import path and call signature against medchem 2.0.5 source before implementation

## Session Continuity

Last session: 2026-03-26
Stopped at: v3.0 roadmap created — Phases 7-14 defined, 39/39 requirements mapped, ROADMAP.md and STATE.md written
Resume file: .planning/ROADMAP.md
