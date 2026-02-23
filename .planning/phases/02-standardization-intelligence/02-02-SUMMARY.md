---
phase: 02-standardization-intelligence
plan: "02"
subsystem: api

tags: [rdkit, standardization, provenance, typescript, react, framer-motion, stereo, aromaticity]

# Dependency graph
requires:
  - phase: 02-01
    provides: ProvenancePipeline, StereoProvenance placeholder, all base provenance schemas

provides:
  - _capture_ring_aromaticity_changes() for ring-level aromaticity diff (STD-05)
  - StereoCenterDetail dataclass with atom_idx, type, before_config, after_config, reason
  - StereoComparison.per_center_detail list populated by extended compare()
  - RingChange Pydantic model and ring_changes field on ProvStageRecord
  - StereoCenterDetailSchema Pydantic model replacing List[dict] on StereoProvenance.per_center
  - ProvenanceTimeline React component (vertical timeline with expandable stage cards)
  - ProvenanceStageCard React component (typed change tables, on-demand MoleculeViewer)
  - All provenance TypeScript interfaces (ChargeChange, BondChange, RingChange, etc.)
  - Conditional provenance rendering in StandardizationResults

affects:
  - frontend standardization UI (ProvenanceTimeline now renders when provenance is present)
  - Any downstream code using StereoProvenance.per_center (now typed, not List[dict])

# Tech tracking
tech-stack:
  added: []
  patterns:
    - On-demand MoleculeViewer rendering via "Show structure" button (not eagerly loaded)
    - Atom count guard for ring aromaticity diff (same guard as charge tracking)
    - Ring aromaticity by ALL-atom test (a ring is aromatic iff all member atoms are aromatic)
    - per_center_detail built from set union of before/after chiral center maps
    - ConditionalProvenance pattern: result.provenance ? ProvenanceTimeline : StepsList

key-files:
  created:
    - backend/tests/test_standardization/test_provenance_m22.py
    - frontend/src/components/standardization/ProvenanceTimeline.tsx
    - frontend/src/components/standardization/ProvenanceStageCard.tsx
  modified:
    - backend/app/services/standardization/stereo_tracker.py
    - backend/app/services/standardization/provenance.py
    - backend/app/schemas/standardization.py
    - frontend/src/types/standardization.ts
    - frontend/src/components/standardization/StandardizationResults.tsx
    - frontend/src/components/standardization/index.ts

key-decisions:
  - "Ring is aromatic iff ALL member atoms report GetIsAromatic()=True — consistent with RDKit's ring-level aromaticity model"
  - "Atom count guard for ring diff (before.GetNumAtoms() == after.GetNumAtoms()) — prevents index mismatch when fragment removal changes atom count"
  - "StereoTracker.compare() extended with optional reason parameter (default: 'standardization') — supports future tautomer_canonicalization and get_parent reason tagging"
  - "StereoProvenance.per_center changed from List[dict] to List[StereoCenterDetailSchema] — eliminates TypeScript any[] risk"
  - "DVAL cross-refs left as empty list with TODO comment — populating requires frontend to pass prior DVAL results; documented as optional in research open questions"
  - "On-demand MoleculeViewer on 'Show structure' button click — avoids RDKit.js rendering cost for all stages by default"
  - "Auto-expand stage cards that have changes on first render; others start collapsed"

# Metrics
duration: 6min
completed: 2026-02-23
---

# Phase 02 Plan 02: Ring Aromaticity Tracking, Stereo Per-Center Detail, and Provenance Timeline UI Summary

**Ring-level aromaticity change tracking (STD-05), per-center stereochemistry normalization detail (STD-06), and full frontend provenance timeline with expandable stage cards, typed change tables, and on-demand molecule rendering**

## Performance

- **Duration:** 6 min
- **Started:** 2026-02-23T13:21:38Z
- **Completed:** 2026-02-23T13:27:12Z
- **Tasks:** 2
- **Files modified:** 9

## Accomplishments

- Ring aromaticity tracking: `_capture_ring_aromaticity_changes()` compares ring membership before/after standardization, identifies rings that changed from kekulized to aromatic (or vice versa), adds `RingChange` records with ring_atoms, ring_size, before_type, after_type to standardizer stage ProvStageRecord
- Stereo per-center detail: `StereoCenterDetail` dataclass and extended `StereoTracker.compare()` with `per_center_detail` list; builds from set union of before/after chiral center and double bond stereo maps
- `StereoCenterDetailSchema` Pydantic model replaces the `List[dict]` placeholder on `StereoProvenance.per_center` — now fully typed with atom_idx, type, before_config, after_config, reason
- Full provenance pipeline integration: stereo_summary in provenance now includes per_center list built from actual before/after stereo comparison
- All TypeScript provenance interfaces added to `standardization.ts`: ChargeChange, BondChange, RadicalChange, RingChange, FragmentRemoval, TautomerProvenance, StereoCenterDetail, StereoProvenance, ProvStageRecord, StandardizationProvenance (no any[] fields)
- `ProvenanceTimeline`: vertical timeline with tautomer summary banner, stereo warning banner, timeline dots, auto-expand stages with changes, applied/total stage count footer
- `ProvenanceStageCard`: expandable with status icon (CheckCircle2/MinusCircle/XCircle), change count badge, before/after SMILES (one-liner when collapsed), typed tables for charge/bond/ring/radical changes, fragment removal table with role badges, DVAL cross-ref list, on-demand MoleculeViewer with atom highlighting
- `StandardizationResults.tsx` updated: renders `ProvenanceTimeline` when `result.provenance` is present; falls back to `StepsList` when provenance is null (backward compatible)
- 12 new M2.2 tests pass (ring tracking, stereo detail, tracker extension); all 77 standardization tests pass

## Task Commits

Each task was committed atomically:

1. **Task 1: Ring aromaticity tracking (STD-05), stereo center detail (STD-06), and tests** - `4fcf73b` (feat)
2. **Task 2: Frontend provenance timeline UI and type extensions** - `9a481cd` (feat)

**Plan metadata:** (docs commit — this summary)

## Files Created/Modified

- `backend/app/services/standardization/stereo_tracker.py` - Added StereoCenterDetail dataclass, extended StereoComparison with per_center_detail, extended compare() with reason param and per-center population logic
- `backend/app/services/standardization/provenance.py` - Added _capture_ring_aromaticity_changes(), integrated ring_changes into standardizer stage, built stereo per_center in stereo_summary
- `backend/app/schemas/standardization.py` - Added RingChange and StereoCenterDetailSchema models; added ring_changes to ProvStageRecord; changed StereoProvenance.per_center from List[dict] to List[StereoCenterDetailSchema]
- `backend/tests/test_standardization/test_provenance_m22.py` - 12 tests: 4 ring tracking, 4 stereo normalization, 4 stereo tracker extension
- `frontend/src/types/standardization.ts` - All provenance TypeScript interfaces; include_provenance on StandardizationOptions; provenance field on StandardizationResult
- `frontend/src/components/standardization/ProvenanceTimeline.tsx` - Vertical timeline component with banners and stage cards
- `frontend/src/components/standardization/ProvenanceStageCard.tsx` - Expandable stage card with typed change tables and on-demand structure rendering
- `frontend/src/components/standardization/StandardizationResults.tsx` - Added ProvenanceTimeline import and conditional render
- `frontend/src/components/standardization/index.ts` - Exported ProvenanceTimeline and ProvenanceStageCard

## Decisions Made

- Ring aromaticity uses ALL-atom test (all ring atoms must be aromatic for ring to be "aromatic")
- Atom count guard prevents ring diff index mismatch on fragment removal stage
- StereoTracker.compare() reason parameter defaults to "standardization" for backward compat
- StereoProvenance.per_center changed from List[dict] to List[StereoCenterDetailSchema]
- DVAL cross-refs remain empty with TODO comment (populating requires frontend DVAL result passthrough — documented as optional in research)
- On-demand MoleculeViewer: renders only when user clicks "Show structure" button

## Deviations from Plan

None — plan executed exactly as written.

## Self-Check: PASSED

- stereo_tracker.py: FOUND
- provenance.py: FOUND
- standardization.py (schemas): FOUND
- test_provenance_m22.py: FOUND
- ProvenanceTimeline.tsx: FOUND
- ProvenanceStageCard.tsx: FOUND
- standardization.ts (types): FOUND
- 4fcf73b (Task 1 commit): FOUND
- 9a481cd (Task 2 commit): FOUND

---
*Phase: 02-standardization-intelligence*
*Completed: 2026-02-23*
