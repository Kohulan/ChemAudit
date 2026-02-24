---
phase: 05-visualizations
plan: 02
status: complete
completed: 2026-02-24
commits:
  - eaec55f  # Task 1: MoleculeComparisonPanel, MoleculePropertyRadar, BatchTimeline
  - 4292378  # Task 2: Comparison flow integration, timeline wiring
requirements_delivered:
  - VIZ-07
  - VIZ-08
  - VIZ-09
---

# Plan 05-02 Summary: Single Molecule Deep View

## What was built

### Task 1: Molecule comparison panel (VIZ-07) with property radar (VIZ-08)

**MoleculePropertyRadar** (`frontend/src/components/batch/MoleculePropertyRadar.tsx`):
- Per-molecule radar chart normalized 0-1 using Recharts RadarChart
- 6 properties: QED, SA Score, Fsp3, Validation Score, Lipinski Violations, Alert Count
- Inverted normalization for "bad" properties (lipinski violations, alert count) so higher = better on chart
- Dataset average overlay in green (rgba(16,185,129,0.4)) with molecule lines in primary/accent colors
- Supports 1-2 molecules overlaid against dataset average
- Custom tooltip showing raw values for each molecule and dataset average
- ResponsiveContainer with 280px height

**MoleculeComparisonPanel** (`frontend/src/components/batch/MoleculeComparisonPanel.tsx`):
- Slide-in drawer from right using Framer Motion AnimatePresence
- Fixed position, 480px wide, max-w-[90vw] for mobile
- Backdrop overlay (semi-transparent black) that closes panel on click
- 2D structure rendering via existing MoleculeViewer component (180x140px per molecule)
- Properties comparison table: Overall Score, QED, SA Score, Fsp3, Lipinski Violations, Alert Count
- Cell highlighting: green for better values, red for worse (respects higher/lower-is-better per property)
- Single column display when only 1 molecule selected, side-by-side for 2
- MoleculePropertyRadar embedded at bottom with all molecules overlaid
- Remove button on each molecule card to deselect individually

**BatchTimeline** (`frontend/src/components/batch/BatchTimeline.tsx`):
- 4-phase horizontal timeline: Upload -> Validation -> Analytics -> Complete
- Status-colored nodes (32px circles) with lucide-react icons
- Connecting segments between nodes colored by completion status
- Responsive: horizontal on lg screens, vertical on mobile
- Staggered entrance animation via Framer Motion (delay: i * 0.1)
- Summary stats below: validation count, processing time, scaffold/outlier counts
- Analytics phase derives status from Record<string, AnalysisStatus>

### Task 2: Comparison flow integration and timeline wiring

**BatchResultsTable modifications**:
- New `onCompare?: () => void` prop
- Floating "Compare (N)" button at bottom center when 1-2 molecules selected
- Button uses Framer Motion scale-in/out animation
- Disabled hint message when >2 molecules selected: "Select up to 2 molecules to compare"
- GitCompare icon from lucide-react

**BatchValidation page modifications**:
- `compareMode` state and `compareMolecules` derived state (first 2 selected)
- `handleCompare`, `handleCloseCompare`, `handleRemoveFromCompare` callbacks
- MoleculeComparisonPanel rendered conditionally when compareMode active
- BatchTimeline rendered above BatchAnalyticsPanel within analytics section
- useBatchAnalytics hook for timeline status and comparison radar dataset stats
- Compare mode resets on "Start New Batch"

## Decisions

- **Radar properties**: Used QED, SA Score, Fsp3, Validation Score, Lipinski Violations, Alert Count (not MW/LogP/TPSA as originally suggested in plan) since these are directly available from BatchResult scoring fields without needing separate property lookups
- **Comparison limit**: Strict 2-molecule maximum enforced at both UI level (button disabled for >2) and panel level (takes first 2 from selection)
- **Timeline placement**: Rendered above the analytics tabbed panel within the same card, providing processing context before chart exploration
- **Floating button**: Used fixed positioning with z-40 to ensure visibility across scroll positions

## Verification

- `npx tsc --noEmit` passes with zero errors
- All 5 files (3 new + 2 modified) compile correctly
- No new npm packages added (all using existing Recharts, Framer Motion, lucide-react)
