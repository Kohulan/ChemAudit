---
phase: 01-deep-validation
plan: "04"
subsystem: frontend-validation
tags:
  - frontend
  - react
  - typescript
  - deep-validation
  - ui-components
dependency_graph:
  requires:
    - frontend/src/types/validation.ts
    - frontend/src/hooks/useLocalStorage.ts
    - frontend/src/pages/SingleValidation.tsx
    - frontend/src/components/validation/ValidationResults.tsx
    - frontend/src/components/validation/IssueCard.tsx
  enables: []
---

## Summary

Built the frontend "Deep Validation" tab for the single validation results page with full component suite.

## Key Files

### Created
- `frontend/src/types/validation.ts` — 17 typed interfaces for deep validation check details + domain grouping constants
- `frontend/src/hooks/useDeepValidationConfig.ts` — localStorage-backed severity configuration hook
- `frontend/src/components/validation/DeepValidationTab.tsx` — Main tab with Category/Severity segmented control and dynamic verdict
- `frontend/src/components/validation/DeepCheckCard.tsx` — Per-check card with structured detail rendering and atom cross-linking
- `frontend/src/components/validation/SeverityConfigPanel.tsx` — Slide-in modal with per-check ERROR/WARNING/INFO toggles
- `frontend/src/components/validation/StereoisomerList.tsx` — Collapsible stereoisomer SMILES list with copy buttons
- `frontend/src/components/validation/FragmentClassificationTable.tsx` — Mini-table with colored classification badges
- `frontend/src/tests/components/DeepValidationTab.test.tsx` — 41 tests across 5 component test suites

### Modified
- `frontend/src/pages/SingleValidation.tsx` — Added "Deep Validation" tab with Microscope icon

## Decisions
- Used Framer Motion AnimatePresence for all expand/collapse and view transitions
- Dynamic verdict computed from effective severities (user overrides) not backend score
- Atom indices rendered as clickable badges that trigger molecule viewer highlighting
- Fragment classification badges: drug=green, salt=blue, solvent=yellow, unknown=gray

## Deviations
- Fixed TypeScript error in severityGroups array (needed `as const` for variant literals after `.filter()`)
- Fixed 3 test failures: selector ambiguity for "severity" buttons, Framer Motion async collapse assertion

## Self-Check: PASSED
- [x] TypeScript compiles without errors
- [x] All 41 frontend tests pass
- [x] All 3 tasks committed atomically
