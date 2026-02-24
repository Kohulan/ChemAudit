# Phase 5: Visualizations - Context

**Gathered:** 2026-02-24
**Status:** Ready for planning

<domain>
## Phase Boundary

9 interactive chart components (Recharts-based, zero new npm packages) making batch analytics and single-molecule scoring data visually explorable. Consumes Phase 3 analytics endpoints (scaffold, chemical space, statistics) and Phase 4 scoring endpoints (drug-likeness, property breakdowns, radars). All are pure presentational components with data fetched by parent containers. No new backend logic — frontend-only phase.

</domain>

<decisions>
## Implementation Decisions

### Batch Dashboard Layout
- **Tabbed sections** for organizing the 6 batch charts — one context at a time, clean focused view
- **Summary badges on tabs** — each tab shows a key metric count (e.g., "3 outliers", "12 scaffolds") so users see what's in each tab without clicking
- **Skeleton placeholders** while analytics data loads — chart-shaped skeleton animations, professional feel
- **PNG download per chart** — each chart gets a small download icon for saving individual chart images

### Chart Interactions
- **Click opens molecule detail panel** — clicking a data point in scatter/histogram opens a side panel or modal showing the molecule's structure, SMILES, and key properties. Chart stays visible (non-disruptive)
- **Brush/lasso selection supported** — users can drag to select a region in scatter plots. Selected molecules highlighted across all charts and exportable
- **Linked highlighting (cross-chart brushing)** — selecting molecules in one chart highlights the same molecules in all other visible charts. Powerful for multi-dimensional exploration
- **PCA default, t-SNE on demand** — chemical space scatter shows PCA by default (always available). Separate "Compute t-SNE" button for batches ≤2000 molecules since it's expensive

### Molecule Comparison Flow
- **Multi-select from batch table + charts** — checkboxes in batch results table AND click-to-add from chart data points. Floating "Compare (N)" button appears when molecules are selected
- **Up to 2 molecules** side-by-side — strict pair comparison for clean, readable layout
- **Comparison shows: structure + properties table + radar** — 2D structure rendering, key properties in a comparison table (MW, LogP, TPSA, etc.), and overlaid radar chart

### Chemical Space Scatter Coloring
- **Dropdown to select color-by property** — users can color scatter plot points by any computed property (MW, LogP, TPSA, score, etc.). Default is overall score. Adds analytical depth

### Scaffold Treemap
- **Categorical colors** — each scaffold family gets a distinct color from a categorical palette. Top scaffolds easily distinguished. "Other" group is gray

### Color & Visual Encoding
- **Match existing app theme** — derive all chart colors from the existing "Laboratory Elegance" palette (crimson primary #c41e3a, warm amber accent #d97706, warm stone neutrals). CSS variable-driven for dark/light mode
- Severity/pass-fail coloring: derive from current theme — existing score tiers (amber gold, orange, red) and status colors (amber success, orange warning, red error) already established

### Claude's Discretion
- Tab grouping strategy (how to distribute 6 charts across 2-3 tabs — use domain analysis: distributions vs structure analysis vs chemical space)
- Chart sizing approach (fixed responsive vs resizable — lean toward fixed responsive for simplicity)
- Placement of analytics section on batch page (below results table vs separate tab)
- Error state handling for failed/partial analytics data
- Radar chart dataset average reference behavior (always visible vs toggle)

</decisions>

<specifics>
## Specific Ideas

### Production-Ready Requirements (User Emphasis)
- Complete, proper implementation — no quick fixes or bandaids
- Must not break existing systems (batch results table, single validation, scoring profiles)
- All components must be fully production-ready with proper error handling

### Research-Informed Technical Patterns
- **Recharts SVG safe for ≤800 points** — scatter plots, histograms, radars, treemaps with <100 cells all fine as SVG
- **Canvas 2D for >800 points** — chemical space scatter (up to 2000+ molecules) and similarity heatmap (100x100) must use Canvas 2D rendering, not Recharts SVG
- **d3-quadtree for canvas hit detection** — O(log n) nearest-neighbor lookup for tooltip hit detection on mousemove over canvas scatter plots
- **Linked brushing via useReducer** — shared `Set<string>` of InChIKeys at page level, each chart receives `selectedKeys` + `onSelectionChange` as props. Wrap charts in `React.memo` to prevent cascade re-renders
- **PNG export approach** — `recharts-to-png` for Recharts-based charts, `canvas.toBlob()` for Canvas charts. Dismiss tooltips before export to avoid foreignObject rendering issues
- **Scaffold treemap cap at 50** — remainder bucketed as "Other". Use Recharts `<Treemap content={...}>` with RDKit.js SVG in `<foreignObject>` for scaffold structure thumbnails
- **Histograms: pre-binned data from backend** — compute bins server-side with NumPy, send `[{bin_start, bin_end, count}]` to frontend. Use `<BarChart>` with `barCategoryGap={0}` and `barGap={0}`
- **100x100 heatmap via Canvas fillRect** — 10,000 cells in <2ms. Never create a DOM element per cell
- **Molecule SVG in tooltips** — call `rdkit.get_mol(smiles).get_svg(150, 120)` in custom tooltip renderer. Debounce 100ms. Always call `.delete()` to avoid WASM memory leaks
- **Accessibility** — wrap each chart in `<figure>` with `aria-label` and `<figcaption>` text summary (Recharts has no built-in ARIA support)

### Existing Codebase Patterns to Follow
- Tab bar pattern from `SingleValidation.tsx`: active = `bg-[var(--color-surface-elevated)]` white pill + primary text + `shadow-sm`, Framer Motion `AnimatePresence mode="wait"` for content transitions
- Card pattern from `ScoringProfilesTab.tsx`: `rounded-2xl p-5 bg-gradient-to-br from-[var(--color-surface-elevated)] to-[var(--color-surface)] border border-[var(--color-border)]`
- Section header pattern: icon in gradient background square + title + subtitle
- Batch table pattern in `BatchResultsTable.tsx`: filter bar in `bg-[var(--color-surface-sunken)]`, click-to-expand bento grid
- Score coloring: `score-excellent` (#fbbf24 amber), `score-good` (#d97706), `score-fair` (#ea580c), `score-poor` (#dc2626)
- All existing Recharts charts use `<ResponsiveContainer width="100%" height={N}>` with `var(--color-border)` grid lines and `var(--color-text-secondary)` axis labels

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 05-visualizations*
*Context gathered: 2026-02-24*
