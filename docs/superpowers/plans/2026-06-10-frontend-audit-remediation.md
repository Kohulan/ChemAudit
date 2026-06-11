# Frontend Audit Remediation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the P1/P2 findings from the 2026-06-10 `/impeccable audit` (score 11/20): WCAG AA contrast and labeling gaps, RDKit cold-start cost, bundle chunking, unresponsive grids, and ~50 hard-coded colors bypassing the CSS-variable token system.

**Architecture:** All changes are in `frontend/`. Work proceeds in four shippable phases: (A) accessibility, (B) performance, (C) responsive, (D) theming. Each task is independently committable. A new shared module `src/lib/chartColors.ts` becomes the single source of truth for score-quality colors in SVG/Recharts fills (where Tailwind classes can't reach); everything else maps onto the existing CSS-variable tokens in `src/index.css` and `tailwind.config.js`.

**Tech Stack:** React 18 + TypeScript, Vite 7, Tailwind 3 (CSS-variable tokens, `darkMode: 'class'`), Framer Motion, Recharts, RDKit.js (CDN WASM), Vitest 4 + Testing Library (jsdom, setup at `src/tests/setup.ts`).

**Key domain facts the engineer must know:**
- The brand deliberately uses a warm hierarchy: **success/excellent = gold/amber, NOT green** (see `PRODUCT.md` and `tailwind.config.js` `status-success`/`score-*` tokens). Replacing green with amber is a fix, not a regression.
- Theme switching toggles a `dark` class on `<html>`; all tokens are CSS variables defined twice in `src/index.css` (`:root` at ~line 18, `.dark` at ~line 97). Charts get `isDark` from `useThemeContext()` (`src/contexts/ThemeContext.tsx`).
- RDKit loads via a CDN `<script>` in `index.html`; `getRDKit()` in `src/hooks/useRDKit.ts` is a lazy singleton — the heavy WASM fetch happens only on the first `initRDKitModule()` call.
- All commands below run from `frontend/` unless stated otherwise. Do NOT add `Co-Authored-By` lines to commits (project rule).

**Setup:**

- [ ] Create a working branch off `development`:

```bash
cd /Volumes/Data_Drive/My_Projects/2026/ChemAudit/ChemAudit
git checkout development && git pull
git checkout -b fix/frontend-audit-remediation
cd frontend && npm install
```

- [ ] Baseline: record current bundle sizes and test status for later comparison:

```bash
npm run test:run 2>&1 | tail -5
npm run build && ls -la dist/assets/ | sort -k5 -n | tail -10
```

Save the size of the largest JS chunk; Task 6 compares against it.

---

## Phase A: Accessibility (WCAG 2.1 AA)

### Task 1: Fix `--color-text-muted` contrast in both themes

The muted text token fails AA in light mode (~3.2:1 on sunken surfaces) and badly in dark mode (~2.1:1).

**Files:**
- Modify: `frontend/src/index.css:26` (light) and `frontend/src/index.css:105` (dark)

- [ ] **Step 1: Verify the replacement values pass 4.5:1 against every surface token**

```bash
node -e '
const L=h=>{const c=[1,3,5].map(i=>parseInt(h.slice(i,i+2),16)/255).map(v=>v<=0.04045?v/12.92:((v+0.055)/1.055)**2.4);return 0.2126*c[0]+0.7152*c[1]+0.0722*c[2]};
const cr=(a,b)=>{const[x,y]=[L(a),L(b)].sort((p,q)=>q-p);return((x+0.05)/(y+0.05)).toFixed(2)};
console.log("light muted #6f655c on sunken #f3f1ed:", cr("#6f655c","#f3f1ed"));
console.log("light muted #6f655c on surface #faf9f7:", cr("#6f655c","#faf9f7"));
console.log("light muted #6f655c on white:", cr("#6f655c","#ffffff"));
console.log("dark muted #8a857e on elevated #1c1917:", cr("#8a857e","#1c1917"));
console.log("dark muted #8a857e on surface #0c0a09:", cr("#8a857e","#0c0a09"));
'
```

Expected: every printed ratio ≥ 4.50 (approximately 5.0, 5.2, 5.7, 4.8, 5.4).

- [ ] **Step 2: Edit the light-mode token**

In `frontend/src/index.css` (inside `:root`, ~line 26):

```css
/* before */
    --color-text-muted: #9c958d;
/* after */
    --color-text-muted: #6f655c;
```

- [ ] **Step 3: Edit the dark-mode token**

In `frontend/src/index.css` (inside `.dark`, ~line 105):

```css
/* before */
    --color-text-muted: #57534e;
/* after */
    --color-text-muted: #8a857e;
```

- [ ] **Step 4: Check for hard-coded duplicates of the old values**

```bash
grep -rn "9c958d\|#57534e" src/ --include="*.tsx" --include="*.ts" --include="*.css"
```

Expected: only matches that are NOT the muted-text token (e.g. `#57534e` appears legitimately as stone-700 in `ScoreChart.tsx` cool-variant fills; leave those, Task 13/14 handles charts). If any other file uses `#9c958d` *as muted text*, replace it with `var(--color-text-muted)`.

- [ ] **Step 5: Visual smoke test**

```bash
npm run dev
```

Open http://localhost:3002, toggle dark mode via the header flask toggle. Hint/subtitle text (e.g. "Choose a format to download batch results" in the export dialog, uppercase tile labels) must be readable in both themes, noticeably more legible in dark mode.

- [ ] **Step 6: Commit**

```bash
git add src/index.css
git commit -m "fix(a11y): raise text-muted contrast to WCAG AA in both themes"
```

### Task 2: ExportDialog dialog semantics and labeled close button

The close button is the app's only unlabeled icon-only button; the dialog also lacks `role="dialog"`.

**Files:**
- Modify: `frontend/src/components/batch/ExportDialog.tsx:138-164`
- Test: `frontend/src/components/batch/ExportDialog.test.tsx` (create)

- [ ] **Step 1: Write the failing test**

Create `frontend/src/components/batch/ExportDialog.test.tsx`:

```tsx
import { render, screen } from '@testing-library/react';
import { describe, it, expect, vi } from 'vitest';
import { ExportDialog } from './ExportDialog';

vi.mock('../../services/api', () => ({
  api: { get: vi.fn(), post: vi.fn() },
}));

describe('ExportDialog accessibility', () => {
  it('exposes dialog role named by its heading', () => {
    render(<ExportDialog jobId="abcd1234" isOpen onClose={() => {}} />);
    expect(
      screen.getByRole('dialog', { name: /export results/i })
    ).toBeInTheDocument();
  });

  it('has an accessible name on the close button', () => {
    render(<ExportDialog jobId="abcd1234" isOpen onClose={() => {}} />);
    expect(
      screen.getByRole('button', { name: /close export dialog/i })
    ).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
npx vitest run src/components/batch/ExportDialog.test.tsx
```

Expected: FAIL, "Unable to find an accessible element with the role 'dialog'".

(If the mock of `../../services/api` fails because the module exports differently, check `src/services/api.ts` exports and mock the actual export names; `ExportDialog.tsx:6` imports `{ api }`.)

- [ ] **Step 3: Add the ARIA attributes**

In `frontend/src/components/batch/ExportDialog.tsx`, three edits:

Backdrop (line 138):

```tsx
// before
    <div className="fixed inset-0 bg-[var(--color-text-primary)]/50 flex items-center justify-center z-50" onClick={onClose}>
// after
    <div className="fixed inset-0 bg-[var(--color-text-primary)]/50 flex items-center justify-center z-50" onClick={onClose} role="presentation">
```

Dialog container (the `motion.div`, lines 139-146) — add three props alongside the existing `onClick`:

```tsx
      <motion.div
        initial={{ opacity: 0, scale: 0.95, y: 10 }}
        animate={{ opacity: 1, scale: 1, y: 0 }}
        exit={{ opacity: 0, scale: 0.95, y: 10 }}
        transition={{ duration: 0.2 }}
        className="bg-[var(--color-surface-elevated)] rounded-2xl shadow-xl max-w-2xl w-full mx-4 max-h-[90vh] overflow-hidden border border-[var(--color-border)] flex flex-col"
        onClick={(e) => e.stopPropagation()}
        role="dialog"
        aria-modal="true"
        aria-labelledby="export-dialog-title"
      >
```

Heading (line 150) and close button (lines 159-164):

```tsx
// heading: add id
            <h2 id="export-dialog-title" className="text-xl font-semibold text-[var(--color-text-primary)] font-display">Export Results</h2>
```

```tsx
// close button: add aria-label
          <button
            onClick={onClose}
            aria-label="Close export dialog"
            className="p-2 rounded-lg hover:bg-[var(--color-surface-sunken)] transition-colors"
          >
            <X className="w-5 h-5 text-[var(--color-text-muted)]" />
          </button>
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
npx vitest run src/components/batch/ExportDialog.test.tsx
```

Expected: 2 passed.

- [ ] **Step 5: Commit**

```bash
git add src/components/batch/ExportDialog.tsx src/components/batch/ExportDialog.test.tsx
git commit -m "fix(a11y): add dialog semantics and close-button label to ExportDialog"
```

### Task 3: Differentiate success and warning badges (color + icon)

`.badge-success` and `.badge-warning` are both amber and visually identical (WCAG 1.4.1 use-of-color). Fix twice over: shift warning to the orange band of the warm palette, AND auto-add a distinguishing icon so color is never the sole signal.

**Files:**
- Modify: `frontend/src/index.css:657-673`
- Modify: `frontend/src/components/ui/Badge.tsx`
- Test: `frontend/src/components/ui/Badge.test.tsx` (create)

- [ ] **Step 1: Write the failing test**

Create `frontend/src/components/ui/Badge.test.tsx`:

```tsx
import { render } from '@testing-library/react';
import { describe, it, expect } from 'vitest';
import { Badge, CountBadge } from './Badge';

describe('Badge status differentiation', () => {
  it('adds a distinguishing icon to success badges', () => {
    const { container } = render(<Badge variant="success">Valid</Badge>);
    expect(container.querySelector('svg')).not.toBeNull();
  });

  it('adds a distinguishing icon to warning badges', () => {
    const { container } = render(<Badge variant="warning">Check</Badge>);
    expect(container.querySelector('svg')).not.toBeNull();
  });

  it('respects an explicit icon override', () => {
    const { container } = render(
      <Badge variant="success" icon={<span data-testid="custom" />}>Valid</Badge>
    );
    expect(container.querySelector('svg')).toBeNull();
  });

  it('keeps numeric CountBadge icon-free', () => {
    const { container } = render(<CountBadge count={5} variant="success" />);
    expect(container.querySelector('svg')).toBeNull();
  });
});
```

- [ ] **Step 2: Run to verify failure**

```bash
npx vitest run src/components/ui/Badge.test.tsx
```

Expected: FAIL on the two "adds a distinguishing icon" cases (no svg rendered).

- [ ] **Step 3: Implement auto-icons in Badge**

In `frontend/src/components/ui/Badge.tsx`, add the import and replace the `Badge` function body:

```tsx
import { cva, type VariantProps } from 'class-variance-authority';
import { CheckCircle2, AlertTriangle } from 'lucide-react';
import { cn } from '../../lib/utils';
```

```tsx
export function Badge({
  children,
  variant,
  size,
  icon,
  dot = false,
  className,
  ...props
}: BadgeProps) {
  // Color alone must not distinguish success from warning (WCAG 1.4.1):
  // both live in the warm amber/orange band, so default a redundant icon.
  // Pass icon={null} to opt out (numeric CountBadge does).
  const resolvedIcon =
    icon !== undefined
      ? icon
      : variant === 'success'
        ? <CheckCircle2 className="w-3 h-3" aria-hidden="true" />
        : variant === 'warning'
          ? <AlertTriangle className="w-3 h-3" aria-hidden="true" />
          : null;

  return (
    <span
      className={cn(badgeVariants({ variant, size }), dot && 'badge-dot', className)}
      {...props}
    >
      {resolvedIcon && <span className="flex-shrink-0">{resolvedIcon}</span>}
      {children}
    </span>
  );
}
```

In `CountBadge` (same file, ~line 69), pass the opt-out:

```tsx
  return (
    <Badge
      variant={variant}
      size="sm"
      icon={null}
      className={cn('min-w-[1.25rem] justify-center tabular-nums', className)}
    >
      {displayCount}
    </Badge>
  );
```

Note: `icon`'s type is `React.ReactNode`, which already admits `null`; the `icon !== undefined` check distinguishes "not passed" from "explicitly none".

- [ ] **Step 4: Shift `.badge-warning` to the orange band**

In `frontend/src/index.css` lines 666-673, replace:

```css
  .badge-warning {
    background: rgba(245, 158, 11, 0.12);
    color: #b45309;
  }
  .dark .badge-warning {
    background: rgba(245, 158, 11, 0.18);
    color: #fbbf24;
  }
```

with:

```css
  .badge-warning {
    background: rgba(234, 88, 12, 0.12);
    color: #c2410c;
  }
  .dark .badge-warning {
    background: rgba(249, 115, 22, 0.18);
    color: #fb923c;
  }
```

(`#c2410c` is orange-700, ~5.2:1 on white; success keeps the gold band per the brand's warm hierarchy.)

- [ ] **Step 5: Run the test to verify it passes, then check existing suites**

```bash
npx vitest run src/components/ui/Badge.test.tsx && npm run test:run 2>&1 | tail -5
```

Expected: Badge tests pass; no new failures elsewhere. If a snapshot or text-content test breaks because badges now contain an svg, update that test, the new markup is intended.

- [ ] **Step 6: Visual check**

In the dev server, validate a molecule (e.g. paste `CCO` on the home page): the "Validated" badge (SingleValidation.tsx:1204) shows a check icon; any warning badges show a triangle and read as orange, not gold.

- [ ] **Step 7: Commit**

```bash
git add src/components/ui/Badge.tsx src/components/ui/Badge.test.tsx src/index.css
git commit -m "fix(a11y): distinguish success/warning badges by icon and hue, not color alone"
```

---

## Phase B: Performance

### Task 4: Route-gate RDKit initialization (kill the 2-3s cold start on non-structure pages)

`AppWithSplash` calls `useRDKit()` unconditionally, so the WASM fetch + 2.6s splash blocks even `/about` and `/privacy`. Gate it by route. `getRDKit()` stays a lazy singleton, so any structure component that mounts later still triggers init on demand.

**Files:**
- Create: `frontend/src/lib/rdkitRoutes.ts`
- Test: `frontend/src/lib/rdkitRoutes.test.ts` (create)
- Modify: `frontend/src/hooks/useRDKit.ts:49-77`
- Modify: `frontend/src/App.tsx:297-348`

- [ ] **Step 1: Write the failing test**

Create `frontend/src/lib/rdkitRoutes.test.ts`:

```ts
import { describe, it, expect } from 'vitest';
import { routeNeedsRDKit } from './rdkitRoutes';

describe('routeNeedsRDKit', () => {
  it.each(['/about', '/privacy', '/diagnostics', '/about/'])(
    'returns false for RDKit-free route %s',
    (path) => {
      expect(routeNeedsRDKit(path)).toBe(false);
    }
  );

  // History and Bookmarks render molecule thumbnails, so they DO need RDKit.
  it.each(['/', '/batch', '/history', '/bookmarks', '/qsar-ready', '/structure-filter', '/dataset-audit', '/report/abc123'])(
    'returns true for structure route %s',
    (path) => {
      expect(routeNeedsRDKit(path)).toBe(true);
    }
  );
});
```

- [ ] **Step 2: Run to verify failure**

```bash
npx vitest run src/lib/rdkitRoutes.test.ts
```

Expected: FAIL, "Cannot find module './rdkitRoutes'" (or equivalent resolve error).

- [ ] **Step 3: Implement the helper**

Create `frontend/src/lib/rdkitRoutes.ts`:

```ts
/**
 * Routes that never render molecular structures and therefore must not
 * pay the RDKit WASM init cost (or the splash wait) on first load.
 *
 * History and Bookmarks are deliberately NOT listed: both render
 * MoleculeViewer thumbnails.
 */
const RDKIT_FREE_ROUTES = new Set(['/about', '/privacy', '/diagnostics']);

export function routeNeedsRDKit(pathname: string): boolean {
  const normalized = pathname.replace(/\/+$/, '') || '/';
  return !RDKIT_FREE_ROUTES.has(normalized);
}
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
npx vitest run src/lib/rdkitRoutes.test.ts
```

Expected: PASS (12 cases).

- [ ] **Step 5: Add an `enabled` flag to `useRDKit`**

In `frontend/src/hooks/useRDKit.ts`, replace the hook (lines 49-77) with:

```ts
export function useRDKit(enabled: boolean = true) {
  const [rdkit, setRDKit] = useState<RDKitModule | null>(null);
  const [loading, setLoading] = useState(enabled);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!enabled) {
      setLoading(false);
      return;
    }

    let mounted = true;
    setLoading(true);

    getRDKit()
      .then((module) => {
        if (mounted) {
          setRDKit(module);
          setLoading(false);
        }
      })
      .catch((err) => {
        if (mounted) {
          setError(err.message);
          setLoading(false);
        }
      });

    return () => {
      mounted = false;
    };
  }, [enabled]);

  return { rdkit, loading, error };
}
```

All existing zero-argument callers keep today's behavior (`enabled` defaults to `true`).

- [ ] **Step 6: Gate the splash in App.tsx**

In `frontend/src/App.tsx`:

Add to the router import (line 2):

```tsx
import { BrowserRouter as Router, Routes, Route, Navigate, useLocation } from 'react-router-dom';
```

Add below the other local imports (after line 12):

```tsx
import { routeNeedsRDKit } from './lib/rdkitRoutes';
```

Replace the top of `AppWithSplash` (lines 297-300):

```tsx
function AppWithSplash() {
  const { pathname } = useLocation();
  const needsRDKit = routeNeedsRDKit(pathname);
  const { loading, error } = useRDKit(needsRDKit);
  // Splash only gates structure routes; /about, /privacy, /diagnostics render immediately.
  const [showSplash, setShowSplash] = useState(needsRDKit);
  const [minTimeElapsed, setMinTimeElapsed] = useState(false);
```

And guard the error branch (line 334):

```tsx
  // Show error state (only meaningful on routes that require RDKit)
  if (needsRDKit && error) {
```

- [ ] **Step 7: Verify behavior in the browser**

```bash
npm run dev
```

- Open http://localhost:3002/about in a fresh tab with DevTools Network open: the page renders with **no splash screen** and **no `.wasm` request** (the small `RDKit_minimal.js` script tag still loads; that's expected).
- Open http://localhost:3002/ : splash appears as before, WASM loads, validation works (paste `CCO`).
- Navigate /about → / via the header: home renders and molecule drawing appears once RDKit finishes loading in the background.

- [ ] **Step 8: Run the full test suite**

```bash
npm run test:run 2>&1 | tail -5
```

Expected: no new failures.

- [ ] **Step 9: Commit**

```bash
git add src/lib/rdkitRoutes.ts src/lib/rdkitRoutes.test.ts src/hooks/useRDKit.ts src/App.tsx
git commit -m "perf: route-gate RDKit WASM init; skip splash on structure-free pages"
```

### Task 5: Pin, integrity-check, and defer the RDKit loader script

The current tag is unpinned (unpkg serves whatever version resolves) and has no Subresource Integrity, so a CDN compromise would execute arbitrary JS in every user's session. Pin to the version in `package.json` (`2024.3.5-1.0.0`), add the SRI hash, and `defer` parsing.

**Files:**
- Modify: `frontend/index.html:24`

- [ ] **Step 1: Verify the SRI hash against the pinned artifact**

The hash below was computed from the live artifact on 2026-06-10; re-verify it yourself rather than trusting the plan:

```bash
curl -sL "https://unpkg.com/@rdkit/rdkit@2024.3.5-1.0.0/dist/RDKit_minimal.js" | openssl dgst -sha384 -binary | openssl base64 -A
```

Expected output: `PxXM29TKNj8nIQ5WQGTV9mEZVT2GEQsMhTdW1ziTUROGZzrdmfAAVQvTX82TAV1q`

If it differs, STOP and investigate before shipping (either the package was republished or the CDN is serving something else).

- [ ] **Step 2: Replace the script tag**

```html
<!-- before -->
    <script src="https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js"></script>
<!-- after -->
    <script
      defer
      src="https://unpkg.com/@rdkit/rdkit@2024.3.5-1.0.0/dist/RDKit_minimal.js"
      integrity="sha384-PxXM29TKNj8nIQ5WQGTV9mEZVT2GEQsMhTdW1ziTUROGZzrdmfAAVQvTX82TAV1q"
      crossorigin="anonymous"
    ></script>
```

Notes:
- `defer` is safe: `getRDKit()` only reads `window.initRDKitModule` inside a React effect, which always runs after deferred scripts execute.
- Pinning the JS also pins the WASM: the loader fetches `RDKit_minimal.wasm` relative to its own (now-versioned) URL. SRI covers the JS loader only; the WASM fetch is not integrity-checked by the browser, which is an accepted residual risk here.
- When bumping `@rdkit/rdkit` in `package.json`, this URL and hash must be updated in the same commit; recompute with the Step 1 command.

- [ ] **Step 3: Verify**

Hard-reload http://localhost:3002/ — splash, then working validation (paste `CCO`, structure renders). No "RDKit.js not loaded" error and no SRI failure in the console (an integrity mismatch logs a distinctive "Failed to find a valid digest" error and blocks the script).

- [ ] **Step 4: Commit**

```bash
git add index.html
git commit -m "perf(security): pin RDKit loader version, add SRI, defer parsing"
```

### Task 6: Vendor chunk splitting in Vite

Recharts, Framer Motion, and GSAP currently land in the main chunk.

**Files:**
- Modify: `frontend/vite.config.ts:23-47`

- [ ] **Step 1: Add `build.rollupOptions.output.manualChunks`**

In the returned config object in `frontend/vite.config.ts`, after the `server` block (line 41) and before `define`, insert:

```ts
    build: {
      rollupOptions: {
        output: {
          manualChunks: {
            recharts: ['recharts'],
            motion: ['framer-motion'],
            gsap: ['gsap'],
            vendor: ['react', 'react-dom', 'react-router-dom'],
          },
        },
      },
    },
```

(Do not add `@rdkit/rdkit`: it loads from CDN, not from the bundle.)

- [ ] **Step 2: Build and compare against the baseline from Setup**

```bash
npm run build && ls -la dist/assets/ | sort -k5 -n | tail -12
```

Expected: new `recharts-*.js`, `motion-*.js`, `gsap-*.js`, `vendor-*.js` files exist, and the main `index-*.js` chunk is substantially smaller than the baseline figure recorded in Setup. Record the new sizes.

- [ ] **Step 3: Smoke-test the production build**

```bash
npm run preview
```

Open the printed URL; home page, batch page, and dark-mode toggle work without console errors.

- [ ] **Step 4: Commit**

```bash
git add vite.config.ts
git commit -m "perf: split recharts, framer-motion, gsap, and react vendor chunks"
```

### Task 7: Lazy-load heavy tab results in SingleValidation

The eagerly-loaded 1,800-line home page imports every tab's result components up front. Defer the two heaviest (`ScoringResults` pulls Recharts; `StandardizationResults` pulls the provenance timeline).

**Files:**
- Modify: `frontend/src/pages/SingleValidation.tsx` (imports at lines 1 and 44-45, plus each render site)

- [ ] **Step 1: Find the render sites**

```bash
grep -n "<ScoringResults\|<StandardizationResults" src/pages/SingleValidation.tsx
```

Record the line numbers (expect 1-3 sites total, inside the tab-content region after line ~1454).

- [ ] **Step 2: Convert the imports**

In `frontend/src/pages/SingleValidation.tsx` line 1, extend the React import:

```tsx
import { lazy, Suspense, useState, useEffect, useLayoutEffect, useRef, useCallback, useMemo } from 'react';
```

Replace lines 44-45:

```tsx
// before
import { ScoringResults } from '../components/scoring/ScoringResults';
import { StandardizationResults } from '../components/standardization/StandardizationResults';
// after
const ScoringResults = lazy(() =>
  import('../components/scoring/ScoringResults').then((m) => ({ default: m.ScoringResults }))
);
const StandardizationResults = lazy(() =>
  import('../components/standardization/StandardizationResults').then((m) => ({ default: m.StandardizationResults }))
);
```

- [ ] **Step 3: Wrap each render site in Suspense**

At every line found in Step 1, wrap the element, keeping its props exactly as they are:

```tsx
<Suspense fallback={null}>
  <ScoringResults /* existing props unchanged */ />
</Suspense>
```

(`fallback={null}` is deliberate: the chunks are small once Recharts is split out in Task 6, and a spinner here would flash.)

- [ ] **Step 4: Type-check, test, and verify in the browser**

```bash
npx tsc --noEmit && npm run test:run 2>&1 | tail -5
```

Expected: clean. Then in the dev server: validate `CCO`, open the scoring section and the standardize tab; results render (a brief blank gap on the very first open is acceptable).

- [ ] **Step 5: Check whether the recharts chunk still loads eagerly**

Hard-reload `/` with the Network tab open and filter "recharts". If `recharts-*.js` loads before any scoring UI is shown, another eagerly-imported component (likely `ScoreTiles` or `AllChecksCard`) imports a chart; note which via `grep -rn "recharts" src/components/validation/`, and record it in the PR description as a follow-up. Do not chase it in this task.

- [ ] **Step 6: Commit**

```bash
git add src/pages/SingleValidation.tsx
git commit -m "perf: lazy-load scoring and standardization result panels on home page"
```

### Task 8: Skip offscreen rendering work in BatchResultsTable

50 rows/page, each with an inline molecule SVG. `content-visibility: auto` lets the browser skip layout/paint for offscreen rows; it degrades gracefully where unsupported. (Full virtualization with react-virtuoso stays a deferred follow-up if profiling still shows jank afterward.)

**Files:**
- Modify: `frontend/src/index.css` (add utility)
- Modify: `frontend/src/components/batch/BatchResultsTable.tsx:458-465`

- [ ] **Step 1: Add the utility class**

At the end of `frontend/src/index.css`, append:

```css
/* ================================================
   RENDER-COST UTILITIES
   ================================================ */
@layer utilities {
  /* Skip layout/paint for offscreen table rows (batch results).
     contain-intrinsic-size approximates the collapsed row height so the
     scrollbar stays stable. No-op in browsers without support. */
  .cv-row {
    content-visibility: auto;
    contain-intrinsic-size: auto 56px;
  }
}
```

- [ ] **Step 2: Apply it to the result row**

In `frontend/src/components/batch/BatchResultsTable.tsx` line 460, add `cv-row` to the row's template-literal class:

```tsx
                    className={`
                      cv-row hover:bg-[var(--color-surface-sunken)] cursor-pointer transition-colors
                      ${result.status === 'error' ? 'bg-red-500/5' : ''}
                      ${expandedRow === result.index ? 'bg-[var(--color-primary)]/5' : ''}
                      ${result.index === focusedMoleculeIndex ? 'ring-2 ring-inset ring-[var(--color-primary)]/40' : ''}
                    `}
```

- [ ] **Step 3: Verify**

Run a batch job (upload a multi-molecule SMILES file on /batch) or open any existing results table. Scroll the 50-row page: scrolling is smooth; expanding a row still works; the focused-row scroll-into-view (keyboard navigation) still lands correctly.

- [ ] **Step 4: Commit**

```bash
git add src/index.css src/components/batch/BatchResultsTable.tsx
git commit -m "perf: skip offscreen row rendering in batch results table"
```

### Task 9: Calm the always-on ambient layer (orbs + footer)

Three full-viewport orbs at `blur-[90-120px]` animate scale+opacity in infinite loops on every page; the footer marquee stacks `backdrop-blur-md` and the aurora runs `blur-[80px]`. Reduce radii, drop the scale animation (opacity-only), and respect reduced-motion.

**Files:**
- Modify: `frontend/src/components/layout/Layout.tsx:1-62`
- Modify: `frontend/src/components/ui/motion-footer.tsx:353,364`

- [ ] **Step 1: Rework the orb configs and animation**

In `frontend/src/components/layout/Layout.tsx`, update the import (line 1):

```tsx
import { motion, useReducedMotion } from 'framer-motion';
```

Replace `orbConfigs` (lines 11-42), dropping the `scale` field and capping blur at 60px:

```tsx
// Refined ambient orb configurations.
// Opacity-only animation: animating scale on a blurred 600px element forces
// a Gaussian re-raster every frame. Blur capped at 60px for the same reason.
const orbConfigs = [
  {
    position: '-top-40 -right-40',
    size: 'w-[600px] h-[600px]',
    color: 'bg-[var(--color-primary)]',
    blur: 'blur-[60px]',
    opacity: [0.06, 0.1, 0.06],
    duration: 10,
    delay: 0,
  },
  {
    position: 'top-1/4 -left-40',
    size: 'w-[500px] h-[500px]',
    color: 'bg-[var(--color-accent)]',
    blur: 'blur-[60px]',
    opacity: [0.04, 0.08, 0.04],
    duration: 12,
    delay: 3,
  },
  {
    position: '-bottom-32 right-1/4',
    size: 'w-[450px] h-[450px]',
    color: 'bg-[var(--color-secondary)]',
    blur: 'blur-[56px]',
    opacity: [0.03, 0.06, 0.03],
    duration: 14,
    delay: 6,
  },
];
```

Then in the `Layout` component, read the reduced-motion preference and render opacity-only (replaces lines 48-62):

```tsx
export function Layout({ children }: LayoutProps) {
  const reduceMotion = useReducedMotion();

  return (
    <div className="relative w-full bg-[var(--color-surface)] min-h-screen">
      <div className="relative z-10 w-full bg-[var(--color-surface)]">
        {/* Ambient gradient orbs */}
        <div className="fixed inset-0 overflow-hidden pointer-events-none z-0">
          {orbConfigs.map((orb, i) => (
            <motion.div
              key={i}
              className={cn('absolute rounded-full', orb.position, orb.size, orb.color, orb.blur)}
              initial={{ opacity: orb.opacity[0] }}
              animate={reduceMotion ? { opacity: orb.opacity[0] } : { opacity: orb.opacity }}
              transition={{ duration: orb.duration, repeat: Infinity, ease: 'easeInOut', delay: orb.delay }}
            />
          ))}
        </div>
```

(The rest of the component is unchanged.)

- [ ] **Step 2: Trim the footer's paint cost**

In `frontend/src/components/ui/motion-footer.tsx`:

Line 353 — aurora blur 80px → 56px:

```tsx
          <div className="footer-aurora absolute left-1/2 top-1/2 h-[60vh] w-[80vw] -translate-x-1/2 -translate-y-1/2 animate-footer-breathe rounded-[50%] blur-[56px] pointer-events-none z-0" />
```

Line 364 — drop `backdrop-blur-md` from the marquee band and compensate with a more opaque background (the band sits over the static footer, so the blur bought nothing):

```tsx
          <div
            className="absolute top-12 left-0 w-full overflow-hidden py-4 z-10 -rotate-2 scale-110 shadow-2xl"
            style={{
              borderTop: '1px solid var(--color-border)',
              borderBottom: '1px solid var(--color-border)',
              backgroundColor: 'color-mix(in oklch, var(--color-surface) 88%, transparent)',
            }}
          >
```

- [ ] **Step 3: Verify visually in both themes**

Dev server: scroll to the footer on `/` in light and dark mode. The ambience reads the same at arm's length (orbs softer-edged is fine); marquee text remains legible over the giant background text. With macOS "Reduce Motion" enabled (System Settings → Accessibility → Display), orbs hold a static opacity.

- [ ] **Step 4: Commit**

```bash
git add src/components/layout/Layout.tsx src/components/ui/motion-footer.tsx
git commit -m "perf: opacity-only ambient orbs, capped blur radii, reduced-motion support"
```

### Task 10: Shrink the 764KB logo asset

**Files:**
- Create: `frontend/public/logo-512.png`, `frontend/public/logo-192.png`
- Modify: every `/logo.png` reference (Header, SplashScreen, ProcessingLogo, About, index.html favicon)

- [ ] **Step 1: Generate resized variants** (sips ships with macOS)

```bash
cd public
sips -Z 512 logo.png --out logo-512.png
sips -Z 192 logo.png --out logo-192.png
ls -la logo*.png
cd ..
```

Expected: `logo-512.png` and `logo-192.png` exist and are dramatically smaller than 764KB (typically under ~120KB and ~30KB). Keep the original `logo.png` in place so nothing breaks if a reference is missed.

- [ ] **Step 2: Repoint the references**

```bash
grep -rn '/logo.png' src/ index.html
```

Expected hits: `index.html:5` (favicon), `src/components/ui/SplashScreen.tsx`, `src/components/layout/Header.tsx`, `src/components/batch/ProcessingLogo.tsx`, `src/pages/About.tsx`. Apply this rule at each hit:

- `index.html` favicon and `Header.tsx` (rendered small) → `/logo-192.png`
- `SplashScreen.tsx`, `ProcessingLogo.tsx`, `About.tsx` (rendered large) → `/logo-512.png`

- [ ] **Step 3: Verify**

```bash
grep -rn '"/logo.png"' src/ index.html
```

Expected: zero matches. Dev server: splash logo, header logo, and About-page logo all render crisply; favicon shows in the tab.

- [ ] **Step 4: Commit**

```bash
git add public/logo-512.png public/logo-192.png src/ index.html
git commit -m "perf: serve resized logo variants instead of 764KB original"
```

---

## Phase C: Responsive

### Task 11: Add mobile breakpoints to the four broken grids

**Files:**
- Modify: `frontend/src/pages/SingleValidation.tsx:1213,1362`
- Modify: `frontend/src/components/batch/BatchProgress.tsx:192`
- Modify: `frontend/src/components/qsar-ready/QSARSingleResult.tsx:26,59`

- [ ] **Step 1: SingleValidation stats row (line 1213)** — five tiles at 375px ≈ 60px each, unusable:

```tsx
// before
                      <div className="grid grid-cols-5 gap-2 text-center">
// after
                      <div className="grid grid-cols-2 sm:grid-cols-3 md:grid-cols-5 gap-2 text-center">
```

- [ ] **Step 2: SingleValidation basic-info grid (line 1362)** — three tiny number tiles fit at 375px but are cramped; tighten the gap on mobile only (deliberate deviation from stacking, three one-line stats full-width would look worse):

```tsx
// before
                      <div className="grid grid-cols-3 gap-3 text-center">
// after
                      <div className="grid grid-cols-3 gap-2 sm:gap-3 text-center">
```

- [ ] **Step 3: BatchProgress stats (line 192)** — Processed/Total/ETA hold longer values, stack on phones:

```tsx
// before
        <div className="grid grid-cols-3 gap-4 py-4 px-2 border-t border-b border-[var(--color-border)] bg-gradient-to-r from-transparent via-[var(--color-surface-sunken)]/50 to-transparent rounded-lg">
// after
        <div className="grid grid-cols-1 sm:grid-cols-3 gap-4 py-4 px-2 border-t border-b border-[var(--color-border)] bg-gradient-to-r from-transparent via-[var(--color-surface-sunken)]/50 to-transparent rounded-lg">
```

- [ ] **Step 4: QSAR before/after comparison (line 59)** — two 200px structure panels cannot fit side-by-side at 375px:

```tsx
// before
      <div className="grid grid-cols-2 gap-6">
// after
      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
```

Also update the doc comment at line 26 so it stays truthful:

```tsx
// before
 * - grid grid-cols-2 gap-6 for before/after structure panels
// after
 * - grid grid-cols-1 md:grid-cols-2 gap-6 for before/after structure panels
```

- [ ] **Step 5: Verify at 375px**

Dev server, DevTools device toolbar at 375×812:
- Home page: validate `CCO`, expand Molecule Info; stat tiles wrap 2-per-row, no horizontal scrollbar.
- /batch during a job: progress stats stack vertically.
- /qsar-ready: process one molecule; before/after panels stack.

- [ ] **Step 6: Commit**

```bash
git add src/pages/SingleValidation.tsx src/components/batch/BatchProgress.tsx src/components/qsar-ready/QSARSingleResult.tsx
git commit -m "fix(responsive): add mobile breakpoints to stat and comparison grids"
```

### Task 12: Bring the hamburger button to a 44px touch target

**Files:**
- Modify: `frontend/src/components/layout/Header.tsx:603-613`

- [ ] **Step 1: Replace the padding-based sizing with explicit minimums**

In the mobile menu `motion.button` (line 606), change the first class string:

```tsx
// before
                'p-2 rounded-xl cursor-pointer',
// after
                'min-w-[44px] min-h-[44px] inline-flex items-center justify-center rounded-xl cursor-pointer',
```

- [ ] **Step 2: Verify**

At 375px width: the hamburger is comfortably tappable, vertically aligned with the theme toggle, and the open/close rotation animation still plays.

- [ ] **Step 3: Commit**

```bash
git add src/components/layout/Header.tsx
git commit -m "fix(responsive): enforce 44px touch target on mobile menu button"
```

---

## Phase D: Theming (token consolidation)

### Task 13: Create the shared chart-color module

Recharts/SVG `fill` props can't use Tailwind classes, which is why charts grew hard-coded hexes. Centralize the score-quality scale.

**Files:**
- Create: `frontend/src/lib/chartColors.ts`
- Test: `frontend/src/lib/chartColors.test.ts` (create)

- [ ] **Step 1: Write the failing test**

Create `frontend/src/lib/chartColors.test.ts`:

```ts
import { describe, it, expect } from 'vitest';
import { scoreFill, scoreBucket, chartTrackFill } from './chartColors';

describe('chartColors', () => {
  it('maps scores to the warm bucket hierarchy', () => {
    expect(scoreBucket(85)).toBe('excellent');
    expect(scoreBucket(80)).toBe('excellent');
    expect(scoreBucket(50)).toBe('good');
    expect(scoreBucket(20)).toBe('moderate');
    expect(scoreBucket(5)).toBe('poor');
  });

  it('returns theme-specific fills (gold for excellent, never green)', () => {
    expect(scoreFill('excellent', false)).toBe('#b45309');
    expect(scoreFill('excellent', true)).toBe('#fcd34d');
    expect(scoreFill('poor', false)).toBe('#dc2626');
  });

  it('uses warm stone for the gauge track', () => {
    expect(chartTrackFill(false)).toBe('#e7e5e4');
    expect(chartTrackFill(true)).toBe('#292524');
  });
});
```

- [ ] **Step 2: Run to verify failure**

```bash
npx vitest run src/lib/chartColors.test.ts
```

Expected: FAIL (module not found).

- [ ] **Step 3: Implement**

Create `frontend/src/lib/chartColors.ts`:

```ts
/**
 * Theme-aware chart color scales.
 *
 * Single source of truth for score-quality colors used in Recharts/SVG
 * fills, where Tailwind classes and CSS variables don't reach gradient
 * stops reliably. Values mirror the warm hierarchy in tailwind.config.js
 * (score-excellent/good/fair/poor): gold = excellent ... red = poor.
 * Success is NEVER green in this app (see PRODUCT.md, "Warm precision").
 */
export type ScoreBucket = 'excellent' | 'good' | 'moderate' | 'poor';

const SCORE_FILL: Record<ScoreBucket, { light: string; dark: string }> = {
  excellent: { light: '#b45309', dark: '#fcd34d' }, // amber-700 / amber-300
  good:      { light: '#d97706', dark: '#fbbf24' }, // amber-600 / amber-400
  moderate:  { light: '#ea580c', dark: '#fb923c' }, // orange-600 / orange-400
  poor:      { light: '#dc2626', dark: '#f87171' }, // red-600 / red-400
};

export function scoreFill(bucket: ScoreBucket, isDark: boolean): string {
  return SCORE_FILL[bucket][isDark ? 'dark' : 'light'];
}

/** Standard score → bucket mapping used by histograms and gauges. */
export function scoreBucket(score: number): ScoreBucket {
  if (score >= 80) return 'excellent';
  if (score >= 50) return 'good';
  if (score >= 20) return 'moderate';
  return 'poor';
}

/** Neutral track behind gauges/radials: warm stone, not cool gray. */
export function chartTrackFill(isDark: boolean): string {
  return isDark ? '#292524' : '#e7e5e4'; // stone-800 / stone-200
}
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
npx vitest run src/lib/chartColors.test.ts
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/lib/chartColors.ts src/lib/chartColors.test.ts
git commit -m "feat(theme): add shared theme-aware chart color module"
```

### Task 14: Migrate the batch score histogram off green/yellow

`ProfileScoreHistogram` in `BatchAnalyticsPanel.tsx` paints "Excellent" green (`#22c55e`), contradicting the gold-is-good brand scale.

**Files:**
- Modify: `frontend/src/components/batch/BatchAnalyticsPanel.tsx:172-186`

- [ ] **Step 1: Wire in the theme and shared scale**

Add the imports at the top of `frontend/src/components/batch/BatchAnalyticsPanel.tsx` (alongside existing imports):

```tsx
import { useThemeContext } from '../../contexts/ThemeContext';
import { scoreFill } from '../../lib/chartColors';
```

In `ProfileScoreHistogram` (line 172), first line of the function body:

```tsx
function ProfileScoreHistogram({ results }: { results: BatchResult[] }) {
  const { isDark } = useThemeContext();
```

- [ ] **Step 2: Replace the four hard-coded fills (lines 183-186)**

```tsx
// before
    { name: 'Excellent', count: buckets.Excellent, fill: '#22c55e' },
    { name: 'Good', count: buckets.Good, fill: '#eab308' },
    { name: 'Moderate', count: buckets.Moderate, fill: '#f97316' },
    { name: 'Poor', count: buckets.Poor, fill: '#ef4444' },
// after
    { name: 'Excellent', count: buckets.Excellent, fill: scoreFill('excellent', isDark) },
    { name: 'Good', count: buckets.Good, fill: scoreFill('good', isDark) },
    { name: 'Moderate', count: buckets.Moderate, fill: scoreFill('moderate', isDark) },
    { name: 'Poor', count: buckets.Poor, fill: scoreFill('poor', isDark) },
```

If `ProfileScoreHistogram` is rendered outside a `ThemeProvider` the hook throws; it isn't (the app root provides it), but if a test renders it directly, wrap with `ThemeProvider`.

- [ ] **Step 3: Verify**

```bash
npx tsc --noEmit
```

Then: run a batch job, open the analytics panel; the histogram reads gold → amber → orange → red in light mode, and brighter variants in dark mode, with no green bar.

- [ ] **Step 4: Commit**

```bash
git add src/components/batch/BatchAnalyticsPanel.tsx
git commit -m "fix(theme): batch score histogram uses warm brand scale, not green"
```

### Task 15: Warm gauge tracks in score charts

`ScoreChart` (and possibly the dataset-audit gauges) use cool grays `#374151`/`#e5e7eb` for the radial track, off-brand against warm stone surfaces.

**Files:**
- Modify: `frontend/src/components/scoring/ScoreChart.tsx:135` (+import)
- Possibly modify: `frontend/src/components/dataset-audit/HealthScoreGauge.tsx`, `frontend/src/components/dataset-audit/SubScoreCards.tsx`

- [ ] **Step 1: Find every cool-gray track**

```bash
grep -rn "#374151\|#e5e7eb\|#d1d5db\|#4b5563" src/components/ --include="*.tsx"
```

- [ ] **Step 2: Swap in `chartTrackFill`**

In `frontend/src/components/scoring/ScoreChart.tsx`, add to imports:

```tsx
import { chartTrackFill } from '../../lib/chartColors';
```

Replace line 135:

```tsx
// before
  const backgroundFill = isDark ? '#374151' : '#e5e7eb';
// after
  const backgroundFill = chartTrackFill(isDark);
```

Apply the identical substitution at every other hit from Step 1 that is a gauge/chart *track or grid* color (these components already have `isDark` from `useThemeContext`, per `grep -n "useThemeContext" <file>`). Leave hits that are not chart fills untouched and list them in the PR description.

- [ ] **Step 2b: Replace off-palette chart reference colors (emerald/green/blue)**

These are outside the brand palette entirely and have no dark-mode variants:

```bash
grep -rn "#059669\|#bbf7d0\|#93c5fd\|#22c55e" src/components/ --include="*.tsx" | grep -v test
```

Expected hits in `src/components/scoring/PropertyDistOverlay.tsx` (lines ~142-159: emerald `stroke="#059669"` / `fill="#059669"` on drug-space reference lines) and `src/components/dataset-audit/ComplexityRadar.tsx` (lines ~120-147: green `fill="#bbf7d0"` reference zone, blue `stroke="#93c5fd"` molecule polygon). Apply:

| Off-palette value | Replacement | Rationale |
|---|---|---|
| `#059669` (emerald reference line/label) | `var(--color-accent)` | amber = the brand's "good/target" signal |
| `#bbf7d0` (green reference zone fill) | `rgba(217, 119, 6, 0.12)` (light) / pass `isDark` and use `rgba(251, 191, 36, 0.12)` (dark) | tinted amber zone, both themes |
| `#93c5fd` (blue molecule polygon stroke) | `var(--color-primary)` | the molecule is the hero; give it the brand primary |

`var(--…)` values work in Recharts color props in this codebase (PropertyDistOverlay already uses `stroke="var(--color-border)"` on its CartesianGrid); use the rgba pairs only where a translucent fill is needed, with `isDark` from `useThemeContext()` as in Task 14.

Note: the score-tier hexes inside `getWarmColor`/`getCoolColor` (ScoreChart lines 48-115) and the amber/red hexes in `HealthScoreGauge.tsx`/`SubScoreCards.tsx` are brand-consistent and theme-switched already; centralizing them is churn without user-visible benefit. Add this comment above `getScoreColor` (line 40-42) instead:

```tsx
/**
 * Get color configuration based on score value.
 * Shared score-bucket fills live in src/lib/chartColors.ts; the gradient
 * stop pairs here are gauge-specific. Keep both in the warm spectrum.
 */
```

- [ ] **Step 3: Verify**

```bash
npx tsc --noEmit && npm run test:run 2>&1 | tail -5
```

Visual: validate a molecule, check the score gauges; the empty track reads warm-neutral in both themes.

- [ ] **Step 4: Commit**

```bash
git add src/components/scoring/ src/components/dataset-audit/ src/lib/chartColors.ts
git commit -m "fix(theme): warm gauge tracks and on-palette chart reference colors"
```

### Task 16: Map AlertGroupCard severity colors onto status tokens

**Files:**
- Modify: `frontend/src/components/safety/AlertGroupCard.tsx:18-29,60-65`

- [ ] **Step 1: Replace `getSeverityColor` (lines 18-29)**

```tsx
/** Severity strip color mapping. */
function getSeverityColor(severity: string): string {
  switch (severity) {
    case 'critical':
      return 'bg-status-error';
    case 'warning':
      return 'bg-status-warning';
    case 'info':
      return 'bg-status-info';
    default:
      return 'bg-chem-dark-400';
  }
}
```

(`status-info` maps to the crimson primary by design, see `tailwind.config.js:100-104`; blue was never a token.)

- [ ] **Step 2: Replace the `NibrBadge` colors (lines 60-65)**, which also lacked dark-mode variants entirely:

```tsx
  const colorClass =
    label === 'Excluded'
      ? 'bg-status-error-light text-status-error-dark dark:bg-status-error/15 dark:text-red-400'
      : label === 'Flag'
        ? 'bg-status-warning-light text-status-warning-dark dark:bg-status-warning/15 dark:text-amber-400'
        : 'bg-chem-dark-100 text-chem-dark-600 dark:bg-chem-dark-800 dark:text-chem-dark-300';
```

- [ ] **Step 3: Verify**

```bash
npx tsc --noEmit
```

Visual: validate a molecule with alerts (example: paste `O=[N+]([O-])c1ccc(Cl)cc1` and open the Alerts tab); severity dots/strips and NIBR badges render correctly in both themes.

- [ ] **Step 4: Commit**

```bash
git add src/components/safety/AlertGroupCard.tsx
git commit -m "fix(theme): AlertGroupCard severity colors use status tokens with dark variants"
```

### Task 17: Sweep `bg-white` / `border-gray-*` / green-success onto tokens

21 files use light-only raw palette classes that clash in dark mode. This is a mechanical mapping applied file-by-file; review each hit's context before substituting (e.g. `bg-white/80` overlay scrims and the decorative theme toggles are exceptions to leave alone).

**Files (audit-identified offenders):** `components/standardization/` (ProvenanceStageCard, ProvenanceTimeline, StandardizationResults, ComparisonView), `components/integrations/` (DatabaseComparisonPanel, DatabaseLookup), `components/batch/` (BatchSummary, BatchProgress, BatchTimeline, BatchResultsTable), plus ClusterMemberGrid, ClusterTable, TaxonomyTable, ScaffoldDisplay, ProfileComparisonView, MoleculeViewer, AlertResults, AlertCard.

**The mapping (memorize, then apply per hit):**

| Raw class | Token replacement |
|---|---|
| `bg-white` (panel/card background) | `bg-[var(--color-surface-elevated)]` |
| `bg-gray-50`, `bg-gray-100` | `bg-[var(--color-surface-sunken)]` |
| `bg-gray-900` (dark-only panel) | `bg-[var(--color-surface-elevated)]` |
| `border-gray-200`, `border-gray-200/80` | `border-[var(--color-border)]` |
| `border-gray-300` | `border-[var(--color-border-strong)]` |
| `border-gray-700` | `border-[var(--color-border)]` |
| `text-gray-500`, `text-gray-600` | `text-[var(--color-text-secondary)]` |
| `text-gray-400` | `text-[var(--color-text-muted)]` |
| `bg-gray-300`, `bg-gray-400` (dots/bars) | `bg-chem-dark-300` / `bg-chem-dark-400` |
| `bg-green-500`, `text-green-*` (success) | `bg-status-success` / `text-status-success-dark` (light) `dark:text-status-success` |
| `bg-emerald-50/20`, `stroke #059669` | `bg-status-success-light/20` / `var(--color-accent)` |

Token classes already render correctly in both themes (CSS variables flip with `.dark`), so most replacements DELETE the need for a `dark:` twin rather than adding one.

- [ ] **Step 1: Enumerate the standardization group**

```bash
grep -rn "bg-white\b\|bg-gray-\|border-gray-\|text-gray-" src/components/standardization/ | grep -v test
```

Apply the mapping to every hit. Then `npx tsc --noEmit` and visually check the Standardize tab (validate `CC(=O)Oc1ccccc1C(=O)O.[Na+].[Cl-]` and open standardization) in both themes.

```bash
git add src/components/standardization/
git commit -m "fix(theme): standardization components use surface/border tokens"
```

- [ ] **Step 2: Enumerate the integrations group**

```bash
grep -rn "bg-white\b\|bg-gray-\|border-gray-\|text-gray-\|emerald" src/components/integrations/ | grep -v test
```

Apply the mapping (keep per-database identity dot colors for now; Task 18 decides those). Visual check: Database tab lookup for `aspirin`. Commit:

```bash
git add src/components/integrations/
git commit -m "fix(theme): integrations panels use surface/border tokens"
```

- [ ] **Step 3: Enumerate batch + remaining components**

```bash
grep -rn "bg-white\b\|bg-gray-\|border-gray-\|bg-green-\|text-green-" src/components/batch/ src/components/scoring/ src/components/profiles/ src/components/molecules/ src/components/alerts/ src/pages/ | grep -v test | grep -v About
```

Apply the mapping (About.tsx is excluded: it is redesigned wholesale in the deferred design pass). Visual check: batch run summary/timeline in both themes. Commit:

```bash
git add src/
git commit -m "fix(theme): batch and shared components use tokens; success is amber not green"
```

- [ ] **Step 4: Final verification sweep**

```bash
grep -rn "bg-white\b\|border-gray-200\|bg-green-500" src/components/ | grep -v test | grep -v ThemeToggle | grep -v About
```

Expected: zero hits (ThemeToggle decorative styling and About.tsx are the allowed exceptions). Run `npm run test:run 2>&1 | tail -5`; no new failures.

### Task 18: Replace decorative purple/indigo accents with brand colors

The brand is crimson/amber/rose/stone. **Decision rule:** per-database identity colors in the integrations comparison UI (ChEMBL=violet, PubChem=blue dots that *distinguish data sources*) are a legitimate full-palette data-viz role and STAY. Decorative purples (callouts, spinners, category chips) move to brand colors.

**Files:**
- Modify: `frontend/src/components/alerts/AlertCard.tsx:45`
- Modify: `frontend/src/components/qsar-ready/QSARMoleculeRow.tsx:26`
- Modify: `frontend/src/components/integrations/DatabaseLookup.tsx:47`
- Modify: `frontend/src/components/integrations/IdentifierResolverCard.tsx:42-43`
- Modify: `frontend/src/components/batch/ProfileSidebar.tsx:87`
- Modify: `frontend/src/pages/About.tsx:404`
- Review: `frontend/src/components/batch/BatchResultsTable.tsx`, `frontend/src/components/scoring/DrugLikenessScore.tsx`, `frontend/src/components/scoring-profiles/LeadFragmentCard.tsx`, `frontend/src/components/validation/MoleculeViewerPanel.tsx`

- [ ] **Step 1: Apply the verified replacements** (quoted classes were grep-confirmed during the audit; confirm each line's context as you edit):

`AlertCard.tsx:45` — "Assay Interference" category chip:

```tsx
// before
'Assay Interference': { bg: 'bg-purple-100', text: 'text-purple-700' },
// after
'Assay Interference': { bg: 'bg-chem-secondary-100 dark:bg-chem-secondary-900/30', text: 'text-chem-secondary-700 dark:text-chem-secondary-300' },
```

`QSARMoleculeRow.tsx:26` — row highlight:

```tsx
// before:  bg-violet-50 text-violet-700
// after:   bg-chem-accent-50 dark:bg-chem-accent-900/20 text-chem-accent-700 dark:text-chem-accent-300
```

`DatabaseLookup.tsx:47` — spinner:

```tsx
// before:  border-indigo-600
// after:   border-[var(--color-primary)]
```

`IdentifierResolverCard.tsx:42-43` — info callout + chip:

```tsx
// before:  bg-indigo-50 border-indigo-200      /  bg-indigo-100 text-indigo-800
// after:   bg-[var(--color-surface-sunken)] border-[var(--color-border)]  /  bg-chem-primary-100 text-chem-primary-800 dark:bg-chem-primary-900/30 dark:text-chem-primary-300
```

`ProfileSidebar.tsx:87` — profile category color:

```tsx
// before
{ bg: 'bg-violet-500/12', text: 'text-violet-500', border: 'border-violet-500/20' }
// after
{ bg: 'bg-chem-secondary-500/12', text: 'text-chem-secondary-500', border: 'border-chem-secondary-500/20' }
```

`About.tsx:404` — off-brand gradient:

```tsx
// before:  from-indigo-500/20 to-blue-500/10
// after:   from-[var(--color-primary)]/15 to-[var(--color-accent)]/10
```

- [ ] **Step 2: Review the remaining purple hits**

```bash
grep -rn "purple-\|violet-\|indigo-\|fuchsia-" src/components/ src/pages/ | grep -v test | grep -v About
```

For each remaining hit, apply the decision rule: database-identity → keep; decorative/status → map to `chem-secondary-*` (rose) or `chem-primary-*` per the Step 1 patterns. `BatchResultsTable`'s `bg-purple-500/10 text-purple-600` and `MoleculeViewerPanel.tsx:58` are salt/fragment indicators — map to `bg-chem-secondary-500/10 text-chem-secondary-600 dark:text-chem-secondary-400` (and the same in `LeadFragmentCard.tsx:56`, `DrugLikenessScore.tsx`).

- [ ] **Step 3: Verify**

```bash
npx tsc --noEmit
grep -rn "indigo-" src/components/ src/pages/ | grep -v test | grep -v About
```

Expected: zero indigo hits outside About/tests. Visual: QSAR page, alerts tab, database tab in both themes; nothing reads as the stock AI purple.

- [ ] **Step 4: Commit**

```bash
git add src/
git commit -m "fix(theme): replace decorative purple/indigo with brand palette (database identity colors kept)"
```

### Task 19: Dark-mode molecule rendering contrast

`MoleculeViewer.tsx:36-41` strips RDKit's white background rect, but RDKit's black bonds/labels stay near-invisible on dark surfaces.

**Files:**
- Modify: `frontend/src/index.css` (molecule-preview block, ~line 713)
- Verify against: `frontend/src/components/molecules/MoleculeViewer.tsx`

- [ ] **Step 1: Confirm the wrapper class the viewer actually uses**

```bash
grep -n "molecule-preview\|className" src/components/molecules/MoleculeViewer.tsx | head -20
```

Expected: the rendered SVG sits inside an element carrying `molecule-preview` (the class is defined at `index.css:714`). If the viewer uses a different wrapper class, target that class in Step 2 instead.

- [ ] **Step 2: Add a dark-mode inversion filter**

In `frontend/src/index.css`, inside the MOLECULE VIEWER `@layer components` block (after the `.molecule-preview > *` rule at ~line 722), add:

```css
  /* RDKit draws black bonds/labels for white backgrounds. In dark mode,
     invert the drawing; hue-rotate(180deg) flips the inverted hues back so
     heteroatom colors (red O, blue N) stay recognizable. */
  .dark .molecule-preview svg {
    filter: invert(0.88) hue-rotate(180deg) saturate(1.35);
  }
```

- [ ] **Step 3: Verify carefully (this one is judgment-sensitive)**

In dark mode, validate molecules covering the common cases:
- `CCO` (plain C/O skeleton): bonds clearly visible, near-white.
- `c1ccncc1` (pyridine, blue N): nitrogen still reads blue-ish.
- A validation result with highlighted atoms (any molecule with issues, e.g. invalid valence example from the Examples list): highlight color must remain visible. If highlights wash out, reduce inversion to `invert(0.82)`.

Light mode must be unchanged.

- [ ] **Step 4: Commit**

```bash
git add src/index.css
git commit -m "fix(theme): legible molecule renderings in dark mode via inversion filter"
```

---

## Wrap-up

- [ ] **Full verification pass**

```bash
npx tsc --noEmit
npm run lint
npm run test:run
npm run build && ls -la dist/assets/ | sort -k5 -n | tail -12
```

Expected: all clean; record final main-chunk size vs the Setup baseline in the PR description.

- [ ] **Manual cross-check of the audit's P1 list** (light + dark, desktop + 375px): muted text legible, export dialog labeled, badges distinguishable, /about loads instantly without WASM, batch table smooth, grids wrap, no green success anywhere, no white panels in dark mode.

- [ ] **Push and open PR against `development`**

```bash
git push -u origin fix/frontend-audit-remediation
```

PR body should list the audit score (11/20), per-dimension fixes, and the deferred items below. End the body with the standard generated-with footer; do NOT add Co-Authored-By.

---

## Explicitly deferred (not in this plan)

| Item | Why deferred | Route |
|---|---|---|
| About.tsx redesign (hero-metric block + 4 identical card grids) | Design work, not mechanical remediation; needs a shape pass | `/impeccable distill About.tsx` |
| ClayButton gradient reduction (6 gradient variants → 2) | Design-system decision affecting every screen | `/impeccable distill` or `extract` |
| react-virtuoso table virtualization | Re-profile after Task 8 (`content-visibility`); only if jank persists | Follow-up perf task |
| `ChemicalSpaceScatter.tsx` `getBoundingClientRect` caching (lines 306-425) | P2, already debounced 50ms; localized | Follow-up perf task |
| Fixed chart heights at 200% zoom (`IssueTreemap`, `PropertyDistOverlay`) | P2 polish | Follow-up |
| ChemThemeToggle teal/slate styling | Decorative, component-scoped, acceptable | None |
| GSAP footer ScrollTrigger/MagneticButton deep rework | Needs profiling evidence before touching working animation code | Follow-up perf task |
| Dark-mode variants for database identity colors (integrations) | Depends on Task 18's keep-decision; small | Follow-up |
