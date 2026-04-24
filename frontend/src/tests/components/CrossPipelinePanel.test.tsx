/**
 * CrossPipelinePanel Tests
 *
 * Covers the fixed-aspect letterboxed thumbnails, the MCS highlight wiring
 * (highlight_atoms / highlight_bonds passed through to MoleculeViewer),
 * pipeline-level error rendering, and auto-trigger behaviour.
 */
import React from 'react';
import { describe, it, expect, vi, beforeEach } from 'vitest';
import { render, screen } from '../setup';
import { CrossPipelinePanel } from '../../components/diagnostics/CrossPipelinePanel';
import type { CrossPipelineResponse, PipelineResult } from '../../types/diagnostics';

// Strip framer-motion-specific props so they don't leak onto real DOM nodes
// (initial/animate/etc. would otherwise trigger React "unknown prop" warnings).
const MOTION_ONLY_PROPS = new Set([
  'initial',
  'animate',
  'exit',
  'transition',
  'variants',
  'whileHover',
  'whileTap',
  'whileFocus',
  'whileInView',
  'layout',
  'layoutId',
]);

function stripMotionProps(props: Record<string, unknown>) {
  const out: Record<string, unknown> = {};
  for (const [k, v] of Object.entries(props)) {
    if (!MOTION_ONLY_PROPS.has(k)) out[k] = v;
  }
  return out;
}

vi.mock('framer-motion', () => {
  const makeComponent =
    (tag: string) =>
    ({ children, ...props }: Record<string, unknown>) => {
      const safe = stripMotionProps(props);
      return React.createElement(tag, safe, children as React.ReactNode);
    };
  // Proxy covers motion.div / motion.button / motion.svg / etc. uniformly.
  const motion = new Proxy({}, { get: (_t, tag: string) => makeComponent(tag) });
  return {
    motion,
    AnimatePresence: ({ children }: { children: React.ReactNode }) => <>{children}</>,
  };
});

// MoleculeViewer requires RDKit WASM — stub it and surface the highlight props
// as data attributes so the tests can assert on them directly.
vi.mock('../../components/molecules/MoleculeViewer', () => ({
  MoleculeViewer: (props: {
    smiles: string | null;
    highlightAtoms?: number[];
    highlightBonds?: number[];
    fit?: string;
    showAtomLabels?: boolean;
  }) => (
    <div
      data-testid="molecule-viewer"
      data-smiles={props.smiles ?? ''}
      data-fit={props.fit ?? 'width'}
      data-show-atom-labels={String(props.showAtomLabels ?? true)}
      data-highlight-atoms={JSON.stringify(props.highlightAtoms ?? [])}
      data-highlight-bonds={JSON.stringify(props.highlightBonds ?? [])}
    />
  ),
}));

function makePipeline(overrides: Partial<PipelineResult> = {}): PipelineResult {
  return {
    name: 'RDKit MolStandardize',
    smiles: 'CCO',
    inchikey: 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N',
    mw: 46.0419,
    formula: 'C2H6O',
    charge: 0,
    stereo_count: 0,
    highlight_atoms: [],
    highlight_bonds: [],
    ...overrides,
  };
}

function makeResponse(overrides: Partial<CrossPipelineResponse> = {}): CrossPipelineResponse {
  const pipelines = overrides.pipelines ?? [
    makePipeline({ name: 'RDKit MolStandardize' }),
    makePipeline({ name: 'ChEMBL Pipeline' }),
    makePipeline({ name: 'Minimal (Sanitize Only)' }),
  ];
  return {
    pipelines,
    disagreements: 0,
    structural_disagreements: 0,
    all_agree: true,
    property_comparison: [],
    ...overrides,
  };
}

describe('CrossPipelinePanel', () => {
  const noopCompare = vi.fn();
  const noopRetry = vi.fn();

  beforeEach(() => {
    noopCompare.mockReset();
    noopRetry.mockReset();
  });

  describe('loading state', () => {
    it('renders 3 skeleton thumbnails while loading', () => {
      const { container } = render(
        <CrossPipelinePanel
          result={null}
          isLoading={true}
          error={null}
          currentSmiles="CCO"
          onComparePipelines={noopCompare}
          onRetry={noopRetry}
        />
      );
      // The skeleton grid uses 3 child divs; assert by counting rounded skeletons
      const rounded = container.querySelectorAll('[class*="rounded"]');
      expect(rounded.length).toBeGreaterThan(0);
      // No MoleculeViewers should be mounted during loading
      expect(screen.queryAllByTestId('molecule-viewer')).toHaveLength(0);
    });
  });

  describe('error state', () => {
    it('renders failure card with retry button', () => {
      render(
        <CrossPipelinePanel
          result={null}
          isLoading={false}
          error="Request failed"
          currentSmiles="CCO"
          onComparePipelines={noopCompare}
          onRetry={noopRetry}
        />
      );
      // Heading + subtext both include the phrase; use getAllByText and
      // assert the retry button separately.
      expect(screen.getAllByText(/Pipeline comparison failed/i).length).toBeGreaterThan(0);
      expect(screen.getByRole('button', { name: /retry/i })).toBeInTheDocument();
    });
  });

  describe('auto-trigger', () => {
    it('calls onComparePipelines once on mount when currentSmiles is set and result is null', () => {
      render(
        <CrossPipelinePanel
          result={null}
          isLoading={false}
          error={null}
          currentSmiles="CCO"
          onComparePipelines={noopCompare}
          onRetry={noopRetry}
        />
      );
      expect(noopCompare).toHaveBeenCalledTimes(1);
      expect(noopCompare).toHaveBeenCalledWith('CCO');
    });

    it('does not auto-trigger when a result is already present', () => {
      render(
        <CrossPipelinePanel
          result={makeResponse()}
          isLoading={false}
          error={null}
          currentSmiles="CCO"
          onComparePipelines={noopCompare}
          onRetry={noopRetry}
        />
      );
      expect(noopCompare).not.toHaveBeenCalled();
    });
  });

  describe('rendered result — all agree', () => {
    it('renders 3 figure thumbnails with captions', () => {
      render(
        <CrossPipelinePanel
          result={makeResponse()}
          isLoading={false}
          error={null}
          currentSmiles="CCO"
          onComparePipelines={noopCompare}
          onRetry={noopRetry}
        />
      );
      const viewers = screen.getAllByTestId('molecule-viewer');
      expect(viewers).toHaveLength(3);
      // Captions live inside <figcaption> — assert there rather than a plain
      // text match, since the truncated label form ("ChEMBL-style") can also
      // appear in the comparison table header.
      const captions = document.querySelectorAll('figcaption');
      expect(captions).toHaveLength(3);
      expect(captions[0].textContent).toContain('RDKit MolStandardize');
      expect(captions[1].textContent).toContain('ChEMBL-style');
      expect(captions[2].textContent).toContain('Minimal Sanitize');
    });

    it('uses letterbox fit and suppresses atom labels for each viewer', () => {
      render(
        <CrossPipelinePanel
          result={makeResponse()}
          isLoading={false}
          error={null}
          currentSmiles="CCO"
          onComparePipelines={noopCompare}
          onRetry={noopRetry}
        />
      );
      for (const viewer of screen.getAllByTestId('molecule-viewer')) {
        expect(viewer.getAttribute('data-fit')).toBe('contain');
        expect(viewer.getAttribute('data-show-atom-labels')).toBe('false');
      }
    });

    it('passes empty highlight arrays when pipelines agree', () => {
      render(
        <CrossPipelinePanel
          result={makeResponse()}
          isLoading={false}
          error={null}
          currentSmiles="CCO"
          onComparePipelines={noopCompare}
          onRetry={noopRetry}
        />
      );
      for (const viewer of screen.getAllByTestId('molecule-viewer')) {
        expect(viewer.getAttribute('data-highlight-atoms')).toBe('[]');
        expect(viewer.getAttribute('data-highlight-bonds')).toBe('[]');
      }
    });
  });

  describe('rendered result — structural disagreement', () => {
    const disagreeingResponse = makeResponse({
      all_agree: false,
      disagreements: 2,
      structural_disagreements: 2,
      pipelines: [
        makePipeline({ name: 'RDKit MolStandardize', smiles: 'CC(=O)O' }),
        makePipeline({
          name: 'ChEMBL Pipeline',
          smiles: 'CC(=O)[O-].[Na+]',
          highlight_atoms: [4],
          highlight_bonds: [],
        }),
        makePipeline({
          name: 'Minimal (Sanitize Only)',
          smiles: 'CC(=O)[O-].[Na+]',
          highlight_atoms: [4],
          highlight_bonds: [],
        }),
      ],
    });

    it('forwards highlight_atoms and highlight_bonds to each diverging viewer', () => {
      render(
        <CrossPipelinePanel
          result={disagreeingResponse}
          isLoading={false}
          error={null}
          currentSmiles="CC(=O)[O-].[Na+]"
          onComparePipelines={noopCompare}
          onRetry={noopRetry}
        />
      );
      const viewers = screen.getAllByTestId('molecule-viewer');
      // RDKit pipeline — empty highlights (agrees with MCS)
      expect(viewers[0].getAttribute('data-highlight-atoms')).toBe('[]');
      // ChEMBL + Minimal — both flag atom index 4 (the Na+ atom)
      expect(viewers[1].getAttribute('data-highlight-atoms')).toBe('[4]');
      expect(viewers[2].getAttribute('data-highlight-atoms')).toBe('[4]');
    });

    it('shows the structural-disagreement badge', () => {
      render(
        <CrossPipelinePanel
          result={disagreeingResponse}
          isLoading={false}
          error={null}
          currentSmiles="CC(=O)[O-].[Na+]"
          onComparePipelines={noopCompare}
          onRetry={noopRetry}
        />
      );
      expect(screen.getByText(/structural disagreement/i)).toBeInTheDocument();
    });
  });

  describe('pipeline-level error', () => {
    it('renders the pipeline error message in place of the molecule viewer', () => {
      const errored = makeResponse({
        all_agree: false,
        disagreements: 6,
        structural_disagreements: 2,
        pipelines: [
          makePipeline({ name: 'RDKit MolStandardize' }),
          makePipeline({ name: 'ChEMBL Pipeline', error: 'Standardization failed' }),
          makePipeline({ name: 'Minimal (Sanitize Only)' }),
        ],
      });
      render(
        <CrossPipelinePanel
          result={errored}
          isLoading={false}
          error={null}
          currentSmiles="CCO"
          onComparePipelines={noopCompare}
          onRetry={noopRetry}
        />
      );
      expect(screen.getByText('Standardization failed')).toBeInTheDocument();
      // Only the 2 healthy pipelines render viewers; the errored one falls back
      // to the error text, so the total viewer count is 2.
      expect(screen.getAllByTestId('molecule-viewer')).toHaveLength(2);
    });
  });
});
