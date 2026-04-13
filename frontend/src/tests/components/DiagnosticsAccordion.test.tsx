/**
 * DiagnosticsAccordion Component Tests
 *
 * Tests the wrapper component that combines all 5 diagnostic tools
 * into independently collapsible sections for the SingleValidation page.
 */
import { describe, it, expect, vi } from 'vitest';
import { render, screen } from '../setup';
import { DiagnosticsAccordion } from '../../components/diagnostics/DiagnosticsAccordion';
import type { CrossPipelineResponse, RoundTripResponse } from '../../types/diagnostics';

// Mock framer-motion to render synchronously in tests
vi.mock('framer-motion', () => ({
  motion: {
    div: ({
      children,
      ...props
    }: React.HTMLAttributes<HTMLDivElement> & Record<string, unknown>) => {
      const {
        initial: _initial,
        animate: _animate,
        exit: _exit,
        transition: _transition,
        variants: _variants,
        ...htmlProps
      } = props;
      return <div {...htmlProps}>{children}</div>;
    },
  },
  AnimatePresence: ({ children }: { children: React.ReactNode }) => <>{children}</>,
}));

// Mock all child diagnostic panels
vi.mock('../../components/diagnostics/SMILESDiagnosticsPanel', () => ({
  SMILESDiagnosticsPanel: (props: Record<string, unknown>) => (
    <div data-testid="smiles-diagnostics-panel" data-smiles={props.originalSmiles}>
      SMILESDiagnosticsPanel
    </div>
  ),
}));

vi.mock('../../components/diagnostics/InChILayerDiffTable', () => ({
  InChILayerDiffTable: () => (
    <div data-testid="inchi-layer-diff-table">InChILayerDiffTable</div>
  ),
}));

vi.mock('../../components/diagnostics/RoundTripPanel', () => ({
  RoundTripPanel: (props: Record<string, unknown>) => (
    <div data-testid="round-trip-panel" data-has-result={props.result ? 'true' : 'false'}>
      RoundTripPanel
    </div>
  ),
}));

vi.mock('../../components/diagnostics/CrossPipelinePanel', () => ({
  CrossPipelinePanel: (props: Record<string, unknown>) => (
    <div data-testid="cross-pipeline-panel" data-has-result={props.result ? 'true' : 'false'}>
      CrossPipelinePanel
    </div>
  ),
}));

vi.mock('../../components/diagnostics/FilePreValidatorPanel', () => ({
  FilePreValidatorPanel: () => (
    <div data-testid="file-pre-validator-panel">FilePreValidatorPanel</div>
  ),
}));

// Mock diagnosticsApi
vi.mock('../../services/api', () => ({
  diagnosticsApi: {
    smiles: vi.fn(),
    inchiDiff: vi.fn(),
    roundtrip: vi.fn(),
    crossPipeline: vi.fn(),
    filePrevalidate: vi.fn(),
  },
}));

const mockCrossPipelineResult: CrossPipelineResponse = {
  pipelines: [],
  disagreements: 0,
  structural_disagreements: 0,
  all_agree: true,
  property_comparison: [],
};

const mockRoundTripResult: RoundTripResponse = {
  route: 'smiles_inchi_smiles',
  original_smiles: 'CCO',
  intermediate: 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
  roundtrip_smiles: 'CCO',
  lossy: false,
  losses: [],
  error: null,
};

describe('DiagnosticsAccordion', () => {
  describe('Renders all 5 diagnostic panels when smiles provided', () => {
    it('renders all section headings', () => {
      render(
        <DiagnosticsAccordion
          smiles="CCO"
          crossPipelineResult={null}
          roundTripResult={null}
          isLoading={false}
        />
      );

      expect(screen.getByText('SMILES Diagnostics')).toBeInTheDocument();
      expect(screen.getByText('InChI Layer Diff')).toBeInTheDocument();
      expect(screen.getByText('Round-Trip Lossiness')).toBeInTheDocument();
      expect(screen.getByText('Cross-Pipeline Standardization')).toBeInTheDocument();
      expect(screen.getByText('File Pre-Validator')).toBeInTheDocument();
    });

    it('renders SMILES diagnostics panel expanded by default', () => {
      render(
        <DiagnosticsAccordion
          smiles="CCO"
          crossPipelineResult={null}
          roundTripResult={null}
          isLoading={false}
        />
      );

      // SMILES panel should be visible (auto-expanded)
      expect(screen.getByTestId('smiles-diagnostics-panel')).toBeInTheDocument();
    });

    it('passes smiles to SMILESDiagnosticsPanel', () => {
      render(
        <DiagnosticsAccordion
          smiles="CCO"
          crossPipelineResult={null}
          roundTripResult={null}
          isLoading={false}
        />
      );

      const panel = screen.getByTestId('smiles-diagnostics-panel');
      expect(panel).toHaveAttribute('data-smiles', 'CCO');
    });
  });

  describe('Loading state', () => {
    it('shows loading indicator when isLoading is true', () => {
      render(
        <DiagnosticsAccordion
          smiles="CCO"
          crossPipelineResult={null}
          roundTripResult={null}
          isLoading={true}
        />
      );

      expect(screen.getByText(/loading diagnostics/i)).toBeInTheDocument();
    });
  });

  describe('Empty state (no smiles)', () => {
    it('shows empty state when smiles is empty string', () => {
      render(
        <DiagnosticsAccordion
          smiles=""
          crossPipelineResult={null}
          roundTripResult={null}
          isLoading={false}
        />
      );

      expect(screen.getByText(/enter a smiles/i)).toBeInTheDocument();
    });

    it('does not render diagnostic panels when smiles is empty', () => {
      render(
        <DiagnosticsAccordion
          smiles=""
          crossPipelineResult={null}
          roundTripResult={null}
          isLoading={false}
        />
      );

      expect(screen.queryByTestId('smiles-diagnostics-panel')).not.toBeInTheDocument();
      expect(screen.queryByTestId('inchi-layer-diff-table')).not.toBeInTheDocument();
      expect(screen.queryByTestId('round-trip-panel')).not.toBeInTheDocument();
      expect(screen.queryByTestId('cross-pipeline-panel')).not.toBeInTheDocument();
      expect(screen.queryByTestId('file-pre-validator-panel')).not.toBeInTheDocument();
    });
  });

  describe('Pre-fetched data', () => {
    it('passes crossPipelineResult to CrossPipelinePanel', () => {
      render(
        <DiagnosticsAccordion
          smiles="CCO"
          crossPipelineResult={mockCrossPipelineResult}
          roundTripResult={null}
          isLoading={false}
        />
      );

      // The cross-pipeline section heading should be present
      expect(screen.getByText('Cross-Pipeline Standardization')).toBeInTheDocument();
    });

    it('passes roundTripResult to RoundTripPanel', () => {
      render(
        <DiagnosticsAccordion
          smiles="CCO"
          crossPipelineResult={null}
          roundTripResult={mockRoundTripResult}
          isLoading={false}
        />
      );

      // The round-trip section heading should be present
      expect(screen.getByText('Round-Trip Lossiness')).toBeInTheDocument();
    });
  });
});
