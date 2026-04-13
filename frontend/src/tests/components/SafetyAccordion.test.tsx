/**
 * SafetyAccordion Component Tests
 *
 * Tests the wrapper component that renders all safety sub-components
 * inside a DrillDownSection accordion for SingleValidation.
 *
 * All child components are mocked to avoid deep dependency trees.
 */
import { describe, it, expect, vi } from 'vitest';
import { render, screen } from '../setup';
import { SafetyAccordion } from '../../components/safety/SafetyAccordion';
import type { AlertScreenResponse, SafetyAssessResponse } from '../../types/safety';

// Mock framer-motion
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

// Mock all safety child components
vi.mock('../../components/safety/SafetySummaryStrip', () => ({
  SafetySummaryStrip: (props: Record<string, unknown>) => (
    <div data-testid="safety-summary-strip">
      SafetySummaryStrip: {String(props.totalAlerts)} alerts
    </div>
  ),
}));

vi.mock('../../components/molecules/MoleculeViewer', () => ({
  MoleculeViewer: (props: Record<string, unknown>) => (
    <div data-testid="molecule-viewer">
      MoleculeViewer: {String(props.smiles)}
    </div>
  ),
}));

vi.mock('../../components/safety/AlertDashboard', () => ({
  AlertDashboard: () => (
    <div data-testid="alert-dashboard">AlertDashboard</div>
  ),
}));

vi.mock('../../components/safety/SafetyFlagsGrid', () => ({
  SafetyFlagsGrid: () => (
    <div data-testid="safety-flags-grid">SafetyFlagsGrid</div>
  ),
}));

vi.mock('../../components/safety/ComplexityRadar', () => ({
  ComplexityRadar: () => (
    <div data-testid="complexity-radar">ComplexityRadar</div>
  ),
}));

/**
 * Build a minimal valid AlertScreenResponse for testing.
 */
function buildMockAlertResult(): AlertScreenResponse {
  return {
    status: 'success',
    molecule_info: {
      smiles: 'CC(=O)Oc1ccccc1C(=O)O',
    },
    alerts: [
      {
        pattern_name: 'acyl_halide',
        description: 'Acyl halide detected',
        severity: 'warning',
        matched_atoms: [0, 1, 2],
        catalog_source: 'PAINS',
      },
    ],
    concern_groups: {
      'Reactive Groups': {
        name: 'Reactive Groups',
        count: 1,
        severity: 'warning',
        alerts: [
          {
            pattern_name: 'acyl_halide',
            description: 'Acyl halide detected',
            severity: 'warning',
            matched_atoms: [0, 1, 2],
            catalog_source: 'PAINS',
          },
        ],
      },
    },
    total_raw: 1,
    total_deduped: 1,
    screened_catalogs: ['PAINS', 'BMS'],
    has_critical: false,
    has_warning: true,
    execution_time_ms: 42,
  };
}

/**
 * Build a minimal valid SafetyAssessResponse for testing.
 */
function buildMockSafetyResult(): SafetyAssessResponse {
  return {
    status: 'success',
    molecule_info: {
      smiles: 'CC(=O)Oc1ccccc1C(=O)O',
    },
    cyp_softspots: {
      sites: [
        {
          site_name: 'CYP3A4-site1',
          reaction_type: 'oxidation',
          matched_atoms: [3, 4],
        },
      ],
      n_sites: 1,
    },
    herg: {
      herg_risk: 'low',
      risk_score: 1,
      max_score: 10,
      flags: [],
      descriptors: { logP: 2.1 },
    },
    bro5: {
      applicable: true,
      passed: true,
      violations: [],
      values: { mw: 250 },
    },
    reos: {
      passed: true,
      violations: [],
      n_violations: 0,
      descriptors: { mw: 250 },
    },
    complexity: {
      properties: {
        'Fsp3': { value: 0.3, p5: 0.1, p95: 0.8, outlier: false, direction: null },
        'NumHBA': { value: 4, p5: 1, p95: 10, outlier: false, direction: null },
      },
      n_outliers: 0,
      outlier_properties: [],
      within_range: true,
    },
    execution_time_ms: 55,
  };
}

describe('SafetyAccordion', () => {
  const defaultSmiles = 'CC(=O)Oc1ccccc1C(=O)O';

  describe('Data loaded state', () => {
    it('renders all safety components when data provided', () => {
      const alertResult = buildMockAlertResult();
      const safetyResult = buildMockSafetyResult();

      render(
        <SafetyAccordion
          smiles={defaultSmiles}
          alertResult={alertResult}
          safetyResult={safetyResult}
          isLoading={false}
          error={null}
        />
      );

      expect(screen.getByTestId('safety-summary-strip')).toBeInTheDocument();
      expect(screen.getByTestId('molecule-viewer')).toBeInTheDocument();
      expect(screen.getByTestId('alert-dashboard')).toBeInTheDocument();
      expect(screen.getByTestId('safety-flags-grid')).toBeInTheDocument();
      expect(screen.getByTestId('complexity-radar')).toBeInTheDocument();
    });

    it('passes correct totalAlerts to SafetySummaryStrip', () => {
      const alertResult = buildMockAlertResult();
      const safetyResult = buildMockSafetyResult();

      render(
        <SafetyAccordion
          smiles={defaultSmiles}
          alertResult={alertResult}
          safetyResult={safetyResult}
          isLoading={false}
          error={null}
        />
      );

      expect(screen.getByTestId('safety-summary-strip')).toHaveTextContent('1 alerts');
    });

    it('passes smiles to MoleculeViewer', () => {
      const alertResult = buildMockAlertResult();
      const safetyResult = buildMockSafetyResult();

      render(
        <SafetyAccordion
          smiles={defaultSmiles}
          alertResult={alertResult}
          safetyResult={safetyResult}
          isLoading={false}
          error={null}
        />
      );

      expect(screen.getByTestId('molecule-viewer')).toHaveTextContent(defaultSmiles);
    });
  });

  describe('Loading state', () => {
    it('shows loading skeleton when isLoading is true', () => {
      const { container } = render(
        <SafetyAccordion
          smiles={defaultSmiles}
          alertResult={null}
          safetyResult={null}
          isLoading={true}
          error={null}
        />
      );

      // Should show skeleton placeholders
      const skeletons = container.querySelectorAll('.animate-pulse');
      expect(skeletons.length).toBeGreaterThan(0);
    });

    it('does not render safety panels when loading', () => {
      render(
        <SafetyAccordion
          smiles={defaultSmiles}
          alertResult={null}
          safetyResult={null}
          isLoading={true}
          error={null}
        />
      );

      expect(screen.queryByTestId('safety-summary-strip')).not.toBeInTheDocument();
      expect(screen.queryByTestId('alert-dashboard')).not.toBeInTheDocument();
      expect(screen.queryByTestId('safety-flags-grid')).not.toBeInTheDocument();
      expect(screen.queryByTestId('complexity-radar')).not.toBeInTheDocument();
    });
  });

  describe('Error state', () => {
    it('shows error message when error is provided', () => {
      render(
        <SafetyAccordion
          smiles={defaultSmiles}
          alertResult={null}
          safetyResult={null}
          isLoading={false}
          error="Safety screening failed"
        />
      );

      expect(screen.getByText('Safety screening failed')).toBeInTheDocument();
    });

    it('does not render safety panels when error present', () => {
      render(
        <SafetyAccordion
          smiles={defaultSmiles}
          alertResult={null}
          safetyResult={null}
          isLoading={false}
          error="Something went wrong"
        />
      );

      expect(screen.queryByTestId('safety-summary-strip')).not.toBeInTheDocument();
      expect(screen.queryByTestId('alert-dashboard')).not.toBeInTheDocument();
    });
  });

  describe('Empty state', () => {
    it('shows empty message when no smiles provided', () => {
      render(
        <SafetyAccordion
          smiles=""
          alertResult={null}
          safetyResult={null}
          isLoading={false}
          error={null}
        />
      );

      expect(
        screen.getByText('Enter a molecule to assess safety')
      ).toBeInTheDocument();
    });

    it('does not render safety panels when no smiles', () => {
      render(
        <SafetyAccordion
          smiles=""
          alertResult={null}
          safetyResult={null}
          isLoading={false}
          error={null}
        />
      );

      expect(screen.queryByTestId('safety-summary-strip')).not.toBeInTheDocument();
      expect(screen.queryByTestId('alert-dashboard')).not.toBeInTheDocument();
      expect(screen.queryByTestId('safety-flags-grid')).not.toBeInTheDocument();
      expect(screen.queryByTestId('complexity-radar')).not.toBeInTheDocument();
    });
  });

  describe('Partial data state', () => {
    it('renders alerts and loading skeleton when alertResult present but safetyResult loading', () => {
      const alertResult = buildMockAlertResult();
      const { container } = render(
        <SafetyAccordion
          smiles={defaultSmiles}
          alertResult={alertResult}
          safetyResult={null}
          isLoading={false}
          error={null}
        />
      );

      // Alert-specific components should render
      expect(screen.getByTestId('alert-dashboard')).toBeInTheDocument();
      expect(screen.getByTestId('molecule-viewer')).toBeInTheDocument();

      // Safety flags and complexity should show skeletons (not full components)
      expect(screen.queryByTestId('safety-flags-grid')).not.toBeInTheDocument();
      expect(screen.queryByTestId('complexity-radar')).not.toBeInTheDocument();

      // Skeleton placeholders for safety flags section
      const skeletons = container.querySelectorAll('.animate-pulse');
      expect(skeletons.length).toBeGreaterThan(0);
    });
  });
});
