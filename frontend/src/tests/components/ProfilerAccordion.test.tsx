/**
 * ProfilerAccordion Component Tests
 *
 * Tests the wrapper component that renders all profiler sub-components
 * inside a DrillDownSection accordion for SingleValidation.
 *
 * All child components are mocked to avoid deep dependency trees.
 */
import { describe, it, expect, vi } from 'vitest';
import { render, screen } from '../setup';
import { ProfilerAccordion } from '../../components/profiler/ProfilerAccordion';
import type { ProfileResponse } from '../../types/profiler';

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

// Mock all profiler child components
vi.mock('../../components/profiler/HeroSection', () => ({
  HeroSection: ({ smiles }: { smiles: string }) => (
    <div data-testid="hero-section">HeroSection: {smiles}</div>
  ),
}));

vi.mock('../../components/profiler/PFIPanel', () => ({
  PFIPanel: () => <div data-testid="pfi-panel">PFIPanel</div>,
}));

vi.mock('../../components/profiler/StarsPanel', () => ({
  StarsPanel: () => <div data-testid="stars-panel">StarsPanel</div>,
}));

vi.mock('../../components/profiler/BioavailabilityPanel', () => ({
  BioavailabilityPanel: () => (
    <div data-testid="bioavailability-panel">BioavailabilityPanel</div>
  ),
}));

vi.mock('../../components/profiler/ConsensusLogPPanel', () => ({
  ConsensusLogPPanel: () => (
    <div data-testid="consensus-logp-panel">ConsensusLogPPanel</div>
  ),
}));

vi.mock('../../components/profiler/SkinPermeationPanel', () => ({
  SkinPermeationPanel: () => (
    <div data-testid="skin-permeation-panel">SkinPermeationPanel</div>
  ),
}));

vi.mock('../../components/profiler/DrugLikenessGrid', () => ({
  DrugLikenessGrid: () => (
    <div data-testid="drug-likeness-grid">DrugLikenessGrid</div>
  ),
}));

vi.mock('../../components/profiler/SAComparisonPanel', () => ({
  SAComparisonPanel: () => (
    <div data-testid="sa-comparison-panel">SAComparisonPanel</div>
  ),
}));

vi.mock('../../components/profiler/Shape3DPanel', () => ({
  Shape3DPanel: () => <div data-testid="shape-3d-panel">Shape3DPanel</div>,
}));

vi.mock('../../components/profiler/LigandEfficiencyPanel', () => ({
  LigandEfficiencyPanel: () => (
    <div data-testid="ligand-efficiency-panel">LigandEfficiencyPanel</div>
  ),
}));

vi.mock('../../components/profiler/MPOPanel', () => ({
  MPOPanel: () => <div data-testid="mpo-panel">MPOPanel</div>,
}));

vi.mock('../../components/profiler/ComparisonBar', () => ({
  ComparisonBar: () => (
    <div data-testid="comparison-bar">ComparisonBar</div>
  ),
}));

// Mock useProfiler hook for lazy compute functions
vi.mock('../../hooks/useProfiler', () => ({
  useProfiler: () => ({
    compute3DShape: vi.fn(),
    computeEfficiency: vi.fn(),
    computeCustomMPO: vi.fn(),
    profileCompound: vi.fn(),
    profile: null,
    isLoading: false,
    error: null,
  }),
}));

/**
 * Build a minimal valid ProfileResponse for testing.
 */
function buildMockProfile(): ProfileResponse {
  return {
    pfi: {
      pfi: 3.2,
      clogp: 2.1,
      aromatic_rings: 1,
      risk: 'low',
    },
    stars: {
      stars: 0,
      details: [],
    },
    abbott: {
      abbott_score: 0.85,
      probability_pct: 85,
      tpsa: 63.6,
      lipinski_violations: 0,
    },
    consensus_logp: {
      consensus_logp: 1.8,
      wildman_crippen: 1.7,
      xlogp3_approx: 1.9,
      xlogp3_is_approximation: false,
    },
    skin_permeation: {
      log_kp: -6.5,
      classification: 'low',
    },
    sa_comparison: {
      sa_score: {
        score: 2.5,
        scale: '1-10',
        classification: 'easy',
        available: true,
      },
      scscore: {
        score: 1.2,
        scale: '1-5',
        classification: 'easy',
        available: true,
      },
      syba: {
        score: 50,
        scale: '-256 to 256',
        classification: 'easy',
        available: true,
      },
      rascore: { available: false, note: 'Not available' },
    },
    cns_mpo: {
      score: 4.5,
      max_score: 6,
      components: {},
    },
    druglikeness: {
      lipinski: { passed: true },
      veber: { passed: true },
    },
  };
}

describe('ProfilerAccordion', () => {
  const defaultSmiles = 'CC(=O)Oc1ccccc1C(=O)O';

  describe('Data loaded state', () => {
    it('renders all profiler panels when profile data is provided', () => {
      const profile = buildMockProfile();

      render(
        <ProfilerAccordion
          smiles={defaultSmiles}
          profile={profile}
          isLoading={false}
          error={null}
        />
      );

      expect(screen.getByTestId('hero-section')).toBeInTheDocument();
      expect(screen.getByTestId('pfi-panel')).toBeInTheDocument();
      expect(screen.getByTestId('stars-panel')).toBeInTheDocument();
      expect(screen.getByTestId('bioavailability-panel')).toBeInTheDocument();
      expect(screen.getByTestId('consensus-logp-panel')).toBeInTheDocument();
      expect(screen.getByTestId('skin-permeation-panel')).toBeInTheDocument();
      expect(screen.getByTestId('drug-likeness-grid')).toBeInTheDocument();
      expect(screen.getByTestId('sa-comparison-panel')).toBeInTheDocument();
      expect(screen.getByTestId('shape-3d-panel')).toBeInTheDocument();
      expect(screen.getByTestId('ligand-efficiency-panel')).toBeInTheDocument();
      expect(screen.getByTestId('mpo-panel')).toBeInTheDocument();
    });

    it('passes smiles to HeroSection', () => {
      const profile = buildMockProfile();

      render(
        <ProfilerAccordion
          smiles={defaultSmiles}
          profile={profile}
          isLoading={false}
          error={null}
        />
      );

      expect(screen.getByTestId('hero-section')).toHaveTextContent(defaultSmiles);
    });
  });

  describe('Loading state', () => {
    it('shows loading skeleton when isLoading is true', () => {
      const { container } = render(
        <ProfilerAccordion
          smiles={defaultSmiles}
          profile={null}
          isLoading={true}
          error={null}
        />
      );

      // Should show skeleton placeholders (animate-pulse divs)
      const skeletons = container.querySelectorAll('.animate-pulse');
      expect(skeletons.length).toBeGreaterThan(0);
    });

    it('does not render profiler panels when loading', () => {
      render(
        <ProfilerAccordion
          smiles={defaultSmiles}
          profile={null}
          isLoading={true}
          error={null}
        />
      );

      expect(screen.queryByTestId('hero-section')).not.toBeInTheDocument();
      expect(screen.queryByTestId('pfi-panel')).not.toBeInTheDocument();
    });
  });

  describe('Error state', () => {
    it('shows error message when error is provided', () => {
      render(
        <ProfilerAccordion
          smiles={defaultSmiles}
          profile={null}
          isLoading={false}
          error="Unable to compute profile"
        />
      );

      expect(screen.getByText('Unable to compute profile')).toBeInTheDocument();
    });

    it('does not render profiler panels when error present', () => {
      render(
        <ProfilerAccordion
          smiles={defaultSmiles}
          profile={null}
          isLoading={false}
          error="Something went wrong"
        />
      );

      expect(screen.queryByTestId('hero-section')).not.toBeInTheDocument();
    });
  });

  describe('Empty state', () => {
    it('shows empty message when no smiles provided', () => {
      render(
        <ProfilerAccordion
          smiles=""
          profile={null}
          isLoading={false}
          error={null}
        />
      );

      expect(
        screen.getByText('Enter a molecule to see its profile')
      ).toBeInTheDocument();
    });

    it('does not render profiler panels when no smiles', () => {
      render(
        <ProfilerAccordion
          smiles=""
          profile={null}
          isLoading={false}
          error={null}
        />
      );

      expect(screen.queryByTestId('hero-section')).not.toBeInTheDocument();
    });
  });

  describe('Drug-likeness conditional rendering', () => {
    it('renders DrugLikenessGrid when druglikeness data is present', () => {
      const profile = buildMockProfile();

      render(
        <ProfilerAccordion
          smiles={defaultSmiles}
          profile={profile}
          isLoading={false}
          error={null}
        />
      );

      expect(screen.getByTestId('drug-likeness-grid')).toBeInTheDocument();
    });

    it('does not render DrugLikenessGrid when druglikeness is missing', () => {
      const profile = buildMockProfile();
      delete profile.druglikeness;

      render(
        <ProfilerAccordion
          smiles={defaultSmiles}
          profile={profile}
          isLoading={false}
          error={null}
        />
      );

      expect(screen.queryByTestId('drug-likeness-grid')).not.toBeInTheDocument();
    });
  });
});
