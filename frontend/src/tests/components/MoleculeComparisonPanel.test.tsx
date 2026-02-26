/**
 * MoleculeComparisonPanel Tests
 *
 * Tests the comparison drawer including the ECFP4 Tanimoto similarity gauge.
 */
import { describe, it, expect, vi, beforeEach } from 'vitest';
import { render, screen, waitFor } from '../setup';
import { MoleculeComparisonPanel } from '../../components/batch/MoleculeComparisonPanel';
import type { BatchResult } from '../../types/batch';

// Mock framer-motion to avoid animation issues in tests
vi.mock('framer-motion', () => ({
  motion: {
    div: ({ children, className, onClick, title }: Record<string, unknown>) => {
      const props: Record<string, unknown> = {};
      if (className) props.className = className;
      if (onClick) props.onClick = onClick;
      if (title) props.title = title;
      return <div {...props}>{children as React.ReactNode}</div>;
    },
  },
  AnimatePresence: ({ children }: { children: React.ReactNode }) => <>{children}</>,
}));

// Mock the API module
const mockGetSimilarity = vi.fn();
vi.mock('../../services/api', () => ({
  validationApi: {
    getSimilarity: (...args: unknown[]) => mockGetSimilarity(...args),
  },
}));

// Mock MoleculeViewer (requires RDKit WASM)
vi.mock('../../components/molecules/MoleculeViewer', () => ({
  MoleculeViewer: ({ smiles }: { smiles: string }) => (
    <div data-testid="molecule-viewer">{smiles}</div>
  ),
}));

// Mock MoleculePropertyRadar (requires recharts)
vi.mock('../../components/batch/MoleculePropertyRadar', () => ({
  MoleculePropertyRadar: () => <div data-testid="property-radar" />,
}));

const createMockBatchResult = (overrides: Partial<BatchResult> = {}): BatchResult => ({
  index: 0,
  smiles: 'CCO',
  original_smiles: 'CCO',
  status: 'success',
  name: 'Ethanol',
  validation: { overall_score: 95, issues: [], all_checks: [] },
  scoring: null,
  alerts: null,
  standardization: null,
  ...overrides,
} as BatchResult);

describe('MoleculeComparisonPanel', () => {
  const mockOnClose = vi.fn();
  const mockOnRemove = vi.fn();

  beforeEach(() => {
    vi.clearAllMocks();
    mockGetSimilarity.mockResolvedValue({
      tanimoto_similarity: 0.78,
      fingerprint_type: 'ECFP4',
      radius: 2,
      n_bits: 2048,
      common_bits: 45,
      bits_a: 52,
      bits_b: 58,
    });
  });

  describe('Basic Rendering', () => {
    it('renders nothing when no molecules provided', () => {
      const { container } = render(
        <MoleculeComparisonPanel
          molecules={[]}
          datasetStats={null}
          onClose={mockOnClose}
          onRemoveMolecule={mockOnRemove}
        />
      );
      expect(container.querySelector('h2')).toBeNull();
    });

    it('renders header with single molecule', () => {
      const mol = createMockBatchResult();
      render(
        <MoleculeComparisonPanel
          molecules={[mol]}
          datasetStats={null}
          onClose={mockOnClose}
          onRemoveMolecule={mockOnRemove}
        />
      );
      expect(screen.getByText('Compare Molecules')).toBeInTheDocument();
    });

    it('renders both molecule cards when two molecules provided', () => {
      const molA = createMockBatchResult({ index: 0, smiles: 'CCO', name: 'Ethanol' });
      const molB = createMockBatchResult({ index: 1, smiles: 'CCCO', name: 'Propanol' });
      render(
        <MoleculeComparisonPanel
          molecules={[molA, molB]}
          datasetStats={null}
          onClose={mockOnClose}
          onRemoveMolecule={mockOnRemove}
        />
      );
      expect(screen.getByText('Molecule A')).toBeInTheDocument();
      expect(screen.getByText('Molecule B')).toBeInTheDocument();
    });
  });

  describe('Similarity Gauge', () => {
    it('does not show similarity gauge with single molecule', () => {
      const mol = createMockBatchResult();
      render(
        <MoleculeComparisonPanel
          molecules={[mol]}
          datasetStats={null}
          onClose={mockOnClose}
          onRemoveMolecule={mockOnRemove}
        />
      );
      expect(screen.queryByText('ECFP4 Tanimoto Similarity')).not.toBeInTheDocument();
    });

    it('fetches and displays similarity when two molecules loaded', async () => {
      const molA = createMockBatchResult({ index: 0, smiles: 'CCO' });
      const molB = createMockBatchResult({ index: 1, smiles: 'CCCO' });
      render(
        <MoleculeComparisonPanel
          molecules={[molA, molB]}
          datasetStats={null}
          onClose={mockOnClose}
          onRemoveMolecule={mockOnRemove}
        />
      );

      await waitFor(() => {
        expect(screen.getByText('ECFP4 Tanimoto Similarity')).toBeInTheDocument();
      });

      expect(mockGetSimilarity).toHaveBeenCalledWith('CCO', 'CCCO');
      expect(screen.getByText('78%')).toBeInTheDocument();
      expect(screen.getByText('Similar')).toBeInTheDocument();
    });

    it('shows loading spinner while fetching', () => {
      // Make the promise never resolve
      mockGetSimilarity.mockReturnValue(new Promise(() => {}));

      const molA = createMockBatchResult({ index: 0, smiles: 'CCO' });
      const molB = createMockBatchResult({ index: 1, smiles: 'CCCO' });
      render(
        <MoleculeComparisonPanel
          molecules={[molA, molB]}
          datasetStats={null}
          onClose={mockOnClose}
          onRemoveMolecule={mockOnRemove}
        />
      );

      expect(screen.getByText('Computing similarity...')).toBeInTheDocument();
    });

    it('shows "Very Similar" for score >= 0.85', async () => {
      mockGetSimilarity.mockResolvedValue({
        tanimoto_similarity: 0.92,
        fingerprint_type: 'ECFP4',
        radius: 2,
        n_bits: 2048,
        common_bits: 50,
        bits_a: 52,
        bits_b: 55,
      });

      const molA = createMockBatchResult({ index: 0, smiles: 'CCO' });
      const molB = createMockBatchResult({ index: 1, smiles: 'CCO' });
      render(
        <MoleculeComparisonPanel
          molecules={[molA, molB]}
          datasetStats={null}
          onClose={mockOnClose}
          onRemoveMolecule={mockOnRemove}
        />
      );

      await waitFor(() => {
        expect(screen.getByText('Very Similar')).toBeInTheDocument();
        expect(screen.getByText('92%')).toBeInTheDocument();
      });
    });

    it('shows "Moderate" for score >= 0.50 and < 0.70', async () => {
      mockGetSimilarity.mockResolvedValue({
        tanimoto_similarity: 0.55,
        fingerprint_type: 'ECFP4',
        radius: 2,
        n_bits: 2048,
        common_bits: 30,
        bits_a: 52,
        bits_b: 58,
      });

      const molA = createMockBatchResult({ index: 0, smiles: 'CCO' });
      const molB = createMockBatchResult({ index: 1, smiles: 'c1ccccc1' });
      render(
        <MoleculeComparisonPanel
          molecules={[molA, molB]}
          datasetStats={null}
          onClose={mockOnClose}
          onRemoveMolecule={mockOnRemove}
        />
      );

      await waitFor(() => {
        expect(screen.getByText('Moderate')).toBeInTheDocument();
        expect(screen.getByText('55%')).toBeInTheDocument();
      });
    });

    it('shows "Dissimilar" for score < 0.50', async () => {
      mockGetSimilarity.mockResolvedValue({
        tanimoto_similarity: 0.23,
        fingerprint_type: 'ECFP4',
        radius: 2,
        n_bits: 2048,
        common_bits: 10,
        bits_a: 52,
        bits_b: 58,
      });

      const molA = createMockBatchResult({ index: 0, smiles: 'CCO' });
      const molB = createMockBatchResult({ index: 1, smiles: 'c1ccccc1' });
      render(
        <MoleculeComparisonPanel
          molecules={[molA, molB]}
          datasetStats={null}
          onClose={mockOnClose}
          onRemoveMolecule={mockOnRemove}
        />
      );

      await waitFor(() => {
        expect(screen.getByText('Dissimilar')).toBeInTheDocument();
        expect(screen.getByText('23%')).toBeInTheDocument();
      });
    });

    it('displays fingerprint bit details', async () => {
      const molA = createMockBatchResult({ index: 0, smiles: 'CCO' });
      const molB = createMockBatchResult({ index: 1, smiles: 'CCCO' });
      render(
        <MoleculeComparisonPanel
          molecules={[molA, molB]}
          datasetStats={null}
          onClose={mockOnClose}
          onRemoveMolecule={mockOnRemove}
        />
      );

      await waitFor(() => {
        expect(screen.getByText('Common: 45 bits')).toBeInTheDocument();
        expect(screen.getByText('A: 52 | B: 58')).toBeInTheDocument();
      });
    });

    it('handles API error gracefully', async () => {
      mockGetSimilarity.mockRejectedValue(new Error('Network error'));

      const molA = createMockBatchResult({ index: 0, smiles: 'CCO' });
      const molB = createMockBatchResult({ index: 1, smiles: 'CCCO' });
      render(
        <MoleculeComparisonPanel
          molecules={[molA, molB]}
          datasetStats={null}
          onClose={mockOnClose}
          onRemoveMolecule={mockOnRemove}
        />
      );

      // Should not crash — gauge simply won't appear
      await waitFor(() => {
        expect(screen.queryByText('Computing similarity...')).not.toBeInTheDocument();
      });
      expect(screen.queryByText('ECFP4 Tanimoto Similarity')).not.toBeInTheDocument();
    });
  });

  describe('Properties Table', () => {
    it('renders property comparison rows', () => {
      const molA = createMockBatchResult({ index: 0, smiles: 'CCO' });
      const molB = createMockBatchResult({ index: 1, smiles: 'CCCO' });
      render(
        <MoleculeComparisonPanel
          molecules={[molA, molB]}
          datasetStats={null}
          onClose={mockOnClose}
          onRemoveMolecule={mockOnRemove}
        />
      );

      expect(screen.getByText('Overall Score')).toBeInTheDocument();
      expect(screen.getByText('QED')).toBeInTheDocument();
      expect(screen.getByText('SA Score')).toBeInTheDocument();
    });
  });
});
