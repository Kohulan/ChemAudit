/**
 * Deep Validation Tab Component Tests
 *
 * Tests for:
 * - DeepValidationTab (main container)
 * - DeepCheckCard (individual check card)
 * - SeverityConfigPanel (severity override modal)
 * - StereoisomerList (collapsible SMILES list)
 * - FragmentClassificationTable (fragment mini-table)
 *
 * RDKit is NOT mocked here — these are pure React component tests.
 * useDeepValidationConfig hook is mocked where needed for control over config state.
 */
import { describe, it, expect, vi, beforeEach } from 'vitest';
import { render, screen, fireEvent } from '../setup';
import { DeepValidationTab } from '../../components/validation/DeepValidationTab';
import { DeepCheckCard } from '../../components/validation/DeepCheckCard';
import { SeverityConfigPanel } from '../../components/validation/SeverityConfigPanel';
import { StereoisomerList } from '../../components/validation/StereoisomerList';
import { FragmentClassificationTable } from '../../components/validation/FragmentClassificationTable';
import { DEEP_CHECK_DOMAINS } from '../../types/validation';
import type { CheckResult, FragmentDetail, DeepValidationConfig } from '../../types/validation';

// ============================================================
// Mock Factory Helpers
// ============================================================

/** Creates a mock CheckResult with sensible defaults */
const createMockDeepCheck = (overrides: Partial<CheckResult> = {}): CheckResult => ({
  check_name: 'stereoisomer_enumeration',
  passed: true,
  severity: 'info',
  message: 'No undefined stereocenters found.',
  affected_atoms: [],
  details: {},
  ...overrides,
});

/** Generates a full set of 16 mock deep check results (one per check name) */
function createAllDeepChecks(): CheckResult[] {
  const allChecks: CheckResult[] = [];
  for (const domain of Object.values(DEEP_CHECK_DOMAINS)) {
    for (const checkName of domain.checks) {
      allChecks.push(
        createMockDeepCheck({
          check_name: checkName,
          passed: true,
          severity: 'info',
          message: `Check ${checkName} passed.`,
        })
      );
    }
  }
  return allChecks;
}

/** Creates a mix of v1 and deep validation checks */
function createMixedChecks(): CheckResult[] {
  const v1Checks: CheckResult[] = [
    createMockDeepCheck({ check_name: 'parsability', passed: true, severity: 'pass' }),
    createMockDeepCheck({ check_name: 'valence', passed: true, severity: 'pass' }),
    createMockDeepCheck({ check_name: 'connectivity', passed: true, severity: 'pass' }),
  ];
  return [...v1Checks, ...createAllDeepChecks()];
}

/** Mock fragment data */
const mockFragments: FragmentDetail[] = [
  {
    smiles: 'CCN',
    molecular_weight: 45.08,
    heavy_atom_count: 3,
    classification: 'drug',
    pattern_name: null,
  },
  {
    smiles: 'Cl',
    molecular_weight: 36.46,
    heavy_atom_count: 1,
    classification: 'salt',
    pattern_name: 'HCl',
  },
  {
    smiles: 'O',
    molecular_weight: 18.02,
    heavy_atom_count: 1,
    classification: 'solvent',
    pattern_name: 'water',
  },
];

// ============================================================
// DeepValidationTab Tests
// ============================================================

describe('DeepValidationTab', () => {
  const onHighlightAtoms = vi.fn();

  beforeEach(() => {
    onHighlightAtoms.mockClear();
  });

  it('renders the tab with segmented control', () => {
    const checks = createAllDeepChecks();
    render(<DeepValidationTab checks={checks} onHighlightAtoms={onHighlightAtoms} />);

    expect(screen.getByRole('button', { name: 'Category' })).toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Severity' })).toBeInTheDocument();
  });

  it('shows three domain sections in category view', () => {
    const checks = createAllDeepChecks();
    render(<DeepValidationTab checks={checks} onHighlightAtoms={onHighlightAtoms} />);

    expect(screen.getByText('Stereo & Tautomers')).toBeInTheDocument();
    expect(screen.getByText('Chemical Composition')).toBeInTheDocument();
    expect(screen.getByText('Structural Complexity')).toBeInTheDocument();
  });

  it('toggles to severity view when Severity button is clicked', () => {
    const checks = createAllDeepChecks();
    render(<DeepValidationTab checks={checks} onHighlightAtoms={onHighlightAtoms} />);

    const severityBtn = screen.getByRole('button', { name: 'Severity' });
    fireEvent.click(severityBtn);

    // After clicking, the Severity button should have the active style classes
    // and Category view sections should no longer be in the document
    expect(severityBtn.className).toContain('shadow-sm');
  });

  it('shows gear icon button for severity config', () => {
    const checks = createAllDeepChecks();
    render(<DeepValidationTab checks={checks} onHighlightAtoms={onHighlightAtoms} />);

    const gearButton = screen.getByLabelText('Open severity configuration');
    expect(gearButton).toBeInTheDocument();
  });

  it('filters to only deep validation checks and excludes v1 checks', () => {
    const mixed = createMixedChecks();
    render(<DeepValidationTab checks={mixed} onHighlightAtoms={onHighlightAtoms} />);

    // Should render domain sections (deep checks present)
    expect(screen.getByText('Stereo & Tautomers')).toBeInTheDocument();

    // Should NOT render v1-only check names as domain sections
    expect(screen.queryByText('Parsability')).not.toBeInTheDocument();
    expect(screen.queryByText('Valence')).not.toBeInTheDocument();
  });

  it('shows PASS verdict when all checks pass', () => {
    const checks = createAllDeepChecks();
    render(<DeepValidationTab checks={checks} onHighlightAtoms={onHighlightAtoms} />);

    expect(screen.getByText(/Deep Validation: PASS/i)).toBeInTheDocument();
  });

  it('shows WARNING verdict when a warning-severity check fails', () => {
    const checks = createAllDeepChecks();
    // Override one check to failing with warning severity
    checks[0] = createMockDeepCheck({
      check_name: checks[0].check_name,
      passed: false,
      severity: 'warning',
      message: 'Warning detected.',
    });
    render(<DeepValidationTab checks={checks} onHighlightAtoms={onHighlightAtoms} />);

    expect(screen.getByText(/Deep Validation: WARNING/i)).toBeInTheDocument();
  });

  it('shows FAIL verdict when an error-severity check fails', () => {
    const checks = createAllDeepChecks();
    checks[0] = createMockDeepCheck({
      check_name: checks[0].check_name,
      passed: false,
      severity: 'error',
      message: 'Error detected.',
    });
    render(<DeepValidationTab checks={checks} onHighlightAtoms={onHighlightAtoms} />);

    expect(screen.getByText(/Deep Validation: FAIL/i)).toBeInTheDocument();
  });

  it('shows empty state when no deep checks are present', () => {
    const v1Only: CheckResult[] = [
      createMockDeepCheck({ check_name: 'parsability', passed: true, severity: 'pass' }),
    ];
    render(<DeepValidationTab checks={v1Only} onHighlightAtoms={onHighlightAtoms} />);

    expect(screen.getByText('No Deep Validation Results')).toBeInTheDocument();
  });

  it('opens severity config panel when gear button is clicked', () => {
    const checks = createAllDeepChecks();
    render(<DeepValidationTab checks={checks} onHighlightAtoms={onHighlightAtoms} />);

    const gearButton = screen.getByLabelText('Open severity configuration');
    fireEvent.click(gearButton);

    // Panel should appear
    expect(screen.getByText('Severity Configuration')).toBeInTheDocument();
  });

  it('shows passed and total check counts', () => {
    const checks = createAllDeepChecks();
    render(<DeepValidationTab checks={checks} onHighlightAtoms={onHighlightAtoms} />);

    // All 16 checks pass → "16/16 passed"
    expect(screen.getByText(`${checks.length}/${checks.length} passed`)).toBeInTheDocument();
  });
});

// ============================================================
// DeepCheckCard Tests
// ============================================================

describe('DeepCheckCard', () => {
  const onHighlightAtoms = vi.fn();

  beforeEach(() => {
    onHighlightAtoms.mockClear();
  });

  it('renders check name in human-readable Title Case', () => {
    const check = createMockDeepCheck({
      check_name: 'stereoisomer_enumeration',
      passed: true,
    });
    render(
      <DeepCheckCard
        check={check}
        effectiveSeverity="info"
        onHighlightAtoms={onHighlightAtoms}
      />
    );

    expect(screen.getByText('Stereoisomer Enumeration')).toBeInTheDocument();
  });

  it('renders the check message', () => {
    const check = createMockDeepCheck({
      message: 'No undefined stereocenters found.',
    });
    render(
      <DeepCheckCard
        check={check}
        effectiveSeverity="info"
        onHighlightAtoms={onHighlightAtoms}
      />
    );

    expect(screen.getByText('No undefined stereocenters found.')).toBeInTheDocument();
  });

  it('renders severity badge with effective severity text', () => {
    const check = createMockDeepCheck({
      passed: false,
      severity: 'info',
    });
    render(
      <DeepCheckCard
        check={check}
        effectiveSeverity="warning"
        onHighlightAtoms={onHighlightAtoms}
      />
    );

    expect(screen.getByText('WARNING')).toBeInTheDocument();
  });

  it('shows original severity when overridden', () => {
    const check = createMockDeepCheck({
      passed: false,
      severity: 'info',
    });
    render(
      <DeepCheckCard
        check={check}
        effectiveSeverity="error"
        onHighlightAtoms={onHighlightAtoms}
      />
    );

    // Should show override label with original
    expect(screen.getByText(/orig: info/i)).toBeInTheDocument();
  });

  it('calls onHighlightAtoms when atom index badge is clicked', () => {
    const check = createMockDeepCheck({
      check_name: 'radical_detection',
      passed: false,
      severity: 'warning',
      message: 'Radical detected.',
      affected_atoms: [0, 1, 2],
    });
    render(
      <DeepCheckCard
        check={check}
        effectiveSeverity="warning"
        onHighlightAtoms={onHighlightAtoms}
      />
    );

    // Click atom badge #0
    const atomBadge = screen.getByText('#0');
    fireEvent.click(atomBadge);

    expect(onHighlightAtoms).toHaveBeenCalledWith([0]);
  });

  it('calls onHighlightAtoms with all atoms when "highlight all" is clicked', () => {
    const check = createMockDeepCheck({
      affected_atoms: [0, 1, 2],
    });
    render(
      <DeepCheckCard
        check={check}
        effectiveSeverity="info"
        onHighlightAtoms={onHighlightAtoms}
      />
    );

    const highlightAllBtn = screen.getByRole('button', { name: /highlight all/i });
    fireEvent.click(highlightAllBtn);

    expect(onHighlightAtoms).toHaveBeenCalledWith([0, 1, 2]);
  });

  it('renders stereoisomer list for stereoisomer_enumeration check', () => {
    const check = createMockDeepCheck({
      check_name: 'stereoisomer_enumeration',
      passed: false,
      severity: 'warning',
      details: {
        undefined_count: 2,
        total_centers: 3,
        atom_indices: [1, 3],
        stereoisomer_smiles: ['C[C@@H](N)C(=O)O', 'C[C@H](N)C(=O)O'],
        enumeration_cap: 128,
        cap_exceeded: false,
      },
    });
    render(
      <DeepCheckCard
        check={check}
        effectiveSeverity="warning"
        onHighlightAtoms={onHighlightAtoms}
      />
    );

    // Show details to reveal stereoisomer list
    const showDetailsBtn = screen.getByRole('button', { name: /show details/i });
    fireEvent.click(showDetailsBtn);

    expect(screen.getByText('2 stereoisomers')).toBeInTheDocument();
  });

  it('renders fragment table for mixture_detection check', () => {
    const check = createMockDeepCheck({
      check_name: 'mixture_detection',
      passed: false,
      severity: 'warning',
      details: {
        num_fragments: 2,
        fragments: [
          {
            smiles: 'CCN',
            molecular_weight: 45.08,
            heavy_atom_count: 3,
            classification: 'drug',
            pattern_name: null,
          },
          {
            smiles: 'Cl',
            molecular_weight: 36.46,
            heavy_atom_count: 1,
            classification: 'salt',
            pattern_name: 'HCl',
          },
        ],
      },
    });
    render(
      <DeepCheckCard
        check={check}
        effectiveSeverity="warning"
        onHighlightAtoms={onHighlightAtoms}
      />
    );

    // Show details
    const showDetailsBtn = screen.getByRole('button', { name: /show details/i });
    fireEvent.click(showDetailsBtn);

    // Fragment table should render
    expect(screen.getByText('CCN')).toBeInTheDocument();
    expect(screen.getByText('Cl')).toBeInTheDocument();
  });

  it('renders PASS badge when check passes', () => {
    const check = createMockDeepCheck({ passed: true });
    render(
      <DeepCheckCard
        check={check}
        effectiveSeverity="info"
        onHighlightAtoms={onHighlightAtoms}
      />
    );

    expect(screen.getByText('PASS')).toBeInTheDocument();
  });

  it('does not show atom badges when no affected atoms', () => {
    const check = createMockDeepCheck({ affected_atoms: [] });
    render(
      <DeepCheckCard
        check={check}
        effectiveSeverity="info"
        onHighlightAtoms={onHighlightAtoms}
      />
    );

    expect(screen.queryByText('Affected atoms:')).not.toBeInTheDocument();
  });
});

// ============================================================
// SeverityConfigPanel Tests
// ============================================================

describe('SeverityConfigPanel', () => {
  const onClose = vi.fn();
  const onSetSeverity = vi.fn();
  const onRemoveOverride = vi.fn();
  const onResetAll = vi.fn();

  const defaultConfig: DeepValidationConfig = {
    severityOverrides: {},
  };

  const defaultChecks = createAllDeepChecks();

  beforeEach(() => {
    onClose.mockClear();
    onSetSeverity.mockClear();
    onRemoveOverride.mockClear();
    onResetAll.mockClear();
  });

  it('renders domain-grouped check list when open', () => {
    render(
      <SeverityConfigPanel
        isOpen={true}
        onClose={onClose}
        checks={defaultChecks}
        config={defaultConfig}
        onSetSeverity={onSetSeverity}
        onRemoveOverride={onRemoveOverride}
        onResetAll={onResetAll}
      />
    );

    // Domain headers
    expect(screen.getByText('Stereo & Tautomers')).toBeInTheDocument();
    expect(screen.getByText('Chemical Composition')).toBeInTheDocument();
    expect(screen.getByText('Structural Complexity')).toBeInTheDocument();

    // Check names should appear as Title Case
    expect(screen.getByText('Stereoisomer Enumeration')).toBeInTheDocument();
    expect(screen.getByText('Mixture Detection')).toBeInTheDocument();
  });

  it('does not render when closed', () => {
    render(
      <SeverityConfigPanel
        isOpen={false}
        onClose={onClose}
        checks={defaultChecks}
        config={defaultConfig}
        onSetSeverity={onSetSeverity}
        onRemoveOverride={onRemoveOverride}
        onResetAll={onResetAll}
      />
    );

    expect(screen.queryByText('Severity Configuration')).not.toBeInTheDocument();
  });

  it('calls onSetSeverity when a severity button is clicked', () => {
    render(
      <SeverityConfigPanel
        isOpen={true}
        onClose={onClose}
        checks={defaultChecks}
        config={defaultConfig}
        onSetSeverity={onSetSeverity}
        onRemoveOverride={onRemoveOverride}
        onResetAll={onResetAll}
      />
    );

    // Click the first ERROR button
    const errorButtons = screen.getAllByRole('button', { name: /^error$/i });
    fireEvent.click(errorButtons[0]);

    expect(onSetSeverity).toHaveBeenCalledWith(
      expect.any(String),
      'error'
    );
  });

  it('calls onResetAll when Reset All Overrides is clicked with active overrides', () => {
    const configWithOverrides: DeepValidationConfig = {
      severityOverrides: { stereoisomer_enumeration: 'error' },
    };
    render(
      <SeverityConfigPanel
        isOpen={true}
        onClose={onClose}
        checks={defaultChecks}
        config={configWithOverrides}
        onSetSeverity={onSetSeverity}
        onRemoveOverride={onRemoveOverride}
        onResetAll={onResetAll}
      />
    );

    const resetBtn = screen.getByRole('button', { name: /Reset All Overrides/i });
    fireEvent.click(resetBtn);

    expect(onResetAll).toHaveBeenCalled();
  });

  it('Reset All button is disabled when there are no overrides', () => {
    render(
      <SeverityConfigPanel
        isOpen={true}
        onClose={onClose}
        checks={defaultChecks}
        config={defaultConfig}
        onSetSeverity={onSetSeverity}
        onRemoveOverride={onRemoveOverride}
        onResetAll={onResetAll}
      />
    );

    const resetBtn = screen.getByRole('button', { name: /Reset All Overrides/i });
    expect(resetBtn).toBeDisabled();
  });

  it('calls onClose when close button is clicked', () => {
    render(
      <SeverityConfigPanel
        isOpen={true}
        onClose={onClose}
        checks={defaultChecks}
        config={defaultConfig}
        onSetSeverity={onSetSeverity}
        onRemoveOverride={onRemoveOverride}
        onResetAll={onResetAll}
      />
    );

    const closeBtn = screen.getByLabelText('Close settings');
    fireEvent.click(closeBtn);

    expect(onClose).toHaveBeenCalled();
  });

  it('shows override indicator dot for checks with active overrides', () => {
    const configWithOverrides: DeepValidationConfig = {
      severityOverrides: { stereoisomer_enumeration: 'error' },
    };
    render(
      <SeverityConfigPanel
        isOpen={true}
        onClose={onClose}
        checks={defaultChecks}
        config={configWithOverrides}
        onSetSeverity={onSetSeverity}
        onRemoveOverride={onRemoveOverride}
        onResetAll={onResetAll}
      />
    );

    // The "overridden" text should appear for the overridden check
    expect(screen.getByText('overridden')).toBeInTheDocument();
  });

  it('shows all 16 deep check names in the panel', () => {
    render(
      <SeverityConfigPanel
        isOpen={true}
        onClose={onClose}
        checks={defaultChecks}
        config={defaultConfig}
        onSetSeverity={onSetSeverity}
        onRemoveOverride={onRemoveOverride}
        onResetAll={onResetAll}
      />
    );

    const allCheckNames = Object.values(DEEP_CHECK_DOMAINS).flatMap((d) => d.checks);
    // Verify all 16 are shown (as Title Case)
    expect(allCheckNames).toHaveLength(16);
    for (const checkName of allCheckNames) {
      const titleCase = checkName
        .split('_')
        .map((w) => w.charAt(0).toUpperCase() + w.slice(1))
        .join(' ');
      expect(screen.getByText(titleCase)).toBeInTheDocument();
    }
  });
});

// ============================================================
// StereoisomerList Tests
// ============================================================

describe('StereoisomerList', () => {
  it('renders collapsed by default showing count only', () => {
    render(
      <StereoisomerList
        smiles={['C[C@@H](N)C(=O)O', 'C[C@H](N)C(=O)O']}
        cap={128}
        capExceeded={false}
      />
    );

    expect(screen.getByText('2 stereoisomers')).toBeInTheDocument();
    // Individual SMILES should NOT be visible when collapsed
    expect(screen.queryByText('C[C@@H](N)C(=O)O')).not.toBeInTheDocument();
  });

  it('expands to show SMILES list when clicked', () => {
    render(
      <StereoisomerList
        smiles={['C[C@@H](N)C(=O)O', 'C[C@H](N)C(=O)O']}
        cap={128}
        capExceeded={false}
      />
    );

    const expandBtn = screen.getByRole('button', { name: /2 stereoisomers/i });
    fireEvent.click(expandBtn);

    expect(screen.getByText('C[C@@H](N)C(=O)O')).toBeInTheDocument();
    expect(screen.getByText('C[C@H](N)C(=O)O')).toBeInTheDocument();
  });

  it('collapses again after second click', () => {
    render(
      <StereoisomerList
        smiles={['C[C@@H](N)C(=O)O']}
        cap={128}
        capExceeded={false}
      />
    );

    const btn = screen.getByRole('button', { name: /1 stereoisomer/i });
    fireEvent.click(btn);
    expect(screen.getByText('C[C@@H](N)C(=O)O')).toBeInTheDocument();

    fireEvent.click(btn);
    // After collapsing, Framer Motion AnimatePresence may keep element briefly; check aria state
    expect(btn).toHaveAttribute('aria-expanded', 'false');
  });

  it('shows cap exceeded warning message when capExceeded is true', () => {
    render(
      <StereoisomerList
        smiles={[]}
        cap={128}
        capExceeded={true}
      />
    );

    expect(screen.getByText(/128 cap exceeded/i)).toBeInTheDocument();
  });

  it('shows singular "stereoisomer" for count of 1', () => {
    render(
      <StereoisomerList
        smiles={['CCO']}
        cap={128}
        capExceeded={false}
      />
    );

    expect(screen.getByText('1 stereoisomer')).toBeInTheDocument();
  });

  it('does not show expand button when cap exceeded', () => {
    render(
      <StereoisomerList
        smiles={[]}
        cap={128}
        capExceeded={true}
      />
    );

    // Should not have an expand toggle button
    expect(screen.queryByRole('button', { name: /stereoisomer/i })).not.toBeInTheDocument();
  });
});

// ============================================================
// FragmentClassificationTable Tests
// ============================================================

describe('FragmentClassificationTable', () => {
  it('renders a table with one row per fragment', () => {
    render(<FragmentClassificationTable fragments={mockFragments} />);

    // Should show SMILES of each fragment
    expect(screen.getByText('CCN')).toBeInTheDocument();
    expect(screen.getByText('Cl')).toBeInTheDocument();
    expect(screen.getByText('O')).toBeInTheDocument();
  });

  it('shows table headers', () => {
    render(<FragmentClassificationTable fragments={mockFragments} />);

    expect(screen.getByText('Fragment')).toBeInTheDocument();
    expect(screen.getByText('MW')).toBeInTheDocument();
    expect(screen.getByText('Classification')).toBeInTheDocument();
    expect(screen.getByText('Pattern')).toBeInTheDocument();
  });

  it('renders classification badges for each fragment type', () => {
    render(<FragmentClassificationTable fragments={mockFragments} />);

    // Each classification should appear as a badge text
    expect(screen.getByText('drug')).toBeInTheDocument();
    expect(screen.getByText('salt')).toBeInTheDocument();
    expect(screen.getByText('solvent')).toBeInTheDocument();
  });

  it('shows pattern names for fragments with patterns', () => {
    render(<FragmentClassificationTable fragments={mockFragments} />);

    expect(screen.getByText('HCl')).toBeInTheDocument();
    expect(screen.getByText('water')).toBeInTheDocument();
    // Null pattern should show dash
    expect(screen.getByText('—')).toBeInTheDocument();
  });

  it('renders empty state message when no fragments provided', () => {
    render(<FragmentClassificationTable fragments={[]} />);

    expect(screen.getByText('No fragment details available.')).toBeInTheDocument();
  });

  it('renders correct number of rows for fragments', () => {
    const twoFragments = mockFragments.slice(0, 2);
    render(<FragmentClassificationTable fragments={twoFragments} />);

    expect(screen.getByText('CCN')).toBeInTheDocument();
    expect(screen.getByText('Cl')).toBeInTheDocument();
    expect(screen.queryByText('O')).not.toBeInTheDocument();
  });
});
