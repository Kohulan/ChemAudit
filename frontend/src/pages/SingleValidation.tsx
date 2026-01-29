import { useState, useEffect } from 'react';
import { useSearchParams } from 'react-router-dom';
import { motion, AnimatePresence } from 'framer-motion';
import {
  Atom,
  Beaker,
  AlertTriangle,
  Database,
  Sparkles,
  RotateCcw,
  Play,
  Search,
  Layers,
  Target,
  Zap,
  CheckCircle2,
  Info,
  Share2,
} from 'lucide-react';
import { StructureInput } from '../components/molecules/StructureInput';
import { MoleculeViewer } from '../components/molecules/MoleculeViewer';
import { IssueCard } from '../components/validation/IssueCard';
import { AlertCard } from '../components/alerts/AlertCard';
import { ScoringResults } from '../components/scoring/ScoringResults';
import { ScoreChart } from '../components/scoring/ScoreChart';
import { StandardizationResults } from '../components/standardization/StandardizationResults';
import { DatabaseLookupResults } from '../components/integrations/DatabaseLookupResults';
import { ClayButton } from '../components/ui/ClayButton';
import { Badge } from '../components/ui/Badge';
import { MoleculeLoader } from '../components/ui/MoleculeLoader';
import { CopyButton } from '../components/ui/CopyButton';
import { useValidation } from '../hooks/useValidation';
import { useMoleculeInfo } from '../hooks/useMoleculeInfo';
import { alertsApi, scoringApi, standardizationApi, integrationsApi } from '../services/api';
import { cn, getScoreLabel } from '../lib/utils';
import type { AlertScreenResponse, AlertError } from '../types/alerts';
import type { ScoringResponse, ScoringError } from '../types/scoring';
import type { StandardizeResponse, StandardizeError } from '../types/standardization';
import type { PubChemResult, ChEMBLResult, COCONUTResult } from '../types/integrations';

const EXAMPLE_MOLECULES = [
  { name: 'Aspirin', smiles: 'CC(=O)Oc1ccccc1C(=O)O' },
  { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' },
  { name: 'L-Alanine (chiral)', smiles: 'C[C@H](N)C(=O)O' },
  { name: 'E-Stilbene (E/Z)', smiles: 'C(/c1ccccc1)=C/c1ccccc1' },
  { name: 'Morphine', smiles: 'CN1CCC23C4=C5C=CC(O)=C4OC2C(O)C=CC3C1C5' },
  { name: 'Rhodanine (PAINS)', smiles: 'O=C1NC(=S)SC1' },
  { name: 'Amine HCl (salt)', smiles: 'CCN.Cl' },
];

type TabType = 'validate' | 'database' | 'alerts' | 'standardize';

interface TabConfig {
  id: TabType;
  label: string;
  icon: React.ReactNode;
  description: string;
}

const TABS: TabConfig[] = [
  {
    id: 'validate',
    label: 'Validate & Score',
    icon: <CheckCircle2 className="w-4 h-4" />,
    description: 'Check structure validity, calculate quality metrics, and assess ML-readiness',
  },
  {
    id: 'database',
    label: 'Database Lookup',
    icon: <Database className="w-4 h-4" />,
    description: 'Search PubChem, ChEMBL, and COCONUT for compound information',
  },
  {
    id: 'alerts',
    label: 'Alerts Screening',
    icon: <AlertTriangle className="w-4 h-4" />,
    description: 'Screen for problematic structural patterns using PAINS, BRENK, NIH, and ZINC filters',
  },
  {
    id: 'standardize',
    label: 'Standardize',
    icon: <Layers className="w-4 h-4" />,
    description: 'Normalize structure representation and remove salts/solvents',
  },
];

export function SingleValidationPage() {
  const [searchParams, setSearchParams] = useSearchParams();
  const [molecule, setMolecule] = useState('');
  const [highlightedAtoms, setHighlightedAtoms] = useState<number[]>([]);
  const [activeTab, setActiveTab] = useState<TabType>('validate');
  const { validate, result, error, isLoading, reset } = useValidation();
  // @ts-expect-error - Used later in event handlers
  const [shareToastVisible, setShareToastVisible] = useState(false);

  // Load molecule from URL on mount
  useEffect(() => {
    const smilesFromUrl = searchParams.get('smiles');
    if (smilesFromUrl) {
      setMolecule(decodeURIComponent(smilesFromUrl));
      // Auto-validate after loading from URL
      setTimeout(() => {
        validate({
          molecule: decodeURIComponent(smilesFromUrl),
          format: 'auto',
          preserve_aromatic: false,
        });
      }, 100);
    }
  }, []); // Only run on mount

  // Get molecule info immediately when molecule is entered
  const { info: moleculeInfo, isLoading: moleculeInfoLoading, error: moleculeInfoError } = useMoleculeInfo(molecule);

  // Check if input has aromatic atoms (lowercase c, n, o, s in SMILES context)
  const hasAromaticInput = /[cnos]/.test(molecule) && molecule.trim().length > 0;

  // Aromatic SMILES preference - show option when aromatic atoms detected
  const [preferAromaticSmiles, setPreferAromaticSmiles] = useState(false);

  // Alert screening state
  const [alertResult, setAlertResult] = useState<AlertScreenResponse | null>(null);
  const [alertError, setAlertError] = useState<AlertError | null>(null);
  const [alertsLoading, setAlertsLoading] = useState(false);
  const [selectedCatalogs, setSelectedCatalogs] = useState<string[]>(['PAINS', 'BRENK']);

  // Scoring state
  const [scoringResult, setScoringResult] = useState<ScoringResponse | null>(null);
  const [scoringError, setScoringError] = useState<ScoringError | null>(null);
  const [scoringLoading, setScoringLoading] = useState(false);

  // Standardization state
  const [standardizationResult, setStandardizationResult] = useState<StandardizeResponse | null>(null);
  const [standardizationError, setStandardizationError] = useState<StandardizeError | null>(null);
  const [standardizationLoading, setStandardizationLoading] = useState(false);
  const [includeTautomer, setIncludeTautomer] = useState(false);

  // SMILES display preference
  const [showKekulized, setShowKekulized] = useState(false);

  // CIP stereochemistry labels toggle
  const [showCIP, setShowCIP] = useState(false);

  // Database lookup state
  const [databaseResults, setDatabaseResults] = useState<{
    pubchem: PubChemResult | null;
    chembl: ChEMBLResult | null;
    coconut: COCONUTResult | null;
  } | null>(null);
  const [databaseLoading, setDatabaseLoading] = useState(false);

  const handleValidate = async () => {
    if (!molecule.trim()) return;
    await validate({
      molecule: molecule.trim(),
      format: 'auto',
      preserve_aromatic: preferAromaticSmiles,
    });
  };

  const handleScreenAlerts = async () => {
    if (!molecule.trim()) return;
    setAlertsLoading(true);
    setAlertError(null);
    try {
      const response = await alertsApi.screenAlerts({
        molecule: molecule.trim(),
        format: 'auto',
        catalogs: selectedCatalogs,
      });
      setAlertResult(response);
    } catch (err) {
      setAlertError(err as AlertError);
      setAlertResult(null);
    } finally {
      setAlertsLoading(false);
    }
  };

  const handleCalculateScores = async () => {
    if (!molecule.trim()) return;
    setScoringLoading(true);
    setScoringError(null);
    try {
      const response = await scoringApi.getScoring(molecule.trim(), 'auto');
      setScoringResult(response);
    } catch (err) {
      setScoringError(err as ScoringError);
      setScoringResult(null);
    } finally {
      setScoringLoading(false);
    }
  };

  const handleStandardize = async () => {
    if (!molecule.trim()) return;
    setStandardizationLoading(true);
    setStandardizationError(null);
    try {
      const response = await standardizationApi.standardize({
        molecule: molecule.trim(),
        format: 'auto',
        options: { include_tautomer: includeTautomer, preserve_stereo: true },
      });
      setStandardizationResult(response);
    } catch (err) {
      setStandardizationError(err as StandardizeError);
      setStandardizationResult(null);
    } finally {
      setStandardizationLoading(false);
    }
  };

  const handleDatabaseLookup = async () => {
    if (!molecule.trim()) return;
    setDatabaseLoading(true);
    try {
      const results = await integrationsApi.lookupAll({ smiles: molecule.trim() });
      setDatabaseResults(results);
    } catch (err) {
      console.error('Database lookup error:', err);
      setDatabaseResults(null);
    } finally {
      setDatabaseLoading(false);
    }
  };

  const handleExampleClick = (smiles: string) => {
    setMolecule(smiles);
    resetAll();
  };

  const resetAll = () => {
    reset();
    setAlertResult(null);
    setAlertError(null);
    setScoringResult(null);
    setScoringError(null);
    setStandardizationResult(null);
    setStandardizationError(null);
    setDatabaseResults(null);
    setHighlightedAtoms([]);
    setShowCIP(false);
  };

  const handleReset = () => {
    setMolecule('');
    resetAll();
    // Clear URL params
    setSearchParams({});
  };

  const handleShare = async () => {
    if (!molecule.trim()) return;

    // Update URL with encoded SMILES
    setSearchParams({ smiles: encodeURIComponent(molecule.trim()) });

    // Copy URL to clipboard
    try {
      await navigator.clipboard.writeText(window.location.href);
      setShareToastVisible(true);
      setTimeout(() => setShareToastVisible(false), 3000);
    } catch (err) {
      console.error('Failed to copy URL:', err);
    }
  };

  const toggleCatalog = (catalog: string) => {
    setSelectedCatalogs((prev) =>
      prev.includes(catalog) ? prev.filter((c) => c !== catalog) : [...prev, catalog]
    );
  };

  const isAnyLoading = isLoading || alertsLoading || scoringLoading || standardizationLoading || databaseLoading;
  const hasError = error || alertError || scoringError || standardizationError;

  // Calculate quality score
  const qualityScore = result?.overall_score ?? scoringResult?.ml_readiness?.score ?? null;
  // ML ready if score >= 70
  const mlReadyScore = scoringResult?.ml_readiness?.score;
  const mlReady = mlReadyScore !== undefined ? mlReadyScore >= 70 : null;

  // Derived state for right column
  const showStandardizationComparison = activeTab === 'standardize' && standardizationResult;
  const hasScores = result?.overall_score !== undefined || scoringResult?.ml_readiness;

  // Get canonical SMILES from any available result
  const canonicalSmiles = result?.molecule_info.canonical_smiles
    || alertResult?.molecule_info.canonical_smiles
    || scoringResult?.molecule_info.canonical_smiles;

  // Current issues to display based on active tab/operation
  const validationIssues = result?.issues || [];
  const alertIssues = alertResult?.alerts || [];

  return (
    <div className="max-w-7xl mx-auto space-y-6 px-4 sm:px-6">
      {/* Header */}
      <motion.div
        className="text-center pt-4"
        initial={{ opacity: 0, y: -20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, ease: [0.25, 0.46, 0.45, 0.94] }}
      >
        <h1 className="text-3xl sm:text-4xl font-bold text-gradient tracking-tight font-display">
          Molecule Validation
        </h1>
        <p className="text-[var(--color-text-secondary)] mt-3 text-base sm:text-lg max-w-2xl mx-auto leading-relaxed">
          Comprehensive validation, scoring, and standardization for chemical structures
        </p>
      </motion.div>

      {/* Example molecules */}
      <motion.div
        className="flex flex-wrap gap-2 justify-center"
        initial={{ opacity: 0 }}
        animate={{ opacity: 1 }}
        transition={{ delay: 0.2, duration: 0.5 }}
      >
        <span className="text-sm text-[var(--color-text-muted)] self-center mr-2">Try:</span>
        {EXAMPLE_MOLECULES.map((example, i) => (
          <motion.button
            key={example.name}
            onClick={() => handleExampleClick(example.smiles)}
            className={cn(
              'px-3 py-1.5 text-sm rounded-full transition-all',
              'bg-[var(--color-surface-elevated)] border border-[var(--color-border)]',
              'text-[var(--color-text-secondary)]',
              'hover:border-[var(--color-primary)] hover:text-[var(--color-primary)]',
              'hover:shadow-[0_0_12px_var(--glow-primary)]'
            )}
            initial={{ opacity: 0, scale: 0.9 }}
            animate={{ opacity: 1, scale: 1 }}
            transition={{ delay: 0.3 + i * 0.05 }}
            whileHover={{ scale: 1.05 }}
            whileTap={{ scale: 0.95 }}
          >
            {example.name}
          </motion.button>
        ))}
      </motion.div>

      {/* Two-Column Layout */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* LEFT COLUMN */}
        <motion.div
          className="space-y-4"
          initial={{ opacity: 0, x: -20 }}
          animate={{ opacity: 1, x: 0 }}
          transition={{ duration: 0.5, delay: 0.2 }}
        >
          {/* Input Section */}
          <div className="card-glass p-5 sm:p-6">
            <div className="flex items-center gap-3 mb-4">
              <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-[var(--color-primary)]/10 to-[var(--color-accent)]/10 flex items-center justify-center text-[var(--color-primary)]">
                <Beaker className="w-5 h-5" />
              </div>
              <div>
                <h4 className="font-semibold text-[var(--color-text-primary)] text-sm tracking-tight">
                  Input Structure
                </h4>
                <p className="text-xs text-[var(--color-text-muted)] mt-0.5">SMILES, InChI, or draw</p>
              </div>
            </div>
            <StructureInput value={molecule} onChange={setMolecule} />

            <div className="mt-4 flex items-center gap-3">
              <ClayButton
                onClick={handleReset}
                disabled={!molecule && !result && !alertResult && !scoringResult && !standardizationResult && !databaseResults}
                leftIcon={<RotateCcw className="w-4 h-4" />}
              >
                Reset
              </ClayButton>
              <ClayButton
                onClick={handleShare}
                disabled={!molecule.trim()}
                leftIcon={<Share2 className="w-4 h-4" />}
              >
                Share
              </ClayButton>
            </div>

            {/* Parsing Failed Message - Suggest trying validation */}
            <AnimatePresence>
              {moleculeInfoError && molecule.trim() && !result && (
                <motion.div
                  initial={{ opacity: 0, height: 0 }}
                  animate={{ opacity: 1, height: 'auto' }}
                  exit={{ opacity: 0, height: 0 }}
                  transition={{ duration: 0.2 }}
                  className="mt-4 pt-4 border-t border-[var(--color-border)]"
                >
                  <div className="flex items-start gap-3 p-3 rounded-xl bg-amber-500/10 border border-amber-500/20">
                    <AlertTriangle className="w-5 h-5 text-amber-500 flex-shrink-0 mt-0.5" />
                    <div className="flex-1">
                      <p className="text-sm font-medium text-amber-600 dark:text-amber-400">
                        Initial parsing failed
                      </p>
                      <p className="text-xs text-[var(--color-text-secondary)] mt-1">
                        The browser-based parser couldn't read this structure, but the server might be able to sanitize and fix it.
                        Try clicking <strong>Validate</strong> to process it with the full RDKit engine.
                      </p>
                    </div>
                  </div>
                </motion.div>
              )}
            </AnimatePresence>

            {/* Aromatic SMILES Preference - Show when aromatic atoms detected (even if parsing failed) */}
            <AnimatePresence>
              {hasAromaticInput && molecule.trim() && (
                <motion.div
                  initial={{ opacity: 0, height: 0 }}
                  animate={{ opacity: 1, height: 'auto' }}
                  exit={{ opacity: 0, height: 0 }}
                  transition={{ duration: 0.2 }}
                  className="mt-3"
                >
                  <label className="flex items-center gap-2 text-xs text-[var(--color-text-secondary)] cursor-pointer p-2 rounded-lg hover:bg-[var(--color-surface-sunken)] transition-colors">
                    <input
                      type="checkbox"
                      checked={preferAromaticSmiles}
                      onChange={(e) => setPreferAromaticSmiles(e.target.checked)}
                      className="rounded border-[var(--color-border-strong)] text-[var(--color-primary)] focus:ring-[var(--color-primary)]/30"
                    />
                    <span>Preserve aromatic notation in output (e.g., <code className="text-[10px] bg-[var(--color-surface-sunken)] px-1 rounded">c1ccccc1</code> instead of <code className="text-[10px] bg-[var(--color-surface-sunken)] px-1 rounded">C1=CC=CC=C1</code>)</span>
                  </label>
                </motion.div>
              )}
            </AnimatePresence>

            {/* Molecule Info - Shows when valid molecule is entered or after validation */}
            <AnimatePresence>
              {(moleculeInfo || result) && !moleculeInfoLoading && (
                <motion.div
                  initial={{ opacity: 0, height: 0 }}
                  animate={{ opacity: 1, height: 'auto' }}
                  exit={{ opacity: 0, height: 0 }}
                  transition={{ duration: 0.2 }}
                  className="mt-4 pt-4 border-t border-[var(--color-border)]"
                >
                  <div className="flex items-center gap-2 mb-3">
                    <Atom className="w-4 h-4 text-[var(--color-primary)]" />
                    <span className="text-xs font-medium text-[var(--color-text-secondary)]">Molecule Info</span>
                    <Badge variant="success" size="sm">{result ? 'Validated' : 'Valid'}</Badge>
                  </div>

                  {/* Show detailed info if validation result available, otherwise basic info */}
                  {result ? (
                    <div className="space-y-3">
                      {/* Stats row */}
                      <div className="grid grid-cols-4 gap-2 text-center">
                        {result.molecule_info.num_atoms && (
                          <div className="bg-[var(--color-surface-sunken)] rounded-lg p-2">
                            <div className="text-lg font-bold text-[var(--color-text-primary)]">{result.molecule_info.num_atoms}</div>
                            <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">Atoms</div>
                          </div>
                        )}
                        {moleculeInfo?.numBonds && (
                          <div className="bg-[var(--color-surface-sunken)] rounded-lg p-2">
                            <div className="text-lg font-bold text-[var(--color-text-primary)]">{moleculeInfo.numBonds}</div>
                            <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">Bonds</div>
                          </div>
                        )}
                        {moleculeInfo?.numRings !== undefined && (
                          <div className="bg-[var(--color-surface-sunken)] rounded-lg p-2">
                            <div className="text-lg font-bold text-[var(--color-text-primary)]">{moleculeInfo.numRings}</div>
                            <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">Rings</div>
                          </div>
                        )}
                        {result.molecule_info.molecular_weight && (
                          <div className="bg-[var(--color-surface-sunken)] rounded-lg p-2">
                            <div className="text-lg font-bold text-[var(--color-text-primary)]">{result.molecule_info.molecular_weight.toFixed(1)}</div>
                            <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">MW</div>
                          </div>
                        )}
                      </div>

                      {/* Detailed info */}
                      <div className="space-y-2 text-sm">
                        {result.molecule_info.molecular_formula && (
                          <div className="flex items-center gap-2">
                            <span className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider w-16">Formula</span>
                            <span className="text-[var(--color-text-primary)] font-medium">{result.molecule_info.molecular_formula}</span>
                          </div>
                        )}
                        {result.molecule_info.canonical_smiles && (
                          <div>
                            <div className="flex items-center justify-between mb-2">
                              <span className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">
                                {showKekulized && moleculeInfo?.kekulizedSmiles ? 'Kekulized SMILES' : 'Canonical SMILES'}
                              </span>
                              <div className="flex items-center gap-2">
                                {moleculeInfo?.kekulizedSmiles && (
                                  <button
                                    onClick={() => setShowKekulized(!showKekulized)}
                                    className={cn(
                                      'flex items-center gap-1.5 text-xs font-medium px-3 py-1.5 rounded-lg transition-all',
                                      'border shadow-sm',
                                      showKekulized
                                        ? 'bg-[var(--color-primary)] text-white border-[var(--color-primary)] shadow-[var(--color-primary)]/20'
                                        : 'bg-[var(--color-surface-elevated)] text-[var(--color-text-secondary)] border-[var(--color-border)] hover:border-[var(--color-primary)] hover:text-[var(--color-primary)]'
                                    )}
                                  >
                                    <svg className="w-3.5 h-3.5" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                                      <path d="M7 16V4m0 0L3 8m4-4l4 4m6 0v12m0 0l4-4m-4 4l-4-4" />
                                    </svg>
                                    {showKekulized ? 'Canonical' : 'Kekulize'}
                                  </button>
                                )}
                                <CopyButton text={showKekulized && moleculeInfo?.kekulizedSmiles ? moleculeInfo.kekulizedSmiles : (result.molecule_info.canonical_smiles || '')} />
                              </div>
                            </div>
                            <code className="text-xs text-[var(--color-text-secondary)] font-mono break-all block bg-[var(--color-surface-sunken)] rounded-lg px-3 py-2">
                              {showKekulized && moleculeInfo?.kekulizedSmiles
                                ? moleculeInfo.kekulizedSmiles
                                : result.molecule_info.canonical_smiles}
                            </code>
                          </div>
                        )}
                        {result.molecule_info.inchi && (
                          <div>
                            <div className="flex items-center justify-between mb-1">
                              <span className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">InChI</span>
                              <CopyButton text={result.molecule_info.inchi} />
                            </div>
                            <code className="text-xs text-[var(--color-text-secondary)] font-mono break-all block bg-[var(--color-surface-sunken)] rounded px-2 py-1.5 max-h-16 overflow-y-auto">
                              {result.molecule_info.inchi}
                            </code>
                          </div>
                        )}
                        {result.molecule_info.inchikey && (
                          <div>
                            <div className="flex items-center justify-between mb-1">
                              <span className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">InChIKey</span>
                              <CopyButton text={result.molecule_info.inchikey} />
                            </div>
                            <code className="text-xs text-[var(--color-text-secondary)] font-mono break-all block bg-[var(--color-surface-sunken)] rounded px-2 py-1.5">
                              {result.molecule_info.inchikey}
                            </code>
                          </div>
                        )}
                      </div>
                    </div>
                  ) : moleculeInfo ? (
                    <>
                      <div className="grid grid-cols-3 gap-3 text-center">
                        <div className="bg-[var(--color-surface-sunken)] rounded-lg p-2">
                          <div className="text-lg font-bold text-[var(--color-text-primary)]">{moleculeInfo.numAtoms}</div>
                          <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">Atoms</div>
                        </div>
                        <div className="bg-[var(--color-surface-sunken)] rounded-lg p-2">
                          <div className="text-lg font-bold text-[var(--color-text-primary)]">{moleculeInfo.numBonds}</div>
                          <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">Bonds</div>
                        </div>
                        <div className="bg-[var(--color-surface-sunken)] rounded-lg p-2">
                          <div className="text-lg font-bold text-[var(--color-text-primary)]">{moleculeInfo.numRings}</div>
                          <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">Rings</div>
                        </div>
                      </div>
                      <div className="mt-3">
                        <div className="flex items-center justify-between mb-2">
                          <span className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">
                            {showKekulized && moleculeInfo.kekulizedSmiles ? 'Kekulized SMILES' : 'Canonical SMILES'}
                          </span>
                          <div className="flex items-center gap-2">
                            {moleculeInfo.kekulizedSmiles && (
                              <button
                                onClick={() => setShowKekulized(!showKekulized)}
                                className={cn(
                                  'flex items-center gap-1.5 text-xs font-medium px-3 py-1.5 rounded-lg transition-all',
                                  'border shadow-sm',
                                  showKekulized
                                    ? 'bg-[var(--color-primary)] text-white border-[var(--color-primary)] shadow-[var(--color-primary)]/20'
                                    : 'bg-[var(--color-surface-elevated)] text-[var(--color-text-secondary)] border-[var(--color-border)] hover:border-[var(--color-primary)] hover:text-[var(--color-primary)]'
                                )}
                              >
                                <svg className="w-3.5 h-3.5" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                                  <path d="M7 16V4m0 0L3 8m4-4l4 4m6 0v12m0 0l4-4m-4 4l-4-4" />
                                </svg>
                                {showKekulized ? 'Canonical' : 'Kekulize'}
                              </button>
                            )}
                            <CopyButton text={showKekulized && moleculeInfo.kekulizedSmiles ? moleculeInfo.kekulizedSmiles : moleculeInfo.canonicalSmiles} />
                          </div>
                        </div>
                        <code className="text-xs text-[var(--color-text-secondary)] font-mono break-all block bg-[var(--color-surface-sunken)] rounded-lg px-3 py-2">
                          {showKekulized && moleculeInfo.kekulizedSmiles
                            ? moleculeInfo.kekulizedSmiles
                            : moleculeInfo.canonicalSmiles}
                        </code>
                      </div>
                    </>
                  ) : null}
                </motion.div>
              )}
            </AnimatePresence>
          </div>

          {/* Combined Tab Bar + Content */}
          <div className="card overflow-hidden">
            {/* Tab Bar */}
            <div className="flex flex-wrap gap-1 p-2 border-b border-[var(--color-border)] bg-[var(--color-surface-sunken)]/50">
              {TABS.map((tab) => (
                <button
                  key={tab.id}
                  onClick={() => setActiveTab(tab.id)}
                  className={cn(
                    'flex items-center gap-2 px-3 py-2 text-sm font-medium rounded-lg transition-all',
                    activeTab === tab.id
                      ? 'bg-[var(--color-surface-elevated)] text-[var(--color-primary)] shadow-sm'
                      : 'text-[var(--color-text-secondary)] hover:bg-[var(--color-surface-elevated)]/50 hover:text-[var(--color-text-primary)]'
                  )}
                >
                  {tab.icon}
                  <span className="hidden sm:inline">{tab.label}</span>
                </button>
              ))}
            </div>

            {/* Tab Content */}
            <AnimatePresence mode="wait">
              <motion.div
                key={activeTab}
                initial={{ opacity: 0, y: 10 }}
                animate={{ opacity: 1, y: 0 }}
                exit={{ opacity: 0, y: -10 }}
                transition={{ duration: 0.2 }}
                className="p-5 sm:p-6"
              >
                {/* Validate & Score Tab */}
                {activeTab === 'validate' && (
                  <div className="space-y-4">
                    <div className="flex items-start gap-3">
                      <div className="w-8 h-8 rounded-lg bg-[var(--color-primary)]/10 flex items-center justify-center text-[var(--color-primary)] flex-shrink-0">
                        <Info className="w-4 h-4" />
                      </div>
                      <p className="text-[var(--color-text-secondary)] text-sm">
                        Validate your chemical structure for correctness and calculate
                        ML-readiness scores for machine learning applications.
                      </p>
                    </div>
                    <div className="flex flex-wrap gap-3">
                      <ClayButton
                        variant="primary"
                        onClick={handleValidate}
                        disabled={!molecule.trim() || isAnyLoading}
                        loading={isLoading}
                        leftIcon={<Play className="w-4 h-4" />}
                      >
                        Validate
                      </ClayButton>
                      <ClayButton
                        variant="accent"
                        onClick={handleCalculateScores}
                        disabled={!molecule.trim() || isAnyLoading}
                        loading={scoringLoading}
                        leftIcon={<Sparkles className="w-4 h-4" />}
                      >
                        Score
                      </ClayButton>
                    </div>
                  </div>
                )}

                {/* Database Lookup Tab */}
                {activeTab === 'database' && (
                  <div className="space-y-4">
                    <div className="flex items-start gap-3">
                      <div className="w-8 h-8 rounded-lg bg-[var(--color-primary)]/10 flex items-center justify-center text-[var(--color-primary)] flex-shrink-0">
                        <Info className="w-4 h-4" />
                      </div>
                      <div>
                        <p className="text-[var(--color-text-secondary)] text-sm mb-3">
                          Search compound databases for additional information:
                        </p>
                        <ul className="list-none space-y-1 text-sm text-[var(--color-text-secondary)]">
                          <li className="flex items-center gap-2">
                            <span className="w-1.5 h-1.5 rounded-full bg-blue-500"></span>
                            <strong className="text-[var(--color-text-primary)]">PubChem</strong> — NIH compound database (100M+ compounds)
                          </li>
                          <li className="flex items-center gap-2">
                            <span className="w-1.5 h-1.5 rounded-full bg-purple-500"></span>
                            <strong className="text-[var(--color-text-primary)]">ChEMBL</strong> — Bioactivity and drug data
                          </li>
                          <li className="flex items-center gap-2">
                            <span className="w-1.5 h-1.5 rounded-full bg-amber-500"></span>
                            <strong className="text-[var(--color-text-primary)]">COCONUT</strong> — Natural products database
                          </li>
                        </ul>
                      </div>
                    </div>
                    <ClayButton
                      variant="primary"
                      onClick={handleDatabaseLookup}
                      disabled={!molecule.trim() || isAnyLoading}
                      loading={databaseLoading}
                      leftIcon={<Search className="w-4 h-4" />}
                    >
                      Look Up
                    </ClayButton>
                  </div>
                )}

                {/* Alerts Screening Tab */}
                {activeTab === 'alerts' && (
                  <div className="space-y-4">
                    <div className="flex items-start gap-3">
                      <div className="w-8 h-8 rounded-lg bg-amber-500/10 flex items-center justify-center text-amber-500 flex-shrink-0">
                        <Info className="w-4 h-4" />
                      </div>
                      <p className="text-[var(--color-text-secondary)] text-sm">
                        Screen for problematic structural patterns that may cause issues
                        in assays or drug development.
                      </p>
                    </div>

                    {/* Catalog selector */}
                    <div>
                      <p className="text-xs text-[var(--color-text-muted)] mb-2">Select catalogs to screen:</p>
                      <div className="flex flex-wrap gap-2">
                        {['PAINS', 'BRENK', 'NIH', 'ZINC'].map((catalog) => (
                          <button
                            key={catalog}
                            onClick={() => toggleCatalog(catalog)}
                            className={cn(
                              'px-3 py-1.5 text-sm rounded-lg transition-all',
                              selectedCatalogs.includes(catalog)
                                ? 'bg-[var(--color-primary)]/15 text-[var(--color-primary)] border border-[var(--color-primary)]/30'
                                : 'bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)] border border-transparent hover:border-[var(--color-border)]'
                            )}
                          >
                            {catalog}
                          </button>
                        ))}
                      </div>
                    </div>

                    <ClayButton
                      variant="primary"
                      onClick={handleScreenAlerts}
                      disabled={!molecule.trim() || isAnyLoading || selectedCatalogs.length === 0}
                      loading={alertsLoading}
                      leftIcon={<AlertTriangle className="w-4 h-4" />}
                    >
                      Screen Alerts
                    </ClayButton>
                  </div>
                )}

                {/* Standardize Tab */}
                {activeTab === 'standardize' && (
                  <div className="space-y-4">
                    <div className="flex items-start gap-3">
                      <div className="w-8 h-8 rounded-lg bg-[var(--color-primary)]/10 flex items-center justify-center text-[var(--color-primary)] flex-shrink-0">
                        <Info className="w-4 h-4" />
                      </div>
                      <div>
                        <p className="text-[var(--color-text-secondary)] text-sm mb-3">
                          Standardize your structure using the ChEMBL pipeline. This includes:
                        </p>
                        <ul className="list-none space-y-1 text-sm text-[var(--color-text-secondary)]">
                          <li className="flex items-center gap-2">
                            <span className="w-1.5 h-1.5 rounded-full bg-[var(--color-primary)]"></span>
                            Salt and solvent removal
                          </li>
                          <li className="flex items-center gap-2">
                            <span className="w-1.5 h-1.5 rounded-full bg-[var(--color-primary)]"></span>
                            Charge neutralization
                          </li>
                          <li className="flex items-center gap-2">
                            <span className="w-1.5 h-1.5 rounded-full bg-[var(--color-primary)]"></span>
                            Stereochemistry normalization
                          </li>
                          <li className="flex items-center gap-2">
                            <span className="w-1.5 h-1.5 rounded-full bg-[var(--color-primary)]"></span>
                            Tautomer canonicalization (optional)
                          </li>
                        </ul>
                      </div>
                    </div>

                    <label className="flex items-center gap-2 text-sm text-[var(--color-text-secondary)] cursor-pointer">
                      <input
                        type="checkbox"
                        checked={includeTautomer}
                        onChange={(e) => setIncludeTautomer(e.target.checked)}
                        className="rounded border-[var(--color-border-strong)] text-[var(--color-primary)] focus:ring-[var(--color-primary)]/30"
                      />
                      Enable tautomer canonicalization
                    </label>

                    <ClayButton
                      variant="primary"
                      onClick={handleStandardize}
                      disabled={!molecule.trim() || isAnyLoading}
                      loading={standardizationLoading}
                      leftIcon={<Layers className="w-4 h-4" />}
                    >
                      Standardize
                    </ClayButton>
                  </div>
                )}
              </motion.div>
            </AnimatePresence>
          </div>

          {/* Database Results - Separate Box */}
          <AnimatePresence>
            {databaseResults && activeTab === 'database' && (
              <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                exit={{ opacity: 0, y: -20 }}
                className="card p-5 sm:p-6"
              >
                <div className="flex items-center gap-3 mb-4">
                  <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-[var(--color-primary)]/10 to-[var(--color-accent)]/10 flex items-center justify-center text-[var(--color-primary)]">
                    <Database className="w-5 h-5" />
                  </div>
                  <div>
                    <h4 className="font-semibold text-[var(--color-text-primary)] text-sm tracking-tight">
                      Database Results
                    </h4>
                    <p className="text-xs text-[var(--color-text-muted)] mt-0.5">Cross-reference results</p>
                  </div>
                </div>
                <DatabaseLookupResults results={databaseResults} />
              </motion.div>
            )}
          </AnimatePresence>

          {/* Loading State */}
          <AnimatePresence>
            {isAnyLoading && (
              <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                exit={{ opacity: 0, y: -20 }}
                className="card-glass p-6"
              >
                <MoleculeLoader
                  size="md"
                  text={
                    isLoading ? 'Running validation checks...' :
                    databaseLoading ? 'Querying external databases...' :
                    scoringLoading ? 'Calculating scores...' :
                    standardizationLoading ? 'Running standardization pipeline...' :
                    'Screening for structural alerts...'
                  }
                />
              </motion.div>
            )}
          </AnimatePresence>

          {/* Error State */}
          <AnimatePresence>
            {hasError && !isAnyLoading && (
              <motion.div
                initial={{ opacity: 0, scale: 0.95 }}
                animate={{ opacity: 1, scale: 1 }}
                exit={{ opacity: 0, scale: 0.95 }}
                className="card p-5 border-red-500/30"
              >
                <div className="flex items-start gap-4">
                  <div className="w-10 h-10 rounded-xl bg-red-500/10 flex items-center justify-center flex-shrink-0">
                    <AlertTriangle className="w-5 h-5 text-red-500" />
                  </div>
                  <div className="flex-1">
                    <h3 className="font-semibold text-red-500 mb-1">
                      {((error?.error || alertError?.error || scoringError?.error || standardizationError?.error) as string)?.includes('parse') ? 'Parse Error' : 'Error'}
                    </h3>
                    <p className="text-sm text-[var(--color-text-secondary)]">
                      {error?.error || alertError?.error || scoringError?.error || standardizationError?.error}
                    </p>
                  </div>
                </div>
              </motion.div>
            )}
          </AnimatePresence>
        </motion.div>

        {/* RIGHT COLUMN */}
        <motion.div
          className="space-y-4"
          initial={{ opacity: 0, x: 20 }}
          animate={{ opacity: 1, x: 0 }}
          transition={{ duration: 0.5, delay: 0.3 }}
        >
          {/* Standardization Comparison View - replaces normal right column */}
          {showStandardizationComparison ? (
            <motion.div
              initial={{ opacity: 0, scale: 0.95 }}
              animate={{ opacity: 1, scale: 1 }}
              className="card p-5 sm:p-6"
            >
              <div className="flex items-center gap-3 mb-4">
                <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-[var(--color-primary)]/10 to-[var(--color-accent)]/10 flex items-center justify-center text-[var(--color-primary)]">
                  <Layers className="w-5 h-5" />
                </div>
                <div className="flex-1">
                  <h4 className="font-semibold text-[var(--color-text-primary)] text-sm tracking-tight">
                    Standardization Results
                  </h4>
                  <p className="text-xs text-[var(--color-text-muted)] mt-0.5">ChEMBL structure pipeline</p>
                </div>
                <Badge variant={standardizationResult.result.steps_applied.filter(s => s.applied).length > 0 ? 'info' : 'success'}>
                  {standardizationResult.result.steps_applied.filter(s => s.applied).length} changes
                </Badge>
              </div>
              <StandardizationResults result={standardizationResult.result} />
              <p className="mt-4 text-xs text-[var(--color-text-muted)] text-right">
                Completed in {standardizationResult.execution_time_ms}ms
              </p>
            </motion.div>
          ) : (
            <>
              {/* Molecule Viewer */}
              <div className="card-glow p-5 sm:p-6">
                <div className="flex items-center gap-3 mb-4">
                  <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-[var(--color-primary)]/10 to-[var(--color-accent)]/10 flex items-center justify-center text-[var(--color-primary)]">
                    <Atom className="w-5 h-5" />
                  </div>
                  <div className="flex-1">
                    <h4 className="font-semibold text-[var(--color-text-primary)] text-sm tracking-tight">
                      Structure Preview
                    </h4>
                    <p className="text-xs text-[var(--color-text-muted)] mt-0.5">
                      {molecule ? 'Rendered with RDKit.js' : 'Enter a SMILES to preview'}
                    </p>
                  </div>
                  {/* CIP Labels Toggle - Show when molecule has stereochemistry */}
                  {moleculeInfo?.hasStereochemistry && (
                    <button
                      onClick={() => setShowCIP(!showCIP)}
                      className={cn(
                        'flex items-center gap-1.5 text-xs font-medium px-3 py-1.5 rounded-lg transition-all',
                        'border shadow-sm',
                        showCIP
                          ? 'bg-[var(--color-primary)] text-white border-[var(--color-primary)] shadow-[var(--color-primary)]/20'
                          : 'bg-[var(--color-surface-elevated)] text-[var(--color-text-secondary)] border-[var(--color-border)] hover:border-[var(--color-primary)] hover:text-[var(--color-primary)]'
                      )}
                      title="Show R/S and E/Z stereochemistry labels"
                    >
                      <svg className="w-3.5 h-3.5" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                        <circle cx="12" cy="12" r="10" />
                        <path d="M12 6v6l4 2" />
                      </svg>
                      {showCIP ? 'Hide CIP' : 'Show CIP'}
                    </button>
                  )}
                </div>
                <div className="molecule-preview rounded-xl min-h-[280px] flex items-center justify-center">
                  <MoleculeViewer
                    smiles={canonicalSmiles || molecule}
                    highlightAtoms={highlightedAtoms}
                    width={380}
                    height={280}
                    showCIP={showCIP}
                  />
                </div>
                {highlightedAtoms.length > 0 && (
                  <motion.p
                    initial={{ opacity: 0 }}
                    animate={{ opacity: 1 }}
                    className="mt-3 text-xs text-center text-amber-500 font-medium"
                  >
                    Highlighting atoms: {highlightedAtoms.join(', ')}
                  </motion.p>
                )}
                {/* Stereochemistry info indicator */}
                {moleculeInfo?.hasStereochemistry && (
                  <motion.div
                    initial={{ opacity: 0 }}
                    animate={{ opacity: 1 }}
                    className="mt-3 flex items-center justify-center gap-2"
                  >
                    <div className="flex items-center gap-1.5 text-xs bg-purple-500/10 text-purple-600 dark:text-purple-400 px-2 py-1 rounded-lg">
                      <svg className="w-3 h-3" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                        <path d="M12 2L2 7l10 5 10-5-10-5zM2 17l10 5 10-5M2 12l10 5 10-5" />
                      </svg>
                      <span>
                        {moleculeInfo.numStereocenters > 0 && `${moleculeInfo.numStereocenters} stereocenter${moleculeInfo.numStereocenters > 1 ? 's' : ''}`}
                        {moleculeInfo.numStereocenters > 0 && /[/\\]/.test(moleculeInfo.canonicalSmiles) && ' + '}
                        {/[/\\]/.test(moleculeInfo.canonicalSmiles) && 'E/Z bonds'}
                      </span>
                    </div>
                    {showCIP && (
                      <span className="text-xs text-[var(--color-text-muted)]">
                        (CIP labels shown)
                      </span>
                    )}
                  </motion.div>
                )}
              </div>

              {/* Score Tiles - Only show after validation/scoring */}
              <AnimatePresence>
                {hasScores && (
                  <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    exit={{ opacity: 0, y: -20 }}
                    className="grid grid-cols-2 gap-4"
                  >
                    {/* Quality Score */}
                    <div className="card-gradient p-4 sm:p-5">
                      <div className="flex items-start justify-between mb-2">
                        <div className="flex items-center gap-2">
                          <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-[var(--color-primary)]/10 to-[var(--color-accent)]/10 flex items-center justify-center text-[var(--color-primary)]">
                            <Target className="w-4 h-4" />
                          </div>
                          <span className="text-xs font-medium text-[var(--color-text-secondary)]">Quality</span>
                        </div>
                        {qualityScore !== null && (
                          <Badge variant={qualityScore >= 70 ? 'success' : qualityScore >= 40 ? 'warning' : 'error'} size="sm">
                            {getScoreLabel(qualityScore)}
                          </Badge>
                        )}
                      </div>
                      <div className="flex items-baseline gap-1">
                        <span className={cn(
                          'text-3xl font-bold tracking-tight',
                          qualityScore !== null && qualityScore >= 70 && 'text-gradient text-glow'
                        )}>
                          {qualityScore !== null ? Math.round(qualityScore) : '--'}
                        </span>
                        <span className="text-xs text-[var(--color-text-muted)]">/100</span>
                      </div>
                    </div>

                    {/* ML Readiness */}
                    <div className="card-accent p-4 sm:p-5">
                      <div className="flex items-center justify-between mb-2">
                        <div className="flex items-center gap-2">
                          <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-[var(--color-accent)]/10 to-[var(--color-primary)]/10 flex items-center justify-center text-[var(--color-accent)]">
                            <Zap className="w-4 h-4" />
                          </div>
                          <span className="text-xs font-medium text-[var(--color-text-secondary)]">ML Ready</span>
                        </div>
                        {mlReady !== null && (
                          <Badge variant={mlReady ? 'success' : 'warning'} dot size="sm">
                            {mlReady ? 'Ready' : 'Review'}
                          </Badge>
                        )}
                      </div>
                      <div className="flex justify-center">
                        {mlReadyScore !== undefined ? (
                          <ScoreChart
                            score={mlReadyScore}
                            label="ML-Readiness"
                            size={100}
                            compact
                          />
                        ) : (
                          <div className="w-[100px] h-[100px] flex items-center justify-center">
                            <span className="text-3xl font-bold text-[var(--color-text-muted)]">--</span>
                          </div>
                        )}
                      </div>
                    </div>
                  </motion.div>
                )}
              </AnimatePresence>

              {/* Issues / Alerts Panel */}
              <AnimatePresence>
                {/* Validation Issues */}
                {result && validationIssues.length > 0 && (
                  <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    exit={{ opacity: 0, y: -20 }}
                    className="card p-5 sm:p-6"
                  >
                    <div className="flex items-center justify-between mb-4">
                      <h4 className="font-semibold text-[var(--color-text-primary)] text-sm">
                        Validation Issues
                      </h4>
                      <Badge variant="warning">{validationIssues.length} found</Badge>
                    </div>
                    <div className="space-y-3 max-h-[400px] overflow-y-auto pr-2">
                      {validationIssues.map((issue, index) => (
                        <IssueCard
                          key={`${issue.check_name}-${index}`}
                          issue={issue}
                          onAtomHover={setHighlightedAtoms}
                        />
                      ))}
                    </div>
                    {result && (
                      <p className="mt-4 text-xs text-[var(--color-text-muted)] text-right">
                        Completed in {result.execution_time_ms.toFixed(0)}ms
                      </p>
                    )}
                  </motion.div>
                )}

                {/* Validation Success - no issues */}
                {result && validationIssues.length === 0 && (
                  <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    exit={{ opacity: 0, y: -20 }}
                    className="rounded-xl p-5 text-center bg-yellow-500/10 border border-yellow-500/20"
                  >
                    <div className="text-4xl mb-2">✓</div>
                    <h3 className="text-lg font-semibold text-amber-600 dark:text-yellow-400 mb-1">
                      No Issues Found
                    </h3>
                    <p className="text-sm text-amber-600/80 dark:text-yellow-400/80">
                      All validation checks passed successfully
                    </p>
                    <p className="mt-3 text-xs text-[var(--color-text-muted)]">
                      Completed in {result.execution_time_ms.toFixed(0)}ms
                    </p>
                  </motion.div>
                )}

                {/* Alert Screening Results */}
                {alertResult && (
                  <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    exit={{ opacity: 0, y: -20 }}
                    className="card p-5 sm:p-6"
                  >
                    <div className="flex items-center justify-between mb-4">
                      <div>
                        <h4 className="font-semibold text-[var(--color-text-primary)] text-sm">
                          Structural Alerts
                        </h4>
                        <p className="text-xs text-[var(--color-text-muted)] mt-0.5">
                          Screened: {alertResult.screened_catalogs.join(', ')}
                        </p>
                      </div>
                      <Badge variant={alertIssues.length === 0 ? 'success' : 'warning'}>
                        {alertIssues.length} alerts
                      </Badge>
                    </div>

                    {alertIssues.length > 0 ? (
                      <div className="space-y-3 max-h-[400px] overflow-y-auto pr-2">
                        {alertIssues.map((alert, index) => (
                          <AlertCard
                            key={`${alert.pattern_name}-${index}`}
                            alert={alert}
                            onAtomHover={setHighlightedAtoms}
                          />
                        ))}
                      </div>
                    ) : (
                      <div className="rounded-xl p-4 text-center bg-yellow-500/10 border border-yellow-500/20">
                        <div className="text-2xl mb-1">✓</div>
                        <p className="text-sm text-amber-600 dark:text-yellow-400">
                          No structural alerts detected
                        </p>
                      </div>
                    )}

                    <p className="mt-4 text-xs text-[var(--color-text-muted)] text-right">
                      Completed in {alertResult.execution_time_ms}ms
                    </p>
                  </motion.div>
                )}
              </AnimatePresence>

            </>
          )}
        </motion.div>
      </div>

      {/* Scoring Results - Full Width Below Grid */}
      <AnimatePresence>
        {scoringResult && activeTab === 'validate' && (
          <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -20 }}
            className="card p-6 sm:p-8"
          >
            <ScoringResults scoringResponse={scoringResult} />
          </motion.div>
        )}
      </AnimatePresence>

      {/* Share URL Toast */}
      <AnimatePresence>
        {shareToastVisible && (
          <motion.div
            initial={{ opacity: 0, y: 50 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: 50 }}
            className="fixed bottom-6 left-1/2 -translate-x-1/2 z-50"
          >
            <div className="flex items-center gap-3 px-4 py-3 rounded-xl bg-[var(--color-surface-elevated)] border border-[var(--color-border-strong)] shadow-2xl">
              <CheckCircle2 className="w-5 h-5 text-green-500" />
              <span className="text-sm font-medium text-[var(--color-text-primary)]">
                Share URL copied to clipboard!
              </span>
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}
