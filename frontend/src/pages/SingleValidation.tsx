import { lazy, Suspense, useState, useEffect, useLayoutEffect, useRef, useCallback, useMemo } from 'react';
import { useSearchParams, useLocation, useNavigate } from 'react-router-dom';
import { useValidationCache, type TabType } from '../contexts/ValidationCacheContext';
import { motion, AnimatePresence } from 'framer-motion';
import {
  Atom,
  Beaker,
  AlertTriangle,
  Database,
  Sparkles,
  RotateCcw,
  Layers,
  CheckCircle2,
  Share2,
  Microscope,
  BarChart3,
  ArrowLeft,
  Hexagon,
  GitMerge,
  Scale,
  FlaskConical,
  Shield,
} from 'lucide-react';
import { StructureInput } from '../components/molecules/StructureInput';
import { ExampleMolecules } from '../components/validation/ExampleMolecules';
import { ErrorPanel, LoadingPanel } from '../components/validation/StatusPanels';
import { ScoreTiles } from '../components/validation/ScoreTiles';
import { ValidationOutcomePanels } from '../components/validation/ValidationOutcomePanels';
import { ValidationIssuesPanel } from '../components/validation/ValidationIssuesPanel';
import { MoleculeViewerPanel } from '../components/validation/MoleculeViewerPanel';
import {
  AllChecksCard,
  type AllChecksBounds,
  type AllChecksPhase,
} from '../components/validation/AllChecksCard';
import {
  ValidateTabPanel,
  DeepValidationTabPanel,
  ScoringProfilesTabPanel,
  CompoundProfileTabPanel,
  AlertsTabPanel,
  StandardizeTabPanel,
} from '../components/validation/TabInputPanels';
const ScoringResults = lazy(() =>
  import('../components/scoring/ScoringResults').then((m) => ({ default: m.ScoringResults }))
);
const StandardizationResults = lazy(() =>
  import('../components/standardization/StandardizationResults').then((m) => ({ default: m.StandardizationResults }))
);
import { DatabaseLookupControls, DatabaseResultsBox } from '../components/integrations/DatabaseLookupTab';
import { ClayButton } from '../components/ui/ClayButton';
import { Badge } from '../components/ui/Badge';
import { CopyButton } from '../components/ui/CopyButton';
import { InfoTooltip } from '../components/ui/Tooltip';
import { TabBar, type TabBarTab, type TabBarResultState } from '../components/ui/TabBar';
import { useValidation } from '../hooks/useValidation';
import { useMoleculeInfo } from '../hooks/useMoleculeInfo';
import { useRecentMolecules } from '../hooks/useRecentMolecules';
import { alertsApi, scoringApi, standardizationApi, integrationsApi } from '../services/api';
import { BookmarkButton } from '../components/bookmarks/BookmarkButton';
import { cn } from '../lib/utils';
import { saveSnapshot, getSnapshot } from '../lib/bookmarkStore';
import { logger } from '../lib/logger';
import type { AlertScreenResponse, AlertError } from '../types/alerts';
import type { ScoringResponse, ScoringError } from '../types/scoring';
import type { StandardizeResponse, StandardizeError } from '../types/standardization';
import type { PubChemResult, ChEMBLResult, COCONUTResult, WikidataResult, ResolvedCompound, ConsistencyResult } from '../types/integrations';
import { IdentifierResolverCard } from '../components/integrations/IdentifierResolverCard';
import { DatabaseComparisonPanel } from '../components/integrations/DatabaseComparisonPanel';
import { useProfiler } from '../hooks/useProfiler';
import { useSafety } from '../hooks/useSafety';

const CHEMBL_CATALOGS = [
  { id: 'CHEMBL_BMS', label: 'BMS HTS Filters' },
  { id: 'CHEMBL_DUNDEE', label: 'Dundee NTD Filters' },
  { id: 'CHEMBL_GLAXO', label: 'Glaxo Hard Filters' },
  { id: 'CHEMBL_INPHARMATICA', label: 'Inpharmatica' },
  { id: 'CHEMBL_LINT', label: 'Lilly MedChem (LINT)' },
  { id: 'CHEMBL_MLSMR', label: 'NIH MLSMR' },
  { id: 'CHEMBL_SURECHEMBL', label: 'SureChEMBL' },
];

type InputType = 'smiles' | 'iupac' | 'identifier' | 'ambiguous';

/** Characters that appear in SMILES but never in IUPAC names or identifiers. */
const SMILES_CHARS = /[[\]()=#@/\\]/;

/** Identifier patterns for database IDs resolved via the backend. */
const IDENTIFIER_PATTERNS: RegExp[] = [
  /^CHEMBL\d+$/i,           // ChEMBL ID
  /^DB\d{5}$/i,             // DrugBank ID
  /^CHEBI:\d+$/i,           // ChEBI ID
  /^\d{2,7}-\d{2}-\d$/,     // CAS number
  /^[A-Z0-9]{14}-[A-Z0-9]{10}-[A-Z0-9]$/, // InChIKey
  /^CID:\d+$/i,             // PubChem CID
];

/**
 * Detect whether input looks like SMILES, IUPAC, or a database identifier.
 * Database identifiers (ChEMBL ID, CAS, DrugBank, etc.) are resolved via the backend.
 */
function detectInputType(value: string): InputType {
  const trimmed = value.trim();
  if (!trimmed) return 'ambiguous';

  // Database identifiers — detected before SMILES/IUPAC
  if (IDENTIFIER_PATTERNS.some((p) => p.test(trimmed))) return 'identifier';
  if (trimmed.includes('wikipedia.org/wiki/')) return 'identifier';

  // Strong SMILES indicators: brackets, ring closures, branch notation
  if (SMILES_CHARS.test(trimmed)) return 'smiles';

  // If starts with InChI marker
  if (trimmed.startsWith('InChI=')) return 'smiles'; // treat InChI as SMILES-like (non-IUPAC)

  // IUPAC suffixes
  const iupacSuffixes = [
    'ane', 'ene', 'yne', 'ol', 'al', 'one', 'oic acid', 'amine', 'amide',
    'ate', 'ester', 'ether', 'ose', 'ide', 'ine', 'ite', 'yl', 'benzene',
    'phenol', 'acetic', 'acid', 'alcohol',
  ];
  const lower = trimmed.toLowerCase();
  const hasIupacSuffix = iupacSuffixes.some((s) => lower.endsWith(s));

  // Contains spaces (typical for IUPAC, never in SMILES)
  if (trimmed.includes(' ') || hasIupacSuffix) return 'iupac';

  // All alphabetic with hyphens and numbers (common IUPAC pattern)
  if (/^[a-zA-Z0-9,'-]+$/.test(trimmed) && trimmed.length > 8 && /[a-z]{3,}/.test(lower)) {
    return 'iupac';
  }

  return 'ambiguous';
}

/** Render a molecular formula like "C8H10N4O2" with the digits subscripted. */
/**
 * Classify a parse / processing error message into a sub-type so the user
 * gets a heading + hint specific to the actual failure mode rather than a
 * generic "Parse Error". String-matches the error against common chemistry
 * failure keywords. Order matters — most specific match first.
 */
type ParseErrorType = 'atom' | 'ring' | 'valence' | 'parse-other' | 'database' | 'generic';

function classifyError(message: string | null | undefined, isDatabase: boolean): ParseErrorType {
  if (isDatabase) return 'database';
  if (!message) return 'generic';
  const lower = message.toLowerCase();
  if (lower.includes('valence')) return 'valence';
  if (lower.includes('ring') || lower.includes('closure')) return 'ring';
  if (lower.includes('atom') || lower.includes('element') || lower.includes('symbol')) return 'atom';
  if (lower.includes('parse') || lower.includes('parser')) return 'parse-other';
  return 'generic';
}

const PARSE_ERROR_HEADINGS: Record<ParseErrorType, string> = {
  atom: 'Unrecognised atom symbol',
  ring: 'Ring closure mismatch',
  valence: 'Invalid valence',
  'parse-other': 'Structure could not be parsed',
  database: 'Database Lookup Failed',
  generic: 'Error',
};

const PARSE_ERROR_HINTS: Record<ParseErrorType, string> = {
  atom: 'Check element symbols are spelled correctly. Two-letter symbols like Cl or Br must be capitalised exactly. Lowercase letters are reserved for aromatic atoms (c, n, o, s).',
  ring: 'Each ring opening digit needs a matching closing digit. Two-digit ring labels use a percent sign (e.g. C1CC%10CC%10CC1).',
  valence: 'Check that nitrogen has no more than 3 bonds unless charged ([N+]), and that sulfur oxidation state is consistent. Hypervalent atoms may need explicit brackets.',
  'parse-other': 'Accepted: SMILES, InChI, IUPAC name, ChEMBL ID, CAS number, or MDL Mol file. Common causes are typos in atom symbols or unbalanced ring closures.',
  database: 'PubChem, ChEMBL, and COCONUT may be temporarily unavailable. Check the SMILES is canonicalisable and try again in a moment.',
  generic: 'Accepted: SMILES, InChI, IUPAC name, ChEMBL ID, CAS number, or MDL Mol file.',
};

/** Compact severity summary tags for the quality card */
// Tabs split across two rows for visual breathing room.
// Row 1: most-frequent core actions (validate, surface safety, normalize, lookup).
// Row 2: deeper / power-user analyses.
const TAB_ROW_1: ReadonlyArray<TabBarTab<TabType>> = [
  {
    id: 'validate',
    label: 'Validate & Score',
    icon: CheckCircle2,
    description: 'Check structure validity, calculate quality metrics, and assess ML-readiness',
  },
  {
    id: 'alerts',
    label: 'Safety',
    icon: Shield,
    description: 'Structural alerts, CYP soft-spots, hERG, bRo5, REOS, and complexity analysis',
  },
  {
    id: 'standardize',
    label: 'Standardize',
    icon: Layers,
    description: 'Normalize structure representation and remove salts/solvents',
  },
  {
    id: 'database',
    label: 'Database Lookup',
    icon: Database,
    description: 'Search PubChem, ChEMBL, and COCONUT for compound information',
  },
];

const TAB_ROW_2: ReadonlyArray<TabBarTab<TabType>> = [
  {
    id: 'deep-validation',
    label: 'Deep Validation',
    icon: Microscope,
    description: 'Advanced structure checks: stereo, tautomers, composition, and complexity analysis',
  },
  {
    id: 'scoring-profiles',
    label: 'Scoring Profiles',
    icon: BarChart3,
    description: 'Consensus drug-likeness, lead-likeness, property breakdowns, bioavailability radar',
  },
  {
    id: 'compound-profile',
    label: 'Compound Profile',
    icon: FlaskConical,
    description: 'Full molecular profiling: PFI, stars, bioavailability, synthesizability, 3D shape',
  },
];

/** Shape of `location.state` when navigating to SingleValidation */
interface LocationState {
  bookmarkId?: number;
  smiles?: string;
  fromBatch?: boolean;
  moleculeName?: string;
  moleculeIndex?: number;
}

export function SingleValidationPage() {
  const [searchParams, setSearchParams] = useSearchParams();
  const location = useLocation();
  const navigate = useNavigate();
  const { getCache, setCache, clearCache } = useValidationCache();
  const [molecule, setMolecule] = useState('');

  // Batch origin context — shows back-to-batch bar when navigated from BatchResultsTable
  const locState = location.state as LocationState | null;
  const batchOrigin = locState?.fromBatch
    ? { moleculeName: locState.moleculeName, moleculeIndex: locState.moleculeIndex }
    : null;
  const [highlightedAtoms, setHighlightedAtoms] = useState<number[]>([]);
  const [highlightLocked, setHighlightLocked] = useState(false);
  const [activeTab, setActiveTab] = useState<TabType>('validate');
  const { validate, result, error, isLoading, reset, restore } = useValidation();
  const [shareToastVisible, setShareToastVisible] = useState(false);
  const { recent, addRecent, removeRecent, clearRecent } = useRecentMolecules();

  // Enrichment: Profiler hook
  const {
    profile: profileResult,
    isLoading: profileLoading,
    error: profileError,
    profileCompound,
  } = useProfiler();

  // Enrichment: Safety hook
  const {
    alertResult: safetyAlertResult,
    safetyResult: safetyAssessResult,
    isLoading: safetyLoading,
    alertError: safetyAlertError,
    safetyError: safetyAssessError,
    screenMolecule: screenSafety,
  } = useSafety();

  // URL deep-link: ?section= param maps to a tab
  const sectionParam = searchParams.get('section');
  useEffect(() => {
    if (sectionParam === 'profile') setActiveTab('compound-profile');
    else if (sectionParam === 'safety') setActiveTab('alerts');
  }, [sectionParam]);

  // Load molecule from URL on mount
  useEffect(() => {
    const smilesFromUrl = searchParams.get('smiles');
    if (smilesFromUrl) {
      // searchParams.get() already decodes URI components — no need for decodeURIComponent
      setMolecule(smilesFromUrl);
      // validate() takes the SMILES as an argument, so we don't need to wait
      // for setMolecule to flush — pass the URL value directly.
      validate({
        molecule: smilesFromUrl,
        format: 'auto',
        preserve_aromatic: false,
      });
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []); // Only run on mount

  // Add to recent molecules after successful validation
  useEffect(() => {
    if (result && !error && molecule.trim()) {
      addRecent(molecule.trim());
    }
  }, [result, error, molecule, addRecent]);

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
  const [chemblExpanded, setChemblExpanded] = useState(false);

  const toggleAllChembl = (enabled: boolean) => {
    const chemblIds = CHEMBL_CATALOGS.map((c) => c.id);
    if (enabled) {
      setSelectedCatalogs((prev) => [...new Set([...prev, ...chemblIds])]);
    } else {
      setSelectedCatalogs((prev) => prev.filter((c) => !chemblIds.includes(c)));
    }
  };

  const allChemblSelected = CHEMBL_CATALOGS.every((c) => selectedCatalogs.includes(c.id));
  const someChemblSelected = CHEMBL_CATALOGS.some((c) => selectedCatalogs.includes(c.id));

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

  // All checks panel — single floating card with sequenced translate-then-grow
  // animation. The card is rendered ONCE as a position:absolute motion.div
  // floating above the page. Two empty anchor divs (one in the right column,
  // one below the grid) reserve layout space and provide measurement targets.
  // The card animates top/left/width/height between the two anchors:
  //   - Expand: top/left first (300ms), then width/height grow (400ms after).
  //   - Collapse: width/height shrink first (400ms), then top/left return (300ms after).
  // Same DOM element throughout. No mount/unmount, no fade.
  const [allChecksPhase, setAllChecksPhase] = useState<AllChecksPhase>('collapsed');
  const allChecksFloatContainerRef = useRef<HTMLDivElement>(null);
  const allChecksCollapsedAnchorRef = useRef<HTMLDivElement>(null);
  const allChecksExpandedAnchorRef = useRef<HTMLDivElement>(null);
  const [allChecksBounds, setAllChecksBounds] = useState<AllChecksBounds | null>(null);

  const measureAllChecksBounds = useCallback(() => {
    const container = allChecksFloatContainerRef.current;
    const collapsed = allChecksCollapsedAnchorRef.current;
    const expanded = allChecksExpandedAnchorRef.current;
    if (!container || !collapsed || !expanded) return;

    // Walk the offsetParent chain from the anchor up to (but not including)
    // the container, accumulating offsetTop/offsetLeft at each step.
    // This gives position relative to the container in the static layout flow,
    // unaffected by CSS transforms (e.g. the right-column slide-in animation).
    // getBoundingClientRect() includes transforms and produces a wrong left
    // during the right-column entry animation, causing the card to jitter.
    function offsetRelativeTo(el: HTMLElement): { top: number; left: number } {
      let top = 0;
      let left = 0;
      let cur: HTMLElement | null = el;
      while (cur && cur !== container) {
        top += cur.offsetTop;
        left += cur.offsetLeft;
        cur = cur.offsetParent as HTMLElement | null;
      }
      return { top, left };
    }

    const cPos = offsetRelativeTo(collapsed);
    const ePos = offsetRelativeTo(expanded);

    setAllChecksBounds({
      collapsed: {
        top: cPos.top,
        left: cPos.left,
        width: collapsed.offsetWidth,
        height: collapsed.offsetHeight,
      },
      expanded: {
        top: ePos.top,
        left: ePos.left,
        width: expanded.offsetWidth,
        height: expanded.offsetHeight,
      },
    });
  }, []);


  // Measure the actual rendered content height with a ResizeObserver so
  // the floating card sizes to its content exactly — no estimation, no
  // phantom whitespace at the bottom of the panel. Callback-ref pattern
  // so the observer attaches the moment the inner content div mounts
  // (which happens AFTER allChecksBounds is set, so a useEffect with
  // stale deps would never fire). Synchronous scrollHeight read on
  // attach gives the floating card the right size on its first paint —
  // no clipping when the user expands the panel.
  const HEADER_STRIP_HEIGHT = 62;
  const allChecksObserverRef = useRef<ResizeObserver | null>(null);
  const [allChecksContentHeight, setAllChecksContentHeight] = useState(0);
  const allChecksContentRef = useCallback((node: HTMLDivElement | null) => {
    if (allChecksObserverRef.current) {
      allChecksObserverRef.current.disconnect();
      allChecksObserverRef.current = null;
    }
    if (!node) return;
    setAllChecksContentHeight(node.scrollHeight);
    if (typeof ResizeObserver === 'undefined') return;
    const ro = new ResizeObserver((entries) => {
      const entry = entries[0];
      if (!entry) return;
      const h = entry.borderBoxSize?.[0]?.blockSize ?? entry.contentRect.height;
      setAllChecksContentHeight(Math.ceil(h));
    });
    ro.observe(node);
    allChecksObserverRef.current = ro;
  }, []);
  const allChecksExpandedHeight = allChecksContentHeight + HEADER_STRIP_HEIGHT;

  useLayoutEffect(() => {
    measureAllChecksBounds();
    window.addEventListener('resize', measureAllChecksBounds);
    return () => window.removeEventListener('resize', measureAllChecksBounds);
  }, [measureAllChecksBounds, result, allChecksPhase, activeTab, allChecksExpandedHeight]);

  // Auto-advance phase when transitions complete (700ms = 300ms move + 400ms grow).
  useEffect(() => {
    if (allChecksPhase === 'expanding') {
      const id = window.setTimeout(() => setAllChecksPhase('expanded'), 750);
      return () => window.clearTimeout(id);
    }
    if (allChecksPhase === 'collapsing') {
      const id = window.setTimeout(() => setAllChecksPhase('collapsed'), 750);
      return () => window.clearTimeout(id);
    }
  }, [allChecksPhase]);

  const toggleAllChecks = useCallback(() => {
    if (allChecksPhase === 'collapsed') setAllChecksPhase('expanding');
    else if (allChecksPhase === 'expanded') setAllChecksPhase('collapsing');
  }, [allChecksPhase]);

  // Molecule preview ref for image download
  const previewRef = useRef<HTMLDivElement>(null);
  const scoringAnchorRef = useRef<HTMLDivElement>(null);
  const comparisonAnchorRef = useRef<HTMLDivElement>(null);

  // Input type auto-detection
  const detectedType = molecule.trim() ? detectInputType(molecule.trim()) : 'ambiguous';
  const [forceInputType, setForceInputType] = useState<'smiles' | 'iupac' | null>(null);
  const effectiveInputType = forceInputType ?? (detectedType === 'ambiguous' || detectedType === 'identifier' ? null : detectedType);

  const isIdentifierInput = detectedType === 'identifier';
  const isIupacInput = effectiveInputType === 'iupac';

  // Database lookup state
  const [databaseResults, setDatabaseResults] = useState<{
    pubchem: PubChemResult | null;
    chembl: ChEMBLResult | null;
    coconut: COCONUTResult | null;
    wikidata: WikidataResult | null;
  } | null>(null);
  const [databaseLoading, setDatabaseLoading] = useState(false);
  const [databaseError, setDatabaseError] = useState<string | null>(null);

  // Identifier resolution state
  const [resolverResult, setResolverResult] = useState<ResolvedCompound | null>(null);

  // Cross-database comparison state
  const [comparisonResult, setComparisonResult] = useState<ConsistencyResult | null>(null);
  const [isComparing, setIsComparing] = useState(false);
  const [autoCompare, setAutoCompare] = useState(true);
  const [dbResultsExpanded, setDbResultsExpanded] = useState(true);

  // Use canonical SMILES from validation result when available (handles IUPAC/common names)
  const resolvedSmiles = result?.molecule_info?.canonical_smiles || molecule.trim();

  const handleValidate = async () => {
    if (!molecule.trim()) return;
    setResolverResult(null);
    setComparisonResult(null);

    // Try identifier resolution for non-SMILES input
    const input = molecule.trim();
    if (!SMILES_CHARS.test(input) && !input.startsWith('InChI=')) {
      try {
        const resolved = await integrationsApi.resolveIdentifier({ identifier: input });
        if (resolved.resolved && resolved.canonical_smiles && resolved.identifier_type_detected !== 'smiles') {
          setResolverResult(resolved);
          // Use resolved SMILES for validation
          await validate({
            molecule: resolved.canonical_smiles,
            format: 'auto',
            preserve_aromatic: preferAromaticSmiles,
            input_type: 'smiles',
          });
          return;
        }
      } catch {
        // Fall through to normal validation
      }
    }

    await validate({
      molecule: input,
      format: 'auto',
      preserve_aromatic: preferAromaticSmiles,
      input_type: effectiveInputType ?? 'auto',
    });
  };

  const handleScreenAlerts = async () => {
    if (!molecule.trim()) return;
    setAlertsLoading(true);
    setAlertError(null);
    try {
      const response = await alertsApi.screenAlerts({
        molecule: resolvedSmiles,
        format: 'smiles',
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
      const response = await scoringApi.getScoring(resolvedSmiles, 'smiles');
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
        molecule: resolvedSmiles,
        format: 'smiles',
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

  /** Run cross-database comparison for the current molecule. */
  const runComparison = async (): Promise<void> => {
    setIsComparing(true);
    try {
      const compareResult = await integrationsApi.compareAcrossDatabases({
        smiles: resolvedSmiles || undefined,
        inchikey: result?.molecule_info?.inchikey || undefined,
      });
      setComparisonResult(compareResult);
    } catch (err) {
      logger.error('Comparison error:', err);
    } finally {
      setIsComparing(false);
    }
  };

  const handleDatabaseLookup = async () => {
    if (!molecule.trim()) return;
    setDatabaseLoading(true);
    setDatabaseError(null);
    setComparisonResult(null);
    try {
      const results = await integrationsApi.lookupAll({ smiles: resolvedSmiles });
      setDatabaseResults(results);
      if (autoCompare) {
        await runComparison();
      }
    } catch (err) {
      logger.error('Database lookup error:', err);
      setDatabaseResults(null);
      const message =
        err instanceof Error
          ? err.message
          : 'Database lookup failed. PubChem, ChEMBL, and COCONUT may be temporarily unavailable.';
      setDatabaseError(message);
    } finally {
      setDatabaseLoading(false);
    }
  };

  const handleExampleClick = (smiles: string) => {
    setMolecule(smiles);
    resetAll();
  };

  const handleSelectRecent = (smiles: string) => {
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
    setDatabaseError(null);
    setResolverResult(null);
    setComparisonResult(null);
    setHighlightedAtoms([]);
    setHighlightLocked(false);
    setShowCIP(false);
  };

  // Restore cached results on mount (from navigation cache or bookmark snapshot)
  const didRestore = useRef(false);
  useEffect(() => {
    if (didRestore.current) return;
    didRestore.current = true;

    // Priority 0: Navigated from batch results — auto-validate the molecule
    if (locState?.fromBatch && locState.smiles) {
      setMolecule(locState.smiles);
      setActiveTab('validate');
      validate({ molecule: locState.smiles, format: 'auto', preserve_aromatic: false });
      return;
    }

    // Priority 1: Restore from bookmark snapshot (navigated from Bookmarks page)
    if (locState?.bookmarkId) {
      getSnapshot(locState.bookmarkId).then((snapshot) => {
        if (snapshot) {
          setMolecule(snapshot.molecule);
          setActiveTab('validate');
          if (snapshot.validationResult) restore(snapshot.validationResult);
          if (snapshot.alertResult) setAlertResult(snapshot.alertResult);
          if (snapshot.scoringResult) setScoringResult(snapshot.scoringResult);
          if (snapshot.standardizationResult) setStandardizationResult(snapshot.standardizationResult);
          if (snapshot.databaseResults) setDatabaseResults({ wikidata: null, ...snapshot.databaseResults });
        } else if (locState.smiles) {
          // No snapshot (different device / cleared) — auto-validate
          setMolecule(locState.smiles);
          validate({ molecule: locState.smiles, preserve_aromatic: false });
        }
      }).catch(() => {
        // IndexedDB error — fallback to auto-validate
        if (locState.smiles) {
          setMolecule(locState.smiles);
          validate({ molecule: locState.smiles, preserve_aromatic: false });
        }
      });
      // Clear location state to prevent re-restore on subsequent navigations
      window.history.replaceState({}, '');
      return;
    }

    // Priority 2: Restore from in-memory navigation cache
    const cached = getCache();
    if (cached && !searchParams.get('smiles')) {
      setMolecule(cached.molecule);
      setActiveTab(cached.activeTab);
      if (cached.result) restore(cached.result);
      if (cached.alertResult) setAlertResult(cached.alertResult);
      if (cached.scoringResult) setScoringResult(cached.scoringResult);
      if (cached.standardizationResult) setStandardizationResult(cached.standardizationResult);
      if (cached.databaseResults) setDatabaseResults({ wikidata: null, ...cached.databaseResults });
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // Persist results to cache whenever they change
  useEffect(() => {
    if (!molecule.trim()) return;
    setCache({
      molecule: molecule.trim(),
      activeTab,
      result,
      alertResult,
      scoringResult,
      standardizationResult,
      databaseResults,
    });
  }, [molecule, activeTab, result, alertResult, scoringResult, standardizationResult, databaseResults, setCache]);

  // Save snapshot to IndexedDB after bookmark is created
  const handleAfterBookmark = useCallback(async (bookmarkId: number) => {
    try {
      await saveSnapshot({
        bookmarkId,
        source: 'single_validation',
        molecule: molecule.trim(),
        savedAt: new Date().toISOString(),
        validationResult: result,
        alertResult,
        scoringResult,
        standardizationResult,
        databaseResults,
      });
    } catch (err) {
      logger.error('Failed to save bookmark snapshot:', err);
    }
  }, [molecule, result, alertResult, scoringResult, standardizationResult, databaseResults]);

  // Handler for locking atom highlighting (persists for SVG download)
  const handleAtomLock = useCallback((atoms: number[]) => {
    if (atoms.length > 0) {
      setHighlightedAtoms(atoms);
      setHighlightLocked(true);
    } else {
      setHighlightedAtoms([]);
      setHighlightLocked(false);
    }
  }, []);

  const handleReset = () => {
    setMolecule('');
    resetAll();
    clearCache();
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
      logger.error('Failed to copy URL:', err);
    }
  };

  const toggleCatalog = (catalog: string) => {
    setSelectedCatalogs((prev) =>
      prev.includes(catalog) ? prev.filter((c) => c !== catalog) : [...prev, catalog]
    );
  };

  const isAnyLoading = isLoading || alertsLoading || scoringLoading || standardizationLoading || databaseLoading;
  const hasError = error || alertError || scoringError || standardizationError || databaseError;

  // Decorate the static tab rows with per-tab "has-result" indicators and a
  // qualified result state (clean / issues / warnings / complete) so a glance
  // at the tab bar communicates not just that an analysis ran, but how it
  // landed. Stays inside the warm-status spectrum per DESIGN.md.
  const tabRowsWithResults = useMemo(() => {
    const validationIssueCount = result?.issues?.length ?? 0;
    const validationCriticalCount = result?.issues?.filter(
      (i) => i.severity === 'critical' || i.severity === 'error',
    ).length ?? 0;
    const alertCount = alertResult?.alerts?.length ?? 0;

    const stateByTab: Partial<Record<TabType, { hasResult: boolean; resultState?: TabBarResultState }>> = {
      'validate': result
        ? {
            hasResult: true,
            resultState:
              validationCriticalCount > 0
                ? 'issues'
                : validationIssueCount > 0
                ? 'warnings'
                : 'clean',
          }
        : { hasResult: false },
      'deep-validation': result ? { hasResult: true, resultState: 'complete' } : { hasResult: false },
      'scoring-profiles': scoringResult ? { hasResult: true, resultState: 'complete' } : { hasResult: false },
      'alerts': (alertResult || safetyAssessResult)
        ? {
            hasResult: true,
            resultState: alertResult ? (alertCount > 0 ? 'issues' : 'clean') : 'complete',
          }
        : { hasResult: false },
      'compound-profile': profileResult ? { hasResult: true, resultState: 'complete' } : { hasResult: false },
      'database': databaseResults ? { hasResult: true, resultState: 'complete' } : { hasResult: false },
      'standardize': standardizationResult ? { hasResult: true, resultState: 'complete' } : { hasResult: false },
    };

    const decorate = (row: ReadonlyArray<TabBarTab<TabType>>) =>
      row.map((tab) => {
        const state = stateByTab[tab.id];
        return {
          ...tab,
          hasResult: !!state?.hasResult,
          resultState: state?.resultState,
        };
      });
    return [decorate(TAB_ROW_1), decorate(TAB_ROW_2)];
  }, [result, scoringResult, alertResult, safetyAssessResult, profileResult, databaseResults, standardizationResult]);

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

  const handleDownloadImage = useCallback(() => {
    const svgEl = previewRef.current?.querySelector('svg');
    if (!svgEl) return;

    const clone = svgEl.cloneNode(true) as SVGSVGElement;
    clone.style.background = '#ffffff';
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');

    const blob = new Blob([new XMLSerializer().serializeToString(clone)], { type: 'image/svg+xml' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `molecule-${(canonicalSmiles || molecule || 'structure').slice(0, 30)}.svg`;
    a.click();
    URL.revokeObjectURL(url);
  }, [canonicalSmiles, molecule]);

  // Auto-scroll to scoring results when they land. Without this the panel
  // appears full-width below the grid and users miss it because their eye
  // is still on the Score button in the left column. Only fires on the
  // null -> set transition so re-renders don't keep triggering.
  const prevScoringRef = useRef<typeof scoringResult>(null);
  useEffect(() => {
    if (!prevScoringRef.current && scoringResult && activeTab === 'validate') {
      scoringAnchorRef.current?.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }
    prevScoringRef.current = scoringResult;
  }, [scoringResult, activeTab]);

  const prevComparisonRef = useRef<typeof comparisonResult>(null);
  useEffect(() => {
    if (!prevComparisonRef.current && comparisonResult && activeTab === 'database') {
      comparisonAnchorRef.current?.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }
    prevComparisonRef.current = comparisonResult;
  }, [comparisonResult, activeTab]);

  // Loading phrase cycling — name the actual chemistry being computed.
  // Phrases rotate on a 1.5s interval while a loading flag is set, so a
  // user looking at a multi-second wait sees what stage the tool is in,
  // not just an opaque "Loading...".
  const [loadingPhase, setLoadingPhase] = useState(0);
  useEffect(() => {
    if (!isAnyLoading) {
      setLoadingPhase(0);
      return;
    }
    const id = setInterval(() => {
      setLoadingPhase((p) => p + 1);
    }, 1500);
    return () => clearInterval(id);
  }, [isAnyLoading]);

  // ── Lazy-load enrichment data when switching to profiler / safety tabs ──
  // Only auto-fires after validation has produced a result. Otherwise the
  // user gets a CTA inside the tab (already rendered) and decides whether
  // to spend the API call. Gating on `result` prevents unsolicited network
  // requests while the user is just exploring the tab bar.
  useEffect(() => {
    if (!result) return;
    if (activeTab === 'compound-profile' && !profileResult && !profileLoading && canonicalSmiles) {
      profileCompound(canonicalSmiles);
    }
    if (activeTab === 'alerts' && !safetyAssessResult && !safetyLoading && canonicalSmiles) {
      screenSafety(canonicalSmiles);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [activeTab, canonicalSmiles, result]);

  function getLoadingText(): string {
    let phrases: ReadonlyArray<string>;
    if (isLoading) {
      phrases = [
        'Sanitizing molecular graph...',
        'Enumerating stereocenters...',
        'Computing InChI layer...',
        'Evaluating valence model...',
      ];
    } else if (databaseLoading) {
      phrases = [
        'Querying PubChem by InChIKey...',
        'Resolving ChEMBL bioactivity links...',
        'Searching COCONUT natural products...',
      ];
    } else if (scoringLoading) {
      phrases = [
        'Computing Lipinski descriptors...',
        'Evaluating ML-readiness criteria...',
        'Calculating bioavailability radar...',
      ];
    } else if (standardizationLoading) {
      phrases = [
        'Applying ChEMBL structure pipeline...',
        'Canonicalizing tautomers...',
        'Removing salts and solvents...',
      ];
    } else {
      phrases = [
        'Screening for structural alerts...',
        'Matching SMARTS patterns...',
        'Checking PAINS and BRENK catalogs...',
      ];
    }
    return phrases[loadingPhase % phrases.length];
  }

  function renderErrorDetail() {
    const activeMessage = (error?.error || alertError?.error || scoringError?.error || standardizationError?.error || databaseError) as string | undefined;
    const errorType = classifyError(activeMessage, !!databaseError);
    return (
      <>
        <h3 className="font-semibold text-red-500 mb-1 font-display">
          {PARSE_ERROR_HEADINGS[errorType]}
        </h3>
        <p className="text-sm text-[var(--color-text-secondary)] break-words">
          {activeMessage}
        </p>
        <p className="mt-2 text-xs text-[var(--color-text-muted)]">
          {PARSE_ERROR_HINTS[errorType]}
        </p>
      </>
    );
  }

  return (
    <div ref={allChecksFloatContainerRef} className="relative max-w-7xl mx-auto space-y-6 px-4 sm:px-6">
      {/* Back to Batch bar — shown when navigated from batch results */}
      {batchOrigin && (
        <motion.div
          initial={{ opacity: 0, y: -12 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.3 }}
          className={cn(
            'sticky top-0 z-30 -mx-4 sm:-mx-6 px-4 sm:px-6 py-2.5',
            'bg-[var(--color-surface-elevated)]/90 backdrop-blur-md',
            'border-b border-[var(--color-primary)]/20',
            'shadow-[0_2px_12px_var(--glow-soft)]'
          )}
        >
          <div className="max-w-7xl mx-auto flex items-center gap-3">
            <button
              onClick={() => navigate('/batch')}
              className={cn(
                'flex items-center gap-2 px-3.5 py-1.5 rounded-lg text-sm font-medium',
                'bg-[var(--color-primary)]/10 text-[var(--color-primary)]',
                'border border-[var(--color-primary)]/20',
                'hover:bg-[var(--color-primary)]/20 hover:border-[var(--color-primary)]/40',
                'transition-all duration-200 cursor-pointer'
              )}
            >
              <ArrowLeft className="w-4 h-4" />
              Back to Batch
            </button>
            <span className="text-xs text-[var(--color-text-muted)]">
              Viewing full analysis
              {batchOrigin.moleculeName && (
                <> for <span className="font-medium text-[var(--color-text-secondary)]">{batchOrigin.moleculeName}</span></>
              )}
              {!batchOrigin.moleculeName && batchOrigin.moleculeIndex !== undefined && (
                <> for molecule <span className="font-medium text-[var(--color-text-secondary)]">#{batchOrigin.moleculeIndex + 1}</span></>
              )}
            </span>
          </div>
        </motion.div>
      )}

      {/* Header */}
      <motion.div
        className="text-center pt-4"
        initial={{ opacity: 0, y: -20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, ease: [0.25, 0.46, 0.45, 0.94] }}
      >
        <h1 className="text-3xl sm:text-4xl font-bold text-[var(--color-text-primary)] tracking-tight font-display">
          Molecule Validation
        </h1>
        <p className="text-[var(--color-text-secondary)] mt-3 text-base sm:text-lg max-w-2xl mx-auto leading-relaxed">
          Comprehensive validation, scoring, and standardization for chemical structures
        </p>
      </motion.div>

      {/* Example molecules + recent */}
      <ExampleMolecules
        onExampleClick={handleExampleClick}
        recent={recent}
        onSelectRecent={handleSelectRecent}
        onRemoveRecent={removeRecent}
        onClearRecent={clearRecent}
      />

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
                <p className="text-xs text-[var(--color-text-muted)] mt-0.5">
                  SMILES, InChI, compound name, ChEMBL ID, CAS, or{' '}
                  <a
                    href="https://app.naturalproducts.net/depict/structuredraw"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-[var(--color-primary)] hover:underline"
                  >
                    draw a structure
                  </a>
                </p>
              </div>
            </div>
            <StructureInput value={molecule} onChange={setMolecule} onSubmit={handleValidate} />

            {/* Input type detection badge */}
            <AnimatePresence>
              {molecule.trim() && detectedType !== 'ambiguous' && (
                <motion.div
                  initial={{ opacity: 0, height: 0 }}
                  animate={{ opacity: 1, height: 'auto' }}
                  exit={{ opacity: 0, height: 0 }}
                  transition={{ duration: 0.2 }}
                  className="mt-3 flex items-center gap-2"
                >
                  <div className={cn(
                    'inline-flex items-center gap-1.5 px-2.5 py-1 rounded-lg text-xs font-medium',
                    // Brand warm-spectrum tints per detected input type. Replaces
                    // the previous off-brand emerald/green/blue palettes.
                    isIdentifierInput && 'bg-[rgba(225,29,72,0.10)] text-[#be123c] dark:text-[#fb7185] border border-[rgba(225,29,72,0.22)]',
                    !isIdentifierInput && isIupacInput && 'bg-[rgba(217,119,6,0.10)] text-[#b45309] dark:text-[#fbbf24] border border-[rgba(217,119,6,0.22)]',
                    !isIdentifierInput && !isIupacInput && 'bg-[rgba(var(--color-primary-rgb),0.10)] text-[var(--color-primary-dark)] dark:text-[var(--color-primary)] border border-[rgba(var(--color-primary-rgb),0.22)]',
                  )}>
                    {isIdentifierInput
                      ? 'Detected as database identifier \u2014 will resolve on validate'
                      : isIupacInput
                        ? 'Detected as IUPAC name'
                        : 'Parsed as SMILES'}
                  </div>
                  {!isIdentifierInput && (
                    <>
                      <button
                        onClick={() => setForceInputType(isIupacInput ? 'smiles' : 'iupac')}
                        className="text-[10px] text-[var(--color-text-muted)] hover:text-[var(--color-primary)] underline"
                      >
                        Force {isIupacInput ? 'SMILES' : 'IUPAC'}
                      </button>
                      {forceInputType && (
                        <button
                          onClick={() => setForceInputType(null)}
                          className="text-[10px] text-[var(--color-text-muted)] hover:text-[var(--color-primary)]"
                        >
                          (auto)
                        </button>
                      )}
                    </>
                  )}
                </motion.div>
              )}
            </AnimatePresence>

            {/* Action button row — uniform default ClayButtons; colour lives in the
                leading icons so each action has a semantic hint while the chrome
                stays grouped. Reset = muted (secondary, undo-style action),
                Share + View Full Profile = primary (deliberate action / go deeper),
                Bookmark = accent amber (warm 'save' positive). */}
            <div className="mt-4 flex flex-wrap items-center gap-3">
              <ClayButton
                onClick={handleReset}
                disabled={!molecule && !result && !alertResult && !scoringResult && !standardizationResult && !databaseResults}
                leftIcon={<RotateCcw className="w-4 h-4 text-[var(--color-text-muted)]" />}
              >
                Reset
              </ClayButton>
              <ClayButton
                onClick={handleShare}
                disabled={!molecule.trim()}
                leftIcon={<Share2 className="w-4 h-4 text-[var(--color-primary)]" />}
              >
                Share
              </ClayButton>
              {/* Bookmark button - show after validation result */}
              {result && molecule.trim() && (
                <BookmarkButton
                  smiles={result.molecule_info.canonical_smiles || molecule.trim()}
                  source="single_validation"
                  onAfterBookmark={handleAfterBookmark}
                />
              )}
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
                      {/* Stats row — distinct icon + warm-spectrum tint per metric so the
                          five tiles read as five different things at a glance, not one
                          repeated card. Stays inside the brand palette per DESIGN.md. */}
                      <div className="grid grid-cols-2 sm:grid-cols-3 md:grid-cols-5 gap-2 text-center">
                        {result.molecule_info.num_atoms != null && (
                          <div className="rounded-lg p-2 bg-[rgba(var(--color-primary-rgb),0.06)] border border-[rgba(var(--color-primary-rgb),0.12)]">
                            <Atom className="w-3.5 h-3.5 mx-auto mb-1 text-[var(--color-primary)]" />
                            <div className="text-lg font-bold text-[var(--color-text-primary)] leading-none">{result.molecule_info.num_atoms}</div>
                            <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider mt-1">Atoms</div>
                          </div>
                        )}
                        {result.molecule_info.num_bonds != null && (
                          <div className="rounded-lg p-2 bg-[rgba(217,119,6,0.06)] border border-[rgba(217,119,6,0.14)]">
                            <GitMerge className="w-3.5 h-3.5 mx-auto mb-1 text-[#d97706] dark:text-[#fbbf24]" />
                            <div className="text-lg font-bold text-[var(--color-text-primary)] leading-none">{result.molecule_info.num_bonds}</div>
                            <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider mt-1">Bonds</div>
                          </div>
                        )}
                        {result.molecule_info.num_rings != null && (
                          <div className="rounded-lg p-2 bg-[rgba(225,29,72,0.06)] border border-[rgba(225,29,72,0.14)]">
                            <Hexagon className="w-3.5 h-3.5 mx-auto mb-1 text-[#e11d48] dark:text-[#fb7185]" />
                            <div className="text-lg font-bold text-[var(--color-text-primary)] leading-none">{result.molecule_info.num_rings}</div>
                            <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider mt-1">Rings</div>
                          </div>
                        )}
                        {result.molecule_info.num_aromatic_rings != null && result.molecule_info.num_aromatic_rings > 0 && (
                          <div className="rounded-lg p-2 bg-[rgba(251,191,36,0.10)] border border-[rgba(251,191,36,0.22)]">
                            <Sparkles className="w-3.5 h-3.5 mx-auto mb-1 text-[#b45309] dark:text-[#fbbf24]" />
                            <div className="text-lg font-bold text-[var(--color-text-primary)] leading-none">{result.molecule_info.num_aromatic_rings}</div>
                            <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider mt-1">Arom.</div>
                          </div>
                        )}
                        {result.molecule_info.molecular_weight && (
                          <div className="rounded-lg p-2 bg-[var(--color-surface-sunken)] border border-[var(--color-border)]">
                            <Scale className="w-3.5 h-3.5 mx-auto mb-1 text-[var(--color-text-secondary)]" />
                            <div className="text-lg font-bold text-[var(--color-text-primary)] leading-none">{result.molecule_info.molecular_weight.toFixed(1)}</div>
                            <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider mt-1">MW</div>
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
                        {resolverResult && resolverResult.resolved && (
                          <IdentifierResolverCard result={resolverResult} />
                        )}
                        {result.molecule_info.canonical_smiles && (
                          <div>
                            <div className="flex items-center justify-between mb-2">
                              <div className="flex items-center gap-1.5">
                                <span className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">
                                  {showKekulized && moleculeInfo?.kekulizedSmiles ? 'Kekulized SMILES' : 'Canonical SMILES'}
                                </span>
                                <InfoTooltip
                                  size="small"
                                  position="right"
                                  title={showKekulized ? 'Kekulized SMILES' : 'Canonical SMILES'}
                                  content={
                                    <div className="space-y-2 text-xs">
                                      <p>
                                        <span className="text-zinc-400">Generated by: </span>
                                        <span className="font-semibold text-white">
                                          {result.molecule_info.canonical_smiles_source || 'RDKit (Python)'}
                                        </span>
                                      </p>
                                      <p>
                                        <span className="text-zinc-400">Notation: </span>
                                        {showKekulized ? 'Kekulized (explicit double bonds)' : 'Canonical (aromatic notation)'}
                                      </p>
                                      <p className="text-zinc-300 leading-relaxed">
                                        Canonical SMILES provides a unique, deterministic string for each molecule.
                                        Different toolkits (RDKit, OpenBabel, CDK, OEChem) use different canonicalization
                                        algorithms and may produce different — but equally valid — canonical forms for the
                                        same molecule.
                                      </p>
                                      <div className="pt-1.5 border-t border-white/15 text-zinc-400 space-y-1">
                                        <p className="italic">
                                          RDKit: Open-Source Cheminformatics.{' '}
                                          <a href="https://www.rdkit.org" target="_blank" rel="noopener noreferrer"
                                            className="text-amber-400 underline underline-offset-2 hover:text-amber-300">
                                            rdkit.org
                                          </a>
                                        </p>
                                        <p>Weininger, D. <em>J. Chem. Inf. Comput. Sci.</em> 1989, 29, 97–101.</p>
                                      </div>
                                    </div>
                                  }
                                />
                                {result.input_interpretation?.detected_as === 'iupac' && (
                                  <span className="text-[10px] text-green-600 dark:text-green-400">(converted from IUPAC name)</span>
                                )}
                              </div>
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
                      <div className="grid grid-cols-3 gap-2 sm:gap-3 text-center">
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
                          <div className="flex items-center gap-1.5">
                            <span className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">
                              {showKekulized && moleculeInfo.kekulizedSmiles ? 'Kekulized SMILES' : 'Canonical SMILES'}
                            </span>
                            <InfoTooltip
                              size="small"
                              position="right"
                              title={showKekulized ? 'Kekulized SMILES' : 'Canonical SMILES'}
                              content={
                                <div className="space-y-2 text-xs">
                                  <p>
                                    <span className="text-zinc-400">Generated by: </span>
                                    <span className="font-semibold text-white">RDKit.js (client-side)</span>
                                  </p>
                                  <p>
                                    <span className="text-zinc-400">Notation: </span>
                                    {showKekulized ? 'Kekulized (explicit double bonds)' : 'Canonical (aromatic notation)'}
                                  </p>
                                  <p className="text-zinc-300 leading-relaxed">
                                    This is a client-side preview generated by RDKit.js (WebAssembly).
                                    After validation, the canonical SMILES will be recomputed by the
                                    server-side RDKit (Python) which may produce a different canonical form.
                                  </p>
                                  <div className="pt-1.5 border-t border-white/15 text-zinc-400 space-y-1">
                                    <p className="italic">
                                      RDKit: Open-Source Cheminformatics.{' '}
                                      <a href="https://www.rdkit.org" target="_blank" rel="noopener noreferrer"
                                        className="text-amber-400 underline underline-offset-2 hover:text-amber-300">
                                        rdkit.org
                                      </a>
                                    </p>
                                    <p>Weininger, D. <em>J. Chem. Inf. Comput. Sci.</em> 1989, 29, 97–101.</p>
                                  </div>
                                </div>
                              }
                            />
                          </div>
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
            <div className="p-3">
              <TabBar
                rows={tabRowsWithResults}
                activeTab={activeTab}
                onTabChange={setActiveTab}
                ariaLabel="Validation analyses"
              />
            </div>

            {/* Accent line below tab bar */}
            <div className="h-px bg-gradient-to-r from-transparent via-[var(--color-border)] to-transparent" />

            {/* Tab Content */}
            <AnimatePresence mode="wait">
              <motion.div
                key={activeTab}
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                exit={{ opacity: 0 }}
                transition={{ duration: 0.2, ease: 'easeInOut' }}
                className="p-5 sm:p-6"
              >
                {/* Validate & Score Tab */}
                {activeTab === 'validate' && (
                  <ValidateTabPanel
                    onValidate={handleValidate}
                    onScore={handleCalculateScores}
                    actionsDisabled={!molecule.trim() || isAnyLoading}
                    hasResult={Boolean(result)}
                    validateLoading={isLoading}
                    scoreLoading={scoringLoading}
                    showScoreHint={!result && Boolean(molecule.trim())}
                  />
                )}

                {/* Deep Validation Tab */}
                {activeTab === 'deep-validation' && (
                  <DeepValidationTabPanel
                    onValidate={handleValidate}
                    actionsDisabled={!molecule.trim() || isAnyLoading}
                    validateLoading={isLoading}
                    checks={result?.all_checks ?? null}
                    onHighlightAtoms={setHighlightedAtoms}
                  />
                )}

                {/* Scoring Profiles Tab */}
                {activeTab === 'scoring-profiles' && (
                  <ScoringProfilesTabPanel
                    onValidate={handleValidate}
                    actionsDisabled={!molecule.trim() || isAnyLoading}
                    validateLoading={isLoading}
                    hasResult={Boolean(result)}
                    smiles={result ? resolvedSmiles : ''}
                  />
                )}

                {/* Database Lookup Tab */}
                {activeTab === 'database' && (
                  <DatabaseLookupControls
                    onLookup={handleDatabaseLookup}
                    lookupDisabled={!molecule.trim() || isAnyLoading}
                    databaseLoading={databaseLoading}
                    isComparing={isComparing}
                    autoCompare={autoCompare}
                    onToggleAutoCompare={() => setAutoCompare(!autoCompare)}
                    showManualCompare={!autoCompare && Boolean(databaseResults) && !comparisonResult}
                    onCompare={runComparison}
                  />
                )}

                {/* Compound Profile Tab */}
                {activeTab === 'compound-profile' && (
                  <CompoundProfileTabPanel
                    onValidate={handleValidate}
                    actionsDisabled={!molecule.trim() || isAnyLoading}
                    validateLoading={isLoading}
                    hasResult={Boolean(result)}
                    smiles={canonicalSmiles || ''}
                    profile={profileResult}
                    profileLoading={profileLoading}
                    profileError={profileError}
                  />
                )}

                {/* Safety Tab (alerts screening + safety assessment) */}
                {activeTab === 'alerts' && (
                  <AlertsTabPanel
                    chemblCatalogs={CHEMBL_CATALOGS}
                    selectedCatalogs={selectedCatalogs}
                    onToggleCatalog={toggleCatalog}
                    chemblExpanded={chemblExpanded}
                    onToggleChemblExpanded={() => setChemblExpanded(!chemblExpanded)}
                    allChemblSelected={allChemblSelected}
                    someChemblSelected={someChemblSelected}
                    onToggleAllChembl={toggleAllChembl}
                    onScreenAlerts={handleScreenAlerts}
                    screenDisabled={!molecule.trim() || isAnyLoading || selectedCatalogs.length === 0}
                    alertsLoading={alertsLoading}
                    safetySmiles={canonicalSmiles || ''}
                    safetyAlertResult={safetyAlertResult}
                    safetyResult={safetyAssessResult}
                    safetyLoading={safetyLoading}
                    safetyError={safetyAlertError || safetyAssessError}
                  />
                )}

                {/* Standardize Tab */}
                {activeTab === 'standardize' && (
                  <StandardizeTabPanel
                    includeTautomer={includeTautomer}
                    onToggleTautomer={setIncludeTautomer}
                    onStandardize={handleStandardize}
                    standardizeDisabled={!molecule.trim() || isAnyLoading}
                    standardizeLoading={standardizationLoading}
                  />
                )}
              </motion.div>
            </AnimatePresence>
          </div>

          {/* Database Results - Collapsible Box (default collapsed) */}
          <AnimatePresence>
            {databaseResults && activeTab === 'database' && (
              <DatabaseResultsBox
                results={databaseResults}
                expanded={dbResultsExpanded}
                onToggleExpanded={() => setDbResultsExpanded(!dbResultsExpanded)}
              />
            )}
          </AnimatePresence>

          {/* Loading State */}
          <LoadingPanel show={isAnyLoading} text={getLoadingText()} />

          {/* Error State */}
          <ErrorPanel
            show={Boolean(hasError) && !isAnyLoading}
            detail={renderErrorDetail()}
            onRetry={handleValidate}
            retryDisabled={!molecule.trim()}
            onDismiss={() => {
              reset();
              setAlertError(null);
              setScoringError(null);
              setStandardizationError(null);
              setDatabaseError(null);
            }}
          />
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
              <Suspense fallback={null}>
                <StandardizationResults result={standardizationResult.result} />
              </Suspense>
              <p className="mt-4 text-xs text-[var(--color-text-muted)] text-right">
                Completed in {standardizationResult.execution_time_ms}ms
              </p>
            </motion.div>
          ) : (
            <>
              {/* Molecule Viewer */}
              <MoleculeViewerPanel
                previewRef={previewRef}
                molecule={molecule}
                canonicalSmiles={canonicalSmiles}
                summary={result?.molecule_info ?? null}
                stereo={moleculeInfo ?? null}
                highlightedAtoms={highlightedAtoms}
                highlightLocked={highlightLocked}
                showCIP={showCIP}
                onToggleCIP={() => setShowCIP(!showCIP)}
                onDownload={handleDownloadImage}
              />

              {/* Validation Issues - Show right after molecule viewer (validate tab only) */}
              <ValidationIssuesPanel
                show={activeTab === 'validate' && Boolean(result) && validationIssues.length > 0}
                issues={validationIssues}
                executionMs={result?.execution_time_ms || 0}
                highlightedAtoms={highlightedAtoms}
                highlightLocked={highlightLocked}
                onAtomHover={setHighlightedAtoms}
                onAtomLock={handleAtomLock}
              />

              {/* Score Tiles - Only show after validation/scoring (validate tab only) */}
              <ScoreTiles
                show={activeTab === 'validate' && Boolean(hasScores)}
                qualityScore={qualityScore}
                mlReadyScore={mlReadyScore}
                mlReady={mlReady}
                issues={result ? result.issues || [] : null}
                totalChecks={result?.all_checks?.length || 0}
              />

              {/* Validation outcome panels (All Clear + structural alerts) */}
              <ValidationOutcomePanels
                showSuccess={activeTab === 'validate' && Boolean(result) && validationIssues.length === 0}
                successExecutionMs={result?.execution_time_ms || 0}
                alertResult={alertResult}
                onAtomHover={setHighlightedAtoms}
              />

              {/* All Checks — COLLAPSED ANCHOR (right column).
                  Empty placeholder that reserves layout space where the
                  floating card visually sits when collapsed. Height
                  animates 62→0 (and back) in sync with the floating card's
                  position animation, so the page reflows smoothly. */}
              {activeTab === 'validate' && result && result.all_checks && result.all_checks.length > 0 && (
                <motion.div
                  ref={allChecksCollapsedAnchorRef}
                  animate={{
                    height:
                      allChecksPhase === 'collapsed' || allChecksPhase === 'collapsing'
                        ? 62
                        : 0,
                  }}
                  initial={false}
                  transition={{
                    duration: 0.35,
                    ease: [0.4, 0, 0.2, 1],
                    delay: allChecksPhase === 'collapsing' ? 0.4 : 0,
                  }}
                  aria-hidden="true"
                />
              )}

            </>
          )}
        </motion.div>
      </div>

      {/* All Checks — EXPANDED ANCHOR (full-width below grid).
          Empty placeholder reserving layout space below the grid for the
          expanded card. Height animates between 0 (collapsed) and the
          card's natural expanded height in sync with the floating card. */}
      {activeTab === 'validate' && result && result.all_checks && result.all_checks.length > 0 && (
        <motion.div
          ref={allChecksExpandedAnchorRef}
          animate={{
            height:
              allChecksPhase === 'expanded' || allChecksPhase === 'expanding'
                ? allChecksExpandedHeight
                : 0,
          }}
          initial={false}
          transition={{
            duration: 0.4,
            ease: [0.4, 0, 0.2, 1],
            delay: allChecksPhase === 'expanding' ? 0.3 : 0,
          }}
          aria-hidden="true"
        />
      )}

      {/* All Checks — FLOATING CARD (the only real card).
          Position absolute, animates top/left/width/height between the
          two anchors above. Sequenced so it MOVES first, then GROWS on
          expand, and SHRINKS first, then MOVES on collapse. Same DOM
          element throughout — no mount/unmount, no fade. */}
      {activeTab === 'validate' && result && result.all_checks && result.all_checks.length > 0 && allChecksBounds && (
        <AllChecksCard
          phase={allChecksPhase}
          bounds={allChecksBounds}
          expandedHeight={allChecksExpandedHeight}
          checks={result.all_checks}
          contentRef={allChecksContentRef}
          onToggle={toggleAllChecks}
        />
      )}

      {/* Cross-Database Comparison — full-width below the grid so the
          comparison surface gets the horizontal room it needs */}
      <AnimatePresence>
        {comparisonResult && activeTab === 'database' && (
          <motion.div
            ref={comparisonAnchorRef}
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -20 }}
          >
            <DatabaseComparisonPanel result={comparisonResult} />
          </motion.div>
        )}
      </AnimatePresence>

      {/* Scoring Results — full-width below the grid; auto-scrolls into
          view when scoring completes so the user is not stranded above
          their own result */}
      <AnimatePresence>
        {scoringResult && activeTab === 'validate' && (
          <motion.div
            ref={scoringAnchorRef}
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -20 }}
            className="card p-6 sm:p-8"
          >
            <Suspense fallback={null}>
              <ScoringResults scoringResponse={scoringResult} />
            </Suspense>
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
