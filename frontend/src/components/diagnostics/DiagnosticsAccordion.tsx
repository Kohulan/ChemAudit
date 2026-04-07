import { useState, useCallback, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown, Stethoscope } from 'lucide-react';
import { SMILESDiagnosticsPanel } from './SMILESDiagnosticsPanel';
import { InChILayerDiffTable } from './InChILayerDiffTable';
import { RoundTripPanel } from './RoundTripPanel';
import { CrossPipelinePanel } from './CrossPipelinePanel';
import { FilePreValidatorPanel } from './FilePreValidatorPanel';
import { diagnosticsApi } from '../../services/api';
import type {
  CrossPipelineResponse,
  RoundTripResponse,
  SMILESDiagnosticsResponse,
  InChIDiffResponse,
  FilePreValidationResponse,
} from '../../types/diagnostics';
import { cn } from '../../lib/utils';

// ─── Props ───────────────────────────────────────────────────────────────────

interface DiagnosticsAccordionProps {
  /** Current molecule SMILES string. */
  smiles: string;
  /** Pre-fetched cross-pipeline result from parent (may be null). */
  crossPipelineResult: CrossPipelineResponse | null;
  /** Pre-fetched round-trip result from parent (may be null). */
  roundTripResult: RoundTripResponse | null;
  /** True while parent is performing initial data fetch. */
  isLoading: boolean;
}

// ─── Collapsible section wrapper ─────────────────────────────────────────────

/**
 * Lightweight collapsible section for each diagnostic tool.
 * Mirrors the DiagnosticSection pattern from the standalone Diagnostics page.
 */
function DiagnosticSection({
  title,
  subtitle,
  children,
  defaultOpen = false,
  disabled = false,
}: {
  title: string;
  subtitle?: string;
  children: React.ReactNode;
  defaultOpen?: boolean;
  disabled?: boolean;
}) {
  const [open, setOpen] = useState(defaultOpen);

  // Auto-open when defaultOpen transitions to true
  useEffect(() => {
    if (defaultOpen) setOpen(true);
  }, [defaultOpen]);

  return (
    <div className="border border-[var(--color-border)] rounded-2xl overflow-hidden">
      <button
        type="button"
        onClick={() => !disabled && setOpen((v) => !v)}
        disabled={disabled}
        aria-expanded={open}
        className={cn(
          'w-full flex items-center justify-between px-5 py-4',
          'text-left transition-colors duration-150',
          disabled
            ? 'cursor-not-allowed opacity-50'
            : 'cursor-pointer hover:bg-[var(--color-surface-sunken)]',
          open ? 'bg-[var(--color-surface-sunken)]' : 'bg-[var(--color-surface-elevated)]',
        )}
      >
        <div className="flex-1 min-w-0">
          <h3 className="text-base font-semibold font-display text-[var(--color-text-primary)]">
            {title}
          </h3>
          {subtitle && (
            <p className="text-xs text-[var(--color-text-muted)] mt-0.5">{subtitle}</p>
          )}
        </div>
        <motion.div
          animate={{ rotate: open ? 180 : 0 }}
          transition={{ duration: 0.25, ease: 'easeOut' }}
          className="shrink-0 ml-3"
        >
          <ChevronDown className="w-5 h-5 text-[var(--color-text-secondary)]" />
        </motion.div>
      </button>

      <AnimatePresence initial={false}>
        {open && (
          <motion.div
            key="content"
            initial={{ height: 0, opacity: 0 }}
            animate={{ height: 'auto', opacity: 1 }}
            exit={{ height: 0, opacity: 0 }}
            transition={{ duration: 0.3, ease: 'easeOut' }}
            className="overflow-hidden"
          >
            <div className="px-5 py-5 border-t border-[var(--color-border)]">{children}</div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}

// ─── Main component ──────────────────────────────────────────────────────────

/**
 * DiagnosticsAccordion — wraps all 5 diagnostic tools into independently
 * collapsible sections suitable for embedding inside a DrillDownSection
 * accordion in SingleValidation.
 *
 * Each tool manages its own loading/error state. Tools that receive
 * pre-fetched data (crossPipeline, roundTrip) show it immediately;
 * others compute on-demand when their section is expanded.
 */
export function DiagnosticsAccordion({
  smiles,
  crossPipelineResult: prefetchedCrossPipeline,
  roundTripResult: prefetchedRoundTrip,
  isLoading,
}: DiagnosticsAccordionProps) {
  // ── Per-tool state ───────────────────────────────────────────────────────

  // SMILES diagnostics (DIAG-01)
  const [smilesResult, setSmilesResult] = useState<SMILESDiagnosticsResponse | null>(null);
  const [smilesLoading, setSmilesLoading] = useState(false);
  const [smilesError, setSmilesError] = useState<string | null>(null);

  // InChI diff (DIAG-02)
  const [inchiDiffResult, setInchiDiffResult] = useState<InChIDiffResponse | null>(null);
  const [inchiDiffLoading, setInchiDiffLoading] = useState(false);
  const [inchiDiffError, setInchiDiffError] = useState<string | null>(null);

  // Round-trip (DIAG-03)
  const [roundtripResult, setRoundtripResult] = useState<RoundTripResponse | null>(
    prefetchedRoundTrip,
  );
  const [roundtripLoading, setRoundtripLoading] = useState(false);
  const [roundtripError, setRoundtripError] = useState<string | null>(null);

  // Cross-pipeline (DIAG-04)
  const [crossPipelineResult, setCrossPipelineResult] = useState<CrossPipelineResponse | null>(
    prefetchedCrossPipeline,
  );
  const [crossPipelineLoading, setCrossPipelineLoading] = useState(false);
  const [crossPipelineError, setCrossPipelineError] = useState<string | null>(null);

  // File pre-validator (DIAG-05)
  const [fileResult, setFileResult] = useState<FilePreValidationResponse | null>(null);
  const [fileLoading, setFileLoading] = useState(false);
  const [fileError, setFileError] = useState<string | null>(null);

  // Sync pre-fetched props when they change
  useEffect(() => {
    if (prefetchedRoundTrip) setRoundtripResult(prefetchedRoundTrip);
  }, [prefetchedRoundTrip]);

  useEffect(() => {
    if (prefetchedCrossPipeline) setCrossPipelineResult(prefetchedCrossPipeline);
  }, [prefetchedCrossPipeline]);

  // ── Action callbacks ─────────────────────────────────────────────────────

  const analyzeMolecule = useCallback(async (inputSmiles: string) => {
    setSmilesLoading(true);
    setSmilesError(null);
    try {
      const result = await diagnosticsApi.smiles(inputSmiles);
      setSmilesResult(result);
    } catch (e: unknown) {
      const err = e as { error?: string; detail?: string };
      setSmilesError(err?.error ?? err?.detail ?? 'SMILES diagnostics failed');
    } finally {
      setSmilesLoading(false);
    }
  }, []);

  const compareInchi = useCallback(async (inchiA: string, inchiB: string) => {
    setInchiDiffLoading(true);
    setInchiDiffError(null);
    try {
      const result = await diagnosticsApi.inchiDiff(inchiA, inchiB);
      setInchiDiffResult(result);
    } catch (e: unknown) {
      const err = e as { error?: string; detail?: string };
      setInchiDiffError(err?.error ?? err?.detail ?? 'InChI comparison failed');
    } finally {
      setInchiDiffLoading(false);
    }
  }, []);

  const checkRoundtrip = useCallback(
    async (inputSmiles: string, route: string = 'smiles_inchi_smiles') => {
      setRoundtripLoading(true);
      setRoundtripError(null);
      try {
        const result = await diagnosticsApi.roundtrip(inputSmiles, route);
        setRoundtripResult(result);
      } catch (e: unknown) {
        const err = e as { error?: string; detail?: string };
        setRoundtripError(err?.error ?? err?.detail ?? 'Round-trip check failed');
      } finally {
        setRoundtripLoading(false);
      }
    },
    [],
  );

  const comparePipelines = useCallback(async (molecule: string) => {
    setCrossPipelineLoading(true);
    setCrossPipelineError(null);
    try {
      const result = await diagnosticsApi.crossPipeline(molecule);
      setCrossPipelineResult(result);
    } catch (e: unknown) {
      const err = e as { error?: string; detail?: string };
      setCrossPipelineError(err?.error ?? err?.detail ?? 'Pipeline comparison failed');
    } finally {
      setCrossPipelineLoading(false);
    }
  }, []);

  const prevalidateFile = useCallback(async (file: File) => {
    setFileLoading(true);
    setFileError(null);
    try {
      const result = await diagnosticsApi.filePrevalidate(file);
      setFileResult(result);
    } catch (e: unknown) {
      const err = e as { error?: string; detail?: string };
      setFileError(err?.error ?? err?.detail ?? 'File pre-validation failed');
    } finally {
      setFileLoading(false);
    }
  }, []);

  const clearFileResult = useCallback(() => {
    setFileResult(null);
    setFileError(null);
  }, []);

  // ── Loading state ────────────────────────────────────────────────────────

  if (isLoading) {
    return (
      <div className="flex items-center gap-3 py-8 justify-center text-[var(--color-text-muted)]">
        <div className="w-5 h-5 border-2 border-current border-t-transparent rounded-full animate-spin" />
        <span className="text-sm">Loading diagnostics...</span>
      </div>
    );
  }

  // ── Empty state (no SMILES) ──────────────────────────────────────────────

  if (!smiles) {
    return (
      <div className="flex flex-col items-center gap-3 py-8 text-center">
        <div className="w-10 h-10 rounded-xl bg-[var(--color-surface-sunken)] flex items-center justify-center">
          <Stethoscope className="w-5 h-5 text-[var(--color-text-muted)]" />
        </div>
        <p className="text-sm text-[var(--color-text-secondary)]">
          Enter a SMILES string to run structure diagnostics.
        </p>
      </div>
    );
  }

  // ── Diagnostic sections ──────────────────────────────────────────────────

  return (
    <div className="space-y-4">
      {/* Section 1: SMILES Diagnostics (DIAG-01) — auto-expanded */}
      <DiagnosticSection
        title="SMILES Diagnostics"
        subtitle="Position-specific parse errors with fix suggestions"
        defaultOpen
      >
        <SMILESDiagnosticsPanel
          result={smilesResult}
          isLoading={smilesLoading}
          error={smilesError}
          originalSmiles={smiles}
          onFixApplied={(corrected) => analyzeMolecule(corrected)}
          onRetry={() => analyzeMolecule(smiles)}
        />
      </DiagnosticSection>

      {/* Section 2: InChI Layer Diff (DIAG-02) */}
      <DiagnosticSection
        title="InChI Layer Diff"
        subtitle="Compare two InChI strings layer-by-layer"
      >
        <InChILayerDiffTable
          result={inchiDiffResult}
          isLoading={inchiDiffLoading}
          error={inchiDiffError}
          onCompare={compareInchi}
          onRetry={() => compareInchi('', '')}
          initialInchiA={undefined}
        />
      </DiagnosticSection>

      {/* Section 3: Round-Trip Lossiness (DIAG-03) */}
      <DiagnosticSection
        title="Round-Trip Lossiness"
        subtitle="Check SMILES round-trip fidelity across formats"
        defaultOpen={!!roundtripResult}
      >
        <RoundTripPanel
          result={roundtripResult}
          isLoading={roundtripLoading}
          error={roundtripError}
          currentSmiles={smiles}
          onCheckRoundtrip={checkRoundtrip}
          onRetry={() => checkRoundtrip(smiles)}
        />
      </DiagnosticSection>

      {/* Section 4: Cross-Pipeline Standardization (DIAG-04) */}
      <DiagnosticSection
        title="Cross-Pipeline Standardization"
        subtitle="Compare RDKit, ChEMBL, and minimal pipeline outputs"
        defaultOpen={!!crossPipelineResult}
      >
        <CrossPipelinePanel
          result={crossPipelineResult}
          isLoading={crossPipelineLoading}
          error={crossPipelineError}
          currentSmiles={smiles}
          onComparePipelines={comparePipelines}
          onRetry={() => comparePipelines(smiles)}
        />
      </DiagnosticSection>

      {/* Section 5: File Pre-Validator (DIAG-05) */}
      <DiagnosticSection
        title="File Pre-Validator"
        subtitle="SDF block integrity and CSV structure checks"
      >
        <FilePreValidatorPanel
          result={fileResult}
          isLoading={fileLoading}
          error={fileError}
          onFileUpload={prevalidateFile}
          onClear={clearFileResult}
        />
      </DiagnosticSection>
    </div>
  );
}
