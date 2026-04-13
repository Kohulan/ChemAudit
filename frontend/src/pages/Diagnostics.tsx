import { useState, useCallback, useEffect, useRef } from 'react';
import { useSearchParams } from 'react-router-dom';
import { motion, AnimatePresence } from 'framer-motion';
import {
  Stethoscope,
  ChevronDown,
} from 'lucide-react';
import { useDiagnostics } from '../hooks/useDiagnostics';
import { DiagnosticsInput } from '../components/diagnostics/DiagnosticsInput';
import { SMILESDiagnosticsPanel } from '../components/diagnostics/SMILESDiagnosticsPanel';
import { InChILayerDiffTable } from '../components/diagnostics/InChILayerDiffTable';
import { RoundTripPanel } from '../components/diagnostics/RoundTripPanel';
import { CrossPipelinePanel } from '../components/diagnostics/CrossPipelinePanel';
import { FilePreValidatorPanel } from '../components/diagnostics/FilePreValidatorPanel';
import { ClayCard } from '../components/ui/ClayCard';

/**
 * Collapsible section container for each diagnostic panel.
 * Collapsed by default; auto-expands when `defaultOpen` becomes true.
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

  // Auto-open when defaultOpen switches from false to true (after first submit)
  useEffect(() => {
    if (defaultOpen) setOpen(true);
  }, [defaultOpen]);

  return (
    <div className="border border-[var(--color-border)] rounded-2xl overflow-hidden">
      {/* Section header button */}
      <button
        type="button"
        onClick={() => !disabled && setOpen(v => !v)}
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
          <h2 className="text-lg font-semibold font-display text-[var(--color-text-primary)]">
            {title}
          </h2>
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

      {/* Collapsible body */}
      <AnimatePresence initial={false}>
        {open && (
          <motion.div
            key="content"
            initial={{ height: 0, opacity: 0 }}
            animate={{ height: 'auto', opacity: 1 }}
            exit={{ height: 0, opacity: 0 }}
            transition={{ duration: 0.5, ease: 'easeOut' }}
            className="overflow-hidden"
          >
            <div className="px-5 py-5 border-t border-[var(--color-border)]">{children}</div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}

/**
 * Empty state shown before any molecule is analyzed.
 */
function EmptyState() {
  return (
    <motion.div
      initial={{ opacity: 0, y: 8 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.3 }}
    >
      <ClayCard variant="flat" size="md" className="text-center py-12">
        <div className="flex flex-col items-center gap-3">
          <div className="w-12 h-12 rounded-2xl bg-[var(--color-surface-sunken)] flex items-center justify-center">
            <Stethoscope className="w-6 h-6 text-[var(--color-text-muted)]" />
          </div>
          <h2 className="text-base font-semibold font-display text-[var(--color-text-primary)]">
            No molecule analyzed yet
          </h2>
          <p className="text-sm text-[var(--color-text-secondary)] max-w-sm">
            Enter a SMILES, InChI, CAS number, or other identifier above to run structure
            diagnostics.
          </p>
        </div>
      </ClayCard>
    </motion.div>
  );
}

/** Utility: join truthy class names into a single string. */
function cn(...classes: (string | undefined | false | null)[]): string {
  return classes.filter(Boolean).join(' ');
}

/**
 * Diagnostics page — /diagnostics
 *
 * Structure Quality Diagnostics (Phase 09):
 * 1. Page heading
 * 2. DiagnosticsInput (molecule + file upload)
 * 3. SMILES Diagnostics panel (D-01 — Plan 04)
 * 4. InChI Layer Diff panel (D-02 — Plan 04)
 * 5. Round-Trip Check panel (D-03 — Plan 04)
 * 6. Cross-Pipeline Comparison panel (D-04 — Plan 05)
 * 7. File Pre-Validator panel (D-05 — Plan 05)
 *
 * URL query param ?smiles= pre-fills and auto-submits (D-02 cross-link).
 */
export function Diagnostics() {
  const [searchParams] = useSearchParams();
  const initialSmiles = searchParams.get('smiles') ?? undefined;

  const {
    smilesResult,
    inchiDiffResult,
    roundtripResult,
    crossPipelineResult,
    fileResult,
    smilesLoading,
    inchiDiffLoading,
    roundtripLoading,
    crossPipelineLoading,
    fileLoading,
    smilesError,
    inchiDiffError,
    roundtripError,
    crossPipelineError,
    fileError,
    currentSmiles,
    analyzeMolecule,
    compareInchi,
    checkRoundtrip,
    comparePipelines,
    prevalidateFile,
    clearFileResult,
  } = useDiagnostics();

  // Track whether any SMILES analysis has been submitted
  const [hasSubmitted, setHasSubmitted] = useState(false);

  // Track last action type so we can reorder results: show the relevant section first
  const [lastAction, setLastAction] = useState<'smiles' | 'file' | null>(null);

  // Ref for scrolling to results after processing completes
  const resultsRef = useRef<HTMLDivElement>(null);
  const pendingScrollRef = useRef(false);

  const handleAnalyze = useCallback(
    (smiles: string) => {
      setHasSubmitted(true);
      setLastAction('smiles');
      pendingScrollRef.current = true;
      analyzeMolecule(smiles);
    },
    [analyzeMolecule],
  );

  const handleFileUpload = useCallback(
    (file: File) => {
      setLastAction('file');
      pendingScrollRef.current = true;
      prevalidateFile(file);
    },
    [prevalidateFile],
  );

  const hasAnyMoleculeResult = !!(
    smilesResult || inchiDiffResult || roundtripResult || crossPipelineResult
  );
  const hasAnyResult = hasAnyMoleculeResult || !!fileResult;
  const showEmptyState = !hasSubmitted && !fileResult;

  // Scroll to results when a user-initiated action completes
  const anyLoading =
    smilesLoading || inchiDiffLoading || roundtripLoading || crossPipelineLoading || fileLoading;
  useEffect(() => {
    if (!pendingScrollRef.current) return;
    if (hasAnyResult && !anyLoading) {
      pendingScrollRef.current = false;
      requestAnimationFrame(() => {
        resultsRef.current?.scrollIntoView({ behavior: 'smooth', block: 'start' });
      });
    }
  }, [hasAnyResult, anyLoading]);

  const requiresMoleculeHint = 'Enter a molecule above to enable this section';

  // Single definition for the File Pre-Validator section, rendered at top or bottom
  // depending on which action the user performed last.
  const filePreValidatorSection = (
    <DiagnosticSection
      title="File Pre-Validator"
      subtitle="SDF block integrity and CSV structure checks"
      defaultOpen={!!fileResult || fileLoading || !!fileError}
    >
      <FilePreValidatorPanel
        result={fileResult}
        isLoading={fileLoading}
        error={fileError}
        onFileUpload={prevalidateFile}
        onClear={clearFileResult}
      />
    </DiagnosticSection>
  );

  return (
    <div className="max-w-[1200px] mx-auto px-4 pt-16 pb-16">
      {/* ── Page heading ── */}
      <div className="mb-6">
        <h1 className="text-2xl font-semibold font-display text-[var(--color-text-primary)]">
          Structure Quality Diagnostics
        </h1>
        <p className="text-sm text-[var(--color-text-secondary)] mt-1">
          Get actionable, position-specific feedback on malformed or ambiguous chemical structures
        </p>
      </div>

      {/* ── Input area ── */}
      <div className="mb-12">
        <DiagnosticsInput
          onAnalyze={handleAnalyze}
          onFileUpload={handleFileUpload}
          isLoading={smilesLoading}
          fileLoading={fileLoading}
          initialSmiles={initialSmiles}
        />
      </div>

      {/* ── Empty state ── */}
      {showEmptyState && (
        <div className="mb-12">
          <EmptyState />
        </div>
      )}

      {/* ── Diagnostic sections ── */}
      {/* The last action's section renders first so the user sees it immediately. */}
      {(hasSubmitted || hasAnyResult) && (
        <div ref={resultsRef} className="space-y-12 scroll-mt-4">

          {/* File Pre-Validator at top when file was the last action */}
          {lastAction === 'file' && filePreValidatorSection}

          {/* SMILES Diagnostics (D-01) */}
          <DiagnosticSection
            title="SMILES Diagnostics"
            subtitle="Position-specific parse errors with fix suggestions"
            defaultOpen={hasSubmitted}
          >
            <SMILESDiagnosticsPanel
              result={smilesResult}
              isLoading={smilesLoading}
              error={smilesError}
              originalSmiles={currentSmiles || ''}
              onFixApplied={(corrected) => analyzeMolecule(corrected)}
              onRetry={() => currentSmiles && analyzeMolecule(currentSmiles)}
            />
          </DiagnosticSection>

          {/* InChI Layer Diff (D-02) */}
          <DiagnosticSection
            title="InChI Layer Diff"
            subtitle="Compare two InChI strings layer-by-layer"
            defaultOpen={!!inchiDiffResult}
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

          {/* Round-Trip Check (D-03) */}
          <DiagnosticSection
            title="Round-Trip Lossiness"
            subtitle="Check SMILES → InChI → SMILES and SMILES → MOL → SMILES fidelity"
            defaultOpen={!!roundtripResult}
          >
            {!currentSmiles && !roundtripResult && !roundtripLoading && !roundtripError && (
              <p className="text-xs text-[var(--color-text-muted)]">{requiresMoleculeHint}</p>
            )}
            {(currentSmiles || roundtripResult || roundtripLoading || roundtripError) && (
              <RoundTripPanel
                result={roundtripResult}
                isLoading={roundtripLoading}
                error={roundtripError}
                currentSmiles={currentSmiles}
                onCheckRoundtrip={checkRoundtrip}
                onRetry={() => currentSmiles ? checkRoundtrip(currentSmiles) : undefined}
              />
            )}
          </DiagnosticSection>

          {/* Cross-Pipeline Comparison (D-04) */}
          <DiagnosticSection
            title="Cross-Pipeline Standardization"
            subtitle="Compare RDKit, ChEMBL, and minimal pipeline outputs"
            defaultOpen={!!crossPipelineResult}
          >
            {!currentSmiles && !crossPipelineResult && !crossPipelineLoading && !crossPipelineError && (
              <p className="text-xs text-[var(--color-text-muted)]">{requiresMoleculeHint}</p>
            )}
            {(currentSmiles || crossPipelineResult || crossPipelineLoading || crossPipelineError) && (
              <CrossPipelinePanel
                result={crossPipelineResult}
                isLoading={crossPipelineLoading}
                error={crossPipelineError}
                currentSmiles={currentSmiles}
                onComparePipelines={comparePipelines}
                onRetry={() => currentSmiles ? comparePipelines(currentSmiles) : undefined}
              />
            )}
          </DiagnosticSection>

          {/* File Pre-Validator at bottom (default position) */}
          {lastAction !== 'file' && filePreValidatorSection}
        </div>
      )}
    </div>
  );
}
