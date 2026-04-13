import { useState, useRef, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { Beaker, RefreshCw } from 'lucide-react';
import { useQSARPipelineConfig } from '../hooks/useQSARPipelineConfig';
import { useQSARReady } from '../hooks/useQSARReady';
import { ClayCard } from '../components/ui/ClayCard';
import { ClayButton } from '../components/ui/ClayButton';
import { Skeleton } from '../components/ui/Skeleton';
import { StructureInput } from '../components/molecules/StructureInput';
import { PipelineConfigPanel } from '../components/qsar-ready/PipelineConfigPanel';
import { QSARSingleResult } from '../components/qsar-ready/QSARSingleResult';
import { QSARBatchInput } from '../components/qsar-ready/QSARBatchInput';
import { QSARBatchSummary } from '../components/qsar-ready/QSARBatchSummary';
import { QSARBatchTable } from '../components/qsar-ready/QSARBatchTable';
import { QSARDownloadPanel } from '../components/qsar-ready/QSARDownloadPanel';
import { QSARProgressBar } from '../components/qsar-ready/QSARProgressBar';
import { cn } from '../lib/utils';

// =============================================================================
// Tab types
// =============================================================================

type QSARTab = 'single' | 'batch';

// =============================================================================
// Page component
// =============================================================================

/**
 * QSAR-Ready Pipeline page — /qsar-ready
 *
 * Phase 10 Plan 04:
 * 1. Page heading + subtitle
 * 2. PipelineConfigPanel (shared above tabs — preset selector + step toggles)
 * 3. Tab switcher: Single Molecule | Batch Upload
 * 4. Single tab: StructureInput + Run Pipeline CTA + result states
 * 5. Batch tab: placeholder (Plan 05)
 *
 * Both hooks are wired here so Plans 04/05 can consume state without
 * prop drilling restructuring.
 */
export function QSARReady() {
  const [activeTab, setActiveTab] = useState<QSARTab>('single');
  const [singleInput, setSingleInput] = useState('');

  // Pipeline config hook — manages preset selection and custom profile persistence
  const pipelineConfig = useQSARPipelineConfig();

  // QSAR processing hook — manages single result, batch WebSocket, pagination
  const qsarState = useQSARReady();

  // Track whether config changed after results were produced
  const lastRunConfigRef = useRef<string | null>(null);
  const [configChanged, setConfigChanged] = useState(false);

  // Track last batch input for re-run
  const lastBatchInputRef = useRef<{ file: File | null; smilesText: string | null }>({ file: null, smilesText: null });

  // Snapshot config when a run completes (single or batch)
  useEffect(() => {
    if (qsarState.singleResult || qsarState.batchResults) {
      if (lastRunConfigRef.current === null) {
        lastRunConfigRef.current = JSON.stringify(pipelineConfig.config);
      }
    }
  }, [qsarState.singleResult, qsarState.batchResults, pipelineConfig.config]);

  // Detect config changes after results exist
  useEffect(() => {
    if (lastRunConfigRef.current === null) return;
    const current = JSON.stringify(pipelineConfig.config);
    setConfigChanged(current !== lastRunConfigRef.current);
  }, [pipelineConfig.config]);

  const isRunning = qsarState.singleLoading || qsarState.batchLoading;

  // Handle single molecule pipeline run
  const handleRunPipeline = () => {
    const trimmed = singleInput.trim();
    if (!trimmed) return;
    lastRunConfigRef.current = JSON.stringify(pipelineConfig.config);
    setConfigChanged(false);
    void qsarState.runSingle(trimmed, pipelineConfig.config);
  };

  // Handle re-run with updated config
  const handleRerun = () => {
    lastRunConfigRef.current = JSON.stringify(pipelineConfig.config);
    setConfigChanged(false);
    if (activeTab === 'single') {
      const trimmed = singleInput.trim();
      if (trimmed) void qsarState.runSingle(trimmed, pipelineConfig.config);
    } else {
      const { file, smilesText } = lastBatchInputRef.current;
      if (file || smilesText) void qsarState.runBatch(file, smilesText, pipelineConfig.config);
    }
  };

  // Wrap runBatch to capture the last batch input
  const handleRunBatch = (file: File | null, smilesText: string | null, config: typeof pipelineConfig.config) => {
    lastBatchInputRef.current = { file, smilesText };
    lastRunConfigRef.current = JSON.stringify(config);
    setConfigChanged(false);
    void qsarState.runBatch(file, smilesText, config);
  };

  return (
    <div className="max-w-[1200px] mx-auto px-4 pt-16 pb-16">

      {/* ── Page heading ── */}
      <div className="mb-6">
        <h1 className="text-2xl font-semibold font-display text-[var(--color-text-primary)]">
          QSAR-Ready Pipeline
        </h1>
        <p className="text-sm text-[var(--color-text-secondary)] mt-1">
          Configure and run a 10-step curation pipeline to produce standardized, ML-ready structures
        </p>
      </div>

      {/* ── Pipeline configuration panel ── */}
      <div className="mb-8">
        <PipelineConfigPanel {...pipelineConfig} />
      </div>

      {/* ── Re-run banner when config changes after results exist ── */}
      <AnimatePresence>
        {configChanged && (
          <motion.div
            initial={{ opacity: 0, y: -10, height: 0 }}
            animate={{ opacity: 1, y: 0, height: 'auto' }}
            exit={{ opacity: 0, y: -10, height: 0 }}
            transition={{ duration: 0.25 }}
            className="mb-6 overflow-hidden"
          >
            <div className="flex items-center justify-between gap-4 px-5 py-3.5 rounded-xl bg-amber-500/10 border border-amber-500/20">
              <div className="flex items-center gap-2.5">
                <RefreshCw className="w-4 h-4 text-amber-600 dark:text-amber-400" />
                <p className="text-sm text-amber-700 dark:text-amber-300">
                  Pipeline configuration changed since last run
                </p>
              </div>
              <ClayButton
                variant="primary"
                size="sm"
                onClick={handleRerun}
                disabled={isRunning}
                leftIcon={<RefreshCw className="w-3.5 h-3.5" />}
              >
                Re-run with new settings
              </ClayButton>
            </div>
          </motion.div>
        )}
      </AnimatePresence>

      {/* ── Tab switcher ── */}
      <div
        role="tablist"
        aria-label="QSAR-Ready processing mode"
        className="flex border-b border-[var(--color-border)] mb-6"
      >
        <button
          type="button"
          role="tab"
          aria-selected={activeTab === 'single'}
          aria-controls="tabpanel-single"
          id="tab-single"
          onClick={() => setActiveTab('single')}
          className={cn(
            'px-5 py-3 text-sm cursor-pointer transition-colors duration-150',
            activeTab === 'single'
              ? 'border-b-2 border-[var(--color-primary)] text-[var(--color-primary)] font-semibold'
              : 'text-[var(--color-text-secondary)] hover:text-[var(--color-text-primary)] border-b-2 border-transparent',
          )}
        >
          Single Molecule
        </button>
        <button
          type="button"
          role="tab"
          aria-selected={activeTab === 'batch'}
          aria-controls="tabpanel-batch"
          id="tab-batch"
          onClick={() => setActiveTab('batch')}
          className={cn(
            'px-5 py-3 text-sm cursor-pointer transition-colors duration-150',
            activeTab === 'batch'
              ? 'border-b-2 border-[var(--color-primary)] text-[var(--color-primary)] font-semibold'
              : 'text-[var(--color-text-secondary)] hover:text-[var(--color-text-primary)] border-b-2 border-transparent',
          )}
        >
          Batch Upload
        </button>
      </div>

      {/* ── Tab content area ── */}
      <AnimatePresence mode="wait">
        {/* ── Single Molecule tab ── */}
        {activeTab === 'single' && (
          <motion.div
            key="single"
            role="tabpanel"
            id="tabpanel-single"
            aria-labelledby="tab-single"
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
            transition={{ duration: 0.2, ease: 'easeOut' }}
            className="space-y-6"
          >
            {/* Input area */}
            <div className="space-y-3">
              <StructureInput
                value={singleInput}
                onChange={setSingleInput}
                onSubmit={handleRunPipeline}
                disabled={qsarState.singleLoading}
                placeholder="Enter a SMILES, InChI, CAS number, or other identifier above, then click Run Pipeline."
              />
              <div className="flex justify-end">
                <ClayButton
                  variant="primary"
                  size="md"
                  onClick={handleRunPipeline}
                  disabled={!singleInput.trim() || qsarState.singleLoading}
                  loading={qsarState.singleLoading}
                >
                  Run Pipeline
                </ClayButton>
              </div>
            </div>

            {/* ── Loading skeleton ── */}
            {qsarState.singleLoading && (
              <div className="space-y-4">
                <div className="grid grid-cols-2 gap-6">
                  <ClayCard variant="flat" size="sm" className="p-4 space-y-3">
                    <Skeleton variant="text" height={12} width="40%" />
                    <Skeleton variant="rectangular" height={200} className="rounded-lg" />
                    <Skeleton variant="text" height={12} width="90%" />
                    <Skeleton variant="text" height={12} width="70%" />
                  </ClayCard>
                  <ClayCard variant="flat" size="sm" className="p-4 space-y-3">
                    <Skeleton variant="text" height={12} width="40%" />
                    <Skeleton variant="rectangular" height={200} className="rounded-lg" />
                    <Skeleton variant="text" height={12} width="90%" />
                    <Skeleton variant="text" height={12} width="70%" />
                  </ClayCard>
                </div>
                <div className="space-y-2">
                  {Array.from({ length: 5 }, (_, i) => (
                    <Skeleton key={i} variant="text" height={40} className="rounded-xl" />
                  ))}
                </div>
              </div>
            )}

            {/* ── Error state ── */}
            {!qsarState.singleLoading && qsarState.singleError && (
              <ClayCard variant="flat" size="sm" className="border border-red-200 bg-red-50 dark:bg-red-900/10 dark:border-red-900/30">
                <p className="text-sm text-red-700 dark:text-red-400">
                  Pipeline failed. {qsarState.singleError}. Check your input structure and try again.
                </p>
              </ClayCard>
            )}

            {/* ── Result ── */}
            {!qsarState.singleLoading && !qsarState.singleError && qsarState.singleResult && (
              <QSARSingleResult result={qsarState.singleResult} />
            )}

            {/* ── Empty state ── */}
            {!qsarState.singleLoading && !qsarState.singleError && !qsarState.singleResult && (
              <ClayCard variant="flat" size="md" className="text-center py-12">
                <div className="flex flex-col items-center gap-3">
                  <div className="w-12 h-12 rounded-2xl bg-[var(--color-surface-sunken)] flex items-center justify-center">
                    <Beaker className="w-6 h-6 text-[var(--color-text-muted)]" />
                  </div>
                  <h2 className="text-base font-semibold font-display text-[var(--color-text-primary)]">
                    No molecule processed yet
                  </h2>
                  <p className="text-sm text-[var(--color-text-secondary)] max-w-sm">
                    Enter a SMILES, InChI, CAS number, or other identifier above, then click Run Pipeline.
                  </p>
                </div>
              </ClayCard>
            )}
          </motion.div>
        )}

        {/* ── Batch Upload tab ── */}
        {activeTab === 'batch' && (
          <motion.div
            key="batch"
            role="tabpanel"
            id="tabpanel-batch"
            aria-labelledby="tab-batch"
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
            transition={{ duration: 0.2, ease: 'easeOut' }}
          >
            <div className="space-y-6">
              {/* ── State: no job yet — empty state + input ── */}
              {!qsarState.batchJobId && !qsarState.batchError && (
                <>
                  <div className="text-center py-4">
                    <h2 className="text-base font-semibold font-display text-[var(--color-text-primary)] mb-1">
                      No batch results yet
                    </h2>
                    <p className="text-sm text-[var(--color-text-secondary)]">
                      Paste SMILES below or drop a CSV or SDF file to begin batch curation.
                    </p>
                  </div>
                  <QSARBatchInput
                    onSubmit={handleRunBatch}
                    config={pipelineConfig.config}
                    loading={qsarState.batchLoading}
                  />
                </>
              )}

              {/* ── State: processing — progress bar ── */}
              {qsarState.batchJobId &&
                qsarState.batchStatus !== null &&
                (qsarState.batchStatus.status === 'pending' ||
                  qsarState.batchStatus.status === 'processing') && (
                <QSARProgressBar status={qsarState.batchStatus} />
              )}

              {/* ── State: results loaded ── */}
              {qsarState.batchResults !== null && (
                <>
                  <QSARBatchSummary summary={qsarState.batchResults.summary} />

                  <QSARDownloadPanel
                    onDownload={qsarState.downloadBatch}
                    disabled={
                      qsarState.batchStatus !== null &&
                      (qsarState.batchStatus.status === 'pending' ||
                        qsarState.batchStatus.status === 'processing')
                    }
                  />

                  <QSARBatchTable
                    results={qsarState.batchResults.results}
                    page={qsarState.batchPage}
                    totalPages={qsarState.batchResults.total_pages}
                    onPageChange={qsarState.fetchBatchPage}
                  />
                </>
              )}

              {/* ── Error state ── */}
              {qsarState.batchError && (
                <div className="rounded-lg border border-red-200 bg-red-50 dark:bg-red-900/10 dark:border-red-900/30 px-4 py-4">
                  <p className="text-sm text-red-700 dark:text-red-400">
                    Batch processing failed. {qsarState.batchError}. Re-upload the file to retry.
                  </p>
                  <button
                    type="button"
                    className="mt-3 text-xs text-red-600 underline"
                    onClick={() => qsarState.clearBatch()}
                  >
                    Clear and try again
                  </button>
                </div>
              )}
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}
