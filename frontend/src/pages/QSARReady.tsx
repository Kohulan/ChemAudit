import { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { Beaker } from 'lucide-react';
import { useQSARPipelineConfig } from '../hooks/useQSARPipelineConfig';
import { useQSARReady } from '../hooks/useQSARReady';
import { ClayCard } from '../components/ui/ClayCard';
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
 * Phase 10 QSAR-Ready Pipeline (Phase 10):
 * 1. Page heading + subtitle
 * 2. Pipeline configuration panel placeholder (Plan 04 will replace)
 * 3. Tab switcher: Single Molecule | Batch Upload
 * 4. Tab content areas (Plan 04 fills Single, Plan 05 fills Batch)
 *
 * Both hooks are wired here so Plan 04/05 can consume state as props
 * without restructuring the component tree.
 */
export function QSARReady() {
  const [activeTab, setActiveTab] = useState<QSARTab>('single');

  // Pipeline config hook — manages preset selection and custom profile persistence
  const pipelineConfig = useQSARPipelineConfig();

  // QSAR processing hook — manages single result, batch WebSocket, pagination
  const qsarState = useQSARReady();

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

      {/* ── Pipeline configuration panel (placeholder for Plan 04) ── */}
      <div id="config-panel-placeholder" className="mb-8">
        <ClayCard variant="default" size="md">
          <div className="flex items-start gap-4">
            <div className="w-10 h-10 rounded-xl bg-[var(--color-surface-sunken)] flex items-center justify-center shrink-0">
              <Beaker className="w-5 h-5 text-[var(--color-text-muted)]" />
            </div>
            <div className="flex-1 min-w-0">
              <h2 className="text-lg font-semibold font-display text-[var(--color-text-primary)]">
                Pipeline Configuration
              </h2>
              <p className="text-sm text-[var(--color-text-secondary)] mt-1">
                Active preset:{' '}
                <span className="font-medium text-[var(--color-text-primary)]">
                  {pipelineConfig.activePreset ?? 'Custom'}
                </span>
                {pipelineConfig.modifiedFrom && (
                  <span className="text-[var(--color-text-muted)]">
                    {' '}(modified from {pipelineConfig.modifiedFrom})
                  </span>
                )}
              </p>
              <p className="text-xs text-[var(--color-text-muted)] mt-2">
                Pipeline step toggles and preset selector will be added in Plan 04.
              </p>
            </div>
          </div>
        </ClayCard>
      </div>

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
          >
            <ClayCard variant="flat" size="md" className="text-center py-12">
              <div className="flex flex-col items-center gap-3">
                <div className="w-12 h-12 rounded-2xl bg-[var(--color-surface-sunken)] flex items-center justify-center">
                  <Beaker className="w-6 h-6 text-[var(--color-text-muted)]" />
                </div>
                <h2 className="text-base font-semibold font-display text-[var(--color-text-primary)]">
                  Single molecule content
                </h2>
                <p className="text-sm text-[var(--color-text-secondary)] max-w-sm">
                  QSARSingleInput + QSARSingleResult components will be added in Plan 04.
                </p>
                {/* Surface the hook state for Plan 04 to verify wiring */}
                {qsarState.singleLoading && (
                  <p className="text-xs text-[var(--color-text-muted)]">Processing…</p>
                )}
                {qsarState.singleError && (
                  <p className="text-xs text-status-error">{qsarState.singleError}</p>
                )}
              </div>
            </ClayCard>
          </motion.div>
        )}

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
            <ClayCard variant="flat" size="md" className="text-center py-12">
              <div className="flex flex-col items-center gap-3">
                <div className="w-12 h-12 rounded-2xl bg-[var(--color-surface-sunken)] flex items-center justify-center">
                  <Beaker className="w-6 h-6 text-[var(--color-text-muted)]" />
                </div>
                <h2 className="text-base font-semibold font-display text-[var(--color-text-primary)]">
                  Batch upload content
                </h2>
                <p className="text-sm text-[var(--color-text-secondary)] max-w-sm">
                  QSARBatchInput + results table components will be added in Plan 05.
                </p>
                {/* Surface the hook state for Plan 05 to verify wiring */}
                {qsarState.batchLoading && qsarState.batchStatus && (
                  <p className="text-xs text-[var(--color-text-muted)]">
                    Processing… {qsarState.batchStatus.progress}%
                  </p>
                )}
                {qsarState.batchError && (
                  <p className="text-xs text-status-error">{qsarState.batchError}</p>
                )}
              </div>
            </ClayCard>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}
