import { useState, useCallback } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { Filter } from 'lucide-react';
import { useStructureFilterConfig } from '../hooks/useStructureFilterConfig';
import { useStructureFilter } from '../hooks/useStructureFilter';
import { StructureFilterInput } from '../components/structure-filter/StructureFilterInput';
import { FunnelChart } from '../components/structure-filter/FunnelChart';
import { StructureFilterConfigPanel } from '../components/structure-filter/StructureFilterConfigPanel';
import { StructureFilterResultsTable } from '../components/structure-filter/StructureFilterResultsTable';
import { StructureFilterDownloadPanel } from '../components/structure-filter/StructureFilterDownloadPanel';
import { StructureFilterProgressBar } from '../components/structure-filter/StructureFilterProgressBar';

/**
 * Structure Filter page — multi-stage funnel pipeline for generative model output.
 *
 * Hooks wired at page top level so Plan 04/05 can add components without
 * prop drilling restructuring (Phase 10 pattern).
 *
 * Layout: D-04 grid (60/40 split: FunnelChart + Config, full-width Results table)
 *
 * State machine rendering (per UI-SPEC):
 * - idle: empty state with filter icon prompt
 * - loading-sync: skeleton pulse in funnel area
 * - loading-async | processing: StructureFilterProgressBar with WebSocket progress
 * - success: FunnelChart + StructureFilterResultsTable + StructureFilterDownloadPanel
 * - error: error card with "Filtering failed" message and retry link
 */
export default function StructureFilterPage() {
  const configHook = useStructureFilterConfig();
  const filterHook = useStructureFilter();
  const [smilesInput, setSmilesInput] = useState('');

  // Derived: count of SMILES lines for async threshold display (used by StructureFilterProgressBar)
  const smilesList = smilesInput
    .split('\n')
    .map((s) => s.trim())
    .filter(Boolean);

  /**
   * Run filter — called by StructureFilterInput with a parsed SMILES list.
   * The input component handles CSV column extraction, SDF detection, etc.
   */
  const handleRunFilter = useCallback((parsedSmilesList: string[]) => {
    if (parsedSmilesList.length === 0) return;
    filterHook.runFilter(parsedSmilesList, configHook.config, configHook.activePreset || undefined);
  }, [configHook.config, configHook.activePreset, filterHook.runFilter]);

  /**
   * Run filter for SDF files — sends the raw file to the backend batch endpoint.
   */
  const handleRunFilterFile = useCallback((file: File) => {
    filterHook.runFilterFile(file, configHook.config, configHook.activePreset || undefined);
  }, [configHook.config, configHook.activePreset, filterHook.runFilterFile]);

  return (
    <div className="max-w-[1200px] mx-auto px-4 pt-16 pb-16">
      {/* Page heading per UI-SPEC */}
      <motion.div
        initial={{ opacity: 0, y: 10 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.3 }}
      >
        <h1 className="text-2xl font-semibold font-display mb-2">
          Structure Filter
        </h1>
        <p className="text-sm text-[var(--color-text-secondary)] mb-8">
          Filter chemical structures through a multi-stage validation funnel.
        </p>
      </motion.div>

      {/* ── Input area: textarea + file drop zone + Run Filter CTA ── */}
      <div className="mb-8">
        <StructureFilterInput
          smilesInput={smilesInput}
          onSmilesChange={setSmilesInput}
          onRunFilter={handleRunFilter}
          onRunFilterFile={handleRunFilterFile}
          isLoading={
            filterHook.state === 'loading-sync' ||
            filterHook.state === 'loading-async' ||
            filterHook.state === 'processing'
          }
        />
      </div>

      {/* Funnel + Config side-by-side — Plan 04 fills this (D-04 layout: 60/40) */}
      <div className="grid grid-cols-1 lg:grid-cols-5 gap-8 mb-8">
        <div className="lg:col-span-3">
          {/* ── Idle state ── */}
          {filterHook.state === 'idle' && (
            <div className="flex flex-col items-center justify-center py-16 text-center">
              <Filter className="w-12 h-12 text-[var(--color-text-secondary)] mb-4" />
              <h2 className="text-base font-semibold font-display mb-2">
                No molecules filtered yet
              </h2>
              <p className="text-sm text-[var(--color-text-secondary)] max-w-md">
                Paste a SMILES list or drop a file above, configure a preset, then click Run Filter
                to see per-stage funnel results.
              </p>
            </div>
          )}

          {/* ── Loading-sync skeleton ── */}
          {filterHook.state === 'loading-sync' && (
            <div className="space-y-3 py-4">
              {/* Skeleton placeholders animate-pulse for sync loading */}
              {Array.from({ length: 5 }).map((_, i) => (
                <div
                  key={i}
                  className="animate-pulse bg-[var(--color-surface-sunken)] rounded-lg"
                  style={{ height: `${48 - i * 4}px`, width: `${100 - i * 8}%` }}
                />
              ))}
              <p className="text-xs text-center text-[var(--color-text-secondary)] mt-2">
                {filterHook.currentStage
                  ? `Filtering: ${filterHook.currentStage}…`
                  : 'Filtering molecules…'}
              </p>
            </div>
          )}

          {/* ── Async progress (loading-async | processing) ── */}
          {(filterHook.state === 'loading-async' || filterHook.state === 'processing') && (
            <div className="flex flex-col items-center justify-center py-8 text-center gap-4">
              <motion.div
                animate={{ rotate: 360 }}
                transition={{ duration: 1.5, repeat: Infinity, ease: 'linear' }}
              >
                <Filter className="w-10 h-10 text-[var(--color-primary)]" />
              </motion.div>
              <p className="text-sm text-[var(--color-text-secondary)]">
                {filterHook.currentStage
                  ? `Filtering: ${filterHook.currentStage}…`
                  : 'Filtering molecules…'}
              </p>
            </div>
          )}

          {/* ── Error state ── */}
          {filterHook.state === 'error' && (
            <div className="flex flex-col items-center justify-center py-16 text-center">
              <p className="text-sm text-red-500">{filterHook.error ?? 'An error occurred'}</p>
              <button
                onClick={filterHook.reset}
                className="mt-3 text-xs text-[var(--color-text-secondary)] hover:text-[var(--color-text-primary)] underline"
              >
                Reset
              </button>
            </div>
          )}

          {/* ── Success: funnel chart + results summary ── */}
          {filterHook.state === 'success' && filterHook.result && (
            <div className="space-y-4">
              {/* Results summary — UI-SPEC copywriting */}
              <p className="text-sm text-[var(--color-text-secondary)]">
                <span className="font-semibold text-[var(--color-text-primary)]">
                  {filterHook.result.output_count.toLocaleString()} molecules
                </span>{' '}
                passed all stages (
                {filterHook.result.input_count > 0
                  ? Math.round(
                      (filterHook.result.output_count / filterHook.result.input_count) * 100,
                    )
                  : 0}
                % of input)
              </p>

              {/* FunnelChart */}
              <FunnelChart
                stages={filterHook.result.stages}
                inputCount={filterHook.result.input_count}
                selectedStage={filterHook.selectedStage}
                onStageClick={filterHook.setSelectedStage}
              />
            </div>
          )}
        </div>

        <div className="lg:col-span-2">
          {/* StructureFilterConfigPanel — wired in Plan 04 */}
          <StructureFilterConfigPanel
            config={configHook.config}
            activePreset={configHook.activePreset}
            modifiedFrom={configHook.modifiedFrom}
            profiles={configHook.profiles}
            onUpdateConfig={configHook.updateConfig}
            onSelectPreset={configHook.selectPreset}
            onSaveProfile={configHook.saveProfile}
            onDeleteProfile={configHook.deleteProfile}
          />
        </div>
      </div>

      {/* ── Progress bar for async batch ── */}
      {(filterHook.state === 'loading-async' || filterHook.state === 'processing') && (
        <div className="mb-6">
          <StructureFilterProgressBar
            progress={filterHook.progress}
            currentStage={filterHook.currentStage}
            totalMolecules={smilesList.length}
          />
        </div>
      )}

      {/* ── Download panel + results table on success ── */}
      {filterHook.state === 'success' && filterHook.result && (
        <AnimatePresence mode="wait">
          <motion.div
            key="results"
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            transition={{ duration: 0.3, ease: [0.4, 0, 0.2, 1] }}
          >
            <StructureFilterDownloadPanel
              molecules={filterHook.result.molecules}
              inputCount={filterHook.result.input_count}
              outputCount={filterHook.result.output_count}
            />
            <div className="mt-6">
              <StructureFilterResultsTable
                molecules={filterHook.result.molecules}
                selectedStage={filterHook.selectedStage}
                onClearFilter={() => filterHook.setSelectedStage(null)}
              />
            </div>
          </motion.div>
        </AnimatePresence>
      )}

      {/* Error is already shown inline above (lines 156-166) */}
    </div>
  );
}
