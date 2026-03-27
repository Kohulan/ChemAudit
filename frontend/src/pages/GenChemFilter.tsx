import { useState, useCallback } from 'react';
import { motion } from 'framer-motion';
import { Filter } from 'lucide-react';
import { useGenChemConfig } from '../hooks/useGenChemConfig';
import { useGenChemFilter } from '../hooks/useGenChemFilter';

/**
 * GenChem Filter page — multi-stage funnel pipeline for generative model output.
 *
 * Hooks wired at page top level so Plan 04/05 can add components without
 * prop drilling restructuring (Phase 10 pattern).
 *
 * Layout: D-04 grid (60/40 split: FunnelChart + Config, full-width Results table)
 */
export default function GenChemFilter() {
  const configHook = useGenChemConfig();
  const filterHook = useGenChemFilter();
  const [smilesInput, setSmilesInput] = useState('');
  const [_uploadedFile, setUploadedFile] = useState<File | null>(null);

  const handleRunFilter = useCallback(() => {
    const smilesList = smilesInput
      .split('\n')
      .map((s) => s.trim())
      .filter(Boolean);
    if (smilesList.length === 0) return;
    filterHook.runFilter(
      smilesList,
      configHook.config,
      configHook.activePreset || undefined,
    );
  }, [smilesInput, configHook.config, configHook.activePreset, filterHook.runFilter]);

  return (
    <div className="max-w-[1200px] mx-auto px-4 pt-16 pb-16">
      {/* Page heading per UI-SPEC */}
      <motion.div
        initial={{ opacity: 0, y: 10 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.3 }}
      >
        <h1 className="text-2xl font-semibold font-display mb-2">
          Generative Chemistry Filter
        </h1>
        <p className="text-sm text-[var(--color-text-secondary)] mb-8">
          Filter generative model output through a multi-stage validation funnel.
        </p>
      </motion.div>

      {/* GenChemInput placeholder — Plan 04 fills this */}
      <div className="mb-8">
        {/* Textarea + file drop zone + Run Filter CTA
            Placeholder: exposes smilesInput / setUploadedFile state for Plan 04 wiring */}
        <div className="card p-4">
          <textarea
            value={smilesInput}
            onChange={(e) => setSmilesInput(e.target.value)}
            placeholder="Paste SMILES (one per line) or drop a file below..."
            className="w-full h-32 bg-transparent text-sm font-mono resize-none focus:outline-none text-[var(--color-text-primary)] placeholder:text-[var(--color-text-secondary)]"
          />
          <div className="flex items-center justify-between mt-3 pt-3 border-t border-[var(--color-border)]">
            <label className="text-xs text-[var(--color-text-secondary)] cursor-pointer hover:text-[var(--color-text-primary)] transition-colors">
              <input
                type="file"
                accept=".txt,.smi,.csv"
                className="hidden"
                onChange={(e) => {
                  const file = e.target.files?.[0] || null;
                  setUploadedFile(file);
                  if (file) {
                    // Read file content and populate textarea
                    const reader = new FileReader();
                    reader.onload = (ev) => {
                      const text = ev.target?.result as string;
                      setSmilesInput(text || '');
                    };
                    reader.readAsText(file);
                  }
                }}
              />
              Upload file (.txt, .smi, .csv)
            </label>
            <button
              onClick={handleRunFilter}
              disabled={filterHook.state === 'loading-sync' || filterHook.state === 'loading-async' || !smilesInput.trim()}
              className="px-4 py-1.5 text-sm font-medium rounded-lg bg-[var(--color-primary)] text-white disabled:opacity-50 disabled:cursor-not-allowed hover:opacity-90 transition-opacity"
            >
              {filterHook.state === 'loading-sync' || filterHook.state === 'loading-async'
                ? 'Filtering...'
                : 'Run Filter'}
            </button>
          </div>
        </div>
      </div>

      {/* Funnel + Config side-by-side — Plan 04 fills this (D-04 layout: 60/40) */}
      <div className="grid grid-cols-1 lg:grid-cols-5 gap-8 mb-8">
        <div className="lg:col-span-3">
          {/* FunnelChart placeholder */}
          {filterHook.state === 'idle' && (
            <div className="flex flex-col items-center justify-center py-16 text-center">
              <Filter className="w-12 h-12 text-[var(--color-text-secondary)] mb-4" />
              <h2 className="text-base font-semibold font-display mb-2">
                No molecules filtered yet
              </h2>
              <p className="text-sm text-[var(--color-text-secondary)] max-w-md">
                Paste a SMILES list or drop a file above, configure a preset, then click
                Run Filter to see per-stage funnel results.
              </p>
            </div>
          )}
          {(filterHook.state === 'loading-sync' || filterHook.state === 'loading-async' || filterHook.state === 'processing') && (
            <div className="flex flex-col items-center justify-center py-16 text-center">
              <motion.div
                animate={{ rotate: 360 }}
                transition={{ duration: 1.5, repeat: Infinity, ease: 'linear' }}
              >
                <Filter className="w-10 h-10 text-[var(--color-primary)]" />
              </motion.div>
              <p className="mt-4 text-sm text-[var(--color-text-secondary)]">
                {filterHook.currentStage
                  ? `Filtering: ${filterHook.currentStage}…`
                  : 'Filtering molecules…'}
              </p>
              {filterHook.progress !== null && (
                <div className="mt-3 w-48 h-1.5 bg-[var(--color-surface-sunken)] rounded-full overflow-hidden">
                  <motion.div
                    className="h-full bg-[var(--color-primary)] rounded-full"
                    initial={{ width: 0 }}
                    animate={{ width: `${filterHook.progress}%` }}
                    transition={{ duration: 0.3 }}
                  />
                </div>
              )}
            </div>
          )}
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
          {filterHook.state === 'success' && filterHook.result && (
            <div className="card p-4">
              {/* FunnelChart will be wired in Plan 04 */}
              <p className="text-sm text-[var(--color-text-secondary)]">
                Filter complete: {filterHook.result.output_count} / {filterHook.result.input_count} molecules passed.
                Funnel chart wired in Plan 04.
              </p>
            </div>
          )}
        </div>

        <div className="lg:col-span-2">
          {/* GenChemConfigPanel placeholder — Plan 04 fills this */}
          <div className="card p-4">
            <p className="text-xs text-[var(--color-text-secondary)]">
              Config panel — active preset:{' '}
              <span className="font-medium text-[var(--color-text-primary)]">
                {configHook.activePreset ?? 'custom'}
              </span>
              . Full panel wired in Plan 04.
            </p>
          </div>
        </div>
      </div>

      {/* Results table + Download — Plan 05 fills this */}
      <div>
        {/* GenChemResultsTable + GenChemDownloadPanel placeholder */}
        {filterHook.state === 'success' && filterHook.result && (
          <div className="card p-4">
            <p className="text-sm text-[var(--color-text-secondary)]">
              Results table and download panel wired in Plan 05.
            </p>
          </div>
        )}
      </div>
    </div>
  );
}
