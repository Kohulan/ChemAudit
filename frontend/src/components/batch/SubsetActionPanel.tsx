/**
 * Slide-over panel for batch subset actions: re-validate, re-score with
 * profile picker, and export selected molecules.
 */

import { useState, useEffect, useCallback } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { X, RefreshCw, BarChart3, Download, ChevronDown, ExternalLink } from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';
import { Badge } from '../ui/Badge';
import { subsetApi, profilesApi } from '../../services/api';
import { EXPORT_FORMATS } from '../../types/workflow';
import { cn } from '../../lib/utils';
import type { ScoringProfile, ExportFormat } from '../../types/workflow';

interface SubsetActionPanelProps {
  jobId: string;
  selectedIndices: Set<number>;
  isOpen: boolean;
  onClose: () => void;
  onNewJob?: (jobId: string) => void;
}

export function SubsetActionPanel({
  jobId,
  selectedIndices,
  isOpen,
  onClose,
  onNewJob,
}: SubsetActionPanelProps) {
  const [profiles, setProfiles] = useState<ScoringProfile[]>([]);
  const [selectedProfileId, setSelectedProfileId] = useState<number | undefined>();
  const [exportFormat, setExportFormat] = useState<ExportFormat>('csv');
  const [isRevalidating, setIsRevalidating] = useState(false);
  const [isRescoring, setIsRescoring] = useState(false);
  const [isExporting, setIsExporting] = useState(false);
  const [result, setResult] = useState<{ action: string; jobId?: string; message: string } | null>(null);

  // Load profiles for the re-score picker
  useEffect(() => {
    if (isOpen) {
      profilesApi.getProfiles().then(setProfiles).catch(() => {});
    }
  }, [isOpen]);

  // Close on Escape
  useEffect(() => {
    const handleKey = (e: KeyboardEvent) => {
      if (e.key === 'Escape') onClose();
    };
    if (isOpen) {
      document.addEventListener('keydown', handleKey);
      return () => document.removeEventListener('keydown', handleKey);
    }
  }, [isOpen, onClose]);

  const indices = Array.from(selectedIndices);

  const handleRevalidate = useCallback(async () => {
    setIsRevalidating(true);
    setResult(null);
    try {
      const res = await subsetApi.revalidateSubset(jobId, indices);
      setResult({
        action: 'revalidate',
        jobId: res.new_job_id,
        message: `Re-validation started for ${res.molecule_count} molecules`,
      });
      onNewJob?.(res.new_job_id);
    } catch (err) {
      setResult({
        action: 'revalidate',
        message: `Failed: ${err instanceof Error ? err.message : 'Unknown error'}`,
      });
    } finally {
      setIsRevalidating(false);
    }
  }, [jobId, indices, onNewJob]);

  const handleRescore = useCallback(async () => {
    setIsRescoring(true);
    setResult(null);
    try {
      const res = await subsetApi.rescoreSubset(jobId, indices, selectedProfileId);
      setResult({
        action: 'rescore',
        jobId: res.new_job_id,
        message: `Re-scoring started for ${res.molecule_count} molecules`,
      });
      onNewJob?.(res.new_job_id);
    } catch (err) {
      setResult({
        action: 'rescore',
        message: `Failed: ${err instanceof Error ? err.message : 'Unknown error'}`,
      });
    } finally {
      setIsRescoring(false);
    }
  }, [jobId, indices, selectedProfileId, onNewJob]);

  const handleExport = useCallback(async () => {
    setIsExporting(true);
    setResult(null);
    try {
      const blob = await subsetApi.exportSubset(jobId, indices, exportFormat);
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `subset_${jobId.slice(0, 8)}.${EXPORT_FORMATS.find((f) => f.value === exportFormat)?.extension ?? 'dat'}`;
      a.click();
      URL.revokeObjectURL(url);
      setResult({ action: 'export', message: 'Download started' });
    } catch (err) {
      setResult({
        action: 'export',
        message: `Failed: ${err instanceof Error ? err.message : 'Unknown error'}`,
      });
    } finally {
      setIsExporting(false);
    }
  }, [jobId, indices, exportFormat]);

  if (!isOpen) return null;

  return (
    <div className="fixed inset-0 z-50 flex justify-end" onClick={onClose}>
      {/* Backdrop */}
      <motion.div
        initial={{ opacity: 0 }}
        animate={{ opacity: 1 }}
        exit={{ opacity: 0 }}
        className="absolute inset-0 bg-black/30"
      />

      {/* Panel */}
      <motion.div
        initial={{ x: '100%' }}
        animate={{ x: 0 }}
        exit={{ x: '100%' }}
        transition={{ type: 'spring', damping: 30, stiffness: 300 }}
        className="relative w-full max-w-md bg-[var(--color-surface-elevated)] border-l border-[var(--color-border)] shadow-2xl flex flex-col h-full"
        onClick={(e) => e.stopPropagation()}
      >
        {/* Header */}
        <div className="flex items-center justify-between px-5 py-4 border-b border-[var(--color-border)]">
          <div>
            <h3 className="text-lg font-semibold text-[var(--color-text-primary)] font-display">
              Subset Actions
            </h3>
            <p className="text-xs text-[var(--color-text-muted)] mt-0.5">
              {selectedIndices.size} molecules selected
            </p>
          </div>
          <button
            onClick={onClose}
            className="p-2 rounded-lg hover:bg-[var(--color-surface-sunken)] transition-colors"
          >
            <X className="w-5 h-5 text-[var(--color-text-muted)]" />
          </button>
        </div>

        {/* Content */}
        <div className="flex-1 overflow-y-auto p-5 space-y-5">
          {/* Selected molecules preview */}
          <div className="bg-[var(--color-surface-sunken)] rounded-xl p-3">
            <div className="text-xs text-[var(--color-text-muted)] mb-1">Selected indices:</div>
            <div className="flex flex-wrap gap-1.5">
              {indices.slice(0, 20).map((idx) => (
                <Badge key={idx} variant="default" size="sm">#{idx}</Badge>
              ))}
              {indices.length > 20 && (
                <Badge variant="default" size="sm">+{indices.length - 20} more</Badge>
              )}
            </div>
          </div>

          {/* Re-validate */}
          <div className="border border-[var(--color-border)] rounded-xl p-4">
            <div className="flex items-center gap-2 mb-2">
              <RefreshCw className="w-4 h-4 text-[var(--color-primary)]" />
              <h4 className="text-sm font-medium text-[var(--color-text-primary)]">Re-validate Selection</h4>
            </div>
            <p className="text-xs text-[var(--color-text-muted)] mb-3">
              Create a new batch job with only the selected molecules.
            </p>
            <ClayButton
              variant="primary"
              size="sm"
              onClick={handleRevalidate}
              disabled={isRevalidating}
              loading={isRevalidating}
              leftIcon={<RefreshCw className="w-3.5 h-3.5" />}
            >
              Re-validate {selectedIndices.size} Molecules
            </ClayButton>
          </div>

          {/* Re-score with profile */}
          <div className="border border-[var(--color-border)] rounded-xl p-4">
            <div className="flex items-center gap-2 mb-2">
              <BarChart3 className="w-4 h-4 text-[var(--color-accent)]" />
              <h4 className="text-sm font-medium text-[var(--color-text-primary)]">Re-score with Profile</h4>
            </div>
            <p className="text-xs text-[var(--color-text-muted)] mb-3">
              Apply a scoring profile to the selected molecules.
            </p>
            {profiles.length > 0 && (
              <div className="relative mb-3">
                <select
                  value={selectedProfileId ?? ''}
                  onChange={(e) => setSelectedProfileId(e.target.value ? Number(e.target.value) : undefined)}
                  className={cn(
                    'w-full px-3 py-2 rounded-lg text-sm appearance-none',
                    'bg-[var(--color-surface-sunken)] border border-[var(--color-border)]',
                    'text-[var(--color-text-primary)]',
                    'focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]/30'
                  )}
                >
                  <option value="">Default scoring</option>
                  {profiles.map((p) => (
                    <option key={p.id} value={p.id}>
                      {p.name} {p.is_preset ? '(Preset)' : ''}
                    </option>
                  ))}
                </select>
                <ChevronDown className="absolute right-3 top-1/2 -translate-y-1/2 w-4 h-4 text-[var(--color-text-muted)] pointer-events-none" />
              </div>
            )}
            <ClayButton
              variant="accent"
              size="sm"
              onClick={handleRescore}
              disabled={isRescoring}
              loading={isRescoring}
              leftIcon={<BarChart3 className="w-3.5 h-3.5" />}
            >
              Re-score Selection
            </ClayButton>
          </div>

          {/* Export subset */}
          <div className="border border-[var(--color-border)] rounded-xl p-4">
            <div className="flex items-center gap-2 mb-2">
              <Download className="w-4 h-4 text-green-500" />
              <h4 className="text-sm font-medium text-[var(--color-text-primary)]">Export Selection</h4>
            </div>
            <p className="text-xs text-[var(--color-text-muted)] mb-3">
              Download only the selected molecules in your chosen format.
            </p>
            <div className="relative mb-3">
              <select
                value={exportFormat}
                onChange={(e) => setExportFormat(e.target.value as ExportFormat)}
                className={cn(
                  'w-full px-3 py-2 rounded-lg text-sm appearance-none',
                  'bg-[var(--color-surface-sunken)] border border-[var(--color-border)]',
                  'text-[var(--color-text-primary)]',
                  'focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]/30'
                )}
              >
                {EXPORT_FORMATS.map((f) => (
                  <option key={f.value} value={f.value}>
                    {f.label} (.{f.extension})
                  </option>
                ))}
              </select>
              <ChevronDown className="absolute right-3 top-1/2 -translate-y-1/2 w-4 h-4 text-[var(--color-text-muted)] pointer-events-none" />
            </div>
            <ClayButton
              size="sm"
              onClick={handleExport}
              disabled={isExporting}
              loading={isExporting}
              leftIcon={<Download className="w-3.5 h-3.5" />}
            >
              Export Selection
            </ClayButton>
          </div>

          {/* Result feedback */}
          <AnimatePresence>
            {result && (
              <motion.div
                initial={{ opacity: 0, y: 10 }}
                animate={{ opacity: 1, y: 0 }}
                exit={{ opacity: 0, y: -10 }}
                className={cn(
                  'p-3 rounded-xl border text-sm',
                  result.message.startsWith('Failed')
                    ? 'bg-red-500/10 border-red-500/20 text-red-600 dark:text-red-400'
                    : 'bg-green-500/10 border-green-500/20 text-green-600 dark:text-green-400'
                )}
              >
                <p>{result.message}</p>
                {result.jobId && (
                  <a
                    href={`/batch?job=${result.jobId}`}
                    className="flex items-center gap-1 text-xs mt-1 underline"
                  >
                    View new job <ExternalLink className="w-3 h-3" />
                  </a>
                )}
              </motion.div>
            )}
          </AnimatePresence>
        </div>
      </motion.div>
    </div>
  );
}
