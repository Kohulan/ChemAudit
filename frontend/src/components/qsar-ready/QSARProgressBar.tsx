import { motion } from 'framer-motion';
import type { QSARBatchStatusResponse } from '../../types/qsar_ready';

interface QSARProgressBarProps {
  /** Current batch job status from the WebSocket / polling hook. */
  status: QSARBatchStatusResponse | null;
}

/**
 * Real-time progress indicator for QSAR batch processing.
 *
 * Driven by WebSocket messages via the useQSARReady hook.
 * Per UI-SPEC:
 * - var(--color-primary) crimson fill, rounded-full h-2
 * - role="progressbar" with aria-valuenow/min/max
 * - "Processing molecule {n} of {total}" in text-xs text-secondary
 * - CSS transition-all 300ms linear for progress fill animation
 * - On complete: shows "Curation complete" banner with fade-in
 */
export function QSARProgressBar({ status }: QSARProgressBarProps) {
  if (!status) {
    return (
      <div className="space-y-3">
        <div className="w-full bg-[var(--color-surface-sunken)] rounded-full h-2 overflow-hidden">
          <div
            role="progressbar"
            aria-valuenow={0}
            aria-valuemin={0}
            aria-valuemax={100}
            className="h-full rounded-full bg-[var(--color-primary)] transition-all duration-300 ease-linear"
            style={{ width: '0%' }}
          />
        </div>
        <p className="text-xs text-[var(--color-text-secondary)]">
          Waiting to start…
        </p>
      </div>
    );
  }

  const progress = Math.max(0, Math.min(100, status.progress));
  const isComplete = status.status === 'complete';
  const isFailed = status.status === 'failed';

  return (
    <div className="space-y-3">
      {/* Progress bar */}
      <div className="w-full bg-[var(--color-surface-sunken)] rounded-full h-2 overflow-hidden">
        <div
          role="progressbar"
          aria-valuenow={status.progress}
          aria-valuemin={0}
          aria-valuemax={100}
          className="h-full rounded-full bg-[var(--color-primary)] transition-all duration-300 ease-linear"
          style={{ width: `${progress}%` }}
        />
      </div>

      {/* Status text */}
      {!isComplete && !isFailed && (
        <p className="text-xs text-[var(--color-text-secondary)]">
          {status.status === 'pending'
            ? 'Preparing batch job…'
            : `Processing molecule ${status.processed} of ${status.total}`}
        </p>
      )}

      {/* Completion banner */}
      {isComplete && (
        <motion.div
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          transition={{ duration: 0.3, ease: 'easeOut' }}
          className="bg-green-50 border border-green-200 rounded-lg px-4 py-3 text-sm text-green-700 font-medium"
        >
          Curation complete &mdash;{' '}
          {status.total > 0 ? status.total : status.processed} molecules curated successfully
        </motion.div>
      )}

      {/* Failure banner */}
      {isFailed && (
        <p className="text-xs text-red-600">
          Batch processing failed. Re-upload the file to retry.
        </p>
      )}
    </div>
  );
}
