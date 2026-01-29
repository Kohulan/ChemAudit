import { useEffect, useState, useRef } from 'react';
import type { BatchProgress as BatchProgressType } from '../../types/batch';

interface BatchProgressProps {
  progress: BatchProgressType | null;
  isConnected: boolean;
  onCancel: () => void;
  onComplete: () => void;
}

/**
 * Progress display component with animated progress bar and ETA.
 */
export function BatchProgress({
  progress,
  isConnected,
  onCancel,
  onComplete,
}: BatchProgressProps) {
  // Smooth progress animation
  const [displayProgress, setDisplayProgress] = useState(0);
  const animationRef = useRef<number | null>(null);

  useEffect(() => {
    const targetProgress = progress?.progress ?? 0;

    // Animate progress smoothly
    const animate = () => {
      setDisplayProgress((current) => {
        const diff = targetProgress - current;
        if (Math.abs(diff) < 0.5) {
          return targetProgress;
        }
        // Ease-out animation
        const step = diff * 0.1;
        animationRef.current = requestAnimationFrame(animate);
        return current + step;
      });
    };

    animationRef.current = requestAnimationFrame(animate);

    return () => {
      if (animationRef.current) {
        cancelAnimationFrame(animationRef.current);
      }
    };
  }, [progress?.progress]);

  // Call onComplete when status changes to complete
  useEffect(() => {
    if (progress?.status === 'complete') {
      onComplete();
    }
  }, [progress?.status, onComplete]);

  const formatETA = (seconds: number | null): string => {
    if (seconds === null || seconds <= 0) return 'Calculating...';

    if (seconds < 60) {
      return `${Math.ceil(seconds)} seconds`;
    }
    if (seconds < 3600) {
      const mins = Math.floor(seconds / 60);
      const secs = Math.ceil(seconds % 60);
      return `${mins}:${secs.toString().padStart(2, '0')}`;
    }
    const hours = Math.floor(seconds / 3600);
    const mins = Math.floor((seconds % 3600) / 60);
    return `${hours}:${mins.toString().padStart(2, '0')}:00`;
  };

  const getStatusColor = () => {
    switch (progress?.status) {
      case 'complete':
        return 'bg-yellow-500';
      case 'failed':
        return 'bg-red-500';
      case 'cancelled':
        return 'bg-amber-500';
      default:
        return 'bg-[var(--color-primary)]';
    }
  };

  const getStatusText = () => {
    switch (progress?.status) {
      case 'pending':
        return 'Waiting to start...';
      case 'processing':
        return `Processing ${progress.processed} of ${progress.total} molecules...`;
      case 'complete':
        return 'Processing complete!';
      case 'failed':
        return `Failed: ${progress.error_message || 'Unknown error'}`;
      case 'cancelled':
        return 'Cancelled';
      default:
        return 'Connecting...';
    }
  };

  return (
    <div className="bg-[var(--color-surface-elevated)] rounded-lg shadow-md border border-[var(--color-border)] p-6 space-y-4">
      {/* Status header */}
      <div className="flex items-center justify-between">
        <h3 className="text-lg font-semibold text-[var(--color-text-primary)]">Batch Processing</h3>
        {/* Connection indicator */}
        <div className="flex items-center space-x-2">
          <span
            className={`w-2 h-2 rounded-full ${
              isConnected ? 'bg-yellow-500' : 'bg-[var(--color-text-muted)]'
            }`}
          />
          <span className="text-xs text-[var(--color-text-muted)]">
            {isConnected ? 'Connected' : 'Disconnected'}
          </span>
        </div>
      </div>

      {/* Progress bar */}
      <div className="space-y-2">
        <div className="h-4 bg-[var(--color-surface-sunken)] rounded-full overflow-hidden">
          <div
            className={`h-full transition-all duration-300 ease-out ${getStatusColor()}`}
            style={{ width: `${Math.max(0, Math.min(100, displayProgress))}%` }}
          />
        </div>

        <div className="flex justify-between text-sm">
          <span className="text-[var(--color-text-secondary)]">{getStatusText()}</span>
          <span className="font-medium text-[var(--color-text-primary)]">
            {Math.round(displayProgress)}%
          </span>
        </div>
      </div>

      {/* Stats row */}
      {progress && progress.status === 'processing' && (
        <div className="grid grid-cols-3 gap-4 py-3 border-t border-b border-[var(--color-border)]">
          <div className="text-center">
            <p className="text-2xl font-semibold text-[var(--color-text-primary)]">
              {progress.processed}
            </p>
            <p className="text-xs text-[var(--color-text-muted)]">Processed</p>
          </div>
          <div className="text-center">
            <p className="text-2xl font-semibold text-[var(--color-text-primary)]">
              {progress.total - progress.processed}
            </p>
            <p className="text-xs text-[var(--color-text-muted)]">Remaining</p>
          </div>
          <div className="text-center">
            <p className="text-2xl font-semibold text-[var(--color-text-primary)]">
              {formatETA(progress.eta_seconds)}
            </p>
            <p className="text-xs text-[var(--color-text-muted)]">ETA</p>
          </div>
        </div>
      )}

      {/* Error message */}
      {progress?.status === 'failed' && progress.error_message && (
        <div className="bg-red-500/10 dark:bg-red-500/20 border border-red-500/30 rounded-lg p-3">
          <p className="text-sm text-red-600 dark:text-red-400">{progress.error_message}</p>
        </div>
      )}

      {/* Cancel button */}
      {progress?.status === 'processing' && (
        <button
          onClick={onCancel}
          className="w-full py-2 px-4 border border-red-500/50 text-red-600 dark:text-red-400 rounded-lg
                     hover:bg-red-500/10 transition-colors duration-200"
        >
          Cancel Processing
        </button>
      )}

      {/* Completion message */}
      {progress?.status === 'complete' && (
        <div className="bg-yellow-500/10 dark:bg-yellow-500/20 border border-yellow-500/30 rounded-lg p-4 text-center">
          <svg
            className="w-8 h-8 text-yellow-500 mx-auto mb-2"
            fill="none"
            strokeLinecap="round"
            strokeLinejoin="round"
            strokeWidth={2}
            viewBox="0 0 24 24"
            stroke="currentColor"
          >
            <path d="M9 12l2 2 4-4m6 2a9 9 0 11-18 0 9 9 0 0118 0z" />
          </svg>
          <p className="text-amber-700 dark:text-yellow-400 font-medium">
            Successfully processed {progress.total} molecules!
          </p>
          <p className="text-amber-600 dark:text-yellow-500 text-sm mt-1">Loading results...</p>
        </div>
      )}
    </div>
  );
}
