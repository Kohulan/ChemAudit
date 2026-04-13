import { useState, useEffect } from 'react';
import { ClayCard } from '../ui/ClayCard';

// =============================================================================
// Types
// =============================================================================

interface StructureFilterProgressBarProps {
  /** Progress percentage (0-100), or null while unknown. */
  progress: number | null;
  /** Current pipeline stage name. */
  currentStage: string | null;
  /** Total number of molecules being filtered. */
  totalMolecules: number;
}

// =============================================================================
// Component
// =============================================================================

/**
 * Async batch progress indicator for Structure Filter.
 *
 * Per UI-SPEC D-22 and state machine specification:
 * - ClayCard wrapper
 * - Heading: "Filtering {N} molecules..."
 * - Progress bar: var(--color-surface) track, var(--color-primary) fill, width: {progress}%
 * - Stage label below bar in text-xs text-secondary
 * - Async threshold notice: amber-50/border-amber-200 banner
 * - 500ms delay before showing to avoid flash on fast jobs
 * - Transition: width 0.3s ease-out on progress fill
 */
export function StructureFilterProgressBar({
  progress,
  currentStage,
  totalMolecules,
}: StructureFilterProgressBarProps) {
  // 500ms delay to prevent flash for fast jobs (UI-SPEC loading state)
  const [visible, setVisible] = useState(false);

  useEffect(() => {
    const timer = setTimeout(() => setVisible(true), 500);
    return () => clearTimeout(timer);
  }, []);

  if (!visible) {
    return null;
  }

  const pct = progress !== null ? Math.max(0, Math.min(100, progress)) : 0;

  return (
    <ClayCard className="p-5 space-y-4">
      {/* Heading */}
      <h3 className="text-sm font-semibold text-[var(--color-text-primary)]">
        Filtering {totalMolecules.toLocaleString()} molecules&hellip;
      </h3>

      {/* Progress bar */}
      <div>
        <div
          className="w-full bg-[var(--color-surface-sunken)] rounded-full h-2 overflow-hidden"
          role="progressbar"
          aria-valuenow={pct}
          aria-valuemin={0}
          aria-valuemax={100}
          aria-label={`Filtering progress: ${pct}%`}
        >
          <div
            className="h-full rounded-full bg-[var(--color-primary)]"
            style={{
              width: `${pct}%`,
              transition: 'width 0.3s ease-out',
            }}
          />
        </div>

        {/* Stage label */}
        <div className="flex items-center justify-between mt-1.5">
          <p className="text-xs text-[var(--color-text-secondary)]">
            {currentStage ? (
              <>
                Stage:{' '}
                <span className="font-medium text-[var(--color-text-primary)]">
                  {currentStage}
                </span>
              </>
            ) : (
              'Preparing…'
            )}
          </p>
          <p className="text-xs text-[var(--color-text-muted)]">{pct}%</p>
        </div>
      </div>

      {/* Async threshold notice banner (amber-50/border-amber-200 per UI-SPEC) */}
      <div className="bg-amber-50 border border-amber-200 rounded-lg px-4 py-3 dark:bg-amber-900/20 dark:border-amber-700/40">
        <p className="text-xs text-amber-700 dark:text-amber-400">
          Filtering {totalMolecules.toLocaleString()} molecules asynchronously. Results will appear
          when complete.
        </p>
      </div>
    </ClayCard>
  );
}
