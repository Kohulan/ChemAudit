import { useState, useEffect } from 'react';
import type { DatasetAuditStatus } from '../../types/dataset_intelligence';
import { ClayCard } from '../ui/ClayCard';

// =============================================================================
// Types
// =============================================================================

interface DatasetProgressBarProps {
  /** Progress percentage (0-100). */
  progress: number;
  /** Current pipeline stage name. */
  currentStage: string | null;
  /** Current audit lifecycle status. */
  status: DatasetAuditStatus;
}

// =============================================================================
// Stage label map
// =============================================================================

/**
 * Human-readable stage labels per UI-SPEC.
 * Maps backend stage identifiers to user-facing copy.
 */
const STAGE_LABELS: Record<string, string> = {
  parsing: 'Parsing molecules...',
  stereo: 'Checking stereochemistry...',
  duplicates: 'Detecting duplicates...',
  alerts: 'Screening structural alerts...',
  standardization: 'Comparing standardization pipelines...',
  scoring: 'Computing health scores...',
  contradictions: 'Detecting contradictory labels...',
};

// =============================================================================
// Component
// =============================================================================

/**
 * Progress bar for dataset audit processing.
 *
 * Follows the StructureFilterProgressBar.tsx pattern exactly:
 * - role="progressbar" with aria attributes
 * - Smooth width transition (0.3s ease-out)
 * - 500ms delay before showing to avoid flash on fast jobs
 * - Stage label displayed below the bar
 * - Returns null when not in processing status
 */
export function DatasetProgressBar({
  progress,
  currentStage,
  status,
}: DatasetProgressBarProps) {
  // 500ms delay to prevent flash for fast jobs (StructureFilterProgressBar pattern)
  const [visible, setVisible] = useState(false);

  useEffect(() => {
    const timer = setTimeout(() => setVisible(true), 500);
    return () => clearTimeout(timer);
  }, []);

  // Only show during processing states
  if (status !== 'processing' && status !== 'uploading') {
    return null;
  }

  if (!visible) {
    return null;
  }

  const pct = Math.max(0, Math.min(100, progress));
  const stageLabel = currentStage
    ? STAGE_LABELS[currentStage] ?? currentStage
    : 'Preparing...';

  return (
    <ClayCard className="p-5 space-y-4">
      {/* Heading */}
      <h3 className="text-sm font-semibold text-[var(--color-text-primary)]">
        Analyzing dataset&hellip;
      </h3>

      {/* Progress bar */}
      <div>
        <div
          className="w-full bg-[var(--color-surface-sunken)] rounded-full h-2 overflow-hidden"
          role="progressbar"
          aria-valuenow={pct}
          aria-valuemin={0}
          aria-valuemax={100}
          aria-label="Dataset analysis progress"
        >
          <div
            className="h-full rounded-full bg-[var(--color-primary)]"
            style={{
              width: `${pct}%`,
              transition: 'width 0.3s ease-out',
            }}
          />
        </div>

        {/* Stage label + percentage */}
        <div className="flex items-center justify-between mt-1.5">
          <p className="text-xs text-[var(--color-text-secondary)]">
            {stageLabel}
          </p>
          <p className="text-xs text-[var(--color-text-muted)]">{pct}%</p>
        </div>
      </div>
    </ClayCard>
  );
}
