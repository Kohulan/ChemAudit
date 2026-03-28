import { ClayCard } from '../ui/ClayCard';

// =============================================================================
// Types
// =============================================================================

interface StdConsistencyPanelProps {
  /** Pipeline comparison data from health audit. */
  comparison: Record<string, unknown>;
  /** Number of molecules sampled for comparison. */
  sampleSize: number;
}

// =============================================================================
// Pipeline configuration
// =============================================================================

const PIPELINE_NAMES = [
  'RDKit MolStandardize',
  'ChEMBL-style',
  'Minimal sanitize',
];

// =============================================================================
// Component
// =============================================================================

/**
 * Standardization pipeline comparison panel.
 *
 * Follows the CrossPipelinePanel.tsx pattern from Phase 9.
 *
 * Shows a table comparing 3 standardization pipelines (agree/disagree counts),
 * with a notice about the sample size used for estimation.
 */
export function StdConsistencyPanel({ comparison, sampleSize }: StdConsistencyPanelProps) {
  // If no data available
  if (sampleSize === 0) {
    return (
      <div className="space-y-3">
        <h3 className="text-base font-semibold font-display text-[var(--color-text-primary)]">
          Standardization Pipeline Comparison
        </h3>
        <div className="flex items-center justify-center py-8 text-sm text-[var(--color-text-muted)]">
          No standardization comparison data available.
        </div>
      </div>
    );
  }

  // Extract pipeline results from comparison data
  const pipelines = PIPELINE_NAMES.map((name) => {
    const key = name.toLowerCase().replace(/\s+/g, '_').replace(/-/g, '_');
    const pipelineData = (comparison[key] ?? comparison[name] ?? {}) as Record<string, unknown>;
    return {
      name,
      agree: typeof pipelineData.agree === 'number' ? pipelineData.agree : 0,
      disagree: typeof pipelineData.disagree === 'number' ? pipelineData.disagree : 0,
    };
  });

  return (
    <div className="space-y-3">
      <h3 className="text-base font-semibold font-display text-[var(--color-text-primary)]">
        Standardization Pipeline Comparison
      </h3>

      <ClayCard variant="flat" size="sm" className="p-0 overflow-hidden">
        <div className="overflow-x-auto">
          <table className="w-full text-sm">
            <thead>
              <tr className="border-b border-[var(--color-border)] bg-[var(--color-surface-sunken)]">
                <th className="text-left px-4 py-2.5 text-xs font-semibold text-[var(--color-text-muted)]">
                  Pipeline
                </th>
                <th className="text-right px-4 py-2.5 text-xs font-semibold text-[var(--color-text-muted)]">
                  Agree
                </th>
                <th className="text-right px-4 py-2.5 text-xs font-semibold text-[var(--color-text-muted)]">
                  Disagree
                </th>
                <th className="text-right px-4 py-2.5 text-xs font-semibold text-[var(--color-text-muted)]">
                  Sample
                </th>
              </tr>
            </thead>
            <tbody>
              {pipelines.map((pipeline) => {
                const total = pipeline.agree + pipeline.disagree;
                const agreePercent = total > 0 ? Math.round((pipeline.agree / total) * 100) : 0;
                return (
                  <tr
                    key={pipeline.name}
                    className="border-b border-[var(--color-border)] last:border-0 hover:bg-[var(--color-surface-sunken)] transition-colors"
                  >
                    <td className="px-4 py-2.5 text-sm text-[var(--color-text-primary)] font-medium">
                      {pipeline.name}
                    </td>
                    <td className="px-4 py-2.5 text-sm text-right text-emerald-600 dark:text-emerald-400 tabular-nums">
                      {pipeline.agree}
                      <span className="text-xs text-[var(--color-text-muted)] ml-1">
                        ({agreePercent}%)
                      </span>
                    </td>
                    <td className="px-4 py-2.5 text-sm text-right text-red-600 dark:text-red-400 tabular-nums">
                      {pipeline.disagree}
                    </td>
                    <td className="px-4 py-2.5 text-sm text-right text-[var(--color-text-secondary)] tabular-nums">
                      {sampleSize}
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>

        {/* Sample notice */}
        <div className="px-4 py-2.5 border-t border-[var(--color-border)] bg-[var(--color-surface-sunken)]">
          <p className="text-xs text-[var(--color-text-muted)] italic">
            Estimated from a random sample of {sampleSize} molecules.
          </p>
        </div>
      </ClayCard>
    </div>
  );
}
