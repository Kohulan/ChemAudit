import { ClayCard } from '../ui/ClayCard';
import { InfoTooltip } from '../ui/Tooltip';

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

const PIPELINES = [
  { key: 'rdkit_molstandardize', label: 'RDKit MolStandardize' },
  { key: 'chembl_style', label: 'ChEMBL-style' },
  { key: 'minimal_sanitize', label: 'Minimal sanitize' },
] as const;

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

  const pipelines = PIPELINES.map(({ key, label }) => {
    const data = (comparison[key] ?? {}) as Record<string, unknown>;
    return {
      name: label,
      agree: typeof data.agree === 'number' ? data.agree : 0,
      disagree: typeof data.disagree === 'number' ? data.disagree : 0,
    };
  });

  return (
    <div className="space-y-3">
      <div className="flex items-center gap-1.5">
        <h3 className="text-base font-semibold font-display text-[var(--color-text-primary)]">
          Standardization Pipeline Comparison
        </h3>
        <InfoTooltip
          title="Standardization Pipeline Comparison"
          content={
            <div className="text-xs space-y-1">
              <p>Compares how three standardization pipelines process each molecule in a random sample.</p>
              <ul className="mt-1 text-white/70 space-y-0.5">
                <li><strong>RDKit MolStandardize</strong> &mdash; Cleanup, largest fragment, uncharger, tautomer canonicalization</li>
                <li><strong>ChEMBL-style</strong> &mdash; ChEMBL structure pipeline rules</li>
                <li><strong>Minimal</strong> &mdash; Parse + sanitize only (baseline)</li>
              </ul>
              <p className="mt-1 text-white/60">&ldquo;Agree&rdquo; means that pipeline&apos;s canonical SMILES matches at least one other pipeline. Disagreements may indicate tautomeric ambiguity or charge handling differences.</p>
            </div>
          }
          size="small"
        />
      </div>

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
