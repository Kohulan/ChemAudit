import { AlertCircle } from 'lucide-react';
import { MetricCard } from './MetricCard';
import type { ConsensusLogPResult } from '../../types/profiler';

interface ConsensusLogPPanelProps {
  data: ConsensusLogPResult;
}

/**
 * Consensus LogP panel.
 *
 * Combines Wildman-Crippen and XLOGP3 (approximation) values to give
 * a consensus estimate. Per Pitfall 7 from RESEARCH.md, the XLOGP3
 * approximation must always be disclosed when xlogp3_is_approximation is true.
 */
export function ConsensusLogPPanel({ data }: ConsensusLogPPanelProps) {
  // Classification per plan spec
  let classification: string;
  let classificationVariant: 'success' | 'warning' | 'error';

  if (data.consensus_logp > 5) {
    classification = 'Very Lipophilic';
    classificationVariant = 'error';
  } else if (data.consensus_logp > 3) {
    classification = 'Lipophilic';
    classificationVariant = 'warning';
  } else if (data.consensus_logp >= 0) {
    classification = 'Moderate';
    classificationVariant = 'success';
  } else {
    // < 0
    classification = 'Very Hydrophilic';
    classificationVariant = 'success';
  }

  return (
    <MetricCard
      title="Consensus LogP"
      score={data.consensus_logp.toFixed(2)}
      classification={classification}
      classificationVariant={classificationVariant}
    >
      <div className="space-y-2 text-sm text-text-secondary">
        <div className="grid grid-cols-2 gap-x-4 gap-y-1">
          <span className="text-text-muted">Wildman-Crippen</span>
          <span className="tabular-nums text-text-primary font-medium">{data.wildman_crippen.toFixed(2)}</span>

          <span className="text-text-muted">
            XLOGP3
            {data.xlogp3_is_approximation && (
              <span className="ml-1 text-xs text-status-warning">(approx.)</span>
            )}
          </span>
          <span className="tabular-nums text-text-primary font-medium">{data.xlogp3_approx.toFixed(2)}</span>

          <span className="text-text-muted">Consensus</span>
          <span className="tabular-nums text-text-primary font-medium">{data.consensus_logp.toFixed(2)}</span>
        </div>

        <p className="text-xs text-text-muted mt-1">
          Consensus = average of both methods
        </p>

        {/* Pitfall 7 approximation disclosure */}
        {data.xlogp3_is_approximation && (
          <div className="mt-3 flex items-start gap-2 p-3 rounded-lg bg-status-warning/10 border border-status-warning/20 text-xs">
            <AlertCircle className="w-3.5 h-3.5 text-status-warning shrink-0 mt-0.5" />
            <p className="text-text-secondary">
              XLOGP3 value is an approximation computed as{' '}
              <code className="font-mono bg-[var(--color-surface-sunken)] px-1 rounded">
                0.92 × WC + 0.12
              </code>
              {' '}(RDKit does not implement XLOGP3 natively).
            </p>
          </div>
        )}

        <div className="mt-3 p-3 rounded-lg bg-[var(--color-surface-sunken)] text-xs space-y-1">
          <p className="font-medium text-text-secondary">Classification</p>
          <p>&lt;0 — very hydrophilic</p>
          <p>0–3 — moderate</p>
          <p>3–5 — lipophilic</p>
          <p>&gt;5 — very lipophilic</p>
        </div>
      </div>
    </MetricCard>
  );
}
