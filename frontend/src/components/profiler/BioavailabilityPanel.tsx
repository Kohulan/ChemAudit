import { MetricCard } from './MetricCard';
import type { AbbottResult } from '../../types/profiler';

interface BioavailabilityPanelProps {
  data: AbbottResult;
}

/**
 * Abbott Bioavailability Score panel.
 *
 * Estimates probability of oral bioavailability via a 4-class model
 * based on TPSA and Lipinski violations.
 * Probability classes: 11%, 17%, 56%, 85%.
 * Reference: Martin J Med Chem 2005.
 */
export function BioavailabilityPanel({ data }: BioavailabilityPanelProps) {
  // Map probability to human-readable classification
  let classification: string;
  let classificationVariant: 'success' | 'warning' | 'error';

  if (data.probability_pct >= 85) {
    classification = 'High';
    classificationVariant = 'success';
  } else if (data.probability_pct >= 56) {
    classification = 'Moderate';
    classificationVariant = 'warning';
  } else {
    classification = 'Low';
    classificationVariant = 'error';
  }

  return (
    <MetricCard
      title="Abbott Bioavailability Score"
      score={`${data.probability_pct}%`}
      classification={classification}
      classificationVariant={classificationVariant}
    >
      <div className="space-y-2 text-sm text-text-secondary">
        <p className="text-xs text-text-muted">
          Based on TPSA and Lipinski violations (Martin J Med Chem 2005)
        </p>

        <div className="grid grid-cols-2 gap-x-4 gap-y-1 mt-2">
          <span className="text-text-muted">TPSA</span>
          <span className="tabular-nums text-text-primary font-medium">{data.tpsa.toFixed(1)} Å²</span>

          <span className="text-text-muted">Lipinski Violations</span>
          <span className="tabular-nums text-text-primary font-medium">{data.lipinski_violations}</span>

          <span className="text-text-muted">Abbott Score</span>
          <span className="tabular-nums text-text-primary font-medium">{data.abbott_score}</span>
        </div>

        <div className="mt-3 p-3 rounded-lg bg-[var(--color-surface-sunken)] text-xs space-y-1">
          <p className="font-medium text-text-secondary">Probability Classes</p>
          <p>11% — lowest (≥2 violations, TPSA &gt;140)</p>
          <p>17% — low (≥2 violations, TPSA ≤140)</p>
          <p>56% — moderate (0–1 violations, TPSA &gt;75)</p>
          <p>85% — high (0–1 violations, TPSA ≤75)</p>
        </div>
      </div>
    </MetricCard>
  );
}
