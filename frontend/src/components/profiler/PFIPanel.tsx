import { MetricCard } from './MetricCard';
import type { PFIResult } from '../../types/profiler';

interface PFIPanelProps {
  data: PFIResult;
}

/**
 * PFI (Property Forecast Index) metric panel.
 *
 * PFI = cLogP + #Aromatic Rings. Lower scores predict better
 * developability. Thresholds: <5 low, 5–7 moderate, >7 high risk.
 * Reference: Young RJ et al., Drug Discov Today 2011.
 */
export function PFIPanel({ data }: PFIPanelProps) {
  const variantMap = {
    low:      'success',
    moderate: 'warning',
    high:     'error',
  } as const;

  const labelMap = {
    low:      'Low Risk',
    moderate: 'Moderate Risk',
    high:     'High Risk',
  } as const;

  return (
    <MetricCard
      title="Property Forecast Index (PFI)"
      score={data.pfi.toFixed(2)}
      classification={labelMap[data.risk]}
      classificationVariant={variantMap[data.risk]}
    >
      <div className="space-y-2 text-sm text-text-secondary">
        <p className="font-medium text-text-primary">
          PFI = cLogP + #Aromatic Rings
        </p>

        <div className="grid grid-cols-2 gap-x-4 gap-y-1 mt-2">
          <span className="text-text-muted">cLogP</span>
          <span className="tabular-nums text-text-primary font-medium">{data.clogp.toFixed(2)}</span>

          <span className="text-text-muted">Aromatic Rings</span>
          <span className="tabular-nums text-text-primary font-medium">{data.aromatic_rings}</span>
        </div>

        <div className="mt-3 p-3 rounded-lg bg-[var(--color-surface-sunken)] text-xs space-y-1">
          <p className="font-medium text-text-secondary">Thresholds</p>
          <p>&lt;5 — low risk</p>
          <p>5–7 — moderate risk</p>
          <p>&gt;7 — high risk</p>
        </div>

        <p className="text-xs text-text-muted mt-2">
          Reference: Young RJ et al., Drug Discov Today 2011
        </p>
      </div>
    </MetricCard>
  );
}
