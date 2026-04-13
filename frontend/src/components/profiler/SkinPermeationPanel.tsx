import { MetricCard } from './MetricCard';
import type { SkinPermeationResult } from '../../types/profiler';

interface SkinPermeationPanelProps {
  data: SkinPermeationResult;
  className?: string;
}

/**
 * Skin Permeation (Potts-Guy) panel.
 *
 * Uses the Potts-Guy model to predict log Kp (dermal penetration rate).
 * Low permeation is generally desirable for most oral/systemic drugs.
 * Reference: Potts & Guy, Pharm Res 1992.
 */
export function SkinPermeationPanel({ data, className }: SkinPermeationPanelProps) {
  const variantMap = {
    low:      'success',
    moderate: 'warning',
    high:     'error',
  } as const;

  const labelMap = {
    low:      'Low',
    moderate: 'Moderate',
    high:     'High',
  } as const;

  return (
    <MetricCard
      title="Skin Permeation (Potts-Guy)"
      score={data.log_kp.toFixed(2)}
      unit="log Kp (cm/s)"
      classification={labelMap[data.classification]}
      classificationVariant={variantMap[data.classification]}
      className={className}
    >
      <div className="grid grid-cols-1 md:grid-cols-3 gap-4 text-xs text-text-secondary">
        <div className="p-3 rounded-lg bg-[var(--color-surface-sunken)] font-mono">
          log Kp = −2.71 + 0.71 × logP − 0.0061 × MW
        </div>

        <div className="p-3 rounded-lg bg-[var(--color-surface-sunken)] space-y-1">
          <p className="font-medium text-text-secondary">Thresholds</p>
          <p>&lt;−6 low · −6 to −4 moderate · &gt;−4 high</p>
        </div>

        <div className="p-3 rounded-lg bg-[var(--color-surface-sunken)] space-y-1">
          <p className="text-text-muted">Low log Kp is favourable for systemic oral drugs.</p>
          <p className="text-text-muted italic">Potts &amp; Guy, Pharm Res 1992</p>
        </div>
      </div>
    </MetricCard>
  );
}
