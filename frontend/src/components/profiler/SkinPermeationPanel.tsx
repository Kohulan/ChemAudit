import { MetricCard } from './MetricCard';
import type { SkinPermeationResult } from '../../types/profiler';

interface SkinPermeationPanelProps {
  data: SkinPermeationResult;
}

/**
 * Skin Permeation (Potts-Guy) panel.
 *
 * Uses the Potts-Guy model to predict log Kp (dermal penetration rate).
 * Low permeation is generally desirable for most oral/systemic drugs.
 * Reference: Potts & Guy, Pharm Res 1992.
 */
export function SkinPermeationPanel({ data }: SkinPermeationPanelProps) {
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
    >
      <div className="space-y-2 text-sm text-text-secondary">
        <div className="p-3 rounded-lg bg-[var(--color-surface-sunken)] text-xs font-mono">
          log Kp = −2.71 + 0.71 × logP − 0.0061 × MW
        </div>

        <div className="mt-3 p-3 rounded-lg bg-[var(--color-surface-sunken)] text-xs space-y-1">
          <p className="font-medium text-text-secondary">Thresholds</p>
          <p>&lt;−6 — low permeation</p>
          <p>−6 to −4 — moderate permeation</p>
          <p>&gt;−4 — high permeation</p>
        </div>

        <p className="text-xs text-text-muted">
          Low permeation (low log Kp) is generally favourable for systemic oral drugs.
        </p>

        <p className="text-xs text-text-muted">
          Reference: Potts &amp; Guy, Pharm Res 1992
        </p>
      </div>
    </MetricCard>
  );
}
