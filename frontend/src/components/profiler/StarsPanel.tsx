import { CheckCircle2, XCircle } from 'lucide-react';
import { MetricCard } from './MetricCard';
import { cn } from '../../lib/utils';
import type { StarsResult } from '../../types/profiler';

interface StarsPanelProps {
  data: StarsResult;
}

/**
 * #Stars panel — QikProp-equivalent property outlier count.
 *
 * Shows all available property ranges with violations highlighted.
 * 0 stars = drug-like, 1–2 = borderline, 3+ = outlier.
 * Per UI-SPEC D-06: all rows shown, violated rows highlighted.
 */
export function StarsPanel({ data }: StarsPanelProps) {
  const starsCount = data.stars;

  let classification: string;
  let classificationVariant: 'success' | 'warning' | 'error';

  if (starsCount === 0) {
    classification = 'Drug-like';
    classificationVariant = 'success';
  } else if (starsCount <= 2) {
    classification = 'Borderline';
    classificationVariant = 'warning';
  } else {
    classification = 'Outlier';
    classificationVariant = 'error';
  }

  return (
    <MetricCard
      title="#Stars (Properties Outside Drug Range)"
      score={starsCount}
      classification={classification}
      classificationVariant={classificationVariant}
    >
      <div className="space-y-1 text-sm">
        {/* Table header */}
        <div className="grid grid-cols-[1fr_auto_auto_auto] gap-x-3 px-2 pb-1 border-b border-[var(--color-border)]">
          <span className="text-xs font-medium text-text-muted uppercase tracking-wider">Property</span>
          <span className="text-xs font-medium text-text-muted uppercase tracking-wider text-right">Value</span>
          <span className="text-xs font-medium text-text-muted uppercase tracking-wider text-right">Range</span>
          <span className="text-xs font-medium text-text-muted uppercase tracking-wider text-center">Status</span>
        </div>

        {/* Property rows */}
        {data.details.map((detail) => (
          <div
            key={detail.property}
            className={cn(
              'grid grid-cols-[1fr_auto_auto_auto] gap-x-3 px-2 py-1.5 rounded-lg',
              detail.violated
                ? 'bg-status-error/10 text-status-error'
                : 'text-text-secondary'
            )}
          >
            <span className="text-xs font-medium truncate" title={detail.property}>
              {detail.property}
            </span>
            <span className="text-xs tabular-nums text-right">
              {typeof detail.value === 'number' ? detail.value.toFixed(2) : detail.value}
            </span>
            <span className="text-xs tabular-nums text-right text-text-muted whitespace-nowrap">
              {detail.range_low}–{detail.range_high}
            </span>
            <span className="flex items-center justify-center">
              {detail.violated ? (
                <XCircle className="w-3.5 h-3.5 text-status-error" aria-label="Violated" />
              ) : (
                <CheckCircle2 className="w-3.5 h-3.5 text-status-success" aria-label="In range" />
              )}
            </span>
          </div>
        ))}

        {data.details.length === 0 && (
          <p className="text-xs text-text-muted px-2 py-2">No property details available.</p>
        )}
      </div>
    </MetricCard>
  );
}
