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
        {/* Property rows — compact 3-col: status icon + name, value, range */}
        {data.details.map((detail) => (
          <div
            key={detail.property}
            className={cn(
              'flex items-center gap-2 px-2 py-1.5 rounded-lg text-xs',
              detail.violated
                ? 'bg-status-error/10 text-status-error'
                : 'text-text-secondary'
            )}
          >
            {detail.violated ? (
              <XCircle className="w-3.5 h-3.5 text-status-error shrink-0" aria-label="Violated" />
            ) : (
              <CheckCircle2 className="w-3.5 h-3.5 text-status-success shrink-0" aria-label="In range" />
            )}
            <span className="font-medium flex-1 min-w-0" title={detail.property}>
              {detail.property}
            </span>
            <span className="tabular-nums shrink-0">
              {typeof detail.value === 'number' ? detail.value.toFixed(1) : detail.value}
            </span>
            <span className="tabular-nums text-text-muted shrink-0 hidden sm:inline">
              ({detail.range_low}–{detail.range_high})
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
