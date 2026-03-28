import { AlertCircle } from 'lucide-react';
import { ClayCard } from '../ui/ClayCard';
import { Badge } from '../ui/Badge';
import type { SAComparisonResult } from '../../types/profiler';
import { cn } from '../../lib/utils';

interface SAComparisonPanelProps {
  data: SAComparisonResult;
}

/**
 * Returns the Badge variant for an SA/SCScore classification.
 */
function classificationVariant(
  classification: string
): 'success' | 'warning' | 'error' {
  if (classification === 'easy') return 'success';
  if (classification === 'moderate') return 'warning';
  return 'error';
}

/**
 * SA Comparison Panel — shows 4 cards: SA Score, SCScore, SYBA, RAscore (slot).
 *
 * Per D-17, D-18, D-19:
 * - SA Score: always available, 1-10 scale
 * - SCScore: may be unavailable (show grey fallback text)
 * - SYBA: may be unavailable (greyed card, not error styling, per D-18)
 * - RAscore: always "coming soon" slot (greyed)
 */
export function SAComparisonPanel({ data }: SAComparisonPanelProps) {
  return (
    <div>
      <h3 className="text-lg font-semibold text-text-primary font-display mb-4">
        Synthesizability Comparison
      </h3>
      {/* D-17 layout: 4 cards per row at md+, 2 per row at sm */}
      <div className="grid grid-cols-2 md:grid-cols-4 gap-3 md:gap-4">

        {/* Card 1 — SA Score (always available) */}
        <ClayCard size="sm" className="flex flex-col gap-2">
          <p className="text-xs font-semibold text-text-secondary uppercase tracking-wide">
            SA Score
          </p>
          <div className="flex items-baseline gap-1">
            <span className="text-2xl font-semibold tabular-nums text-text-primary font-display">
              {data.sa_score.score.toFixed(2)}
            </span>
            <span className="text-xs text-text-muted">(1–10)</span>
          </div>
          <Badge
            variant={classificationVariant(data.sa_score.classification)}
            size="sm"
            className="self-start"
          >
            {data.sa_score.classification}
          </Badge>
        </ClayCard>

        {/* Card 2 — SCScore */}
        {data.scscore.available ? (
          <ClayCard size="sm" className="flex flex-col gap-2">
            <p className="text-xs font-semibold text-text-secondary uppercase tracking-wide">
              SCScore
            </p>
            <div className="flex items-baseline gap-1">
              <span className="text-2xl font-semibold tabular-nums text-text-primary font-display">
                {data.scscore.score.toFixed(2)}
              </span>
              <span className="text-xs text-text-muted">(1–5)</span>
            </div>
            <Badge
              variant={classificationVariant(data.scscore.classification)}
              size="sm"
              className="self-start"
            >
              {data.scscore.classification}
            </Badge>
          </ClayCard>
        ) : (
          <ClayCard
            size="sm"
            className={cn('flex flex-col gap-2 bg-surface-sunken')}
          >
            <p className="text-xs font-semibold text-text-muted uppercase tracking-wide">
              SCScore
            </p>
            <p className="text-sm text-text-muted">
              SCScore unavailable
            </p>
            {data.scscore.error && (
              <p className="text-xs text-text-muted">{data.scscore.error}</p>
            )}
          </ClayCard>
        )}

        {/* Card 3 — SYBA (may be unavailable per D-18) */}
        {data.syba.available ? (
          <ClayCard size="sm" className="flex flex-col gap-2">
            <p className="text-xs font-semibold text-text-secondary uppercase tracking-wide">
              SYBA
            </p>
            <div className="flex items-baseline gap-1">
              <span className="text-2xl font-semibold tabular-nums text-text-primary font-display">
                {data.syba.score.toFixed(1)}
              </span>
              <span className="text-xs text-text-muted">(-166 to +275)</span>
            </div>
            <Badge
              variant={classificationVariant(data.syba.classification)}
              size="sm"
              className="self-start"
            >
              {data.syba.classification}
            </Badge>
          </ClayCard>
        ) : (
          /* D-18: greyed card, not error styling */
          <ClayCard
            size="sm"
            className={cn('flex flex-col gap-2 bg-surface-sunken')}
          >
            <p className="text-xs font-semibold text-text-muted uppercase tracking-wide">
              SYBA
            </p>
            <div className="flex items-start gap-1.5 text-text-muted">
              <AlertCircle className="w-3.5 h-3.5 mt-0.5 flex-shrink-0" />
              <p className="text-xs leading-relaxed">
                SYBA score unavailable (subprocess isolation required). SA Score and SCScore are shown.
              </p>
            </div>
          </ClayCard>
        )}

        {/* Card 4 — RAscore (always greyed, coming soon per D-19) */}
        <ClayCard
          size="sm"
          className={cn('flex flex-col gap-2 bg-surface-sunken')}
        >
          <p className="text-xs font-semibold text-text-muted uppercase tracking-wide">
            RAscore
          </p>
          <p className="text-sm text-text-muted italic">
            Coming in a future release
          </p>
        </ClayCard>

      </div>
    </div>
  );
}
