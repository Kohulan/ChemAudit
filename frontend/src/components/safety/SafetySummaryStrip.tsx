import { StatusDot, CountBadge } from '../ui/Badge';
import { cn } from '../../lib/utils';

interface SafetySummaryStripProps {
  totalAlerts: number;
  concernGroupCount: number;
  hasCritical: boolean;
  cypStatus: 'default' | 'success' | 'warning' | 'error';
  hergStatus: 'success' | 'warning' | 'error';
  bro5Status: 'success' | 'error' | 'default';
  reosStatus: 'success' | 'warning' | 'error';
  complexityOutliers: number;
}

type DotStatus = 'success' | 'warning' | 'error' | 'neutral';

/**
 * Map safety flag status strings to StatusDot status type.
 * 'default' maps to 'neutral' (grey dot, per UI-SPEC for N/A states).
 */
function toStatusDot(status: 'default' | 'success' | 'warning' | 'error'): DotStatus {
  if (status === 'default') return 'neutral';
  return status;
}

/**
 * Sticky summary strip shown after alert screening completes.
 *
 * Layout (per D-06):
 * - Left: alert count badge + "{N} structural alert(s) · {M} concern group(s)"
 * - Center: 4 traffic-light dots for CYP / hERG / bRo5 / REOS
 * - Right: "{N} complexity outlier(s)" or "Complexity within normal range"
 *
 * Colors per UI-SPEC:
 * - success → status-success (warm amber #fbbf24, NOT green)
 * - warning → status-warning (#f59e0b)
 * - error   → status-error (#ef4444)
 * - default → --color-text-muted (#9c958d)
 */
export function SafetySummaryStrip({
  totalAlerts,
  concernGroupCount,
  hasCritical,
  cypStatus,
  hergStatus,
  bro5Status,
  reosStatus,
  complexityOutliers,
}: SafetySummaryStripProps) {
  const alertVariant = hasCritical ? 'error' : totalAlerts > 0 ? 'warning' : 'success';

  const flags = [
    { key: 'CYP', label: 'CYP soft-spots', status: toStatusDot(cypStatus) },
    { key: 'hERG', label: 'hERG liability', status: toStatusDot(hergStatus) },
    { key: 'bRo5', label: 'beyond Rule-of-5', status: toStatusDot(bro5Status) },
    { key: 'REOS', label: 'REOS filter', status: toStatusDot(reosStatus) },
  ] as const;

  return (
    <div
      className={cn(
        'sticky top-[76px] z-30',
        'flex items-center justify-between gap-4',
        'px-5 py-3',
        'bg-[var(--color-surface-elevated)]',
        'border-b border-[var(--color-border)]',
        'shadow-sm',
        'rounded-2xl',
        'mb-6',
      )}
    >
      {/* Left: alert count */}
      <div className="flex items-center gap-2 min-w-0">
        <CountBadge
          count={totalAlerts}
          variant={alertVariant}
        />
        <span className="text-sm text-[var(--color-text-secondary)] truncate">
          {totalAlerts === 0
            ? 'No structural alerts detected'
            : `${totalAlerts} structural alert${totalAlerts !== 1 ? 's' : ''} · ${concernGroupCount} concern group${concernGroupCount !== 1 ? 's' : ''}`}
        </span>
      </div>

      {/* Center: 4 traffic-light dots */}
      <div className="flex items-center gap-4 shrink-0">
        {flags.map(({ key, label, status }) => (
          <div key={key} className="flex items-center gap-1.5">
            <StatusDot status={status} />
            <span
              className="text-xs font-semibold text-[var(--color-text-muted)] font-display"
              aria-hidden="true"
            >
              {key}
            </span>
            {/* Visually hidden accessible label (per D-06 + WCAG) */}
            <span className="sr-only">{label}: {status}</span>
          </div>
        ))}
      </div>

      {/* Right: complexity summary */}
      <div className="flex items-center gap-1.5 shrink-0">
        <span
          className={cn(
            'text-sm font-medium',
            complexityOutliers > 0
              ? 'text-status-warning'
              : 'text-[var(--color-text-secondary)]',
          )}
        >
          {complexityOutliers > 0
            ? `${complexityOutliers} complexity outlier${complexityOutliers !== 1 ? 's' : ''}`
            : 'Complexity within normal range'}
        </span>
      </div>
    </div>
  );
}
