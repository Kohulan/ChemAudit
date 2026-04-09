import { useState, useId } from 'react';
import { AnimatePresence, motion } from 'framer-motion';
import { ChevronDown } from 'lucide-react';
import { ClayCard } from '../ui/ClayCard';
import { CountBadge } from '../ui/Badge';
import type { AlertResult } from '../../types/safety';

interface AlertGroupCardProps {
  groupName: string;
  alerts: AlertResult[];
  severity: string;
  defaultExpanded?: boolean;
  onHighlightAtoms: (atoms: number[]) => void;
  highlightedAlert: string | null;
}

/** Severity strip color mapping. */
function getSeverityColor(severity: string): string {
  switch (severity) {
    case 'critical':
      return 'bg-red-500';
    case 'warning':
      return 'bg-amber-500';
    case 'info':
      return 'bg-blue-400';
    default:
      return 'bg-gray-400';
  }
}

/** CountBadge variant matching severity. */
function getSeverityBadgeVariant(severity: string): 'error' | 'warning' | 'default' {
  switch (severity) {
    case 'critical':
      return 'error';
    case 'warning':
      return 'warning';
    default:
      return 'default';
  }
}

/**
 * Return a human-readable NIBR severity badge label from catalog_source / filter_set.
 * Returns null when this is not an NIBR-origin alert.
 */
function getNibrLabel(alert: AlertResult): string | null {
  const source = (alert.catalog_source ?? '').toLowerCase();
  const filterSet = (alert.filter_set ?? '').toLowerCase();

  if (!source.includes('nibr') && !filterSet.includes('nibr')) return null;

  if (filterSet.includes('exclu') || source.includes('exclu')) return 'Excluded';
  if (filterSet.includes('flag') || source.includes('flag')) return 'Flag';
  return 'Annotation';
}

/** Inline badge for NIBR severity coding (shown alongside the alert row). */
function NibrBadge({ label }: { label: string }) {
  const colorClass =
    label === 'Excluded'
      ? 'bg-red-100 text-red-700'
      : label === 'Flag'
        ? 'bg-amber-100 text-amber-700'
        : 'bg-gray-100 text-gray-600';

  return (
    <span className={`px-1.5 py-0.5 text-[10px] font-semibold rounded ${colorClass}`}>
      {label}
    </span>
  );
}

/**
 * AlertGroupCard — expandable concern-group or catalog-source group row.
 *
 * Renders a header button with a severity colour strip, group name, count badge,
 * and a chevron. The expandable section lists individual alert rows as buttons.
 * Clicking an alert row fires `onHighlightAtoms`; re-clicking the same alert clears.
 */
export function AlertGroupCard({
  groupName,
  alerts,
  severity,
  defaultExpanded = true,
  onHighlightAtoms,
  highlightedAlert,
}: AlertGroupCardProps) {
  const [expanded, setExpanded] = useState(defaultExpanded);
  const panelId = useId();

  const stripColor = getSeverityColor(severity);
  const badgeVariant = getSeverityBadgeVariant(severity);

  const handleAlertClick = (alert: AlertResult) => {
    if (alert.pattern_name === highlightedAlert) {
      // Re-click: clear highlight
      onHighlightAtoms([]);
    } else {
      onHighlightAtoms(alert.matched_atoms);
    }
  };

  return (
    <ClayCard size="sm" className="overflow-hidden p-0">
      {/* Header toggle button */}
      <button
        type="button"
        onClick={() => setExpanded(!expanded)}
        aria-expanded={expanded}
        aria-controls={panelId}
        className="w-full flex items-center gap-3 px-4 py-3 text-left hover:bg-[var(--color-surface-elevated)] transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-chem-primary-500 focus-visible:ring-inset"
      >
        {/* Severity strip */}
        <div className={`w-1 self-stretch rounded-full shrink-0 ${stripColor}`} />

        {/* Group name */}
        <span className="flex-1 text-sm font-semibold text-[var(--color-text-primary)] font-sans">
          {groupName}
        </span>

        {/* Count badge */}
        <CountBadge count={alerts.length} variant={badgeVariant} />

        {/* Chevron */}
        <motion.span
          animate={{ rotate: expanded ? 180 : 0 }}
          transition={{ duration: 0.2, ease: 'easeOut' }}
          className="shrink-0 text-[var(--color-text-muted)]"
        >
          <ChevronDown className="w-4 h-4" />
        </motion.span>
      </button>

      {/* Expandable alert list */}
      <AnimatePresence initial={false}>
        {expanded && (
          <motion.div
            id={panelId}
            key="panel"
            initial={{ height: 0, opacity: 0 }}
            animate={{ height: 'auto', opacity: 1 }}
            exit={{ height: 0, opacity: 0 }}
            transition={{ duration: 0.3, ease: 'easeOut' }}
            className="overflow-hidden"
          >
            <div className="divide-y divide-[var(--color-border)]">
              {alerts.map((alert) => {
                const isHighlighted = alert.pattern_name === highlightedAlert;
                const nibrLabel = getNibrLabel(alert);

                return (
                  <button
                    key={alert.pattern_name}
                    type="button"
                    onClick={() => handleAlertClick(alert)}
                    aria-label={`${alert.pattern_name} — click to highlight matched atoms`}
                    className={[
                      'w-full flex items-start gap-3 px-4 py-2.5 text-left transition-colors',
                      'hover:bg-[var(--color-surface-elevated)]',
                      'focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-chem-primary-500 focus-visible:ring-inset',
                      'min-h-[44px]', // WCAG 2.5.5 touch target
                      isHighlighted ? 'border-l-2 border-amber-400 pl-[14px]' : '',
                    ]
                      .filter(Boolean)
                      .join(' ')}
                  >
                    {/* Pattern name + NIBR badge */}
                    <div className="flex-1 min-w-0">
                      <div className="flex items-center gap-1.5 flex-wrap">
                        <span className="text-sm text-[var(--color-text-primary)]">
                          {alert.pattern_name}
                        </span>
                        {nibrLabel && <NibrBadge label={nibrLabel} />}
                      </div>
                      {/* Source attribution */}
                      <p className="text-xs text-[var(--color-text-muted)] mt-0.5">
                        {alert.catalog_source}
                      </p>
                    </div>
                  </button>
                );
              })}
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </ClayCard>
  );
}
