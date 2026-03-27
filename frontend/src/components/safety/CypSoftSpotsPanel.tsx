import { cn } from '../../lib/utils';
import type { CypSite } from '../../types/safety';

interface CypSoftSpotsPanelProps {
  sites: CypSite[];
  onHighlightAtoms: (atoms: number[]) => void;
  highlightedAlert: string | null;
}

/**
 * Expandable CYP soft-spots detail panel.
 *
 * Renders a list of CYP metabolism sites. Each site is a clickable button
 * that highlights matched atoms on the molecule viewer. Clicking the same
 * site again clears the highlight (toggle behavior).
 *
 * Used as the children slot of the CYP MetricCard in SafetyFlagsGrid.
 * Per D-14: shares the same onHighlightAtoms callback as AlertDashboard.
 */
export function CypSoftSpotsPanel({
  sites,
  onHighlightAtoms,
  highlightedAlert,
}: CypSoftSpotsPanelProps) {
  if (sites.length === 0) {
    return (
      <p className="text-sm text-text-muted">No CYP soft-spots detected</p>
    );
  }

  return (
    <div className="space-y-1.5">
      {sites.map((site) => {
        const isActive = highlightedAlert === site.site_name;
        return (
          <button
            key={site.site_name}
            type="button"
            onClick={() => {
              if (isActive) {
                onHighlightAtoms([]);
              } else {
                onHighlightAtoms(site.matched_atoms);
              }
            }}
            aria-label={`${site.site_name} (${site.reaction_type}) -- click to highlight matched atoms`}
            className={cn(
              'w-full text-left px-3 py-2 rounded-lg border transition-colors',
              'hover:bg-[var(--color-surface-sunken)]',
              isActive
                ? 'border-l-4 border-l-amber-400 border-[var(--color-border)] bg-amber-50/30 dark:bg-amber-900/10'
                : 'border border-[var(--color-border)]'
            )}
          >
            <span className="text-sm font-medium text-text-primary">{site.site_name}</span>
            <span className="block text-xs text-text-muted mt-0.5">{site.reaction_type}</span>
          </button>
        );
      })}
    </div>
  );
}
