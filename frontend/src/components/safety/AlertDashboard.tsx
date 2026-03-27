import { useState } from 'react';
import { Skeleton } from '../ui/Skeleton';
import { AlertGroupCard } from './AlertGroupCard';
import type { AlertScreenResponse, AlertResult } from '../../types/safety';

interface AlertDashboardProps {
  alertResult: AlertScreenResponse;
  onHighlightAtoms: (atoms: number[]) => void;
  highlightedAlert: string | null;
  isLoading?: boolean;
}

type ViewMode = 'concern' | 'catalog';

/** Determine worst severity among a list of alerts. */
function worstSeverity(alerts: AlertResult[]): string {
  if (alerts.some((a) => a.severity === 'critical')) return 'critical';
  if (alerts.some((a) => a.severity === 'warning')) return 'warning';
  return 'info';
}

/** Group a flat alert list by catalog_source. */
function groupByCatalogSource(
  alerts: AlertResult[],
  enabledCatalogs: Set<string>
): Array<{ catalog: string; alerts: AlertResult[]; severity: string }> {
  const map = new Map<string, AlertResult[]>();
  for (const alert of alerts) {
    if (!enabledCatalogs.has(alert.catalog_source)) continue;
    if (!map.has(alert.catalog_source)) map.set(alert.catalog_source, []);
    map.get(alert.catalog_source)!.push(alert);
  }
  return Array.from(map.entries()).map(([catalog, items]) => ({
    catalog,
    alerts: items,
    severity: worstSeverity(items),
  }));
}

/** Loading skeleton: 3 placeholder rows. */
function DashboardSkeleton() {
  return (
    <div className="space-y-3" aria-busy="true" aria-label="Loading alert groups">
      {[0, 1, 2].map((i) => (
        <div key={i} className="clay-card-sm p-4 flex items-center gap-3 rounded-xl">
          <Skeleton variant="rounded" width={4} height={32} />
          <Skeleton variant="text" height={16} className="flex-1" />
          <Skeleton variant="rounded" width={28} height={20} />
          <Skeleton variant="circular" width={16} height={16} />
        </div>
      ))}
    </div>
  );
}

/**
 * AlertDashboard — dual-view alert visualization with concern-group (default)
 * and catalog-source views, client-side catalog filter toggles, and
 * click-to-highlight atom callback chain.
 *
 * Interaction contracts D-08 through D-11 from 08-UI-SPEC.md.
 */
export function AlertDashboard({
  alertResult,
  onHighlightAtoms,
  highlightedAlert,
  isLoading = false,
}: AlertDashboardProps) {
  const [viewMode, setViewMode] = useState<ViewMode>('concern');
  const [enabledCatalogs, setEnabledCatalogs] = useState<Set<string>>(
    () => new Set(alertResult.screened_catalogs)
  );

  if (isLoading) {
    return (
      <section aria-label="Structural Alerts">
        <h2 className="text-lg font-semibold font-display text-[var(--color-text-primary)] mb-4">
          Structural Alerts
        </h2>
        <DashboardSkeleton />
      </section>
    );
  }

  const allFiltersEnabled =
    enabledCatalogs.size === alertResult.screened_catalogs.length;

  const toggleCatalog = (catalog: string) => {
    setEnabledCatalogs((prev) => {
      const next = new Set(prev);
      if (next.has(catalog)) {
        next.delete(catalog);
      } else {
        next.add(catalog);
      }
      return next;
    });
  };

  const resetFilters = () => {
    setEnabledCatalogs(new Set(alertResult.screened_catalogs));
  };

  // ---- Concern Group View ----
  const concernGroups = Object.entries(alertResult.concern_groups)
    .map(([name, group]) => ({
      name,
      severity: group.severity,
      alerts: group.alerts.filter((a) => enabledCatalogs.has(a.catalog_source)),
    }))
    .filter((g) => g.alerts.length > 0);

  // ---- Catalog Source View ----
  const catalogGroups = groupByCatalogSource(alertResult.alerts, enabledCatalogs);

  const isEmpty =
    (viewMode === 'concern' && concernGroups.length === 0) ||
    (viewMode === 'catalog' && catalogGroups.length === 0);

  return (
    <section aria-label="Structural Alerts">
      {/* Section heading */}
      <h2 className="text-lg font-semibold font-display text-[var(--color-text-primary)] mb-4">
        Structural Alerts
      </h2>

      {/* View toggle — D-08 */}
      <div
        role="radiogroup"
        aria-label="Alert view"
        className="inline-flex rounded-lg overflow-hidden border border-[var(--color-border)] mb-4"
      >
        <button
          type="button"
          role="radio"
          aria-checked={viewMode === 'concern'}
          onClick={() => setViewMode('concern')}
          className={[
            'px-4 py-2 text-sm font-medium transition-colors',
            viewMode === 'concern'
              ? 'bg-chem-primary-600 text-white'
              : 'bg-[var(--color-surface)] text-[var(--color-text-secondary)] border-r border-[var(--color-border)] hover:bg-[var(--color-surface-elevated)]',
          ].join(' ')}
        >
          Concern Group
        </button>
        <button
          type="button"
          role="radio"
          aria-checked={viewMode === 'catalog'}
          onClick={() => setViewMode('catalog')}
          className={[
            'px-4 py-2 text-sm font-medium transition-colors',
            viewMode === 'catalog'
              ? 'bg-chem-primary-600 text-white'
              : 'bg-[var(--color-surface)] text-[var(--color-text-secondary)] hover:bg-[var(--color-surface-elevated)]',
          ].join(' ')}
        >
          Catalog Source
        </button>
      </div>

      {/* Catalog filter toggles — D-11 */}
      {alertResult.screened_catalogs.length > 0 && (
        <div className="flex flex-wrap items-center gap-2 mb-4">
          {alertResult.screened_catalogs.map((catalog) => {
            const enabled = enabledCatalogs.has(catalog);
            return (
              <button
                key={catalog}
                type="button"
                role="checkbox"
                aria-checked={enabled}
                onClick={() => toggleCatalog(catalog)}
                style={{ minHeight: '36px' }}
                className={[
                  'px-3 py-1.5 text-xs font-medium rounded-lg border transition-colors',
                  enabled
                    ? 'bg-[var(--color-surface-elevated)] border-chem-primary-600/40 text-[var(--color-text-primary)]'
                    : 'bg-[var(--color-surface-sunken)] border-[var(--color-border)] text-[var(--color-text-muted)] line-through',
                ].join(' ')}
              >
                {catalog}
              </button>
            );
          })}

          {/* Reset filters button — only when any filter is off */}
          {!allFiltersEnabled && (
            <button
              type="button"
              onClick={resetFilters}
              className="px-3 py-1.5 text-xs text-chem-primary-600 hover:underline transition-colors"
              style={{ minHeight: '36px' }}
            >
              Reset filters
            </button>
          )}
        </div>
      )}

      {/* Empty filtered state */}
      {isEmpty && (
        <div className="text-sm text-[var(--color-text-muted)] py-6 text-center">
          All alerts filtered — toggle catalogs above to show results
        </div>
      )}

      {/* Concern Group View */}
      {!isEmpty && viewMode === 'concern' && (
        <div className="space-y-3">
          {concernGroups.map((group) => (
            <AlertGroupCard
              key={group.name}
              groupName={group.name}
              alerts={group.alerts}
              severity={group.severity}
              defaultExpanded={true}
              onHighlightAtoms={onHighlightAtoms}
              highlightedAlert={highlightedAlert}
            />
          ))}
        </div>
      )}

      {/* Catalog Source View */}
      {!isEmpty && viewMode === 'catalog' && (
        <div className="space-y-3">
          {catalogGroups.map((group) => (
            <AlertGroupCard
              key={group.catalog}
              groupName={group.catalog}
              alerts={group.alerts}
              severity={group.severity}
              defaultExpanded={true}
              onHighlightAtoms={onHighlightAtoms}
              highlightedAlert={highlightedAlert}
            />
          ))}
        </div>
      )}
    </section>
  );
}
