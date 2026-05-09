import type { ComponentType } from 'react';
import { cn } from '../../lib/utils';

export type TabBarResultState = 'clean' | 'issues' | 'warnings' | 'complete';

export interface TabBarTab<T extends string = string> {
  id: T;
  label: string;
  icon: ComponentType<{ className?: string }>;
  description?: string;
  /** When true, render a small dot on the inactive tab to signal it has results. */
  hasResult?: boolean;
  /**
   * When `hasResult` is true, this further qualifies the dot colour so a glance
   * at the tab bar communicates whether the result is clean, has issues, has
   * warnings, or is just complete (no chemistry-valence distinction applies).
   * All four colours stay inside the warm-status spectrum per DESIGN.md.
   */
  resultState?: TabBarResultState;
}

const RESULT_STATE_DOT_CLASS: Record<TabBarResultState, string> = {
  clean: 'bg-[#fbbf24]', // score-excellent gold-amber
  issues: 'bg-[#ea580c]', // score-fair orange-flame
  warnings: 'bg-[#d97706]', // score-good amber
  complete: 'bg-[var(--color-primary)]',
};

const RESULT_STATE_SR_LABEL: Record<TabBarResultState, string> = {
  clean: '(all clear)',
  issues: '(issues found)',
  warnings: '(warnings)',
  complete: '(has results)',
};

export interface TabBarProps<T extends string = string> {
  /** Single row of tabs. Mutually exclusive with `rows`. */
  tabs?: ReadonlyArray<TabBarTab<T>>;
  /** Multiple rows of tabs rendered inside one tray. Mutually exclusive with `tabs`. */
  rows?: ReadonlyArray<ReadonlyArray<TabBarTab<T>>>;
  activeTab: T;
  onTabChange: (id: T) => void;
  ariaLabel: string;
  idPrefix?: string;
  layout?: 'equal' | 'auto';
  className?: string;
}

export function TabBar<T extends string = string>({
  tabs,
  rows,
  activeTab,
  onTabChange,
  ariaLabel,
  idPrefix = 'tab',
  layout = 'equal',
  className,
}: TabBarProps<T>) {
  const tabRows: ReadonlyArray<ReadonlyArray<TabBarTab<T>>> = rows ?? (tabs ? [tabs] : []);

  return (
    <div
      className={cn('p-1.5 rounded-2xl space-y-1', className)}
      style={{
        backgroundColor: 'var(--color-surface-sunken)',
        boxShadow: `
          2px 4px 8px 0 rgba(var(--color-primary-rgb), 0.08),
          inset 2px 2px 6px rgba(255, 255, 255, 0.5),
          inset -2px -2px 6px rgba(0, 0, 0, 0.06)
        `,
      }}
      role="tablist"
      aria-label={ariaLabel}
    >
      {tabRows.map((row, rowIdx) => (
        <div
          key={rowIdx}
          className={cn('flex gap-1', layout === 'auto' && 'flex-wrap')}
        >
          {row.map((tab) => {
            const isActive = activeTab === tab.id;
            return (
              <button
                key={tab.id}
                type="button"
                role="tab"
                id={`${idPrefix}-${tab.id}`}
                aria-selected={isActive}
                aria-controls={`${idPrefix}panel-${tab.id}`}
                onClick={() => onTabChange(tab.id)}
                title={tab.description}
                className={cn(
                  'relative px-3 py-2.5 text-sm font-medium rounded-xl',
                  'flex items-center justify-center gap-2',
                  'transition-all duration-400 ease-out',
                  'focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-[var(--color-primary)] focus-visible:ring-offset-2',
                  layout === 'equal' && 'flex-1',
                  isActive
                    ? 'text-[var(--color-primary-dark)] font-semibold'
                    : 'text-[var(--color-text-secondary)] hover:text-[var(--color-text-primary)]',
                )}
                style={
                  isActive
                    ? {
                        backgroundColor: 'var(--color-surface-elevated)',
                        boxShadow: `
                          1px 2px 4px 0 rgba(var(--color-primary-rgb), 0.12),
                          2px 4px 10px 0 rgba(var(--color-primary-rgb), 0.08),
                          inset 1px 1px 4px rgba(255, 255, 255, 0.6),
                          inset -1px -1px 4px rgba(var(--color-primary-rgb), 0.06)
                        `,
                      }
                    : { backgroundColor: 'transparent', boxShadow: 'none' }
                }
              >
                <tab.icon
                  className={cn(
                    'w-4 h-4 transition-transform duration-300',
                    isActive && 'scale-110',
                  )}
                />
                <span className="font-display">{tab.label}</span>
                {tab.hasResult && !isActive && (
                  <>
                    <span
                      aria-hidden="true"
                      className={cn(
                        'absolute top-1.5 right-1.5 w-1.5 h-1.5 rounded-full',
                        RESULT_STATE_DOT_CLASS[tab.resultState ?? 'complete'],
                      )}
                    />
                    <span className="sr-only">
                      {RESULT_STATE_SR_LABEL[tab.resultState ?? 'complete']}
                    </span>
                  </>
                )}
              </button>
            );
          })}
        </div>
      ))}
    </div>
  );
}
