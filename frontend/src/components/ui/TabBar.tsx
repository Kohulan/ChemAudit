import type { ComponentType } from 'react';
import { cn } from '../../lib/utils';

export interface TabBarTab<T extends string = string> {
  id: T;
  label: string;
  icon: ComponentType<{ className?: string }>;
  description?: string;
}

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
              </button>
            );
          })}
        </div>
      ))}
    </div>
  );
}
