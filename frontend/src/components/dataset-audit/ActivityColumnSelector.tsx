import type { NumericColumnInfo } from '../../types/dataset_intelligence';

// =============================================================================
// Types
// =============================================================================

interface ActivityColumnSelectorProps {
  /** Detected numeric columns with priority levels. */
  columns: NumericColumnInfo[];
  /** Currently selected column name. */
  selectedColumn: string;
  /** Callback when column selection changes. */
  onColumnChange: (column: string) => void;
}

// =============================================================================
// Component
// =============================================================================

/**
 * Dropdown for selecting the activity column used in contradictory label detection.
 *
 * Groups columns by priority:
 * - Priority 1: Likely activity columns (exact name matches like "IC50", "pIC50")
 * - Priority 2: Possible activity columns (substring matches)
 * - Priority 3: Other numeric columns
 *
 * Pre-selects the first priority-1 column, falling back to priority-2, then priority-3.
 * Shows a message when no numeric columns are detected.
 */
export function ActivityColumnSelector({
  columns,
  selectedColumn,
  onColumnChange,
}: ActivityColumnSelectorProps) {
  if (columns.length === 0) {
    return (
      <div className="text-sm text-[var(--color-text-muted)]">
        No numeric activity columns detected. Upload a CSV with numeric activity
        data to detect contradictions.
      </div>
    );
  }

  // Group columns by priority
  const priority1 = columns.filter((c) => c.priority === 1);
  const priority2 = columns.filter((c) => c.priority === 2);
  const priority3 = columns.filter((c) => c.priority >= 3);

  return (
    <div className="flex flex-col gap-1">
      <label
        htmlFor="activity-column-select"
        className="text-xs font-medium text-[var(--color-text-secondary)]"
      >
        Activity Column
      </label>
      <select
        id="activity-column-select"
        value={selectedColumn}
        onChange={(e) => onColumnChange(e.target.value)}
        className={[
          'px-3 py-2 text-sm rounded-lg',
          'bg-[var(--color-surface)] text-[var(--color-text-primary)]',
          'border border-[var(--color-border)]',
          'focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)] focus:border-transparent',
          'cursor-pointer',
        ].join(' ')}
      >
        {priority1.length > 0 && (
          <optgroup label="Likely activity columns">
            {priority1.map((col) => (
              <option key={col.name} value={col.name}>
                {col.name}
              </option>
            ))}
          </optgroup>
        )}
        {priority2.length > 0 && (
          <optgroup label="Possible activity columns">
            {priority2.map((col) => (
              <option key={col.name} value={col.name}>
                {col.name}
              </option>
            ))}
          </optgroup>
        )}
        {priority3.length > 0 && (
          <optgroup label="Other numeric columns">
            {priority3.map((col) => (
              <option key={col.name} value={col.name}>
                {col.name}
              </option>
            ))}
          </optgroup>
        )}
      </select>
    </div>
  );
}
