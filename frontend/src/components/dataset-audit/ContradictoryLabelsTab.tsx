import { useState, useMemo, useCallback } from 'react';
import { ChevronLeft, ChevronRight } from 'lucide-react';
import type {
  ContradictoryLabelResult,
  NumericColumnInfo,
} from '../../types/dataset_intelligence';
import {
  FOLD_THRESHOLD_OPTIONS,
  DEFAULT_FOLD_THRESHOLD,
} from '../../types/dataset_intelligence';
import { ActivityColumnSelector } from './ActivityColumnSelector';
import { ContradictionCard } from './ContradictionCard';
import { ClayButton } from '../ui/ClayButton';

// =============================================================================
// Constants
// =============================================================================

const CARDS_PER_PAGE = 20;

// =============================================================================
// Types
// =============================================================================

interface ContradictoryLabelsTabProps {
  /** All contradictory label results from the audit. */
  contradictions: ContradictoryLabelResult[];
  /** Detected numeric columns for activity selection. */
  numericColumns: NumericColumnInfo[];
}

// =============================================================================
// Helpers
// =============================================================================

/**
 * Get the default selected column from the numeric columns list.
 * Pre-selects: first priority-1, then priority-2, then priority-3.
 */
function getDefaultColumn(columns: NumericColumnInfo[]): string {
  if (columns.length === 0) return '';
  const sorted = [...columns].sort((a, b) => a.priority - b.priority);
  return sorted[0].name;
}

// =============================================================================
// Component
// =============================================================================

/**
 * Contradictory Labels tab content for the Dataset Audit page.
 *
 * Shows an activity column selector (auto-detected with priority heuristic),
 * a fold-threshold dropdown (3x/5x/10x/50x/100x), a count summary, and
 * a paginated list of ContradictionCards sorted by fold-difference descending.
 *
 * Client-side filtering: the backend computes contradictions for the default
 * threshold (10x). For thresholds < 10x, all contradictions are available.
 * For thresholds > 10x, results are filtered client-side.
 */
export function ContradictoryLabelsTab({
  contradictions,
  numericColumns,
}: ContradictoryLabelsTabProps) {
  // State
  const [selectedColumn, setSelectedColumn] = useState<string>(() =>
    getDefaultColumn(numericColumns),
  );
  const [foldThreshold, setFoldThreshold] = useState<number>(
    DEFAULT_FOLD_THRESHOLD,
  );
  const [expandedInChIKey, setExpandedInChIKey] = useState<string | null>(null);
  const [currentPage, setCurrentPage] = useState<number>(0);

  // Filter contradictions by fold-difference threshold
  const filtered = useMemo(() => {
    return contradictions
      .filter((c) => c.fold_difference >= foldThreshold)
      .sort((a, b) => b.fold_difference - a.fold_difference);
  }, [contradictions, foldThreshold]);

  // Pagination
  const totalPages = Math.max(1, Math.ceil(filtered.length / CARDS_PER_PAGE));
  const pageItems = useMemo(() => {
    const start = currentPage * CARDS_PER_PAGE;
    return filtered.slice(start, start + CARDS_PER_PAGE);
  }, [filtered, currentPage]);

  // Handlers
  const handleColumnChange = useCallback((column: string) => {
    setSelectedColumn(column);
    setCurrentPage(0);
    setExpandedInChIKey(null);
  }, []);

  const handleThresholdChange = useCallback(
    (e: React.ChangeEvent<HTMLSelectElement>) => {
      setFoldThreshold(Number(e.target.value));
      setCurrentPage(0);
      setExpandedInChIKey(null);
    },
    [],
  );

  const handleToggle = useCallback(
    (inchikey: string) => {
      setExpandedInChIKey((prev) => (prev === inchikey ? null : inchikey));
    },
    [],
  );

  // Handle no-activity-column case
  if (numericColumns.length === 0) {
    return (
      <div className="min-h-[200px] flex items-center justify-center">
        <p className="text-sm text-[var(--color-text-muted)]">
          No numeric activity columns detected. Upload a CSV with numeric
          activity data to detect contradictions.
        </p>
      </div>
    );
  }

  return (
    <div className="min-h-[200px] space-y-6">
      {/* Top row: Activity column selector + Fold threshold dropdown */}
      <div className="flex flex-wrap items-end gap-4">
        <ActivityColumnSelector
          columns={numericColumns}
          selectedColumn={selectedColumn}
          onColumnChange={handleColumnChange}
        />

        {/* Fold-Difference Threshold dropdown */}
        <div className="flex flex-col gap-1">
          <label
            htmlFor="fold-threshold-select"
            className="text-xs font-medium text-[var(--color-text-secondary)]"
          >
            Fold-Difference Threshold
          </label>
          <select
            id="fold-threshold-select"
            value={foldThreshold}
            onChange={handleThresholdChange}
            className={[
              'px-3 py-2 text-sm rounded-lg',
              'bg-[var(--color-surface)] text-[var(--color-text-primary)]',
              'border border-[var(--color-border)]',
              'focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)] focus:border-transparent',
              'cursor-pointer',
            ].join(' ')}
          >
            {FOLD_THRESHOLD_OPTIONS.map((opt) => (
              <option key={opt} value={opt}>
                {opt}x
              </option>
            ))}
          </select>
        </div>
      </div>

      {/* Count summary */}
      <p className="text-sm text-[var(--color-text-secondary)]">
        {filtered.length > 0
          ? `${filtered.length} contradictions found at ${foldThreshold}x threshold`
          : `No contradictory labels detected at ${foldThreshold}x threshold.`}
      </p>

      {/* Contradiction card list */}
      {pageItems.length > 0 && (
        <div className="space-y-3">
          {pageItems.map((c) => (
            <ContradictionCard
              key={c.inchikey}
              contradiction={c}
              isExpanded={expandedInChIKey === c.inchikey}
              onToggle={() => handleToggle(c.inchikey)}
            />
          ))}
        </div>
      )}

      {/* Pagination controls */}
      {totalPages > 1 && (
        <div className="flex items-center justify-between pt-2">
          <span className="text-xs text-[var(--color-text-muted)]">
            {CARDS_PER_PAGE} per page
          </span>

          <div className="flex items-center gap-3">
            <ClayButton
              variant="ghost"
              size="icon"
              onClick={() => setCurrentPage((p) => Math.max(0, p - 1))}
              disabled={currentPage <= 0}
              aria-label="Previous page"
            >
              <ChevronLeft className="w-4 h-4" />
            </ClayButton>

            <span className="text-xs text-[var(--color-text-secondary)]">
              Page {currentPage + 1} of {totalPages}
            </span>

            <ClayButton
              variant="ghost"
              size="icon"
              onClick={() =>
                setCurrentPage((p) => Math.min(totalPages - 1, p + 1))
              }
              disabled={currentPage >= totalPages - 1}
              aria-label="Next page"
            >
              <ChevronRight className="w-4 h-4" />
            </ClayButton>
          </div>
        </div>
      )}
    </div>
  );
}
