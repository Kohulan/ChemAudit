import { useState, useMemo, useCallback } from 'react';
import { motion } from 'framer-motion';
import { ChevronLeft, ChevronRight, AlertTriangle, TrendingUp, Layers } from 'lucide-react';
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
import { InfoTooltip } from '../ui/Tooltip';

// =============================================================================
// Constants
// =============================================================================

const CARDS_PER_PAGE = 20;
const SEVERITY_LABEL = '\u2265'; // ≥ symbol

// =============================================================================
// Summary stat card
// =============================================================================

function SummaryStatCard({ icon, label, value, color, bgColor, delay = 0 }: {
  icon: React.ReactNode;
  label: string;
  value: number | string;
  color: string;
  bgColor: string;
  delay?: number;
}) {
  return (
    <motion.div
      initial={{ opacity: 0, y: 12, scale: 0.95 }}
      animate={{ opacity: 1, y: 0, scale: 1 }}
      transition={{ delay, duration: 0.35, ease: [0.4, 0, 0.2, 1] }}
      className={`relative flex items-center gap-3 px-4 py-3.5 rounded-2xl overflow-hidden ${bgColor}`}
      style={{
        boxShadow: `
          2px 4px 8px 0 rgba(0, 0, 0, 0.06),
          4px 8px 12px 0 rgba(0, 0, 0, 0.04),
          inset 2px 2px 5px rgba(255, 255, 255, 0.6),
          inset -2px -2px 5px rgba(0, 0, 0, 0.04)
        `,
      }}
    >
      <div className={`relative ${color}`}>{icon}</div>
      <div className="relative">
        <p className={`text-lg font-bold tabular-nums ${color}`}>{value}</p>
        <p className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wide font-medium">{label}</p>
      </div>
    </motion.div>
  );
}

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
      {/* Info header */}
      <div className="flex items-start gap-2 px-4 py-3 rounded-xl bg-[var(--color-surface-sunken)] border border-[var(--color-border)]">
        <InfoTooltip
          title="Contradictory Labels Detection"
          content={
            <div className="text-xs space-y-1">
              <p>Identifies molecules that appear multiple times (same InChIKey) but with significantly different activity values.</p>
              <p className="mt-1"><strong>How it works:</strong></p>
              <ul className="text-white/70 space-y-0.5">
                <li>1. Groups molecules by InChIKey (structural identity)</li>
                <li>2. For each group with 2+ entries, computes the fold-difference (max / min activity)</li>
                <li>3. Flags groups where fold-difference exceeds the threshold</li>
              </ul>
              <p className="mt-1 text-white/60">Select an activity column (e.g., IC50, pIC50, EC50) &mdash; not a molecular descriptor. Adjust the fold-threshold to control sensitivity: 3x catches subtle discrepancies, 100x catches only gross errors.</p>
            </div>
          }
          size="small"
        />
        <p className="text-xs text-[var(--color-text-secondary)] leading-relaxed">
          Detects same-structure molecules (by InChIKey) with conflicting activity values.
          Select a <strong>numeric activity column</strong> (e.g., IC50, pIC50, EC50) and a fold-difference threshold to control sensitivity.
        </p>
      </div>

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

      {/* Summary stats */}
      {filtered.length > 0 ? (
        <motion.div
          initial={{ opacity: 0, y: 8 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.3 }}
          className="grid grid-cols-2 lg:grid-cols-4 gap-3"
        >
          <SummaryStatCard
            icon={<AlertTriangle className="w-4 h-4" />}
            label="Contradictions"
            value={filtered.length}
            color="text-[var(--color-primary)]"
            bgColor="bg-[var(--color-surface-elevated)]"
            delay={0}
          />
          <SummaryStatCard
            icon={<TrendingUp className="w-4 h-4" />}
            label="Max Fold"
            value={`${filtered[0].fold_difference.toFixed(0)}x`}
            color="text-[var(--color-accent-dark)]"
            bgColor="bg-[var(--color-surface-elevated)]"
            delay={0.07}
          />
          <SummaryStatCard
            icon={<Layers className="w-4 h-4" />}
            label="Total Entries"
            value={filtered.reduce((sum, c) => sum + c.entry_count, 0)}
            color="text-[var(--color-text-primary)]"
            bgColor="bg-[var(--color-surface-elevated)]"
            delay={0.14}
          />
          <SummaryStatCard
            icon={<AlertTriangle className="w-4 h-4" />}
            label={`Critical (${SEVERITY_LABEL}100x)`}
            value={filtered.filter((c) => c.fold_difference >= 100).length}
            color="text-[var(--color-primary-dark)]"
            bgColor="bg-[var(--color-surface-elevated)]"
            delay={0.21}
          />
        </motion.div>
      ) : (
        <motion.div
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          className="flex flex-col items-center justify-center py-10 text-center"
        >
          <div className="w-12 h-12 rounded-2xl bg-emerald-50 dark:bg-emerald-950/20 flex items-center justify-center mb-3">
            <AlertTriangle className="w-6 h-6 text-emerald-500" />
          </div>
          <p className="text-sm font-medium text-[var(--color-text-primary)]">
            No contradictory labels detected
          </p>
          <p className="text-xs text-[var(--color-text-muted)] mt-1">
            at {foldThreshold}x fold-difference threshold
          </p>
        </motion.div>
      )}

      {/* Contradiction card list */}
      {pageItems.length > 0 && (
        <div className="space-y-3">
          {pageItems.map((c, i) => (
            <ContradictionCard
              key={c.inchikey}
              contradiction={c}
              isExpanded={expandedInChIKey === c.inchikey}
              onToggle={() => handleToggle(c.inchikey)}
              index={i}
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
