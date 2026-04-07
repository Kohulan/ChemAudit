import { useState, useMemo, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import type { IssueEntry } from '../../types/dataset_intelligence';
import { ClayButton } from '../ui/ClayButton';

// =============================================================================
// Types
// =============================================================================

interface TreemapDrillDownProps {
  /** Currently selected issue type (null = collapsed). */
  issueType: string | null;
  /** Full list of issues from health audit. */
  issues: IssueEntry[];
  /** Whether the panel is open. */
  isOpen: boolean;
}

// =============================================================================
// Constants
// =============================================================================

const PAGE_SIZE = 50;

// =============================================================================
// Component
// =============================================================================

/**
 * Expandable drill-down panel showing affected SMILES for a selected issue type.
 *
 * Per UI-SPEC animation contract:
 * - Expand: height 0 -> auto, 0.3s ease-out
 * - Collapse: height auto -> 0, 0.25s ease-in
 * - Re-clicking same cell collapses
 * - Clicking different cell swaps content (instant swap, no collapse/re-expand)
 *
 * Paginated at 50 rows per page with Previous/Next controls.
 */
export function TreemapDrillDown({ issueType, issues, isOpen }: TreemapDrillDownProps) {
  const [page, setPage] = useState(0);

  // Reset pagination when issue type changes
  useEffect(() => {
    setPage(0);
  }, [issueType]);

  // Filter issues by selected type
  const filtered = useMemo(() => {
    if (!issueType) return [];
    return issues.filter((i) => i.issue_type === issueType);
  }, [issues, issueType]);

  const totalPages = Math.ceil(filtered.length / PAGE_SIZE);
  const start = page * PAGE_SIZE;
  const end = Math.min(start + PAGE_SIZE, filtered.length);
  const pageItems = filtered.slice(start, end);

  return (
    <AnimatePresence mode="wait">
      {isOpen && issueType && filtered.length > 0 && (
        <motion.div
          key={issueType}
          initial={{ height: 0, opacity: 0 }}
          animate={{ height: 'auto', opacity: 1 }}
          exit={{ height: 0, opacity: 0 }}
          transition={{
            height: { duration: 0.3, ease: 'easeOut' },
            opacity: { duration: 0.2 },
          }}
          className="overflow-hidden"
        >
          <div className="border border-[var(--color-border)] rounded-xl bg-[var(--color-surface)] p-4 mt-4">
            {/* Heading */}
            <h4 className="text-sm font-semibold font-display text-[var(--color-text-primary)] mb-3">
              Affected Molecules: {issueType}
            </h4>

            {/* Announcer */}
            <div aria-live="polite" className="sr-only">
              Showing {filtered.length} molecules with {issueType}
            </div>

            {/* Table */}
            <div className="overflow-x-auto">
              <table className="w-full text-sm">
                <thead>
                  <tr className="border-b border-[var(--color-border)]">
                    <th className="text-left px-3 py-2 text-xs font-semibold text-[var(--color-text-muted)] w-[60px]">
                      Row
                    </th>
                    <th className="text-left px-3 py-2 text-xs font-semibold text-[var(--color-text-muted)]">
                      SMILES
                    </th>
                    <th className="text-left px-3 py-2 text-xs font-semibold text-[var(--color-text-muted)] w-[100px]">
                      Severity
                    </th>
                    <th className="text-left px-3 py-2 text-xs font-semibold text-[var(--color-text-muted)]">
                      Detail
                    </th>
                  </tr>
                </thead>
                <tbody>
                  {pageItems.map((issue, idx) => (
                    <tr
                      key={`${issue.row_index}-${idx}`}
                      className="border-b border-[var(--color-border)] last:border-0 hover:bg-[var(--color-surface-sunken)] transition-colors"
                    >
                      <td className="px-3 py-2 text-xs text-[var(--color-text-secondary)] tabular-nums">
                        {issue.row_index}
                      </td>
                      <td className="px-3 py-2">
                        <span
                          className="text-xs font-mono text-[var(--color-text-primary)] block max-w-[300px] truncate"
                          title={issue.smiles}
                        >
                          {issue.smiles}
                        </span>
                      </td>
                      <td className="px-3 py-2">
                        <SeverityBadge severity={issue.severity} />
                      </td>
                      <td className="px-3 py-2 text-xs text-[var(--color-text-secondary)]">
                        {issue.description}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>

            {/* Pagination */}
            {totalPages > 1 && (
              <div className="flex items-center justify-between mt-3 pt-3 border-t border-[var(--color-border)]">
                <span className="text-xs text-[var(--color-text-secondary)]">
                  Showing {start + 1}-{end} of {filtered.length} affected molecules
                </span>
                <div className="flex gap-2">
                  <ClayButton
                    variant="ghost"
                    size="sm"
                    disabled={page === 0}
                    onClick={() => setPage((p) => Math.max(0, p - 1))}
                    className="text-xs"
                  >
                    Previous
                  </ClayButton>
                  <ClayButton
                    variant="ghost"
                    size="sm"
                    disabled={page >= totalPages - 1}
                    onClick={() => setPage((p) => Math.min(totalPages - 1, p + 1))}
                    className="text-xs"
                  >
                    Next
                  </ClayButton>
                </div>
              </div>
            )}
          </div>
        </motion.div>
      )}
    </AnimatePresence>
  );
}

// =============================================================================
// Severity badge helper
// =============================================================================

function SeverityBadge({ severity }: { severity: string }) {
  const colorMap: Record<string, string> = {
    critical: 'bg-red-500/15 text-red-600 dark:text-red-400',
    error: 'bg-red-500/10 text-red-600 dark:text-red-400',
    warning: 'bg-amber-500/10 text-amber-600 dark:text-amber-400',
    info: 'bg-blue-500/10 text-blue-600 dark:text-blue-400',
  };
  const classes = colorMap[severity.toLowerCase()] ?? 'bg-gray-500/10 text-gray-600 dark:text-gray-400';

  return (
    <span className={`inline-flex px-2 py-0.5 rounded-full text-xs font-medium ${classes}`}>
      {severity}
    </span>
  );
}
