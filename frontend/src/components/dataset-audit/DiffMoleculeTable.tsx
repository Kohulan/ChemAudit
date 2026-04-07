import { useState, useMemo, useCallback } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown, ChevronLeft, ChevronRight } from 'lucide-react';
import type { DiffMolecule, DatasetDiffResults } from '../../types/dataset_intelligence';
import { ClayButton } from '../ui/ClayButton';

// =============================================================================
// Constants
// =============================================================================

const ROWS_PER_PAGE = 50;

type DiffCategory = 'all' | 'added' | 'removed' | 'modified';

const CATEGORY_FILTERS: Array<{ value: DiffCategory; label: string }> = [
  { value: 'all', label: 'All' },
  { value: 'added', label: 'Added' },
  { value: 'removed', label: 'Removed' },
  { value: 'modified', label: 'Modified' },
];

// =============================================================================
// Types
// =============================================================================

interface DiffMoleculeTableProps {
  /** All diff molecules (will be filtered by category). */
  molecules: DiffMolecule[];
  /** Currently active category filter. */
  category: DiffCategory;
  /** Full diff results for column mismatch notice. */
  diffResults: DatasetDiffResults;
}

/** Internal type combining molecule with its category. */
interface TaggedMolecule extends DiffMolecule {
  _category: 'added' | 'removed' | 'modified';
}

// =============================================================================
// Helpers
// =============================================================================

function getCategoryBadge(category: 'added' | 'removed' | 'modified'): {
  label: string;
  classes: string;
} {
  switch (category) {
    case 'added':
      return {
        label: 'Added',
        classes: 'bg-green-100 text-green-700 dark:bg-green-900/30 dark:text-green-400',
      };
    case 'removed':
      return {
        label: 'Removed',
        classes: 'bg-red-100 text-red-700 dark:bg-red-900/30 dark:text-red-400',
      };
    case 'modified':
      return {
        label: 'Modified',
        classes: 'bg-amber-100 text-amber-700 dark:bg-amber-900/30 dark:text-amber-400',
      };
  }
}

/**
 * Truncate a string in the middle for display.
 */
function truncateMiddle(s: string, maxLen: number): string {
  if (s.length <= maxLen) return s;
  const half = Math.floor((maxLen - 3) / 2);
  return `${s.slice(0, half)}...${s.slice(-half)}`;
}

// =============================================================================
// Component
// =============================================================================

/**
 * Filterable, paginated table of diff molecules with expandable modified rows.
 *
 * Per UI-SPEC D-13:
 * - Category filter tabs: All | Added | Removed | Modified (radio group)
 * - Table columns: InChIKey (mono, truncated), SMILES (mono, truncated), Row #, Category badge
 * - Modified rows expandable to show per-property changes (Column | Old | New)
 * - Paginated at 50 rows per page
 * - Column mismatch notice when unique columns exist in either dataset
 */
export function DiffMoleculeTable({
  molecules: _molecules,
  category: _initialCategory,
  diffResults,
}: DiffMoleculeTableProps) {
  const [category, setCategory] = useState<DiffCategory>(_initialCategory);
  const [currentPage, setCurrentPage] = useState<number>(0);
  const [expandedKey, setExpandedKey] = useState<string | null>(null);

  // Build tagged molecule list combining all categories
  const allTagged = useMemo<TaggedMolecule[]>(() => {
    const tagged: TaggedMolecule[] = [];
    for (const m of diffResults.added) {
      tagged.push({ ...m, _category: 'added' });
    }
    for (const m of diffResults.removed) {
      tagged.push({ ...m, _category: 'removed' });
    }
    for (const m of diffResults.modified) {
      tagged.push({ ...m, _category: 'modified' });
    }
    return tagged;
  }, [diffResults]);

  // Filter by selected category
  const filtered = useMemo(() => {
    if (category === 'all') return allTagged;
    return allTagged.filter((m) => m._category === category);
  }, [allTagged, category]);

  // Pagination
  const totalPages = Math.max(1, Math.ceil(filtered.length / ROWS_PER_PAGE));
  const pageItems = useMemo(() => {
    const start = currentPage * ROWS_PER_PAGE;
    return filtered.slice(start, start + ROWS_PER_PAGE);
  }, [filtered, currentPage]);

  const handleCategoryChange = useCallback((cat: DiffCategory) => {
    setCategory(cat);
    setCurrentPage(0);
    setExpandedKey(null);
  }, []);

  const handleToggleExpand = useCallback((key: string) => {
    setExpandedKey((prev) => (prev === key ? null : key));
  }, []);

  const uniqueColumnsTotal =
    diffResults.unique_columns_primary + diffResults.unique_columns_comparison;

  return (
    <div className="space-y-4">
      {/* Column mismatch notice */}
      {uniqueColumnsTotal > 0 && (
        <div className="px-4 py-3 rounded-lg bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800/40">
          <p className="text-xs text-blue-700 dark:text-blue-300">
            Note: Only columns present in both datasets are compared.{' '}
            {uniqueColumnsTotal} columns are unique to each file.
          </p>
        </div>
      )}

      {/* Category filter radio group */}
      <div
        className="flex gap-1 p-1 rounded-lg bg-[var(--color-surface-sunken)]"
        role="radiogroup"
        aria-label="Filter diff results by category"
      >
        {CATEGORY_FILTERS.map((f) => {
          const isActive = category === f.value;
          return (
            <button
              key={f.value}
              role="radio"
              aria-checked={isActive}
              onClick={() => handleCategoryChange(f.value)}
              className={[
                'px-3 py-1.5 text-xs font-medium rounded-md transition-colors duration-150',
                isActive
                  ? 'bg-[var(--color-surface-elevated)] text-[var(--color-text-primary)] shadow-sm'
                  : 'text-[var(--color-text-secondary)] hover:text-[var(--color-text-primary)]',
              ].join(' ')}
            >
              {f.label}
            </button>
          );
        })}
      </div>

      {/* Table */}
      <div className="overflow-hidden rounded-lg border border-[var(--color-border)]">
        <div className="overflow-x-auto">
          <table className="w-full text-sm">
            <thead>
              <tr className="bg-[var(--color-surface-sunken)] border-b border-[var(--color-border)]">
                <th className="px-4 py-3 text-left text-xs font-semibold text-[var(--color-text-secondary)] w-12">
                  #
                </th>
                <th className="px-4 py-3 text-left text-xs font-semibold text-[var(--color-text-secondary)]">
                  InChIKey
                </th>
                <th className="px-4 py-3 text-left text-xs font-semibold text-[var(--color-text-secondary)]">
                  SMILES
                </th>
                <th className="px-4 py-3 text-left text-xs font-semibold text-[var(--color-text-secondary)]">
                  Category
                </th>
                {/* Expand chevron column */}
                <th className="w-8" />
              </tr>
            </thead>
            <tbody>
              {pageItems.map((mol) => {
                const rowKey = `${mol._category}-${mol.inchikey}-${mol.row_index}`;
                const isExpanded = expandedKey === rowKey;
                const badge = getCategoryBadge(mol._category);
                const isModified = mol._category === 'modified';

                return (
                  <MoleculeRow
                    key={rowKey}
                    mol={mol}
                    rowKey={rowKey}
                    isExpanded={isExpanded}
                    isModified={isModified}
                    badge={badge}
                    onToggle={() => handleToggleExpand(rowKey)}
                  />
                );
              })}
              {pageItems.length === 0 && (
                <tr>
                  <td
                    colSpan={5}
                    className="px-4 py-8 text-center text-sm text-[var(--color-text-muted)]"
                  >
                    No molecules in this category.
                  </td>
                </tr>
              )}
            </tbody>
          </table>
        </div>

        {/* Pagination footer */}
        {totalPages > 1 && (
          <div className="flex items-center justify-between px-4 py-3 border-t border-[var(--color-border)] bg-[var(--color-surface-sunken)]/50">
            <span className="text-xs text-[var(--color-text-muted)]">
              {ROWS_PER_PAGE} per page
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
    </div>
  );
}

// =============================================================================
// MoleculeRow sub-component
// =============================================================================

function MoleculeRow({
  mol,
  rowKey,
  isExpanded,
  isModified,
  badge,
  onToggle,
}: {
  mol: TaggedMolecule;
  rowKey: string;
  isExpanded: boolean;
  isModified: boolean;
  badge: { label: string; classes: string };
  onToggle: () => void;
}) {
  return (
    <>
      <tr
        className={[
          'border-t border-[var(--color-border)]/50',
          isModified ? 'cursor-pointer hover:bg-[var(--color-surface-sunken)]/50' : '',
        ].join(' ')}
        onClick={isModified ? onToggle : undefined}
        aria-expanded={isModified ? isExpanded : undefined}
      >
        <td className="px-4 py-3 text-xs text-[var(--color-text-muted)] tabular-nums">
          {mol.row_index}
        </td>
        <td
          className="px-4 py-3 text-xs font-mono text-[var(--color-text-primary)] truncate max-w-[180px]"
          title={mol.inchikey}
        >
          {truncateMiddle(mol.inchikey, 20)}
        </td>
        <td
          className="px-4 py-3 text-xs font-mono text-[var(--color-text-primary)] truncate max-w-[200px]"
          title={mol.smiles}
        >
          {truncateMiddle(mol.smiles, 30)}
        </td>
        <td className="px-4 py-3">
          <span
            className={`inline-flex items-center px-2 py-0.5 text-xs font-medium rounded-full ${badge.classes}`}
          >
            {badge.label}
          </span>
        </td>
        <td className="px-2 py-3">
          {isModified && (
            <motion.span
              animate={{ rotate: isExpanded ? 180 : 0 }}
              transition={{ duration: 0.2 }}
            >
              <ChevronDown className="w-4 h-4 text-[var(--color-text-muted)]" />
            </motion.span>
          )}
        </td>
      </tr>

      {/* Expandable property changes for modified rows */}
      <AnimatePresence initial={false}>
        {isExpanded && isModified && mol.changes.length > 0 && (
          <tr>
            <td colSpan={5} className="p-0">
              <motion.div
                key={`${rowKey}-detail`}
                initial={{ height: 0, opacity: 0 }}
                animate={{ height: 'auto', opacity: 1 }}
                exit={{ height: 0, opacity: 0 }}
                transition={{ duration: 0.2, ease: 'easeOut' }}
                className="overflow-hidden"
              >
                <div className="px-8 py-3 bg-[var(--color-surface-sunken)]/30">
                  <table className="w-full text-xs">
                    <thead>
                      <tr className="text-[var(--color-text-secondary)]">
                        <th className="text-left py-1 pr-4 font-medium">
                          Column Name
                        </th>
                        <th className="text-left py-1 pr-4 font-medium">
                          Old Value
                        </th>
                        <th className="text-left py-1 font-medium">
                          New Value
                        </th>
                      </tr>
                    </thead>
                    <tbody>
                      {mol.changes.map((change) => (
                        <tr
                          key={change.column}
                          className="border-t border-[var(--color-border)]/30"
                        >
                          <td className="py-1.5 pr-4 font-medium text-[var(--color-text-primary)]">
                            {change.column}
                          </td>
                          <td className="py-1.5 pr-4 font-mono text-red-600 dark:text-red-400">
                            {String(change.old_value ?? '')}
                          </td>
                          <td className="py-1.5 font-mono text-green-600 dark:text-green-400">
                            {String(change.new_value ?? '')}
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </motion.div>
            </td>
          </tr>
        )}
      </AnimatePresence>
    </>
  );
}
