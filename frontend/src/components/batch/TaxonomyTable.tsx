/**
 * TaxonomyTable
 *
 * Paginated table of ALL taxonomy categories (not just top 20).
 * Filterable by search query, sorted by count descending.
 * Shows category name, count, example SMILES, and description.
 */

import React, { useMemo, useState } from 'react';
import { ChevronLeft, ChevronRight } from 'lucide-react';
import type { TaxonomyResult } from '../../types/analytics';

interface TaxonomyTableProps {
  /** Full taxonomy result with per_molecule and category_counts */
  taxonomyResult: TaxonomyResult;
  /** Search query for filtering categories (case-insensitive) */
  searchQuery: string;
}

const PAGE_SIZE = 50;

interface CategoryRow {
  name: string;
  count: number;
  exampleSmiles: string;
  description: string;
}

export function TaxonomyTable({ taxonomyResult, searchQuery }: TaxonomyTableProps) {
  const [currentPage, setCurrentPage] = useState(0);

  // Build rows from category_counts, enriched with per_molecule data
  const allRows: CategoryRow[] = useMemo(() => {
    const rows: CategoryRow[] = [];

    for (const [categoryName, count] of Object.entries(taxonomyResult.category_counts)) {
      // Find first molecule that has this category
      let exampleSmiles = '';
      let description = '';

      for (const mol of taxonomyResult.per_molecule) {
        const match = mol.categories.find(
          (c) => c.name === categoryName || c.category === categoryName,
        );
        if (match) {
          exampleSmiles = mol.smiles;
          description = match.description;
          break;
        }
      }

      rows.push({ name: categoryName, count, exampleSmiles, description });
    }

    // Sort by count descending
    rows.sort((a, b) => b.count - a.count);
    return rows;
  }, [taxonomyResult]);

  // Filter by search query
  const filteredRows = useMemo(() => {
    if (!searchQuery.trim()) return allRows;
    const q = searchQuery.toLowerCase();
    return allRows.filter((r) => r.name.toLowerCase().includes(q));
  }, [allRows, searchQuery]);

  // Reset page when search changes
  React.useEffect(() => {
    setCurrentPage(0);
  }, [searchQuery]);

  const totalPages = Math.max(1, Math.ceil(filteredRows.length / PAGE_SIZE));
  const pageRows = filteredRows.slice(
    currentPage * PAGE_SIZE,
    (currentPage + 1) * PAGE_SIZE,
  );

  if (filteredRows.length === 0 && searchQuery.trim()) {
    return (
      <p className="text-sm text-[var(--color-text-muted)] py-4 text-center">
        No categories match &apos;{searchQuery}&apos;.
      </p>
    );
  }

  return (
    <div className="space-y-3">
      <div className="overflow-x-auto">
        <table className="w-full border-collapse text-left">
          <thead>
            <tr className="border-b border-[var(--color-border)]">
              <th className="px-3 py-2 text-xs font-medium text-[var(--color-text-muted)] uppercase">
                Category
              </th>
              <th className="px-3 py-2 text-xs font-medium text-[var(--color-text-muted)] uppercase text-right">
                Count
              </th>
              <th className="px-3 py-2 text-xs font-medium text-[var(--color-text-muted)] uppercase">
                Example
              </th>
              <th className="px-3 py-2 text-xs font-medium text-[var(--color-text-muted)] uppercase">
                Description
              </th>
            </tr>
          </thead>
          <tbody className="divide-y divide-[var(--color-border)]">
            {pageRows.map((row) => (
              <tr
                key={row.name}
                className="hover:bg-[var(--color-surface-sunken)] transition-colors"
              >
                <td className="px-3 py-2 text-sm text-[var(--color-text-primary)]">
                  {row.name}
                </td>
                <td className="px-3 py-2 text-xs font-mono tabular-nums text-[var(--color-text-secondary)] text-right">
                  {row.count}
                </td>
                <td className="px-3 py-2 text-xs font-mono text-[var(--color-text-muted)] max-w-[200px] truncate" title={row.exampleSmiles}>
                  {row.exampleSmiles || '-'}
                </td>
                <td className="px-3 py-2 text-sm text-[var(--color-text-secondary)] max-w-[300px]">
                  {row.description || '-'}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {/* Pagination controls */}
      {totalPages > 1 && (
        <div className="flex items-center justify-between px-2">
          <span className="text-xs text-[var(--color-text-muted)]">
            {filteredRows.length} categories
          </span>
          <div className="flex items-center gap-2">
            <button
              onClick={() => setCurrentPage((p) => Math.max(0, p - 1))}
              disabled={currentPage === 0}
              className="p-1.5 rounded-lg hover:bg-[var(--color-surface-sunken)] disabled:opacity-40 disabled:cursor-not-allowed transition-colors text-[var(--color-text-muted)]"
              aria-label="Previous page"
            >
              <ChevronLeft className="w-4 h-4" />
            </button>
            <span className="text-xs text-[var(--color-text-muted)]">
              {currentPage + 1} / {totalPages}
            </span>
            <button
              onClick={() => setCurrentPage((p) => Math.min(totalPages - 1, p + 1))}
              disabled={currentPage >= totalPages - 1}
              className="p-1.5 rounded-lg hover:bg-[var(--color-surface-sunken)] disabled:opacity-40 disabled:cursor-not-allowed transition-colors text-[var(--color-text-muted)]"
              aria-label="Next page"
            >
              <ChevronRight className="w-4 h-4" />
            </button>
          </div>
        </div>
      )}
    </div>
  );
}
