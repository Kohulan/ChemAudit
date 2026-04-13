import { useState, useEffect } from 'react';
import { X, ChevronLeft, ChevronRight } from 'lucide-react';
import { AnimatePresence, motion } from 'framer-motion';
import { ClayButton } from '../ui/ClayButton';
import type { MoleculeResult } from '../../types/structure_filter';

// =============================================================================
// Types
// =============================================================================

interface StructureFilterResultsTableProps {
  molecules: MoleculeResult[];
  selectedStage: string | null;
  onClearFilter: () => void;
}

// =============================================================================
// Constants
// =============================================================================

/** Rows per page — per UI-SPEC pagination spec. */
const PAGE_SIZE = 100;

// =============================================================================
// Status chip helpers
// =============================================================================

function getStatusChipClass(status: MoleculeResult['status']): string {
  switch (status) {
    case 'passed':
      return 'bg-green-100 text-green-700 dark:bg-green-900/30 dark:text-green-400';
    case 'rejected':
      return 'bg-red-100 text-red-700 dark:bg-red-900/30 dark:text-red-400';
    case 'duplicate':
      return 'bg-amber-100 text-amber-700 dark:bg-amber-900/30 dark:text-amber-400';
    case 'error':
    default:
      return 'bg-stone-100 text-stone-600 dark:bg-stone-800/50 dark:text-stone-400';
  }
}

function statusLabel(status: MoleculeResult['status']): string {
  switch (status) {
    case 'passed':
      return 'Passed';
    case 'rejected':
      return 'Rejected';
    case 'duplicate':
      return 'Duplicate';
    case 'error':
      return 'Error';
  }
}

// =============================================================================
// Component
// =============================================================================

/**
 * Filterable, paginated table for Structure Filter results.
 *
 * Per UI-SPEC D-09 and results table spec:
 * - Columns: SMILES (mono, truncated 40 chars) | Status chip | Failed At | Rejection Reason
 * - Active filter chip above table when selectedStage is set
 * - Pagination at 100 rows per page
 * - AnimatePresence on table content (opacity 0→1, 0.3s ease-out)
 * - Accessibility: aria-live on caption region and filter chip
 */
export function StructureFilterResultsTable({
  molecules,
  selectedStage,
  onClearFilter,
}: StructureFilterResultsTableProps) {
  const [currentPage, setCurrentPage] = useState(1);

  // Reset pagination when filter changes
  useEffect(() => {
    setCurrentPage(1);
  }, [selectedStage]);

  // Apply stage filter
  const filteredMolecules = selectedStage
    ? molecules.filter((m) => m.failed_at === selectedStage)
    : molecules;

  const totalPages = Math.max(1, Math.ceil(filteredMolecules.length / PAGE_SIZE));
  const pageStart = (currentPage - 1) * PAGE_SIZE;
  const pageEnd = currentPage * PAGE_SIZE;
  const pageRows = filteredMolecules.slice(pageStart, pageEnd);

  const totalCount = filteredMolecules.length;
  const captionText = selectedStage
    ? `Filter results — ${totalCount} molecule${totalCount !== 1 ? 's' : ''}, filtered by ${selectedStage}`
    : `Filter results — ${totalCount} molecule${totalCount !== 1 ? 's' : ''}`;

  return (
    <div className="space-y-3">
      {/* Active filter chip */}
      <AnimatePresence>
        {selectedStage && (
          <motion.div
            key="filter-chip"
            initial={{ opacity: 0, y: -4 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -4 }}
            transition={{ duration: 0.2 }}
            aria-live="polite"
            className="flex items-center gap-2"
          >
            <span className="inline-flex items-center gap-1.5 px-3 py-1 rounded-full text-xs font-medium bg-[var(--color-primary)]/10 text-[var(--color-primary)] border border-[var(--color-primary)]/20">
              Showing: rejected at{' '}
              <span className="font-semibold">{selectedStage}</span>
              <button
                onClick={onClearFilter}
                aria-label={`Clear filter: rejected at ${selectedStage}`}
                className="ml-0.5 hover:opacity-70 transition-opacity"
              >
                <X className="w-3 h-3" />
              </button>
            </span>
          </motion.div>
        )}
      </AnimatePresence>

      {/* Table container */}
      <div className="overflow-hidden rounded-lg border border-[var(--color-border)]">
        <div className="overflow-x-auto">
          <table className="w-full text-sm">
            {/* Caption with aria-live for screen readers */}
            <caption
              aria-live="polite"
              className="sr-only"
            >
              {captionText}
            </caption>

            <thead>
              <tr className="bg-[var(--color-surface-sunken)] border-b border-[var(--color-border)]">
                <th className="px-4 py-3 text-left text-xs font-semibold text-[var(--color-text-secondary)]">
                  SMILES
                </th>
                <th className="px-4 py-3 text-left text-xs font-semibold text-[var(--color-text-secondary)] w-28">
                  Status
                </th>
                <th className="px-4 py-3 text-left text-xs font-semibold text-[var(--color-text-secondary)] w-36">
                  Failed At
                </th>
                <th className="px-4 py-3 text-left text-xs font-semibold text-[var(--color-text-secondary)]">
                  Rejection Reason
                </th>
              </tr>
            </thead>

            <tbody>
              {pageRows.map((molecule, idx) => {
                const displaySmiles = molecule.smiles.length > 40
                  ? `${molecule.smiles.slice(0, 40)}…`
                  : molecule.smiles;

                return (
                  <tr
                    key={`${pageStart + idx}-${molecule.smiles}`}
                    className="border-b border-[var(--color-border)]/50 hover:bg-[var(--color-surface-sunken)]/50 transition-colors"
                  >
                    {/* SMILES — monospace, truncated at 40 chars with title for full */}
                    <td className="px-4 py-3">
                      <span
                        className="font-mono text-sm text-[var(--color-text-primary)]"
                        title={molecule.smiles}
                      >
                        {displaySmiles}
                      </span>
                    </td>

                    {/* Status chip */}
                    <td className="px-4 py-3">
                      <span
                        className={`inline-flex items-center px-2 py-0.5 rounded-full text-xs font-medium ${getStatusChipClass(molecule.status)}`}
                        aria-label={`${molecule.smiles}: ${molecule.status}`}
                      >
                        {statusLabel(molecule.status)}
                      </span>
                    </td>

                    {/* Failed At */}
                    <td className="px-4 py-3">
                      <span className="text-xs text-[var(--color-text-secondary)]">
                        {molecule.failed_at ?? '—'}
                      </span>
                    </td>

                    {/* Rejection Reason */}
                    <td className="px-4 py-3">
                      <span className="text-xs text-[var(--color-text-secondary)]">
                        {molecule.rejection_reason ?? '—'}
                      </span>
                    </td>
                  </tr>
                );
              })}

              {pageRows.length === 0 && (
                <tr>
                  <td
                    colSpan={4}
                    className="px-4 py-8 text-center text-sm text-[var(--color-text-muted)]"
                  >
                    {selectedStage
                      ? `No molecules rejected at stage "${selectedStage}".`
                      : 'No results to display.'}
                  </td>
                </tr>
              )}
            </tbody>
          </table>
        </div>

        {/* Pagination footer */}
        <div className="flex items-center justify-between px-4 py-3 border-t border-[var(--color-border)] bg-[var(--color-surface-sunken)]/50">
          <span className="text-xs text-[var(--color-text-muted)]">
            {PAGE_SIZE} per page &bull; {filteredMolecules.length} total
          </span>

          <div className="flex items-center gap-3">
            <ClayButton
              variant="ghost"
              size="icon"
              onClick={() => setCurrentPage((p) => p - 1)}
              disabled={currentPage <= 1}
              aria-label="Previous page"
            >
              <ChevronLeft className="w-4 h-4" />
            </ClayButton>

            <span className="text-xs text-[var(--color-text-secondary)]">
              Page {currentPage} of {totalPages}
            </span>

            <ClayButton
              variant="ghost"
              size="icon"
              onClick={() => setCurrentPage((p) => p + 1)}
              disabled={currentPage >= totalPages}
              aria-label="Next page"
            >
              <ChevronRight className="w-4 h-4" />
            </ClayButton>
          </div>
        </div>
      </div>
    </div>
  );
}
