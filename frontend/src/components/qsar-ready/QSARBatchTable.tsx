import { ChevronLeft, ChevronRight } from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';
import { QSARMoleculeRow } from './QSARMoleculeRow';
import type { QSARReadyResult } from '../../types/qsar_ready';

interface QSARBatchTableProps {
  results: QSARReadyResult[];
  page: number;
  totalPages: number;
  onPageChange: (page: number) => void;
}

/**
 * Paginated table of QSAR-Ready batch results.
 *
 * Per UI-SPEC batch result layout:
 * - Full-width table with no max-height (pagination handles size)
 * - Table header: index, Original SMILES, Curated SMILES, InChIKey, Status
 * - Rows: QSARMoleculeRow per result (expandable with mini step chips)
 * - Min row height: 48px
 * - Pagination footer: "Page {page} of {totalPages}", prev/next ClayButtons
 * - "50 per page" label
 */
export function QSARBatchTable({
  results,
  page,
  totalPages,
  onPageChange,
}: QSARBatchTableProps) {
  return (
    <div className="overflow-hidden rounded-lg border border-[var(--color-border)]">
      {/* ── Table ── */}
      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="bg-[var(--color-surface-sunken)] border-b border-[var(--color-border)]">
              <th className="px-4 py-3 text-left text-xs font-semibold text-[var(--color-text-secondary)] w-12">
                #
              </th>
              <th className="px-4 py-3 text-left text-xs font-semibold text-[var(--color-text-secondary)]">
                Original SMILES
              </th>
              <th className="px-4 py-3 text-left text-xs font-semibold text-[var(--color-text-secondary)]">
                Curated SMILES
              </th>
              <th className="px-4 py-3 text-left text-xs font-semibold text-[var(--color-text-secondary)]">
                InChIKey
              </th>
              <th className="px-4 py-3 text-left text-xs font-semibold text-[var(--color-text-secondary)]">
                Status
              </th>
              {/* Expand chevron column — no header */}
              <th className="w-8" />
            </tr>
          </thead>
          <tbody>
            {results.map((result, idx) => (
              <QSARMoleculeRow
                key={`${page}-${idx}`}
                result={result}
                index={idx}
              />
            ))}
            {results.length === 0 && (
              <tr>
                <td
                  colSpan={6}
                  className="px-4 py-8 text-center text-sm text-[var(--color-text-muted)]"
                >
                  No results on this page.
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>

      {/* ── Pagination footer ── */}
      <div className="flex items-center justify-between px-4 py-3 border-t border-[var(--color-border)] bg-[var(--color-surface-sunken)]/50">
        <span className="text-xs text-[var(--color-text-muted)]">50 per page</span>

        <div className="flex items-center gap-3">
          <ClayButton
            variant="ghost"
            size="icon"
            onClick={() => onPageChange(page - 1)}
            disabled={page <= 1}
            aria-label="Previous page"
          >
            <ChevronLeft className="w-4 h-4" />
          </ClayButton>

          <span className="text-xs text-[var(--color-text-secondary)]">
            Page {page} of {totalPages}
          </span>

          <ClayButton
            variant="ghost"
            size="icon"
            onClick={() => onPageChange(page + 1)}
            disabled={page >= totalPages}
            aria-label="Next page"
          >
            <ChevronRight className="w-4 h-4" />
          </ClayButton>
        </div>
      </div>
    </div>
  );
}
