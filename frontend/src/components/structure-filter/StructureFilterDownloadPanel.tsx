import { Download } from 'lucide-react';
import { ClayCard } from '../ui/ClayCard';
import { ClayButton } from '../ui/ClayButton';
import type { MoleculeResult } from '../../types/structure_filter';

// =============================================================================
// Types
// =============================================================================

interface StructureFilterDownloadPanelProps {
  molecules: MoleculeResult[];
  inputCount: number;
  outputCount: number;
}

// =============================================================================
// Helpers
// =============================================================================

/** Escape a field value for CSV (wrap in double-quotes if it contains comma/quote/newline). */
function csvEscape(value: string | null): string {
  if (value === null) return '';
  // Escape internal double-quotes by doubling them, then wrap if necessary
  if (value.includes(',') || value.includes('"') || value.includes('\n')) {
    return `"${value.replace(/"/g, '""')}"`;
  }
  return value;
}

// =============================================================================
// Component
// =============================================================================

/**
 * Download panel for Structure Filter results.
 *
 * Per D-10:
 * - "Download Passed SMILES (.txt)": passed molecules only, one SMILES per line
 * - "Download Full Results (.csv)": all molecules, columns smiles/status/failed_at/rejection_reason
 * - Summary text above buttons: "{N} molecules passed all stages ({pct}% of input)"
 * - Accessibility: aria-label on each download button
 */
export function StructureFilterDownloadPanel({
  molecules,
  inputCount,
  outputCount,
}: StructureFilterDownloadPanelProps) {
  const pct = inputCount > 0 ? Math.round((outputCount / inputCount) * 100) : 0;

  // ── Download: Passed SMILES (.txt) ──────────────────────────────────────────
  const handleDownloadTxt = () => {
    const passedSmiles = molecules
      .filter((m) => m.status === 'passed')
      .map((m) => m.smiles);

    const blob = new Blob([passedSmiles.join('\n')], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'structure_filter_passed.txt';
    a.click();
    URL.revokeObjectURL(url);
  };

  // ── Download: Full Results (.csv) ───────────────────────────────────────────
  const handleDownloadCsv = () => {
    const header = 'smiles,status,failed_at,rejection_reason';
    const rows = molecules.map(
      (m) =>
        [
          csvEscape(m.smiles),
          csvEscape(m.status),
          csvEscape(m.failed_at),
          csvEscape(m.rejection_reason),
        ].join(','),
    );
    const csvContent = [header, ...rows].join('\n');

    const blob = new Blob([csvContent], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'structure_filter_results.csv';
    a.click();
    URL.revokeObjectURL(url);
  };

  return (
    <ClayCard className="p-5">
      <div className="flex flex-col sm:flex-row items-start sm:items-center justify-between gap-4">
        {/* Summary text */}
        <p className="text-sm text-[var(--color-text-secondary)]">
          <span className="font-semibold text-[var(--color-text-primary)]">{outputCount}</span>{' '}
          molecule{outputCount !== 1 ? 's' : ''} passed all stages{' '}
          <span className="text-[var(--color-text-muted)]">({pct}% of input)</span>
        </p>

        {/* Download buttons */}
        <div className="flex items-center gap-2 flex-wrap">
          <ClayButton
            variant="primary"
            size="md"
            leftIcon={<Download className="w-4 h-4" />}
            onClick={handleDownloadTxt}
            aria-label="Download passed SMILES as text file"
          >
            Download Passed SMILES (.txt)
          </ClayButton>

          <ClayButton
            variant="default"
            size="md"
            leftIcon={<Download className="w-4 h-4" />}
            onClick={handleDownloadCsv}
            aria-label="Download full filter results as CSV"
          >
            Download Full Results (.csv)
          </ClayButton>
        </div>
      </div>
    </ClayCard>
  );
}
