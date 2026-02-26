/**
 * FragmentClassificationTable Component
 *
 * Mini-table displaying mixture fragment classifications for mixture detection results.
 * Shows SMILES, molecular weight, classification badge, and pattern name for each fragment.
 * Designed per user decision: "a proper mini-table, not just a summary line".
 */
import { Badge } from '../ui/Badge';
import type { FragmentDetail } from '../../types/validation';

interface FragmentClassificationTableProps {
  fragments: FragmentDetail[];
}

/** Returns the Badge variant for a fragment classification */
function getClassificationVariant(
  classification: FragmentDetail['classification']
): 'success' | 'info' | 'warning' | 'default' {
  switch (classification) {
    case 'drug':
      return 'success';
    case 'salt':
      return 'info';
    case 'solvent':
      return 'warning';
    default:
      return 'default';
  }
}

/**
 * Renders a mini-table of fragment classifications for mixture detection results.
 */
export function FragmentClassificationTable({ fragments }: FragmentClassificationTableProps) {
  if (fragments.length === 0) {
    return (
      <p className="text-xs text-[var(--color-text-muted)] mt-2">
        No fragment details available.
      </p>
    );
  }

  return (
    <div className="mt-3 overflow-x-auto rounded-lg border border-[var(--color-border)]">
      <table className="w-full text-xs min-w-[400px]">
        <thead>
          <tr className="border-b border-[var(--color-border)] bg-[var(--color-surface-sunken)]">
            <th className="text-left px-3 py-2 font-medium text-[var(--color-text-muted)] uppercase tracking-wider text-[10px]">
              Fragment
            </th>
            <th className="text-right px-3 py-2 font-medium text-[var(--color-text-muted)] uppercase tracking-wider text-[10px]">
              MW
            </th>
            <th className="text-center px-3 py-2 font-medium text-[var(--color-text-muted)] uppercase tracking-wider text-[10px]">
              Classification
            </th>
            <th className="text-left px-3 py-2 font-medium text-[var(--color-text-muted)] uppercase tracking-wider text-[10px]">
              Pattern
            </th>
          </tr>
        </thead>
        <tbody>
          {fragments.map((fragment, idx) => (
            <tr
              key={idx}
              className="border-b border-[var(--color-border)] last:border-0 hover:bg-[var(--color-surface-sunken)] transition-colors"
            >
              <td className="px-3 py-2">
                <code className="font-mono text-[var(--color-text-secondary)] break-all">
                  {fragment.smiles}
                </code>
              </td>
              <td className="px-3 py-2 text-right text-[var(--color-text-secondary)] tabular-nums whitespace-nowrap">
                {fragment.molecular_weight.toFixed(1)}
              </td>
              <td className="px-3 py-2 text-center">
                <Badge
                  variant={getClassificationVariant(fragment.classification)}
                  size="sm"
                >
                  {fragment.classification}
                </Badge>
              </td>
              <td className="px-3 py-2 text-[var(--color-text-muted)]">
                {fragment.pattern_name ?? '—'}
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
