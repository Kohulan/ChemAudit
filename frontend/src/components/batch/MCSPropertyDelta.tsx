/**
 * MCSPropertyDelta
 *
 * Table comparing molecular properties between two molecules with
 * delta column showing favorable/unfavorable differences.
 */

import type { MCSPropertyDelta as MCSPropertyDeltaType } from '../../types/analytics';

interface MCSPropertyDeltaProps {
  deltas: MCSPropertyDeltaType[];
}

/** Properties where higher values are generally favorable. */
const FAVORABLE_HIGHER = new Set(['QED']);
/** Properties where lower values are generally favorable. */
const FAVORABLE_LOWER = new Set<string>();

/**
 * Format property values with appropriate decimal places.
 */
function formatValue(property: string, value: number): string {
  switch (property) {
    case 'MW':
      return value.toFixed(1);
    case 'LogP':
      return value.toFixed(2);
    case 'TPSA':
      return value.toFixed(1);
    case 'QED':
      return value.toFixed(3);
    case 'HBA':
    case 'HBD':
    case 'RotBonds':
    case 'RingCount':
      return Math.round(value).toString();
    default:
      return value.toFixed(2);
  }
}

/**
 * Get color class for delta value based on whether higher/lower is favorable.
 */
function getDeltaColor(property: string, delta: number): string {
  if (delta === 0) return 'text-[var(--color-text-muted)]';
  if (FAVORABLE_HIGHER.has(property)) {
    return delta > 0
      ? 'text-emerald-600 dark:text-emerald-400'
      : 'text-red-600 dark:text-red-400';
  }
  if (FAVORABLE_LOWER.has(property)) {
    return delta < 0
      ? 'text-emerald-600 dark:text-emerald-400'
      : 'text-red-600 dark:text-red-400';
  }
  // Neutral properties (MW, LogP, TPSA, HBA, HBD, RotBonds, RingCount)
  return 'text-[var(--color-text-muted)]';
}

function formatDelta(delta: number, property: string): string {
  const formatted = formatValue(property, Math.abs(delta));
  if (delta > 0) return `+${formatted}`;
  if (delta < 0) return `-${formatted}`;
  return '0';
}

export function MCSPropertyDelta({ deltas }: MCSPropertyDeltaProps) {
  return (
    <div
      className="rounded-xl border border-[var(--color-border)] overflow-hidden"
      role="table"
      aria-label="Property comparison between Molecule A and Molecule B"
    >
      <table className="w-full text-xs">
        <thead>
          <tr className="bg-[var(--color-surface-sunken)]">
            <th className="px-3 py-2 text-left text-[var(--color-text-muted)] font-medium">
              Property
            </th>
            <th className="px-3 py-2 text-right text-[var(--color-text-muted)] font-medium">
              Mol A
            </th>
            <th className="px-3 py-2 text-right text-[var(--color-text-muted)] font-medium">
              Mol B
            </th>
            <th className="px-3 py-2 text-right text-[var(--color-text-muted)] font-medium">
              Delta
            </th>
          </tr>
        </thead>
        <tbody className="divide-y divide-[var(--color-border)]">
          {deltas.map((row) => (
            <tr key={row.property}>
              <td className="px-3 py-2 text-[var(--color-text-secondary)] font-medium">
                {row.property}
              </td>
              <td className="px-3 py-2 text-right font-mono tabular-nums text-[var(--color-text-primary)]">
                {formatValue(row.property, row.mol_a)}
              </td>
              <td className="px-3 py-2 text-right font-mono tabular-nums text-[var(--color-text-primary)]">
                {formatValue(row.property, row.mol_b)}
              </td>
              <td
                className={`px-3 py-2 text-right font-mono tabular-nums ${getDeltaColor(row.property, row.delta)}`}
              >
                {formatDelta(row.delta, row.property)}
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
