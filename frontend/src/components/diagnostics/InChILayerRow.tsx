import { Check, X } from 'lucide-react';

interface InChILayerRowProps {
  /** Layer name as returned by the API (e.g. "connections", "double_bond_stereo"). */
  layer: string;
  /** Value from InChI-A for this layer, or null if absent. */
  valueA: string | null;
  /** Value from InChI-B for this layer, or null if absent. */
  valueB: string | null;
  /** True when both values are present and identical. */
  match: boolean;
}

/**
 * Single row in the InChI layer diff table (D-07).
 *
 * Columns: Layer name | InChI-A value | InChI-B value | Match indicator.
 * Rows are non-interactive (display only).
 */
export function InChILayerRow({ layer, valueA, valueB, match }: InChILayerRowProps) {
  // Humanise the layer name: replace underscores, capitalise first letter.
  const humanLayer = layer
    .replace(/_/g, ' ')
    .replace(/^\w/, c => c.toUpperCase());

  const bothPresent = valueA !== null && valueB !== null;

  return (
    <tr className="border-b border-[var(--color-border)] last:border-0">
      {/* Layer name */}
      <td className="px-4 py-2 text-sm font-semibold font-display text-[var(--color-text-primary)] whitespace-nowrap">
        {humanLayer}
      </td>

      {/* InChI-A value */}
      <td className="px-4 py-2 font-mono text-xs text-[var(--color-text-secondary)] break-all max-w-[200px]">
        {valueA !== null ? (
          valueA
        ) : (
          <span className="text-[var(--color-text-muted)]">--</span>
        )}
      </td>

      {/* InChI-B value */}
      <td className="px-4 py-2 font-mono text-xs text-[var(--color-text-secondary)] break-all max-w-[200px]">
        {valueB !== null ? (
          valueB
        ) : (
          <span className="text-[var(--color-text-muted)]">--</span>
        )}
      </td>

      {/* Match indicator */}
      <td className="px-4 py-2 text-center">
        {!bothPresent ? (
          <span
            className="text-[var(--color-text-muted)] text-sm"
            aria-label="Layer not present"
          >
            --
          </span>
        ) : match ? (
          <Check
            className="w-4 h-4 text-status-success inline-block"
            aria-label="Match"
          />
        ) : (
          <X
            className="w-4 h-4 text-status-error inline-block"
            aria-label="Mismatch"
          />
        )}
      </td>
    </tr>
  );
}
