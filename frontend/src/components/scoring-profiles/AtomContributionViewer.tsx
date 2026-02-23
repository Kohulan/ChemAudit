import { useMemo } from 'react';
import type { AtomContribution } from '../../types/scoring';
import { cn } from '../../lib/utils';

interface AtomContributionViewerProps {
  contributions: AtomContribution[];
  label: string;
  unit?: string;
}

/**
 * Displays per-atom contributions as a colored heatmap grid.
 * Red = high positive contribution, blue = high negative, gray = near zero.
 */
export function AtomContributionViewer({ contributions, label, unit }: AtomContributionViewerProps) {
  const { maxAbs, coloredAtoms } = useMemo(() => {
    const values = contributions.map(c => Math.abs(c.contribution));
    const max = Math.max(...values, 0.001);

    const colored = contributions
      .filter(c => Math.abs(c.contribution) > 0.001)
      .sort((a, b) => Math.abs(b.contribution) - Math.abs(a.contribution));

    return { maxAbs: max, coloredAtoms: colored };
  }, [contributions]);

  if (coloredAtoms.length === 0) {
    return (
      <p className="text-xs text-[var(--color-text-muted)]">
        No significant per-atom contributions found.
      </p>
    );
  }

  function getColor(value: number): string {
    const ratio = value / maxAbs;
    if (ratio > 0.5) return 'bg-red-500/30 text-red-700 dark:text-red-300';
    if (ratio > 0.2) return 'bg-orange-500/20 text-orange-700 dark:text-orange-300';
    if (ratio > 0) return 'bg-amber-500/15 text-amber-700 dark:text-amber-300';
    if (ratio < -0.5) return 'bg-blue-500/30 text-blue-700 dark:text-blue-300';
    if (ratio < -0.2) return 'bg-sky-500/20 text-sky-700 dark:text-sky-300';
    if (ratio < 0) return 'bg-cyan-500/15 text-cyan-700 dark:text-cyan-300';
    return 'bg-gray-500/10 text-gray-600 dark:text-gray-400';
  }

  return (
    <div>
      <p className="text-xs font-medium text-[var(--color-text-muted)] mb-2">
        {label} per-atom contributions {unit && <span className="text-[var(--color-text-muted)]">({unit})</span>}
      </p>
      <div className="flex flex-wrap gap-1">
        {coloredAtoms.map((c) => (
          <div
            key={c.atom_index}
            className={cn(
              'inline-flex flex-col items-center justify-center w-10 h-10 rounded-lg text-xs font-mono',
              getColor(c.contribution)
            )}
            title={`${c.symbol}${c.atom_index}: ${c.contribution > 0 ? '+' : ''}${c.contribution.toFixed(3)}`}
          >
            <span className="font-semibold text-xs leading-none">{c.symbol}</span>
            <span className="text-[8px] leading-none mt-0.5 opacity-70">
              {c.contribution > 0 ? '+' : ''}{c.contribution.toFixed(1)}
            </span>
          </div>
        ))}
      </div>
      <div className="flex gap-3 mt-2 text-xs text-[var(--color-text-muted)]">
        <span className="flex items-center gap-1">
          <span className="w-2.5 h-2.5 rounded bg-red-500/30" /> high +
        </span>
        <span className="flex items-center gap-1">
          <span className="w-2.5 h-2.5 rounded bg-amber-500/15" /> low +
        </span>
        <span className="flex items-center gap-1">
          <span className="w-2.5 h-2.5 rounded bg-blue-500/30" /> high -
        </span>
        <span className="flex items-center gap-1">
          <span className="w-2.5 h-2.5 rounded bg-cyan-500/15" /> low -
        </span>
      </div>
    </div>
  );
}
