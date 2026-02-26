/**
 * ScaffoldTreemap (VIZ-05)
 *
 * Scaffold frequency treemap showing top 50 scaffolds + "Other" bucket.
 * Click on a scaffold cell selects all molecule indices in that scaffold group.
 */

import React, { useMemo, useCallback } from 'react';
import { Treemap, ResponsiveContainer, Tooltip } from 'recharts';
import type { ScaffoldResult } from '../../../types/analytics';

interface ScaffoldTreemapProps {
  data: ScaffoldResult | null;
  selectedIndices: Set<number>;
  onSelectionChange: (indices: Set<number>) => void;
}

const SCAFFOLD_COLORS = [
  '#c41e3a', '#d97706', '#ea580c', '#059669', '#2563eb', '#7c3aed',
  '#db2777', '#0891b2', '#65a30d', '#dc2626', '#f59e0b', '#06b6d4',
  '#8b5cf6', '#ec4899', '#14b8a6', '#f97316',
];
const OTHER_COLOR = '#6b7280';

interface TreemapEntry {
  name: string;
  value: number;
  color: string;
  moleculeIndices: number[];
  fullSmiles: string;
  [key: string]: unknown;
}

/**
 * Custom content renderer for scaffold treemap cells.
 */
function CustomContent(props: Record<string, unknown>) {
  const { x, y, width, height, name, value, color } = props as {
    x: number;
    y: number;
    width: number;
    height: number;
    name: string;
    value: number;
    color: string;
  };

  if (width < 4 || height < 4) return null;

  const maxChars = Math.floor(width / 7);
  const label = name === 'Other' ? 'Other' : (name.length > maxChars ? name.slice(0, maxChars) + '...' : name);

  return (
    <g>
      <rect
        x={x}
        y={y}
        width={width}
        height={height}
        fill={color || OTHER_COLOR}
        fillOpacity={0.85}
        stroke="var(--color-surface)"
        strokeWidth={2}
        rx={4}
        style={{ cursor: 'pointer' }}
      />
      {width > 50 && height > 30 && (
        <text
          x={x + width / 2}
          y={y + height / 2 - (height > 50 ? 6 : 0)}
          textAnchor="middle"
          dominantBaseline="middle"
          fill="#fff"
          fontSize={10}
          fontWeight={600}
        >
          {label}
        </text>
      )}
      {width > 40 && height > 50 && (
        <text
          x={x + width / 2}
          y={y + height / 2 + 12}
          textAnchor="middle"
          dominantBaseline="middle"
          fill="rgba(255,255,255,0.8)"
          fontSize={10}
        >
          {value}
        </text>
      )}
    </g>
  );
}

export const ScaffoldTreemap = React.memo(function ScaffoldTreemap({
  data,
  onSelectionChange,
}: ScaffoldTreemapProps) {
  const treemapData = useMemo(() => {
    if (!data || !data.scaffolds || data.scaffolds.length === 0) return [];

    // Sort by count descending
    const sorted = [...data.scaffolds].sort((a, b) => b.count - a.count);

    // Take top 50, rest as "Other"
    const top50 = sorted.slice(0, 50);
    const rest = sorted.slice(50);

    const entries: TreemapEntry[] = top50.map((s, i) => ({
      name: s.scaffold_smiles || '(acyclic)',
      value: s.count,
      color: SCAFFOLD_COLORS[i % SCAFFOLD_COLORS.length],
      moleculeIndices: s.molecule_indices,
      fullSmiles: s.scaffold_smiles,
    }));

    if (rest.length > 0) {
      const otherIndices = rest.flatMap((s) => s.molecule_indices);
      const otherCount = rest.reduce((sum, s) => sum + s.count, 0);
      entries.push({
        name: 'Other',
        value: otherCount,
        color: OTHER_COLOR,
        moleculeIndices: otherIndices,
        fullSmiles: `${rest.length} scaffolds`,
      });
    }

    return entries;
  }, [data]);

  const handleClick = useCallback(
    (node: any) => {
      if (node?.moleculeIndices) {
        onSelectionChange(new Set(node.moleculeIndices as number[]));
      }
    },
    [onSelectionChange]
  );

  if (!data) {
    return (
      <div className="flex items-center justify-center h-[320px] text-sm text-[var(--color-text-muted)]">
        Computing scaffolds...
      </div>
    );
  }

  if (treemapData.length === 0) {
    return (
      <div className="flex items-center justify-center h-[320px] text-sm text-[var(--color-text-muted)]">
        No scaffolds found
      </div>
    );
  }

  const total = treemapData.reduce((sum, d) => sum + d.value, 0);

  return (
    <ResponsiveContainer width="100%" height={320}>
      <Treemap
        data={treemapData}
        dataKey="value"
        aspectRatio={4 / 3}
        content={<CustomContent />}
        onClick={handleClick}
      >
        <Tooltip
          content={({ active, payload }) => {
            if (!active || !payload?.[0]) return null;
            const entry = payload[0].payload as TreemapEntry;
            const pct = total > 0 ? ((entry.value / total) * 100).toFixed(1) : '0';
            return (
              <div className="rounded-lg bg-[var(--color-surface-elevated)] border border-[var(--color-border)] shadow-lg p-2 text-xs max-w-[300px]">
                <p className="font-mono font-semibold text-[var(--color-text-primary)] break-all">
                  {entry.fullSmiles || '(acyclic)'}
                </p>
                <p className="text-[var(--color-text-secondary)]">
                  Count: {entry.value} ({pct}%)
                </p>
              </div>
            );
          }}
        />
      </Treemap>
    </ResponsiveContainer>
  );
});
