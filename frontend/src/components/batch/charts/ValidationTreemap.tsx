/**
 * ValidationTreemap (VIZ-04)
 *
 * Treemap of validation issue counts with categorical colors.
 * Groups issues by category prefix when possible.
 */

import React, { useMemo } from 'react';
import { Treemap, ResponsiveContainer, Tooltip } from 'recharts';

interface ValidationTreemapProps {
  data: Record<string, number>;
  selectedIndices: Set<number>;
  onSelectionChange: (indices: Set<number>) => void;
}

const CATEGORY_COLORS = [
  '#c41e3a', // primary
  '#d97706', // amber
  '#ea580c', // orange
  '#059669', // emerald
  '#2563eb', // blue
  '#7c3aed', // purple
  '#db2777', // pink
  '#6b7280', // gray
];

interface TreemapNode {
  name: string;
  value: number;
  color: string;
  children?: TreemapNode[];
}

/**
 * Custom content renderer for treemap cells.
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

  return (
    <g>
      <rect
        x={x}
        y={y}
        width={width}
        height={height}
        fill={color || '#6b7280'}
        fillOpacity={0.85}
        stroke="var(--color-surface)"
        strokeWidth={2}
        rx={4}
      />
      {width > 60 && height > 30 && (
        <text
          x={x + width / 2}
          y={y + height / 2 - (height > 50 ? 6 : 0)}
          textAnchor="middle"
          dominantBaseline="middle"
          fill="#fff"
          fontSize={11}
          fontWeight={600}
        >
          {name && name.length > Math.floor(width / 7)
            ? name.slice(0, Math.floor(width / 7)) + '...'
            : name}
        </text>
      )}
      {width > 60 && height > 50 && (
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

export const ValidationTreemap = React.memo(function ValidationTreemap({
  data,
}: ValidationTreemapProps) {
  const treemapData = useMemo(() => {
    const entries = Object.entries(data);
    if (entries.length === 0) return [];

    // Try to group by prefix (e.g., "stereo_", "composition_")
    const groups: Record<string, Array<{ name: string; value: number }>> = {};

    for (const [name, value] of entries) {
      const underscoreIdx = name.indexOf('_');
      const prefix =
        underscoreIdx > 0 && underscoreIdx < name.length - 1
          ? name.slice(0, underscoreIdx)
          : 'other';
      if (!groups[prefix]) groups[prefix] = [];
      groups[prefix].push({ name, value });
    }

    // If most items fall into "other", use flat treemap
    const groupKeys = Object.keys(groups);
    const otherPct =
      (groups['other']?.length ?? 0) / entries.length;

    if (groupKeys.length <= 2 || otherPct > 0.6) {
      // Flat treemap
      return entries.map(([name, value], i) => ({
        name,
        value,
        color: CATEGORY_COLORS[i % CATEGORY_COLORS.length],
      }));
    }

    // Hierarchical treemap
    const result: TreemapNode[] = [];
    let colorIdx = 0;
    for (const [, items] of Object.entries(groups)) {
      const color = CATEGORY_COLORS[colorIdx % CATEGORY_COLORS.length];
      colorIdx++;
      for (const item of items) {
        result.push({
          name: item.name,
          value: item.value,
          color,
        });
      }
    }
    return result;
  }, [data]);

  if (treemapData.length === 0) {
    return (
      <div className="flex items-center justify-center h-[280px] text-sm text-[var(--color-text-muted)]">
        No validation issues
      </div>
    );
  }

  return (
    <ResponsiveContainer width="100%" height={280}>
      <Treemap
        data={treemapData}
        dataKey="value"
        aspectRatio={4 / 3}
        content={<CustomContent />}
      >
        <Tooltip
          content={({ active, payload }) => {
            if (!active || !payload?.[0]) return null;
            const entry = payload[0].payload;
            return (
              <div className="rounded-lg bg-[var(--color-surface-elevated)] border border-[var(--color-border)] shadow-lg p-2 text-xs">
                <p className="font-semibold text-[var(--color-text-primary)]">{entry.name}</p>
                <p className="text-[var(--color-text-secondary)]">
                  Count: {entry.value}
                </p>
              </div>
            );
          }}
        />
      </Treemap>
    </ResponsiveContainer>
  );
});
