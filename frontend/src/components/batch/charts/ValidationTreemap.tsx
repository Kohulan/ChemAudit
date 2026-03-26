/**
 * ValidationTreemap (VIZ-04)
 *
 * Treemap of validation issue counts with categorical colors.
 * Cells are clickable to filter the results table by issue type.
 * Tooltip shows check description on hover.
 */

import React, { useMemo, useCallback } from 'react';
import { Treemap, ResponsiveContainer, Tooltip } from 'recharts';
import { CHECK_DESCRIPTIONS, formatCheckName } from '../../../constants/checkDescriptions';

interface ValidationTreemapProps {
  data: Record<string, number>;
  onIssueClick?: (checkName: string) => void;
  activeIssueFilter?: string | null;
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
  const { x, y, width, height, name, value, color, activeIssueFilter, onIssueClick } = props as {
    x: number;
    y: number;
    width: number;
    height: number;
    name: string;
    value: number;
    color: string;
    activeIssueFilter?: string | null;
    onIssueClick?: (name: string) => void;
  };

  if (width < 4 || height < 4) return null;

  const isActive = activeIssueFilter === name;
  const isDimmed = activeIssueFilter && !isActive;

  return (
    <g
      onClick={(e) => {
        e.stopPropagation();
        onIssueClick?.(name);
      }}
      style={{ cursor: onIssueClick ? 'pointer' : 'default' }}
    >
      <rect
        x={x}
        y={y}
        width={width}
        height={height}
        fill={color || '#6b7280'}
        fillOpacity={isDimmed ? 0.4 : 0.85}
        stroke={isActive ? '#fff' : 'var(--color-surface)'}
        strokeWidth={isActive ? 3 : 2}
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
          style={{ pointerEvents: 'none' }}
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
          style={{ pointerEvents: 'none' }}
        >
          {value}
        </text>
      )}
    </g>
  );
}

export const ValidationTreemap = React.memo(function ValidationTreemap({
  data,
  onIssueClick,
  activeIssueFilter,
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

  const handleCellClick = useCallback(
    (name: string) => {
      if (!onIssueClick) return;
      // Toggle: clicking active filter clears it
      onIssueClick(activeIssueFilter === name ? '' : name);
    },
    [onIssueClick, activeIssueFilter]
  );

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
        content={
          <CustomContent
            activeIssueFilter={activeIssueFilter}
            onIssueClick={handleCellClick}
          />
        }
      >
        <Tooltip
          content={({ active, payload }) => {
            if (!active || !payload?.[0]) return null;
            const entry = payload[0].payload;
            const description = CHECK_DESCRIPTIONS[entry.name];
            return (
              <div className="rounded-lg bg-[var(--color-surface-elevated)] border border-[var(--color-border)] shadow-lg p-3 text-xs max-w-[280px]">
                <p className="font-semibold text-[var(--color-text-primary)] mb-1">
                  {formatCheckName(entry.name)}
                </p>
                {description && (
                  <p className="text-[var(--color-text-secondary)] mb-2 leading-relaxed">
                    {description}
                  </p>
                )}
                <p className="text-[var(--color-text-muted)]">
                  {entry.value} molecule{entry.value !== 1 ? 's' : ''} affected
                </p>
                {onIssueClick && (
                  <p className="text-[var(--color-primary)] mt-1.5 font-medium">
                    Click to filter results
                  </p>
                )}
              </div>
            );
          }}
        />
      </Treemap>
    </ResponsiveContainer>
  );
});
