/**
 * TaxonomyBarChart
 *
 * Horizontal bar chart (Recharts layout="vertical") showing top 20 chemotype
 * categories sorted by count descending. Supports click-to-filter with
 * toggle behavior (clicking active bar deselects).
 */

import React, { useMemo, useCallback } from 'react';
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  Tooltip,
  ResponsiveContainer,
  Cell,
} from 'recharts';

interface TaxonomyBarChartProps {
  /** Map of category name -> molecule count */
  categoryCounts: Record<string, number>;
  /** Currently selected category (null = none) */
  activeCategory: string | null;
  /** Toggle filter for a category */
  onCategoryClick: (category: string) => void;
}

/** Truncate long category names for Y-axis display */
function truncateLabel(label: string, maxLen = 28): string {
  return label.length > maxLen ? label.slice(0, maxLen - 1) + '\u2026' : label;
}

export function TaxonomyBarChart({
  categoryCounts,
  activeCategory,
  onCategoryClick,
}: TaxonomyBarChartProps) {
  // Sort descending, take top 20 (D-10)
  const chartData = useMemo(() => {
    return Object.entries(categoryCounts)
      .sort((a, b) => b[1] - a[1])
      .slice(0, 20)
      .map(([name, count]) => ({ name, count }));
  }, [categoryCounts]);

  const handleBarClick = useCallback(
    (data: { name?: string }) => {
      if (data.name) onCategoryClick(data.name);
    },
    [onCategoryClick],
  );

  const handleKeyDown = useCallback(
    (e: React.KeyboardEvent, name: string) => {
      if (e.key === 'Enter' || e.key === ' ') {
        e.preventDefault();
        onCategoryClick(name);
      }
    },
    [onCategoryClick],
  );

  const barHeight = Math.max(44, 30); // minimum 44px touch target
  const chartHeight = Math.max(300, chartData.length * barHeight + 40);

  return (
    <div>
      <ResponsiveContainer width="100%" height={chartHeight}>
        <BarChart
          data={chartData}
          layout="vertical"
          margin={{ top: 5, right: 30, bottom: 5, left: 140 }}
        >
          <XAxis
            type="number"
            allowDecimals={false}
            tick={{ fill: 'var(--color-text-muted)', fontSize: 12 }}
          />
          <YAxis
            type="category"
            dataKey="name"
            tick={{ fill: 'var(--color-text-muted)', fontSize: 14 }}
            tickFormatter={(v: string) => truncateLabel(v)}
            width={135}
          />
          <Tooltip
            contentStyle={{
              backgroundColor: 'var(--color-surface-elevated)',
              border: '1px solid var(--color-border)',
              borderRadius: '8px',
              color: 'var(--color-text-primary)',
            }}
            formatter={(value: number | undefined) => [`${value ?? 0} molecules`, 'Count']}
          />
          <Bar
            dataKey="count"
            radius={[0, 4, 4, 0]}
            onClick={handleBarClick}
            cursor="pointer"
            minPointSize={2}
          >
            {chartData.map((entry) => {
              const isActive = activeCategory === entry.name;
              const isDimmed = activeCategory !== null && !isActive;

              let fill: string;
              let opacity: number;

              if (isActive) {
                fill = 'rgba(196, 30, 58, 1.0)';
                opacity = 1;
              } else if (isDimmed) {
                fill = 'var(--color-text-muted)';
                opacity = 0.2;
              } else {
                fill = 'rgba(196, 30, 58, 0.6)';
                opacity = 1;
              }

              return (
                <Cell
                  key={entry.name}
                  fill={fill}
                  opacity={opacity}
                  role="button"
                  aria-label={`${entry.name}: ${entry.count} molecules`}
                  tabIndex={0}
                  onKeyDown={(e: React.KeyboardEvent) => handleKeyDown(e, entry.name)}
                  style={{ transition: 'opacity 0.2s ease-out, fill 0.15s ease-out' }}
                />
              );
            })}
          </Bar>
        </BarChart>
      </ResponsiveContainer>
      <p className="text-xs text-[var(--color-text-muted)] text-center mt-2">
        Click a bar to filter the batch table
      </p>
    </div>
  );
}
