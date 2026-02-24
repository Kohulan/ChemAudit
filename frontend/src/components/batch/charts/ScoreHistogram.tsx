/**
 * ScoreHistogram (VIZ-01)
 *
 * Score distribution histogram showing excellent/good/moderate/poor categories.
 * Click on a bar selects all molecule indices in that score range.
 */

import React, { useCallback, useMemo } from 'react';
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  Cell,
} from 'recharts';
import type { BatchResult } from '../../../types/batch';

interface ScoreHistogramProps {
  data: { excellent: number; good: number; moderate: number; poor: number };
  results: BatchResult[];
  selectedIndices: Set<number>;
  onSelectionChange: (indices: Set<number>) => void;
}

const CATEGORIES = [
  { key: 'excellent', label: 'Excellent (80-100)', color: '#fbbf24', min: 80, max: 101 },
  { key: 'good', label: 'Good (60-80)', color: '#d97706', min: 60, max: 80 },
  { key: 'moderate', label: 'Moderate (40-60)', color: '#ea580c', min: 40, max: 60 },
  { key: 'poor', label: 'Poor (0-40)', color: '#dc2626', min: 0, max: 40 },
] as const;

export const ScoreHistogram = React.memo(function ScoreHistogram({
  data,
  results,
  selectedIndices,
  onSelectionChange,
}: ScoreHistogramProps) {
  const total = data.excellent + data.good + data.moderate + data.poor;

  const chartData = useMemo(
    () =>
      CATEGORIES.map((cat) => ({
        name: cat.label,
        count: data[cat.key as keyof typeof data],
        color: cat.color,
        min: cat.min,
        max: cat.max,
      })),
    [data]
  );

  // Build index maps per score range for click-to-select
  const rangeIndices = useMemo(() => {
    const map: Record<string, number[]> = {};
    for (const cat of CATEGORIES) {
      map[cat.label] = [];
    }
    for (const r of results) {
      if (r.status !== 'success' || !r.validation) continue;
      const score = r.validation.overall_score;
      for (const cat of CATEGORIES) {
        if (score >= cat.min && score < cat.max) {
          map[cat.label].push(r.index);
          break;
        }
      }
    }
    return map;
  }, [results]);

  const handleClick = useCallback(
    (entry: { name: string }) => {
      const indices = rangeIndices[entry.name];
      if (indices) {
        onSelectionChange(new Set(indices));
      }
    },
    [rangeIndices, onSelectionChange]
  );

  if (total === 0) {
    return (
      <div className="flex items-center justify-center h-[240px] text-sm text-[var(--color-text-muted)]">
        No score data available
      </div>
    );
  }

  const hasSelection = selectedIndices.size > 0;

  return (
    <ResponsiveContainer width="100%" height={240}>
      <BarChart
        data={chartData}
        barCategoryGap={0}
        barGap={0}
        margin={{ top: 5, right: 10, left: 0, bottom: 5 }}
      >
        <CartesianGrid strokeDasharray="3 3" stroke="var(--color-border)" />
        <XAxis
          dataKey="name"
          tick={{ fontSize: 10, fill: 'var(--color-text-secondary)' }}
          tickLine={false}
        />
        <YAxis
          tick={{ fontSize: 10, fill: 'var(--color-text-secondary)' }}
          tickLine={false}
          axisLine={false}
        />
        <Tooltip
          content={({ active, payload }) => {
            if (!active || !payload?.[0]) return null;
            const entry = payload[0].payload;
            const pct = total > 0 ? ((entry.count / total) * 100).toFixed(1) : '0';
            return (
              <div className="rounded-lg bg-[var(--color-surface-elevated)] border border-[var(--color-border)] shadow-lg p-2 text-xs">
                <p className="font-semibold text-[var(--color-text-primary)]">{entry.name}</p>
                <p className="text-[var(--color-text-secondary)]">
                  {entry.count} molecule{entry.count !== 1 ? 's' : ''} ({pct}%)
                </p>
              </div>
            );
          }}
        />
        <Bar
          dataKey="count"
          cursor="pointer"
          onClick={(_data: unknown, index: number) => handleClick(chartData[index])}
          radius={[4, 4, 0, 0]}
        >
          {chartData.map((entry, index) => {
            // Determine if any selected molecules are in this bar's range
            const barIndices = rangeIndices[entry.name] || [];
            const hasSelectedInBar =
              hasSelection && barIndices.some((i) => selectedIndices.has(i));
            const opacity = hasSelection ? (hasSelectedInBar ? 1 : 0.4) : 1;

            return <Cell key={index} fill={entry.color} fillOpacity={opacity} />;
          })}
        </Bar>
      </BarChart>
    </ResponsiveContainer>
  );
});
