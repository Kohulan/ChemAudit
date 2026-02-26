/**
 * AlertFrequencyChart (VIZ-03)
 *
 * Horizontal bar chart of alert type frequencies, sorted descending.
 */

import React, { useMemo } from 'react';
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
} from 'recharts';

interface AlertFrequencyChartProps {
  data: Record<string, number>;
  selectedIndices: Set<number>;
  onSelectionChange: (indices: Set<number>) => void;
}

export const AlertFrequencyChart = React.memo(function AlertFrequencyChart({
  data,
}: AlertFrequencyChartProps) {
  const chartData = useMemo(() => {
    return Object.entries(data)
      .map(([name, count]) => ({ name, count }))
      .sort((a, b) => b.count - a.count);
  }, [data]);

  if (chartData.length === 0) {
    return (
      <div className="flex items-center justify-center h-[200px] text-sm text-[var(--color-text-muted)]">
        No alerts detected
      </div>
    );
  }

  const chartHeight = Math.max(200, chartData.length * 32);

  return (
    <ResponsiveContainer width="100%" height={chartHeight}>
      <BarChart
        layout="vertical"
        data={chartData}
        margin={{ top: 5, right: 20, left: 0, bottom: 5 }}
      >
        <CartesianGrid strokeDasharray="3 3" stroke="var(--color-border)" horizontal={false} />
        <XAxis
          type="number"
          tick={{ fontSize: 10, fill: 'var(--color-text-secondary)' }}
          tickLine={false}
          axisLine={false}
        />
        <YAxis
          type="category"
          dataKey="name"
          width={120}
          tick={{ fontSize: 10, fill: 'var(--color-text-secondary)' }}
          tickLine={false}
          axisLine={false}
          tickFormatter={(value: string) =>
            value.length > 18 ? value.slice(0, 18) + '...' : value
          }
        />
        <Tooltip
          content={({ active, payload }) => {
            if (!active || !payload?.[0]) return null;
            const entry = payload[0].payload;
            return (
              <div className="rounded-lg bg-[var(--color-surface-elevated)] border border-[var(--color-border)] shadow-lg p-2 text-xs">
                <p className="font-semibold text-[var(--color-text-primary)]">{entry.name}</p>
                <p className="text-[var(--color-text-secondary)]">
                  Count: {entry.count}
                </p>
              </div>
            );
          }}
        />
        <Bar
          dataKey="count"
          fill="var(--color-primary)"
          radius={[0, 4, 4, 0]}
          barSize={20}
        />
      </BarChart>
    </ResponsiveContainer>
  );
});
