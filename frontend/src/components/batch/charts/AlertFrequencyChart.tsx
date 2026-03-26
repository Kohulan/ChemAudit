/**
 * AlertFrequencyChart (VIZ-03)
 *
 * Horizontal bar chart of alert type frequencies, sorted descending.
 * Click on a bar to filter the results table to molecules with that alert catalog.
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

interface AlertFrequencyChartProps {
  data: Record<string, number>;
  onAlertClick?: (catalogName: string) => void;
  activeAlertFilter?: string | null;
}

export const AlertFrequencyChart = React.memo(function AlertFrequencyChart({
  data,
  onAlertClick,
  activeAlertFilter,
}: AlertFrequencyChartProps) {
  const chartData = useMemo(() => {
    return Object.entries(data)
      .map(([name, count]) => ({ name, count }))
      .sort((a, b) => b.count - a.count);
  }, [data]);

  const handleClick = useCallback(
    (entry: { name: string }) => {
      if (!onAlertClick) return;
      // Toggle: clicking active filter clears it
      onAlertClick(activeAlertFilter === entry.name ? '' : entry.name);
    },
    [onAlertClick, activeAlertFilter]
  );

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
                  {entry.count} molecule{entry.count !== 1 ? 's' : ''} flagged
                </p>
                {onAlertClick && (
                  <p className="text-[var(--color-primary)] mt-1 font-medium">
                    Click to filter results
                  </p>
                )}
              </div>
            );
          }}
        />
        <Bar
          dataKey="count"
          cursor={onAlertClick ? 'pointer' : 'default'}
          onClick={(_data: unknown, index: number) => handleClick(chartData[index])}
          radius={[0, 4, 4, 0]}
          barSize={20}
        >
          {chartData.map((entry, index) => {
            const isActive = activeAlertFilter === entry.name;
            const hasFilter = !!activeAlertFilter;
            return (
              <Cell
                key={index}
                fill="var(--color-primary)"
                fillOpacity={hasFilter ? (isActive ? 1 : 0.3) : 0.85}
                stroke={isActive ? 'var(--color-text-primary)' : 'none'}
                strokeWidth={isActive ? 1.5 : 0}
              />
            );
          })}
        </Bar>
      </BarChart>
    </ResponsiveContainer>
  );
});
