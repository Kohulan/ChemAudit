/**
 * PropertyScatterPlot (VIZ-02)
 *
 * Property vs property scatter plot with configurable X, Y, and color-by axes.
 * Click on a point toggles it in selectedIndices.
 */

import React, { useMemo, useCallback } from 'react';
import {
  ScatterChart,
  Scatter,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
} from 'recharts';
import type { BatchResult } from '../../../types/batch';

interface PropertyScatterPlotProps {
  results: BatchResult[];
  xProperty: string;
  yProperty: string;
  colorByProperty: string;
  selectedIndices: Set<number>;
  onSelectionChange: (indices: Set<number>) => void;
  onPropertyChange: (axis: 'x' | 'y' | 'color', property: string) => void;
}

const PROPERTY_OPTIONS = [
  'MW', 'LogP', 'TPSA', 'QED', 'overall_score', 'SA_score', 'Fsp3',
];

/**
 * Extract a property value from a BatchResult.
 */
function getProperty(r: BatchResult, prop: string): number | null {
  if (r.status !== 'success') return null;

  switch (prop) {
    case 'MW':
      // MW not directly stored; approximate from validation or use null
      return null;
    case 'LogP':
      return null; // Not directly available in BatchResult
    case 'TPSA':
      return null;
    case 'QED':
      return r.scoring?.druglikeness?.qed_score ?? null;
    case 'overall_score':
      return r.validation?.overall_score ?? null;
    case 'SA_score':
      return r.scoring?.admet?.sa_score ?? null;
    case 'Fsp3':
      return r.scoring?.admet?.fsp3 ?? null;
    default:
      return null;
  }
}

function valueToColor(value: number, min: number, max: number): string {
  const t = max > min ? (value - min) / (max - min) : 0.5;
  // Green (good) -> Yellow -> Red (poor)
  const r = Math.round(255 * Math.min(1, t * 2));
  const g = Math.round(255 * Math.min(1, (1 - t) * 2));
  return `rgb(${r}, ${g}, 60)`;
}

export const PropertyScatterPlot = React.memo(function PropertyScatterPlot({
  results,
  xProperty,
  yProperty,
  colorByProperty,
  selectedIndices,
  onSelectionChange,
  onPropertyChange,
}: PropertyScatterPlotProps) {
  // Build data points
  const { points, colorMin, colorMax } = useMemo(() => {
    const pts: Array<{
      x: number;
      y: number;
      colorValue: number;
      index: number;
      smiles: string;
    }> = [];
    let cMin = Infinity;
    let cMax = -Infinity;

    for (const r of results) {
      const xVal = getProperty(r, xProperty);
      const yVal = getProperty(r, yProperty);
      const cVal = getProperty(r, colorByProperty);
      if (xVal === null || yVal === null) continue;
      const cv = cVal ?? 0;
      if (cv < cMin) cMin = cv;
      if (cv > cMax) cMax = cv;
      pts.push({ x: xVal, y: yVal, colorValue: cv, index: r.index, smiles: r.smiles });
    }
    return { points: pts, colorMin: cMin === Infinity ? 0 : cMin, colorMax: cMax === -Infinity ? 1 : cMax };
  }, [results, xProperty, yProperty, colorByProperty]);

  const hasSelection = selectedIndices.size > 0;

  const handlePointClick = useCallback(
    (entry: { index: number }) => {
      const next = new Set(selectedIndices);
      if (next.has(entry.index)) {
        next.delete(entry.index);
      } else {
        next.add(entry.index);
      }
      onSelectionChange(next);
    },
    [selectedIndices, onSelectionChange]
  );

  if (points.length === 0) {
    return (
      <div className="space-y-3">
        <PropertyDropdowns
          xProperty={xProperty}
          yProperty={yProperty}
          colorByProperty={colorByProperty}
          onPropertyChange={onPropertyChange}
        />
        <div className="flex items-center justify-center h-[280px] text-sm text-[var(--color-text-muted)]">
          No data available for selected properties
        </div>
      </div>
    );
  }

  return (
    <div className="space-y-3">
      <PropertyDropdowns
        xProperty={xProperty}
        yProperty={yProperty}
        colorByProperty={colorByProperty}
        onPropertyChange={onPropertyChange}
      />
      <ResponsiveContainer width="100%" height={320}>
        <ScatterChart margin={{ top: 10, right: 10, left: 0, bottom: 10 }}>
          <CartesianGrid strokeDasharray="3 3" stroke="var(--color-border)" />
          <XAxis
            type="number"
            dataKey="x"
            name={xProperty}
            tick={{ fontSize: 10, fill: 'var(--color-text-secondary)' }}
            label={{
              value: xProperty,
              position: 'bottom',
              fontSize: 11,
              fill: 'var(--color-text-muted)',
              offset: -5,
            }}
          />
          <YAxis
            type="number"
            dataKey="y"
            name={yProperty}
            tick={{ fontSize: 10, fill: 'var(--color-text-secondary)' }}
            label={{
              value: yProperty,
              angle: -90,
              position: 'insideLeft',
              fontSize: 11,
              fill: 'var(--color-text-muted)',
            }}
          />
          <Tooltip
            content={({ active, payload }) => {
              if (!active || !payload?.[0]) return null;
              const pt = payload[0].payload;
              const smiles =
                pt.smiles.length > 40 ? pt.smiles.slice(0, 40) + '...' : pt.smiles;
              return (
                <div className="rounded-lg bg-[var(--color-surface-elevated)] border border-[var(--color-border)] shadow-lg p-2 text-xs max-w-[300px]">
                  <p className="font-mono text-[var(--color-text-primary)] truncate">{smiles}</p>
                  <p className="text-[var(--color-text-secondary)]">
                    {xProperty}: {pt.x?.toFixed(2)} | {yProperty}: {pt.y?.toFixed(2)}
                  </p>
                  <p className="text-[var(--color-text-muted)]">
                    {colorByProperty}: {pt.colorValue?.toFixed(2)}
                  </p>
                </div>
              );
            }}
          />
          <Scatter
            data={points}
            cursor="pointer"
            onClick={(_data: unknown, _index: number, event: unknown) => {
              // Recharts scatter passes the data entry differently
              const e = event as { payload?: { index: number } };
              if (e?.payload) handlePointClick(e.payload);
            }}
            shape={(props: unknown) => {
              const p = props as {
                cx: number;
                cy: number;
                payload: { colorValue: number; index: number };
              };
              const isSelected = selectedIndices.has(p.payload.index);
              const opacity = hasSelection ? (isSelected ? 1 : 0.4) : 0.8;
              const r = isSelected ? 6 : 4;
              const fill = valueToColor(p.payload.colorValue, colorMin, colorMax);

              return (
                <circle
                  cx={p.cx}
                  cy={p.cy}
                  r={r}
                  fill={fill}
                  fillOpacity={opacity}
                  stroke={isSelected ? '#c41e3a' : 'none'}
                  strokeWidth={isSelected ? 2 : 0}
                />
              );
            }}
          />
        </ScatterChart>
      </ResponsiveContainer>
    </div>
  );
});

interface PropertyDropdownsProps {
  xProperty: string;
  yProperty: string;
  colorByProperty: string;
  onPropertyChange: (axis: 'x' | 'y' | 'color', property: string) => void;
}

function PropertyDropdowns({
  xProperty,
  yProperty,
  colorByProperty,
  onPropertyChange,
}: PropertyDropdownsProps) {
  return (
    <div className="flex flex-wrap gap-3">
      <label className="flex items-center gap-1.5 text-xs text-[var(--color-text-secondary)]">
        X:
        <select
          value={xProperty}
          onChange={(e) => onPropertyChange('x', e.target.value)}
          className="text-xs px-2 py-1 rounded-md bg-[var(--color-surface-sunken)] border border-[var(--color-border)] text-[var(--color-text-primary)]"
        >
          {PROPERTY_OPTIONS.map((p) => (
            <option key={p} value={p}>
              {p}
            </option>
          ))}
        </select>
      </label>
      <label className="flex items-center gap-1.5 text-xs text-[var(--color-text-secondary)]">
        Y:
        <select
          value={yProperty}
          onChange={(e) => onPropertyChange('y', e.target.value)}
          className="text-xs px-2 py-1 rounded-md bg-[var(--color-surface-sunken)] border border-[var(--color-border)] text-[var(--color-text-primary)]"
        >
          {PROPERTY_OPTIONS.map((p) => (
            <option key={p} value={p}>
              {p}
            </option>
          ))}
        </select>
      </label>
      <label className="flex items-center gap-1.5 text-xs text-[var(--color-text-secondary)]">
        Color:
        <select
          value={colorByProperty}
          onChange={(e) => onPropertyChange('color', e.target.value)}
          className="text-xs px-2 py-1 rounded-md bg-[var(--color-surface-sunken)] border border-[var(--color-border)] text-[var(--color-text-primary)]"
        >
          {PROPERTY_OPTIONS.map((p) => (
            <option key={p} value={p}>
              {p}
            </option>
          ))}
        </select>
      </label>
    </div>
  );
}
