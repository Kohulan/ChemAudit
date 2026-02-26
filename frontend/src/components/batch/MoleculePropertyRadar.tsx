/**
 * MoleculePropertyRadar (VIZ-08)
 *
 * Per-molecule radar chart normalized 0-1, with optional dataset average overlay.
 * Properties: MW, LogP, TPSA, QED, SA Score, Fsp3.
 */

import React, { useMemo } from 'react';
import {
  RadarChart,
  PolarGrid,
  PolarAngleAxis,
  PolarRadiusAxis,
  Radar,
  ResponsiveContainer,
  Tooltip,
  Legend,
} from 'recharts';
import type { BatchResult } from '../../types/batch';
import type { PropertyStats } from '../../types/analytics';

interface MoleculePropertyRadarProps {
  molecules: BatchResult[];
  datasetStats: PropertyStats[] | null;
}

interface RadarProperty {
  name: string;
  key: string;
  min: number;
  max: number;
}

const RADAR_PROPERTIES: RadarProperty[] = [
  { name: 'QED', key: 'qed', min: 0, max: 1 },
  { name: 'SA Score', key: 'sa_score', min: 1, max: 10 },
  { name: 'Fsp3', key: 'fsp3', min: 0, max: 1 },
  { name: 'Val. Score', key: 'overall_score', min: 0, max: 100 },
  { name: 'Lipinski Viol.', key: 'lipinski_violations', min: 0, max: 5 },
  { name: 'Alert Count', key: 'alert_count', min: 0, max: 10 },
];

function normalize(value: number, min: number, max: number): number {
  if (max <= min) return 0.5;
  return Math.max(0, Math.min(1, (value - min) / (max - min)));
}

function getMoleculeValue(r: BatchResult, key: string): number | null {
  switch (key) {
    case 'qed':
      return r.scoring?.druglikeness?.qed_score ?? null;
    case 'sa_score':
      return r.scoring?.admet?.sa_score ?? null;
    case 'fsp3':
      return r.scoring?.admet?.fsp3 ?? null;
    case 'overall_score':
      return r.validation?.overall_score ?? null;
    case 'lipinski_violations':
      return r.scoring?.druglikeness?.lipinski_violations ?? null;
    case 'alert_count':
      return r.alerts?.alert_count ?? null;
    default:
      return null;
  }
}

function getDatasetMean(stats: PropertyStats[] | null, key: string): number | null {
  if (!stats) return null;
  // Map radar keys to backend property_name values
  const mapping: Record<string, string> = {
    qed: 'qed_score',
    sa_score: 'sa_score',
    fsp3: 'fsp3',
    overall_score: 'overall_score',
    lipinski_violations: 'lipinski_violations',
    alert_count: 'alert_count',
  };
  const propName = mapping[key] || key;
  const stat = stats.find((s) => s.property_name === propName);
  return stat?.mean ?? null;
}

const MOLECULE_COLORS = [
  { stroke: 'rgba(196, 30, 58, 0.9)', fill: 'rgba(196, 30, 58, 0.2)', name: 'Molecule A' },
  { stroke: 'rgba(217, 119, 6, 0.9)', fill: 'rgba(217, 119, 6, 0.2)', name: 'Molecule B' },
];

const AVERAGE_STYLE = {
  stroke: 'rgba(16, 185, 129, 0.4)',
  fill: 'rgba(16, 185, 129, 0.08)',
};

export const MoleculePropertyRadar = React.memo(function MoleculePropertyRadar({
  molecules,
  datasetStats,
}: MoleculePropertyRadarProps) {
  const chartData = useMemo(() => {
    return RADAR_PROPERTIES.map((prop) => {
      const entry: Record<string, unknown> = {
        name: prop.name,
        fullMark: 1,
      };

      // Dataset average
      const mean = getDatasetMean(datasetStats, prop.key);
      if (mean !== null) {
        // For "bad" properties (lipinski violations, alert count), invert normalization
        const isInverted = prop.key === 'lipinski_violations' || prop.key === 'alert_count';
        entry.average = isInverted
          ? 1 - normalize(mean, prop.min, prop.max)
          : normalize(mean, prop.min, prop.max);
        entry.averageRaw = mean;
      }

      // Each molecule
      molecules.forEach((mol, i) => {
        const val = getMoleculeValue(mol, prop.key);
        if (val !== null) {
          const isInverted = prop.key === 'lipinski_violations' || prop.key === 'alert_count';
          entry[`mol${i}`] = isInverted
            ? 1 - normalize(val, prop.min, prop.max)
            : normalize(val, prop.min, prop.max);
          entry[`mol${i}Raw`] = val;
        } else {
          entry[`mol${i}`] = 0;
          entry[`mol${i}Raw`] = null;
        }
      });

      return entry;
    });
  }, [molecules, datasetStats]);

  return (
    <div className="rounded-xl bg-[var(--color-surface-sunken)] p-4">
      <ResponsiveContainer width="100%" height={280}>
        <RadarChart data={chartData} outerRadius="70%">
          <PolarGrid stroke="var(--color-border)" />
          <PolarAngleAxis
            dataKey="name"
            tick={{ fontSize: 10, fill: 'var(--color-text-secondary)' }}
          />
          <PolarRadiusAxis
            angle={90}
            domain={[0, 1]}
            tick={{ fontSize: 9, fill: 'var(--color-text-muted)' }}
            tickCount={5}
          />

          {/* Dataset Average */}
          {datasetStats && (
            <Radar
              name="Dataset Average"
              dataKey="average"
              stroke={AVERAGE_STYLE.stroke}
              fill={AVERAGE_STYLE.fill}
              fillOpacity={1}
            />
          )}

          {/* Molecule lines */}
          {molecules.map((_, i) => (
            <Radar
              key={i}
              name={molecules.length > 1 ? MOLECULE_COLORS[i].name : 'This Molecule'}
              dataKey={`mol${i}`}
              stroke={MOLECULE_COLORS[i].stroke}
              fill={MOLECULE_COLORS[i].fill}
              fillOpacity={1}
              dot
            />
          ))}

          <Tooltip
            content={({ active, payload }) => {
              if (!active || !payload?.[0]) return null;
              const entry = payload[0].payload;
              return (
                <div className="rounded-lg bg-[var(--color-surface-elevated)] border border-[var(--color-border)] shadow-lg p-2 text-xs">
                  <p className="font-semibold text-[var(--color-text-primary)] mb-1">
                    {entry.name}
                  </p>
                  {molecules.map((_, i) => {
                    const raw = entry[`mol${i}Raw`];
                    return (
                      <p key={i} className="text-[var(--color-text-secondary)]">
                        {molecules.length > 1 ? `Mol ${String.fromCharCode(65 + i)}: ` : 'Value: '}
                        {raw !== null && raw !== undefined ? Number(raw).toFixed(2) : 'N/A'}
                      </p>
                    );
                  })}
                  {entry.averageRaw !== undefined && entry.averageRaw !== null && (
                    <p className="text-emerald-600 dark:text-emerald-400">
                      Dataset avg: {Number(entry.averageRaw).toFixed(2)}
                    </p>
                  )}
                </div>
              );
            }}
          />
          <Legend wrapperStyle={{ fontSize: 11 }} />
        </RadarChart>
      </ResponsiveContainer>
    </div>
  );
});
