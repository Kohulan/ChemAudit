import {
  RadarChart,
  PolarGrid,
  PolarAngleAxis,
  PolarRadiusAxis,
  Radar,
  ResponsiveContainer,
  Tooltip,
  ScatterChart,
  XAxis,
  YAxis,
  CartesianGrid,
  Scatter,
} from 'recharts';
import type { BioavailabilityRadar, BoiledEgg } from '../../types/scoring';
import { cn } from '../../lib/utils';
import { InfoTooltip } from '../ui/Tooltip';

interface BioavailabilityCardProps {
  radar: BioavailabilityRadar | null;
  boiledEgg: BoiledEgg | null;
}

function getRegionColor(region: string) {
  switch (region) {
    case 'yolk': return 'text-amber-500';
    case 'white': return 'text-emerald-500';
    default: return 'text-gray-500';
  }
}

function getRegionBg(region: string) {
  switch (region) {
    case 'yolk': return 'bg-amber-500/10 border-amber-500/20';
    case 'white': return 'bg-emerald-500/10 border-emerald-500/20';
    default: return 'bg-gray-500/10 border-gray-500/20';
  }
}

export function BioavailabilityCard({ radar, boiledEgg }: BioavailabilityCardProps) {
  return (
    <div className="space-y-6">
      {/* Radar Chart */}
      {radar && (
        <div>
          <div className="flex items-center justify-between mb-3">
            <div className="flex items-center gap-1.5">
              <h4 className="text-sm font-semibold text-[var(--color-text-primary)]">
                Bioavailability Radar
              </h4>
              <InfoTooltip
                title="Bioavailability Radar"
                content={
                  <div className="text-xs space-y-2">
                    <p>
                      6-axis normalized assessment of oral bioavailability based on key
                      physicochemical properties. Each axis evaluates whether a property
                      falls within the optimal drug-like range.
                    </p>
                    <table className="w-full text-[10px]">
                      <thead>
                        <tr className="border-b border-white/20">
                          <th className="text-left py-0.5 text-zinc-400">Axis</th>
                          <th className="text-left py-0.5 text-zinc-400">Property</th>
                          <th className="text-left py-0.5 text-zinc-400">Range</th>
                        </tr>
                      </thead>
                      <tbody>
                        <tr><td>LIPO</td><td>WLOGP</td><td>-0.7 to 5.0</td></tr>
                        <tr><td>SIZE</td><td>MW</td><td>150–500 Da</td></tr>
                        <tr><td>POLAR</td><td>TPSA</td><td>20–130 A^2</td></tr>
                        <tr><td>INSOLU</td><td>LogS</td><td>-6 to 0</td></tr>
                        <tr><td>INSATU</td><td>Fsp3</td><td>{'\u2265'} 0.25</td></tr>
                        <tr><td>FLEX</td><td>RotBonds</td><td>{'\u2264'} 9</td></tr>
                      </tbody>
                    </table>
                    <p className="text-zinc-400 italic">
                      Daina, Michielin & Zoete. Sci. Rep. 2017, 7, 42717
                    </p>
                  </div>
                }
              />
            </div>
            <span className={cn(
              'text-xs px-2 py-0.5 rounded-full font-medium',
              radar.overall_in_range_count >= 5
                ? 'bg-emerald-500/10 text-emerald-600 dark:text-emerald-400'
                : radar.overall_in_range_count >= 3
                  ? 'bg-amber-500/10 text-amber-600 dark:text-amber-400'
                  : 'bg-red-500/10 text-red-600 dark:text-red-400'
            )}>
              {radar.overall_in_range_count}/6 in range
            </span>
          </div>

          <div className="rounded-xl bg-[var(--color-surface-sunken)] p-4">
            <ResponsiveContainer width="100%" height={280}>
              <RadarChart
                data={radar.axes.map(a => ({
                  name: a.name,
                  value: a.normalized,
                  optimal: 1.0,
                  fullMark: 1.0,
                }))}
                outerRadius="70%"
              >
                <PolarGrid stroke="var(--color-border)" />
                <PolarAngleAxis
                  dataKey="name"
                  tick={{ fontSize: 11, fill: 'var(--color-text-secondary)' }}
                />
                <PolarRadiusAxis
                  angle={90}
                  domain={[0, 1]}
                  tick={{ fontSize: 9, fill: 'var(--color-text-muted)' }}
                  tickCount={5}
                />
                {/* Optimal zone (reference) */}
                <Radar
                  name="Optimal"
                  dataKey="optimal"
                  stroke="rgba(16, 185, 129, 0.4)"
                  fill="rgba(16, 185, 129, 0.08)"
                  fillOpacity={1}
                />
                {/* Actual values */}
                <Radar
                  name="Molecule"
                  dataKey="value"
                  stroke="rgba(99, 102, 241, 0.9)"
                  fill="rgba(99, 102, 241, 0.2)"
                  fillOpacity={1}
                  dot
                />
                <Tooltip
                  content={({ active, payload }) => {
                    if (!active || !payload?.[0]) return null;
                    const axis = radar.axes.find(
                      a => a.name === payload[0]?.payload?.name
                    );
                    if (!axis) return null;
                    return (
                      <div className="rounded-lg bg-[var(--color-surface-elevated)] border border-[var(--color-border)] shadow-lg p-2 text-xs">
                        <p className="font-semibold text-[var(--color-text-primary)]">
                          {axis.name} ({axis.property_name})
                        </p>
                        <p className="text-[var(--color-text-secondary)]">
                          Value: {axis.actual_value.toFixed(2)} {axis.unit}
                        </p>
                        <p className="text-[var(--color-text-muted)]">
                          Range: {axis.optimal_min} - {axis.optimal_max}
                        </p>
                        <p className={axis.in_range ? 'text-emerald-500' : 'text-red-500'}>
                          {axis.in_range ? 'In range' : 'Out of range'}
                        </p>
                      </div>
                    );
                  }}
                />
              </RadarChart>
            </ResponsiveContainer>
          </div>

          {/* Axis detail table */}
          <div className="mt-3 grid grid-cols-2 gap-2">
            {radar.axes.map((a) => (
              <div
                key={a.name}
                className={cn(
                  'flex items-center justify-between px-3 py-2 rounded-lg text-xs',
                  a.in_range
                    ? 'bg-emerald-500/5 border border-emerald-500/10'
                    : 'bg-red-500/5 border border-red-500/10'
                )}
              >
                <span className="font-medium text-[var(--color-text-primary)]">{a.name}</span>
                <span className={cn(
                  'font-mono',
                  a.in_range ? 'text-emerald-600 dark:text-emerald-400' : 'text-red-600 dark:text-red-400'
                )}>
                  {a.actual_value.toFixed(1)}
                </span>
              </div>
            ))}
          </div>

          <p className="text-xs text-[var(--color-text-muted)] mt-2">{radar.interpretation}</p>
        </div>
      )}

      {/* BOILED-Egg Plot */}
      {boiledEgg && (
        <div>
          <div className="flex items-center justify-between mb-3">
            <div className="flex items-center gap-1.5">
              <h4 className="text-sm font-semibold text-[var(--color-text-primary)]">
                BOILED-Egg Classification
              </h4>
              <InfoTooltip
                title="BOILED-Egg Model"
                content={
                  <div className="text-xs space-y-2">
                    <p>
                      <strong>B</strong>rain <strong>O</strong>r <strong>I</strong>ntestina<strong>L</strong>{' '}
                      <strong>E</strong>stimate<strong>D</strong> permeation method. Predicts passive
                      gastrointestinal (GI) absorption and blood-brain barrier (BBB) penetration
                      using an elliptical model on TPSA vs WLOGP.
                    </p>
                    <div className="space-y-1">
                      <div className="flex items-center gap-1.5">
                        <span className="w-2.5 h-2.5 rounded-full bg-emerald-500 inline-block" />
                        <span><strong>White</strong> — GI absorbed (93% accuracy)</span>
                      </div>
                      <div className="flex items-center gap-1.5">
                        <span className="w-2.5 h-2.5 rounded-full bg-amber-500 inline-block" />
                        <span><strong>Yolk</strong> — BBB permeant (90% accuracy)</span>
                      </div>
                      <div className="flex items-center gap-1.5">
                        <span className="w-2.5 h-2.5 rounded-full bg-gray-400 inline-block" />
                        <span><strong>Grey</strong> — Neither absorbed nor permeant</span>
                      </div>
                    </div>
                    <p className="text-zinc-400">
                      Validated on 46 FDA NCEs with 83% accuracy.
                    </p>
                    <p className="text-zinc-400 italic">
                      Daina & Zoete. ChemMedChem 2016, 11, 1117–1121
                    </p>
                  </div>
                }
              />
            </div>
            <span className={cn(
              'text-xs px-2 py-0.5 rounded-full font-medium capitalize',
              getRegionColor(boiledEgg.region).replace('text-', 'bg-').replace('500', '500/10'),
              getRegionColor(boiledEgg.region)
            )}>
              {boiledEgg.region}
            </span>
          </div>

          <div className="rounded-xl bg-[var(--color-surface-sunken)] p-4">
            <ResponsiveContainer width="100%" height={280}>
              <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 20 }}>
                <CartesianGrid stroke="var(--color-border)" strokeDasharray="3 3" />
                <XAxis
                  type="number"
                  dataKey="tpsa"
                  name="TPSA"
                  domain={[0, 250]}
                  label={{ value: 'TPSA (A^2)', position: 'bottom', fontSize: 11, fill: 'var(--color-text-muted)' }}
                  tick={{ fontSize: 10, fill: 'var(--color-text-muted)' }}
                />
                <YAxis
                  type="number"
                  dataKey="wlogp"
                  name="WLOGP"
                  domain={[-3, 8]}
                  label={{ value: 'WLOGP', angle: -90, position: 'insideLeft', fontSize: 11, fill: 'var(--color-text-muted)' }}
                  tick={{ fontSize: 10, fill: 'var(--color-text-muted)' }}
                />
                {/* SVG overlay for ellipses */}
                {boiledEgg.gi_ellipse && (
                  <svg>
                    <defs>
                      <clipPath id="chartClip">
                        <rect x={0} y={0} width="100%" height="100%" />
                      </clipPath>
                    </defs>
                  </svg>
                )}
                <Scatter
                  name="Molecule"
                  data={[{ tpsa: boiledEgg.tpsa, wlogp: boiledEgg.wlogp }]}
                  fill={
                    boiledEgg.region === 'yolk'
                      ? '#f59e0b'
                      : boiledEgg.region === 'white'
                        ? '#10b981'
                        : '#6b7280'
                  }
                  shape="circle"
                  r={8}
                />
                <Tooltip
                  content={({ active }) => {
                    if (!active) return null;
                    return (
                      <div className="rounded-lg bg-[var(--color-surface-elevated)] border border-[var(--color-border)] shadow-lg p-2 text-xs">
                        <p className="font-semibold text-[var(--color-text-primary)]">Molecule</p>
                        <p className="text-[var(--color-text-secondary)]">TPSA: {boiledEgg.tpsa.toFixed(1)} A^2</p>
                        <p className="text-[var(--color-text-secondary)]">WLOGP: {boiledEgg.wlogp.toFixed(2)}</p>
                        <p className={cn('font-medium capitalize', getRegionColor(boiledEgg.region))}>
                          Region: {boiledEgg.region}
                        </p>
                      </div>
                    );
                  }}
                />
              </ScatterChart>
            </ResponsiveContainer>
          </div>

          {/* Region classification */}
          <div className={cn(
            'mt-3 p-3 rounded-xl border',
            getRegionBg(boiledEgg.region)
          )}>
            <div className="flex items-center gap-3">
              <div className="flex gap-2">
                <div className="flex items-center gap-1">
                  <span className={cn(
                    'w-3 h-3 rounded-full',
                    boiledEgg.gi_absorbed ? 'bg-emerald-500' : 'bg-gray-300 dark:bg-gray-600'
                  )} />
                  <span className="text-xs text-[var(--color-text-secondary)]">GI Absorbed</span>
                </div>
                <div className="flex items-center gap-1">
                  <span className={cn(
                    'w-3 h-3 rounded-full',
                    boiledEgg.bbb_permeant ? 'bg-amber-500' : 'bg-gray-300 dark:bg-gray-600'
                  )} />
                  <span className="text-xs text-[var(--color-text-secondary)]">BBB Permeant</span>
                </div>
              </div>
            </div>
            <p className="text-xs text-[var(--color-text-muted)] mt-2">{boiledEgg.interpretation}</p>
          </div>
        </div>
      )}
    </div>
  );
}
