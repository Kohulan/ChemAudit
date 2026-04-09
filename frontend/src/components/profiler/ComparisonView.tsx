import { motion } from 'framer-motion';
import { X } from 'lucide-react';
import {
  RadarChart,
  PolarGrid,
  PolarAngleAxis,
  Radar,
  Legend,
  ResponsiveContainer,
} from 'recharts';
import { MoleculeViewer } from '../molecules/MoleculeViewer';
import { ClayButton } from '../ui/ClayButton';
import type { PinnedMolecule } from './ComparisonBar';

interface ComparisonViewProps {
  molecules: PinnedMolecule[];
  onClose: () => void;
}

/**
 * Drug-like property ranges for radar normalization (0-1).
 */
const RADAR_RANGES: Record<string, { min: number; max: number }> = {
  MW:       { min: 0,  max: 500 },
  LogP:     { min: -2, max: 5 },
  TPSA:     { min: 0,  max: 140 },
  HBD:      { min: 0,  max: 5 },
  HBA:      { min: 0,  max: 10 },
  RotBonds: { min: 0,  max: 10 },
};

const MOLECULE_COLORS = [
  'var(--color-primary)',
  '#3b82f6',
  '#10b981',
  '#f59e0b',
  '#8b5cf6',
];

function normalize(value: number, min: number, max: number): number {
  return Math.min(1, Math.max(0, (value - min) / (max - min)));
}

function buildRadarData(mol: PinnedMolecule) {
  const profile = mol.profile;
  const detailMap: Record<string, number> = {};
  for (const d of profile.stars.details) {
    detailMap[d.property.toLowerCase()] = d.value;
  }

  const mw = detailMap['mw'] ?? detailMap['mol weight'] ?? detailMap['molecular weight'] ?? 0;
  const logp = profile.consensus_logp?.consensus_logp ?? profile.pfi?.clogp ?? 0;
  const tpsa = detailMap['tpsa'] ?? 0;
  const hbd = detailMap['hbd'] ?? detailMap['h-bond donors'] ?? detailMap['hbond donors'] ?? 0;
  const hba = detailMap['hba'] ?? detailMap['h-bond acceptors'] ?? detailMap['hbond acceptors'] ?? 0;
  const rotbonds = detailMap['rotatable bonds'] ?? detailMap['rotbonds'] ?? detailMap['rot bonds'] ?? 0;

  return {
    MW: normalize(mw, RADAR_RANGES.MW.min, RADAR_RANGES.MW.max),
    LogP: normalize(logp, RADAR_RANGES.LogP.min, RADAR_RANGES.LogP.max),
    TPSA: normalize(tpsa, RADAR_RANGES.TPSA.min, RADAR_RANGES.TPSA.max),
    HBD: normalize(hbd, RADAR_RANGES.HBD.min, RADAR_RANGES.HBD.max),
    HBA: normalize(hba, RADAR_RANGES.HBA.min, RADAR_RANGES.HBA.max),
    RotBonds: normalize(rotbonds, RADAR_RANGES.RotBonds.min, RADAR_RANGES.RotBonds.max),
  };
}

function buildOverlayRadarData(molecules: PinnedMolecule[]) {
  const axes = Object.keys(RADAR_RANGES);
  const perMolecule = molecules.map(buildRadarData);

  return axes.map((axis) => {
    const entry: Record<string, string | number> = { axis };
    molecules.forEach((mol, i) => {
      entry[mol.label] = perMolecule[i][axis as keyof ReturnType<typeof buildRadarData>] ?? 0;
    });
    return entry;
  });
}

interface DeltaProperty {
  key: string;
  label: string;
  getValue: (mol: PinnedMolecule) => number | null;
  decimals: number;
}

const DELTA_PROPERTIES: DeltaProperty[] = [
  {
    key: 'mw', label: 'MW', decimals: 1,
    getValue: (mol) => {
      const d = mol.profile.stars.details.find((x) =>
        x.property.toLowerCase().includes('mw') || x.property.toLowerCase().includes('mol weight')
      );
      return d?.value ?? null;
    },
  },
  {
    key: 'logp', label: 'LogP', decimals: 2,
    getValue: (mol) => mol.profile.consensus_logp?.consensus_logp ?? mol.profile.pfi?.clogp ?? null,
  },
  {
    key: 'tpsa', label: 'TPSA', decimals: 1,
    getValue: (mol) => {
      const d = mol.profile.stars.details.find((x) => x.property.toLowerCase().includes('tpsa'));
      return d?.value ?? null;
    },
  },
  {
    key: 'hbd', label: 'HBD', decimals: 0,
    getValue: (mol) => {
      const d = mol.profile.stars.details.find((x) =>
        x.property.toLowerCase().includes('hbd') || x.property.toLowerCase().includes('h-bond donor')
      );
      return d?.value ?? null;
    },
  },
  {
    key: 'hba', label: 'HBA', decimals: 0,
    getValue: (mol) => {
      const d = mol.profile.stars.details.find((x) =>
        x.property.toLowerCase().includes('hba') || x.property.toLowerCase().includes('h-bond acceptor')
      );
      return d?.value ?? null;
    },
  },
  {
    key: 'rotbonds', label: 'Rot. Bonds', decimals: 0,
    getValue: (mol) => {
      const d = mol.profile.stars.details.find((x) => x.property.toLowerCase().includes('rotat'));
      return d?.value ?? null;
    },
  },
  {
    key: 'pfi', label: 'PFI', decimals: 2,
    getValue: (mol) => mol.profile.pfi?.pfi ?? null,
  },
  {
    key: 'abbott', label: 'Abbott %', decimals: 0,
    getValue: (mol) => mol.profile.abbott?.probability_pct ?? null,
  },
];

function formatDelta(delta: number, decimals: number): string {
  const rounded = delta.toFixed(decimals);
  return delta >= 0 ? `+${rounded}` : rounded;
}

function getDruglikenessPassCount(mol: PinnedMolecule): { pass: number; total: number } {
  const dl = mol.profile.druglikeness as Record<string, { passed: boolean }> | undefined;
  if (!dl) return { pass: 0, total: 0 };
  const entries = Object.values(dl).filter((v) => v && typeof v === 'object' && 'passed' in v);
  const pass = entries.filter((v) => v.passed).length;
  return { pass, total: entries.length };
}

/**
 * ComparisonView — right sidebar panel for side-by-side molecule comparison.
 *
 * Slides in from the right with a semi-transparent backdrop. Panel width
 * grows with molecule count (min 480px for 2, up to full width for 5).
 */
export function ComparisonView({ molecules, onClose }: ComparisonViewProps) {
  const n = molecules.length;
  const radarData = buildOverlayRadarData(molecules);

  // Panel width scales with molecule count
  const panelWidth = Math.min(n * 320 + 80, 1400);

  return (
    <div className="fixed inset-x-0 bottom-0 top-[76px] z-40 flex justify-end">
      {/* Backdrop — fades in smoothly */}
      <motion.div
        className="absolute inset-0 bg-black/50"
        initial={{ opacity: 0 }}
        animate={{ opacity: 1 }}
        exit={{ opacity: 0 }}
        transition={{ duration: 0.4, ease: 'easeOut' }}
        onClick={onClose}
      />

      {/* Sidebar panel — smooth slide from right with spring physics */}
      <motion.div
        className="relative h-full overflow-y-auto bg-[var(--color-surface-base)] border-l border-[var(--color-border)] shadow-[-8px_0_30px_-12px_rgba(0,0,0,0.25)]"
        style={{ width: `min(${panelWidth}px, 92vw)` }}
        initial={{ x: '100%', opacity: 0.5 }}
        animate={{ x: 0, opacity: 1 }}
        exit={{ x: '100%', opacity: 0 }}
        transition={{
          x: { type: 'spring', stiffness: 300, damping: 30, mass: 0.8 },
          opacity: { duration: 0.25, ease: 'easeOut' },
        }}
      >
        {/* Header — sticky at top with gradient background */}
        <div className="sticky top-0 z-10 px-6 py-4 border-b border-zinc-700 bg-gradient-to-r from-zinc-800 to-zinc-700">
          <div className="flex items-center justify-between">
            <h2 className="text-lg font-semibold font-display text-white">
              Comparison ({n} molecules)
            </h2>
            <ClayButton
              variant="primary"
              size="sm"
              leftIcon={<X className="w-4 h-4" />}
              onClick={onClose}
            >
              Close
            </ClayButton>
          </div>
        </div>

        {/* Scrollable content */}
        <div className="p-6 space-y-6">
          {/* Per-molecule cards */}
          <div
            className="grid gap-4"
            style={{ gridTemplateColumns: `repeat(${n}, 1fr)` }}
          >
            {molecules.map((mol, i) => {
              const { pass, total } = getDruglikenessPassCount(mol);
              return (
                <div
                  key={mol.smiles}
                  className="rounded-xl border border-[var(--color-border)] bg-[var(--color-surface-elevated)] p-4 space-y-3"
                >
                  {/* Molecule number + label */}
                  <div className="flex items-center gap-2">
                    <span
                      className="w-5 h-5 rounded-full flex items-center justify-center text-white text-xs font-semibold shrink-0"
                      style={{ backgroundColor: MOLECULE_COLORS[i] ?? '#6b7280' }}
                    >
                      {i + 1}
                    </span>
                    <span
                      className="text-xs font-mono text-[var(--color-text-secondary)] truncate"
                      title={mol.smiles}
                    >
                      {mol.smiles}
                    </span>
                  </div>

                  {/* 2D structure — white background for proper contrast */}
                  <div className="w-full rounded-lg overflow-hidden bg-white border border-[var(--color-border)]">
                    <MoleculeViewer smiles={mol.smiles} width={280} height={160} />
                  </div>

                  {/* Key scores */}
                  <div className="space-y-1.5">
                    <div className="flex justify-between text-xs">
                      <span className="text-[var(--color-text-muted)]">PFI</span>
                      <span className="font-semibold tabular-nums text-[var(--color-text-primary)]">
                        {mol.profile.pfi?.pfi?.toFixed(2) ?? '—'}
                        {mol.profile.pfi?.risk && (
                          <span className={`ml-1 text-[10px] ${
                            mol.profile.pfi.risk === 'low' ? 'text-green-600' :
                            mol.profile.pfi.risk === 'moderate' ? 'text-amber-600' : 'text-red-600'
                          }`}>
                            {mol.profile.pfi.risk}
                          </span>
                        )}
                      </span>
                    </div>
                    <div className="flex justify-between text-xs">
                      <span className="text-[var(--color-text-muted)]">Abbott %</span>
                      <span className="font-semibold tabular-nums text-[var(--color-text-primary)]">
                        {mol.profile.abbott?.probability_pct != null ? `${mol.profile.abbott.probability_pct}%` : '—'}
                      </span>
                    </div>
                    <div className="flex justify-between text-xs">
                      <span className="text-[var(--color-text-muted)]">SA Score</span>
                      <span className="font-semibold tabular-nums text-[var(--color-text-primary)]">
                        {mol.profile.sa_comparison?.sa_score?.score != null
                          ? mol.profile.sa_comparison.sa_score.score.toFixed(2) : '—'}
                      </span>
                    </div>
                    {total > 0 && (
                      <div className="flex justify-between text-xs">
                        <span className="text-[var(--color-text-muted)]">Drug-likeness</span>
                        <span className={`font-semibold tabular-nums ${
                          pass === total ? 'text-green-600' : pass >= total / 2 ? 'text-amber-600' : 'text-red-600'
                        }`}>
                          {pass}/{total} rules
                        </span>
                      </div>
                    )}
                  </div>
                </div>
              );
            })}
          </div>

          {/* Overlaid radar chart */}
          <div className="rounded-xl border border-[var(--color-border)] bg-[var(--color-surface-elevated)] p-5">
            <h3 className="text-sm font-semibold font-display text-[var(--color-text-primary)] mb-4">
              Property Radar (normalized)
            </h3>
            <ResponsiveContainer width="100%" height={280} minWidth={200} minHeight={200}>
              <RadarChart data={radarData}>
                <PolarGrid stroke="var(--color-border)" />
                <PolarAngleAxis
                  dataKey="axis"
                  tick={{ fontSize: 11, fill: 'var(--color-text-secondary)' }}
                />
                {molecules.map((mol, i) => (
                  <Radar
                    key={mol.smiles}
                    name={mol.label}
                    dataKey={mol.label}
                    stroke={MOLECULE_COLORS[i] ?? '#6b7280'}
                    fill={MOLECULE_COLORS[i] ?? '#6b7280'}
                    fillOpacity={0.12}
                    strokeWidth={2}
                  />
                ))}
                <Legend wrapperStyle={{ fontSize: '11px' }} />
              </RadarChart>
            </ResponsiveContainer>
          </div>

          {/* Delta table */}
          <div className="rounded-xl border border-[var(--color-border)] bg-[var(--color-surface-elevated)] p-5 overflow-x-auto">
            <h3 className="text-sm font-semibold font-display text-[var(--color-text-primary)] mb-1">
              Delta from molecule 1
            </h3>
            <p className="text-xs text-[var(--color-text-muted)] mb-4">
              Differences relative to {molecules[0]?.smiles?.substring(0, 20)}.
            </p>
            <table className="w-full text-xs">
              <thead>
                <tr className="border-b border-[var(--color-border)]">
                  <th className="text-left py-2 pr-3 text-[var(--color-text-muted)] font-medium">Property</th>
                  {molecules.map((mol, i) => (
                    <th
                      key={mol.smiles}
                      className="text-right py-2 px-2 font-medium"
                      style={{ color: MOLECULE_COLORS[i] ?? '#6b7280' }}
                    >
                      {mol.smiles.substring(0, 15)}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {DELTA_PROPERTIES.map((prop) => {
                  const values = molecules.map((mol) => prop.getValue(mol));
                  const base = values[0];
                  return (
                    <tr key={prop.key} className="border-b border-[var(--color-border)]/50">
                      <td className="py-2 pr-3 text-[var(--color-text-secondary)]">{prop.label}</td>
                      {values.map((val, i) => {
                        if (val === null) {
                          return <td key={i} className="text-right py-2 px-2 text-[var(--color-text-muted)]">—</td>;
                        }
                        if (i === 0) {
                          return (
                            <td key={i} className="text-right py-2 px-2 text-[var(--color-text-primary)] tabular-nums font-medium">
                              {val.toFixed(prop.decimals)}
                            </td>
                          );
                        }
                        const delta = base !== null ? val - base : null;
                        return (
                          <td
                            key={i}
                            className="text-right py-2 px-2 tabular-nums"
                            style={{
                              color: delta === null ? 'var(--color-text-muted)'
                                : delta === 0 ? 'var(--color-text-secondary)'
                                : delta > 0 ? '#16a34a' : '#dc2626',
                            }}
                          >
                            {delta !== null ? formatDelta(delta, prop.decimals) : '—'}
                          </td>
                        );
                      })}
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>

          {/* Bottom padding to clear ComparisonBar */}
          <div className="h-16" />
        </div>
      </motion.div>
    </div>
  );
}
