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
import { ClayCard } from '../ui/ClayCard';
import type { PinnedMolecule } from './ComparisonBar';

interface ComparisonViewProps {
  molecules: PinnedMolecule[];
  onClose: () => void;
}

/**
 * Drug-like property ranges for radar normalization (0-1).
 * Same axes as HeroSection.
 */
const RADAR_RANGES: Record<string, { min: number; max: number; label: string }> = {
  MW:       { min: 0,  max: 500,  label: 'MW' },
  LogP:     { min: -2, max: 5,    label: 'LogP' },
  TPSA:     { min: 0,  max: 140,  label: 'TPSA' },
  HBD:      { min: 0,  max: 5,    label: 'HBD' },
  HBA:      { min: 0,  max: 10,   label: 'HBA' },
  RotBonds: { min: 0,  max: 10,   label: 'Rot. Bonds' },
};

/**
 * Per-molecule colors: first molecule uses primary crimson.
 */
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

/**
 * Build radar data point for a single molecule.
 */
function buildRadarData(mol: PinnedMolecule) {
  const profile = mol.profile;
  const detailMap: Record<string, number> = {};
  for (const d of profile.stars.details) {
    detailMap[d.property.toLowerCase()] = d.value;
  }

  const mw = detailMap['mw'] ?? detailMap['mol weight'] ?? detailMap['molecular weight'] ?? 0;
  const logp = profile.consensus_logp?.consensus_logp ?? profile.pfi?.clogp ?? 0;
  const tpsa = detailMap['tpsa'] ?? 0;
  const hbd = detailMap['hbd'] ?? detailMap['h-bond donors'] ?? 0;
  const hba = detailMap['hba'] ?? detailMap['h-bond acceptors'] ?? 0;
  const rotbonds = detailMap['rotatable bonds'] ?? detailMap['rotbonds'] ?? 0;

  return {
    MW: normalize(mw, RADAR_RANGES.MW.min, RADAR_RANGES.MW.max),
    LogP: normalize(logp, RADAR_RANGES.LogP.min, RADAR_RANGES.LogP.max),
    TPSA: normalize(tpsa, RADAR_RANGES.TPSA.min, RADAR_RANGES.TPSA.max),
    HBD: normalize(hbd, RADAR_RANGES.HBD.min, RADAR_RANGES.HBD.max),
    HBA: normalize(hba, RADAR_RANGES.HBA.min, RADAR_RANGES.HBA.max),
    RotBonds: normalize(rotbonds, RADAR_RANGES.RotBonds.min, RADAR_RANGES.RotBonds.max),
  };
}

/**
 * Build overlay radar chart data (one entry per radar axis).
 * Each entry has values for all molecules.
 */
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

/**
 * Key properties for the delta table.
 */
interface DeltaProperty {
  key: string;
  label: string;
  getValue: (mol: PinnedMolecule) => number | null;
  decimals: number;
}

const DELTA_PROPERTIES: DeltaProperty[] = [
  {
    key: 'mw',
    label: 'MW',
    getValue: (mol) => {
      const d = mol.profile.stars.details.find((x) =>
        x.property.toLowerCase().includes('mw') ||
        x.property.toLowerCase().includes('mol weight') ||
        x.property.toLowerCase().includes('molecular weight')
      );
      return d?.value ?? null;
    },
    decimals: 1,
  },
  {
    key: 'logp',
    label: 'LogP',
    getValue: (mol) => mol.profile.consensus_logp?.consensus_logp ?? mol.profile.pfi?.clogp ?? null,
    decimals: 2,
  },
  {
    key: 'tpsa',
    label: 'TPSA',
    getValue: (mol) => {
      const d = mol.profile.stars.details.find((x) => x.property.toLowerCase().includes('tpsa'));
      return d?.value ?? null;
    },
    decimals: 1,
  },
  {
    key: 'hbd',
    label: 'HBD',
    getValue: (mol) => {
      const d = mol.profile.stars.details.find((x) =>
        x.property.toLowerCase().includes('hbd') ||
        x.property.toLowerCase().includes('h-bond donor')
      );
      return d?.value ?? null;
    },
    decimals: 0,
  },
  {
    key: 'hba',
    label: 'HBA',
    getValue: (mol) => {
      const d = mol.profile.stars.details.find((x) =>
        x.property.toLowerCase().includes('hba') ||
        x.property.toLowerCase().includes('h-bond acceptor')
      );
      return d?.value ?? null;
    },
    decimals: 0,
  },
  {
    key: 'rotbonds',
    label: 'Rot. Bonds',
    getValue: (mol) => {
      const d = mol.profile.stars.details.find((x) =>
        x.property.toLowerCase().includes('rotat')
      );
      return d?.value ?? null;
    },
    decimals: 0,
  },
  {
    key: 'pfi',
    label: 'PFI',
    getValue: (mol) => mol.profile.pfi?.pfi ?? null,
    decimals: 2,
  },
  {
    key: 'abbott',
    label: 'Abbott %',
    getValue: (mol) => mol.profile.abbott?.probability_pct ?? null,
    decimals: 0,
  },
];

function formatDelta(delta: number, decimals: number): string {
  const rounded = delta.toFixed(decimals);
  return delta >= 0 ? `+${rounded}` : rounded;
}

/**
 * Drug-likeness pass count for a molecule (simple count from druglikeness field).
 */
function getDruglikenessPassCount(mol: PinnedMolecule): { pass: number; total: number } {
  const dl = mol.profile.druglikeness as Record<string, { passed: boolean }> | undefined;
  if (!dl) return { pass: 0, total: 0 };
  const entries = Object.values(dl).filter((v) => v && typeof v === 'object' && 'passed' in v);
  const pass = entries.filter((v) => v.passed).length;
  return { pass, total: entries.length };
}

/**
 * ComparisonView — side-by-side molecule comparison.
 *
 * Per D-21, D-22: columns per molecule, overlaid radar chart, delta table.
 */
export function ComparisonView({ molecules, onClose }: ComparisonViewProps) {
  const n = molecules.length;
  const radarData = buildOverlayRadarData(molecules);

  return (
    <motion.div
      className="mt-8"
      initial={{ opacity: 0, y: 16 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.4, ease: 'easeOut' }}
    >
      {/* Header */}
      <div className="flex items-center justify-between mb-6">
        <h2 className="text-2xl font-semibold font-display text-text-primary">
          Comparison ({n} molecules)
        </h2>
        <ClayButton
          variant="ghost"
          size="sm"
          leftIcon={<X className="w-4 h-4" />}
          onClick={onClose}
        >
          Exit comparison
        </ClayButton>
      </div>

      {/* Per-molecule columns */}
      <div
        className={`grid gap-4 overflow-x-auto`}
        style={{ gridTemplateColumns: `repeat(${n}, minmax(200px, 1fr))` }}
      >
        {molecules.map((mol, i) => {
          const { pass, total } = getDruglikenessPassCount(mol);
          return (
            <ClayCard key={mol.smiles} size="sm" className="min-w-[200px]">
              {/* Molecule number + label */}
              <div className="flex items-center gap-2 mb-3">
                <span
                  className="w-5 h-5 rounded-full flex items-center justify-center text-white text-xs font-semibold flex-shrink-0"
                  style={{ backgroundColor: MOLECULE_COLORS[i] ?? '#6b7280' }}
                >
                  {i + 1}
                </span>
                <span
                  className="text-xs font-mono text-text-secondary truncate"
                  title={mol.smiles}
                >
                  {mol.label}
                </span>
              </div>

              {/* 2D structure */}
              <div className="w-full h-32 mb-3 rounded-lg overflow-hidden bg-surface-sunken">
                <MoleculeViewer smiles={mol.smiles} width={220} height={128} />
              </div>

              {/* Key scores */}
              <div className="space-y-2">
                <div className="flex justify-between text-xs">
                  <span className="text-text-muted">PFI</span>
                  <span className="font-semibold tabular-nums text-text-primary">
                    {mol.profile.pfi?.pfi?.toFixed(2) ?? '—'}
                    <span
                      className="ml-1 text-[10px]"
                      style={{
                        color:
                          mol.profile.pfi?.risk === 'low'
                            ? 'var(--color-score-good, #d97706)'
                            : mol.profile.pfi?.risk === 'moderate'
                            ? '#ea580c'
                            : '#dc2626',
                      }}
                    >
                      {mol.profile.pfi?.risk}
                    </span>
                  </span>
                </div>
                <div className="flex justify-between text-xs">
                  <span className="text-text-muted">Abbott %</span>
                  <span className="font-semibold tabular-nums text-text-primary">
                    {mol.profile.abbott?.probability_pct != null
                      ? `${mol.profile.abbott.probability_pct}%`
                      : '—'}
                  </span>
                </div>
                <div className="flex justify-between text-xs">
                  <span className="text-text-muted">SA Score</span>
                  <span className="font-semibold tabular-nums text-text-primary">
                    {mol.profile.sa_comparison?.sa_score?.score != null
                      ? mol.profile.sa_comparison.sa_score.score.toFixed(2)
                      : '—'}
                  </span>
                </div>
                {total > 0 && (
                  <div className="flex justify-between text-xs">
                    <span className="text-text-muted">Drug-likeness</span>
                    <span
                      className="font-semibold tabular-nums"
                      style={{
                        color: pass === total ? '#16a34a' : pass >= total / 2 ? '#d97706' : '#dc2626',
                      }}
                    >
                      {pass}/{total} rules
                    </span>
                  </div>
                )}
              </div>
            </ClayCard>
          );
        })}
      </div>

      {/* Overlaid radar chart */}
      <ClayCard className="mt-6">
        <h3 className="text-base font-semibold font-display text-text-primary mb-4">
          Property Radar (normalized)
        </h3>
        <ResponsiveContainer width="100%" height={300}>
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
            <Legend
              wrapperStyle={{ fontSize: '11px', color: 'var(--color-text-secondary)' }}
            />
          </RadarChart>
        </ResponsiveContainer>
      </ClayCard>

      {/* Delta table */}
      <ClayCard className="mt-6 overflow-x-auto">
        <h3 className="text-base font-semibold font-display text-text-primary mb-1">
          Delta from molecule 1
        </h3>
        <p className="text-xs text-text-muted mb-4">
          Differences relative to {molecules[0]?.label}.
        </p>
        <table className="w-full text-xs min-w-[400px]">
          <thead>
            <tr className="border-b border-border">
              <th className="text-left py-2 pr-3 text-text-muted font-medium">Property</th>
              {molecules.map((mol, i) => (
                <th
                  key={mol.smiles}
                  className="text-right py-2 px-2 font-medium"
                  style={{ color: MOLECULE_COLORS[i] ?? '#6b7280' }}
                >
                  {mol.label}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {DELTA_PROPERTIES.map((prop) => {
              const values = molecules.map((mol) => prop.getValue(mol));
              const base = values[0];
              return (
                <tr key={prop.key} className="border-b border-border/50 hover:bg-surface-sunken/50">
                  <td className="py-2 pr-3 text-text-secondary">{prop.label}</td>
                  {values.map((val, i) => {
                    if (val === null) {
                      return (
                        <td key={i} className="text-right py-2 px-2 text-text-muted">
                          —
                        </td>
                      );
                    }
                    if (i === 0) {
                      return (
                        <td key={i} className="text-right py-2 px-2 text-text-primary tabular-nums font-medium">
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
                          color:
                            delta === null
                              ? 'var(--color-text-muted)'
                              : delta === 0
                              ? 'var(--color-text-secondary)'
                              : delta > 0
                              ? '#16a34a'
                              : '#dc2626',
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
      </ClayCard>
    </motion.div>
  );
}
