import { motion } from 'framer-motion';
import { Pin } from 'lucide-react';
import {
  RadarChart,
  PolarGrid,
  PolarAngleAxis,
  PolarRadiusAxis,
  Radar,
  ResponsiveContainer,
} from 'recharts';
import { MoleculeViewer } from '../molecules/MoleculeViewer';
import { ClayButton } from '../ui/ClayButton';
import { ClayCard } from '../ui/ClayCard';
import type { ProfileResponse } from '../../types/profiler';

interface HeroSectionProps {
  smiles: string;
  profile: ProfileResponse;
  onPin?: () => void;
}

/**
 * Drug-like property ranges for radar normalization.
 * Each axis is normalized to [0, 1] within these bounds.
 */
const RADAR_RANGES: Record<string, { min: number; max: number; label: string }> = {
  MW:        { min: 0,   max: 500,  label: 'MW' },
  LogP:      { min: -2,  max: 5,    label: 'LogP' },
  TPSA:      { min: 0,   max: 140,  label: 'TPSA' },
  HBD:       { min: 0,   max: 5,    label: 'HBD' },
  HBA:       { min: 0,   max: 10,   label: 'HBA' },
  RotBonds:  { min: 0,   max: 10,   label: 'Rot. Bonds' },
};

/** Normalize a value to the [0, 1] range, clamped. */
function normalize(value: number, min: number, max: number): number {
  return Math.min(1, Math.max(0, (value - min) / (max - min)));
}

/** Extract radar axis values from the profile response. */
function buildRadarData(profile: ProfileResponse) {
  // Extract property values from the stars details when available,
  // falling back to profiler fields and consensus LogP.
  const detailMap: Record<string, number> = {};
  for (const d of profile.stars.details) {
    detailMap[d.property.toLowerCase()] = d.value;
  }

  const mw = detailMap['mw'] ?? detailMap['mol weight'] ?? detailMap['molecular weight'] ?? 0;
  const logp = profile.consensus_logp?.consensus_logp ?? profile.pfi?.clogp ?? 0;
  const tpsa = detailMap['tpsa'] ?? profile.abbott?.tpsa ?? 0;
  const hbd = detailMap['hbd'] ?? detailMap['h bond donors'] ?? detailMap['hbond donors'] ?? 0;
  const hba = detailMap['hba'] ?? detailMap['h bond acceptors'] ?? detailMap['hbond acceptors'] ?? 0;
  const rotbonds = detailMap['rotatable bonds'] ?? detailMap['rot bonds'] ?? detailMap['rotbonds'] ?? 0;

  return [
    { axis: 'MW',       label: 'MW',          value: normalize(mw,       RADAR_RANGES.MW.min,       RADAR_RANGES.MW.max) },
    { axis: 'LogP',     label: 'LogP',        value: normalize(logp,     RADAR_RANGES.LogP.min,     RADAR_RANGES.LogP.max) },
    { axis: 'TPSA',     label: 'TPSA',        value: normalize(tpsa,     RADAR_RANGES.TPSA.min,     RADAR_RANGES.TPSA.max) },
    { axis: 'HBD',      label: 'HBD',         value: normalize(hbd,      RADAR_RANGES.HBD.min,      RADAR_RANGES.HBD.max) },
    { axis: 'HBA',      label: 'HBA',         value: normalize(hba,      RADAR_RANGES.HBA.min,      RADAR_RANGES.HBA.max) },
    { axis: 'RotBonds', label: 'Rot. Bonds',  value: normalize(rotbonds, RADAR_RANGES.RotBonds.min, RADAR_RANGES.RotBonds.max) },
  ];
}

/**
 * HeroSection — primary visual anchor for the Compound Profiler page.
 *
 * Left column: 2D molecule structure with pin button.
 * Right column: 6-axis property radar chart (MW, LogP, TPSA, HBD, HBA, RotBonds).
 *
 * Layout per UI-SPEC D-02:
 * - lg+: grid grid-cols-[2fr_3fr] gap-8 (side-by-side)
 * - md/sm: stacked (grid-cols-1)
 */
export function HeroSection({ smiles, profile, onPin }: HeroSectionProps) {
  const radarData = buildRadarData(profile);

  // Derive a short display name for the pin aria-label
  const displayName = smiles.length > 20 ? `${smiles.slice(0, 20)}…` : smiles;

  return (
    <motion.div
      className="grid grid-cols-1 lg:grid-cols-[2fr_3fr] gap-4 lg:gap-8"
      initial={{ opacity: 0, y: 20 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.6, ease: 'easeOut' }}
    >
      {/* Left column — 2D structure */}
      <ClayCard size="md" className="relative">
        {/* Pin button — top-right corner (D-20) */}
        {onPin && (
          <div className="absolute top-3 right-3 z-10">
            <ClayButton
              variant="ghost"
              size="icon"
              onClick={onPin}
              aria-label={`Pin ${displayName} for comparison`}
            >
              <Pin className="w-4 h-4" />
            </ClayButton>
          </div>
        )}

        <MoleculeViewer
          smiles={smiles}
          width={400}
          height={280}
          className="w-full"
        />
      </ClayCard>

      {/* Right column — 6-axis property radar */}
      <ClayCard size="md">
        <p className="text-sm font-semibold text-text-secondary font-display mb-3">
          Property Profile
        </p>
        <p className="text-xs text-text-muted mb-4">
          Normalized to drug-like ranges: MW (0–500), LogP (−2 to 5), TPSA (0–140), HBD (0–5), HBA (0–10), RotBonds (0–10)
        </p>

        <div className="w-full" style={{ height: 260 }}>
          <ResponsiveContainer width="100%" height="100%">
            <RadarChart
              data={radarData}
              cx="50%"
              cy="50%"
              outerRadius="70%"
            >
              <PolarGrid
                stroke="var(--color-border)"
                strokeOpacity={0.4}
              />
              <PolarAngleAxis
                dataKey="label"
                tick={{
                  fontSize: 12,
                  fill: 'var(--color-text-secondary)',
                  fontFamily: 'Outfit, sans-serif',
                }}
              />
              <PolarRadiusAxis
                domain={[0, 1]}
                tick={false}
                axisLine={false}
                tickCount={5}
              />
              <Radar
                name="Properties"
                dataKey="value"
                stroke="var(--color-primary)"
                fill="var(--color-primary)"
                fillOpacity={0.3}
                strokeWidth={2}
                dot={{ r: 3, fill: 'var(--color-primary)', strokeWidth: 0 }}
              />
            </RadarChart>
          </ResponsiveContainer>
        </div>

        <p className="text-xs text-text-muted text-center mt-2">
          Outer edge = drug-like upper bound. Values outside ideal range are visible beyond the boundary.
        </p>
      </ClayCard>
    </motion.div>
  );
}
