import { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import {
  Shield, FlaskConical, Puzzle, Binary,
  ChevronDown, AlertTriangle, Sparkles, Check, X, Info,
} from 'lucide-react';
import type { MLReadinessResult, MLDimension, MLDimensionItem } from '../../types/scoring';
import { ScoreChart } from './ScoreChart';
import { InfoTooltip } from '../ui/Tooltip';
import { cn } from '../../lib/utils';

/** Format a score: integers display without decimals, fractional values get 1 decimal place. */
function fmt(n: number): string {
  return Number.isInteger(n) ? String(n) : n.toFixed(1);
}

interface MLReadinessScoreProps {
  result: MLReadinessResult;
  /** When true, hides the header/score chart and only shows breakdown */
  breakdownOnly?: boolean;
}

/** Color config by score thresholds */
function getScoreColor(score: number, maxScore: number) {
  const pct = maxScore > 0 ? (score / maxScore) * 100 : 0;
  if (pct >= 85) return {
    gradient: 'from-emerald-500 to-green-400',
    bg: 'bg-emerald-500/10',
    bgRaw: 'emerald',
    text: 'text-emerald-600 dark:text-emerald-400',
    border: 'border-emerald-500/20',
    ring: 'ring-emerald-500/20',
    stroke: '#10b981',
    track: 'rgba(16, 185, 129, 0.15)',
  };
  if (pct >= 70) return {
    gradient: 'from-teal-500 to-cyan-400',
    bg: 'bg-teal-500/10',
    bgRaw: 'teal',
    text: 'text-teal-600 dark:text-teal-400',
    border: 'border-teal-500/20',
    ring: 'ring-teal-500/20',
    stroke: '#14b8a6',
    track: 'rgba(20, 184, 166, 0.15)',
  };
  if (pct >= 50) return {
    gradient: 'from-amber-500 to-yellow-400',
    bg: 'bg-amber-500/10',
    bgRaw: 'amber',
    text: 'text-amber-600 dark:text-amber-400',
    border: 'border-amber-500/20',
    ring: 'ring-amber-500/20',
    stroke: '#f59e0b',
    track: 'rgba(245, 158, 11, 0.15)',
  };
  if (pct >= 30) return {
    gradient: 'from-orange-500 to-orange-400',
    bg: 'bg-orange-500/10',
    bgRaw: 'orange',
    text: 'text-orange-600 dark:text-orange-400',
    border: 'border-orange-500/20',
    ring: 'ring-orange-500/20',
    stroke: '#f97316',
    track: 'rgba(249, 115, 22, 0.15)',
  };
  return {
    gradient: 'from-red-500 to-red-400',
    bg: 'bg-red-500/10',
    bgRaw: 'red',
    text: 'text-red-600 dark:text-red-400',
    border: 'border-red-500/20',
    ring: 'ring-red-500/20',
    stroke: '#ef4444',
    track: 'rgba(239, 68, 68, 0.15)',
  };
}

/** Tier-based guidance — visible at a glance, no tooltip needed */
const TIER_GUIDANCE: Record<string, string> = {
  Excellent: 'Suitable for most ML workflows without modification.',
  Good: 'Suitable for ML with minor preprocessing.',
  Moderate: 'Usable, but may need filtering or augmentation in ML pipelines.',
  Limited: 'Likely to cause issues in ML models — consider alternatives.',
  Poor: 'Not recommended for ML use without significant preprocessing.',
};

/** Short labels for dimension tags */
const DIM_SHORT: Record<string, string> = {
  structural_quality: 'Structure',
  property_profile: 'Properties',
  complexity_feasibility: 'Complexity',
  representation_quality: 'Representation',
};

/** Dimension icon mapping */
const DIMENSION_ICONS: Record<string, React.ReactNode> = {
  structural_quality: <Shield className="w-5 h-5" />,
  property_profile: <FlaskConical className="w-5 h-5" />,
  complexity_feasibility: <Puzzle className="w-5 h-5" />,
  representation_quality: <Binary className="w-5 h-5" />,
};

/** Dimension tooltips with What/How/Why/Citation */
const DIMENSION_TOOLTIPS: Record<string, { title: string; content: React.ReactNode }> = {
  structural_quality: {
    title: 'Structural Quality (20 pts)',
    content: (
      <div className="text-xs space-y-2">
        <p><strong>What:</strong> Whether the molecule is structurally clean for ML pipelines.</p>
        <p><strong>How:</strong> Binary checks: single component, organic elements, no radicals, reasonable charge, no dummy atoms.</p>
        <p><strong>Why:</strong> Multi-component mixtures, metals, and radicals cause descriptor calculation failures and model bias.</p>
        <p className="text-white/50">Lipinski et al. (1997) Adv. Drug Deliv. Rev. 23:3-25</p>
      </div>
    ),
  },
  property_profile: {
    title: 'Property Profile (35 pts)',
    content: (
      <div className="text-xs space-y-2">
        <p><strong>What:</strong> Physicochemical properties scored against QED-derived ideal ranges.</p>
        <p><strong>How:</strong> Desirability function (0-1) for MW, LogP, TPSA, HBD, HBA, rotatable bonds, aromatic rings.</p>
        <p><strong>Why:</strong> ML models trained on drug-like chemical space perform best when test molecules fall within typical property ranges.</p>
        <p className="text-white/50">Bickerton et al. (2012) Nature Chemistry 4:90-98</p>
      </div>
    ),
  },
  complexity_feasibility: {
    title: 'Complexity & Feasibility (25 pts)',
    content: (
      <div className="text-xs space-y-2">
        <p><strong>What:</strong> Drug-likeness quality, synthetic feasibility, 3D character, stereocenter count.</p>
        <p><strong>How:</strong> QED (0-1 scaled to 8pts), SA Score (inverse, 8pts), Fsp3 (desirability, 4pts), stereocenters (5pts with undefined penalty).</p>
        <p><strong>Why:</strong> Overly complex or synthetically intractable molecules are poor ML training data.</p>
        <p className="text-white/50">Ertl &amp; Schuffenhauer (2009) J. Cheminf. 1:8; Lovering et al. (2009) J. Med. Chem. 52:6752</p>
      </div>
    ),
  },
  representation_quality: {
    title: 'Representation Quality (20 pts)',
    content: (
      <div className="text-xs space-y-2">
        <p><strong>What:</strong> How well the molecule can be represented computationally.</p>
        <p><strong>How:</strong> Descriptor completeness (451 RDKit), fingerprint generation (7 types), bit density (1-30% ideal), conformer generation (ETKDGv3).</p>
        <p><strong>Why:</strong> Molecules that fail descriptor/FP calculations produce NaN values that break ML models. Conformer generation tests 3D feasibility.</p>
        <p className="text-white/50">Rogers &amp; Hahn (2010) J. Chem. Inf. Model. 50:742-754</p>
      </div>
    ),
  },
};

/** Mini circular score ring for dimension cards */
function MiniScoreRing({
  score,
  maxScore,
  color,
  delay = 0,
}: {
  score: number;
  maxScore: number;
  color: ReturnType<typeof getScoreColor>;
  delay?: number;
}) {
  const pct = maxScore > 0 ? score / maxScore : 0;
  const radius = 22;
  const circumference = 2 * Math.PI * radius;
  const strokeDashoffset = circumference * (1 - pct);

  return (
    <div className="relative w-14 h-14 flex-shrink-0">
      <svg className="w-14 h-14 -rotate-90" viewBox="0 0 52 52">
        <circle
          cx="26" cy="26" r={radius}
          fill="none"
          stroke={color.track}
          strokeWidth="4"
        />
        <motion.circle
          cx="26" cy="26" r={radius}
          fill="none"
          stroke={color.stroke}
          strokeWidth="4"
          strokeLinecap="round"
          strokeDasharray={circumference}
          initial={{ strokeDashoffset: circumference }}
          animate={{ strokeDashoffset }}
          transition={{ duration: 1, delay: delay + 0.3, ease: 'easeOut' }}
        />
      </svg>
      <div className="absolute inset-0 flex flex-col items-center justify-center">
        <span className={cn('text-sm font-bold leading-none', color.text)}>
          {fmt(score)}
        </span>
        <span className="text-[9px] text-[var(--color-text-muted)] leading-tight">
          /{fmt(maxScore)}
        </span>
      </div>
    </div>
  );
}

/** Expandable dimension card */
function DimensionCard({
  dimKey,
  dimension,
  delay = 0,
}: {
  dimKey: string;
  dimension: MLDimension;
  delay?: number;
}) {
  const [expanded, setExpanded] = useState(false);
  const color = getScoreColor(dimension.score, dimension.max_score);
  const icon = DIMENSION_ICONS[dimKey];
  const tooltip = DIMENSION_TOOLTIPS[dimKey];
  const percentage = dimension.max_score > 0 ? (dimension.score / dimension.max_score) * 100 : 0;

  return (
    <motion.div
      initial={{ opacity: 0, y: 20 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.4, delay }}
      className={cn(
        'group relative rounded-2xl overflow-hidden',
        'bg-gradient-to-br from-[var(--color-surface-elevated)] to-[var(--color-surface)]',
        'border border-[var(--color-border)]',
        'hover:border-[var(--color-primary)]/30',
        'hover:shadow-lg hover:shadow-[var(--color-primary)]/5',
        'transition-all duration-200'
      )}
    >
      {/* Dual background glow */}
      <div className="absolute inset-0 overflow-hidden rounded-2xl pointer-events-none">
        <div className={cn(
          'absolute -top-16 -right-16 w-40 h-40 rounded-full blur-3xl opacity-[0.12]',
          'transition-opacity duration-300 group-hover:opacity-[0.2]',
          color.bg,
        )} />
        <div className={cn(
          'absolute -bottom-8 -left-8 w-24 h-24 rounded-full blur-2xl opacity-[0.06]',
          color.bg,
        )} />
      </div>

      {/* Header — always visible */}
      <button
        onClick={() => setExpanded(!expanded)}
        className="relative w-full p-4 text-left cursor-pointer"
      >
        <div className="flex items-center gap-3">
          {/* Icon */}
          <div className={cn(
            'w-10 h-10 rounded-xl flex items-center justify-center',
            'transition-all duration-200',
            color.bg, color.text,
            'group-hover:shadow-md',
          )}>
            {icon}
          </div>

          {/* Name + meta */}
          <div className="flex-1 min-w-0">
            <div className="flex items-center gap-1.5">
              <h4 className="font-semibold text-[var(--color-text-primary)] text-sm truncate">
                {dimension.name}
              </h4>
              {tooltip && <InfoTooltip title={tooltip.title} content={tooltip.content} asSpan />}
            </div>
            <p className="text-[11px] text-[var(--color-text-muted)] mt-0.5">
              {dimension.items.length} check{dimension.items.length !== 1 ? 's' : ''}
              <span className="mx-1.5 opacity-30">|</span>
              <span className={cn('font-medium', color.text)}>
                {Math.round(percentage)}%
              </span>
            </p>
          </div>

          {/* Mini score ring */}
          <MiniScoreRing score={dimension.score} maxScore={dimension.max_score} color={color} delay={delay} />

          {/* Chevron */}
          <motion.div
            animate={{ rotate: expanded ? 180 : 0 }}
            transition={{ duration: 0.2 }}
          >
            <ChevronDown className="w-4 h-4 text-[var(--color-text-muted)]" />
          </motion.div>
        </div>

        {/* Progress bar */}
        <div className="mt-3 h-1.5 bg-[var(--color-surface-sunken)] rounded-full overflow-hidden">
          <motion.div
            initial={{ width: 0 }}
            animate={{ width: `${percentage}%` }}
            transition={{ duration: 0.8, delay: delay + 0.2, ease: 'easeOut' }}
            className={cn('h-full rounded-full bg-gradient-to-r', color.gradient)}
          />
        </div>
      </button>

      {/* Expandable items */}
      <AnimatePresence>
        {expanded && dimension.items.length > 0 && (
          <motion.div
            initial={{ height: 0, opacity: 0 }}
            animate={{ height: 'auto', opacity: 1 }}
            exit={{ height: 0, opacity: 0 }}
            transition={{ duration: 0.25, ease: 'easeOut' }}
            className="overflow-hidden"
          >
            <div className="px-4 pb-4 space-y-0.5">
              <div className="border-t border-[var(--color-border)]/60 mb-2" />
              {dimension.items.map((item, idx) => (
                <motion.div
                  key={idx}
                  initial={{ opacity: 0, x: -8 }}
                  animate={{ opacity: 1, x: 0 }}
                  transition={{ duration: 0.2, delay: idx * 0.04 }}
                >
                  <ItemRow item={item} />
                </motion.div>
              ))}
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </motion.div>
  );
}

/** Single scored item row */
function ItemRow({ item }: { item: MLDimensionItem }) {
  const pct = item.max_score > 0 ? (item.score / item.max_score) * 100 : 0;
  const color = getScoreColor(item.score, item.max_score);

  return (
    <div className={cn(
      'flex items-center justify-between py-2 px-2.5 rounded-xl',
      'hover:bg-[var(--color-surface-sunken)]/60 transition-colors duration-150',
    )}>
      <div className="flex items-center gap-2.5 min-w-0">
        {/* Status indicator */}
        {item.passed !== undefined && item.passed !== null ? (
          <div className={cn(
            'w-5 h-5 rounded-md flex items-center justify-center flex-shrink-0',
            item.passed
              ? 'bg-emerald-500/15 text-emerald-500'
              : 'bg-red-500/15 text-red-500',
          )}>
            {item.passed
              ? <Check className="w-3 h-3" strokeWidth={3} />
              : <X className="w-3 h-3" strokeWidth={3} />
            }
          </div>
        ) : (
          <div className="w-5 h-5 rounded-md flex items-center justify-center flex-shrink-0 bg-[var(--color-surface-sunken)]">
            <Info className="w-3 h-3 text-[var(--color-text-muted)]" />
          </div>
        )}

        {/* Name + subtitle */}
        <div className="min-w-0">
          <div className="flex items-center gap-1">
            <span className="text-xs font-medium text-[var(--color-text-secondary)]">{item.name}</span>
            {item.tooltip && (
              <InfoTooltip title={item.name} content={<p className="text-xs">{item.tooltip}</p>} />
            )}
          </div>
          {item.subtitle && (
            <p className="text-[10px] text-[var(--color-text-muted)] truncate font-mono mt-0.5">
              {item.subtitle}
            </p>
          )}
        </div>
      </div>

      {/* Score + bar */}
      <div className="flex items-center gap-2.5 flex-shrink-0">
        <div className="w-20 h-1.5 bg-[var(--color-surface-sunken)] rounded-full overflow-hidden">
          <div
            className={cn('h-full rounded-full bg-gradient-to-r transition-all duration-500', color.gradient)}
            style={{ width: `${pct}%` }}
          />
        </div>
        <span className={cn('text-xs font-semibold tabular-nums min-w-[3.5rem] text-right', color.text)}>
          {fmt(item.score)}/{fmt(item.max_score)}
        </span>
      </div>
    </div>
  );
}

/** Caveats banner */
function CaveatsBanner({ caveats }: { caveats: string[] }) {
  if (caveats.length === 0) return null;

  return (
    <motion.div
      initial={{ opacity: 0, y: -10 }}
      animate={{ opacity: 1, y: 0 }}
      className={cn(
        'relative rounded-xl overflow-hidden',
        'bg-amber-500/5 border border-amber-500/15',
      )}
    >
      {/* Left accent bar */}
      <div className="absolute left-0 top-0 bottom-0 w-1 bg-gradient-to-b from-amber-400 to-amber-600 rounded-l-xl" />

      <div className="flex items-start gap-2.5 p-3 pl-4">
        <div className="w-6 h-6 rounded-lg bg-amber-500/15 flex items-center justify-center flex-shrink-0 mt-0.5">
          <AlertTriangle className="w-3.5 h-3.5 text-amber-500" />
        </div>
        <div className="space-y-1.5">
          <p className="text-xs font-semibold text-amber-600 dark:text-amber-400 tracking-wide uppercase">
            Caveats
          </p>
          <ul className="text-xs text-[var(--color-text-secondary)] space-y-1">
            {caveats.map((caveat, i) => (
              <li key={i} className="flex items-start gap-1.5">
                <span className="w-1 h-1 rounded-full bg-amber-400 mt-1.5 flex-shrink-0" />
                {caveat}
              </li>
            ))}
          </ul>
        </div>
      </div>
    </motion.div>
  );
}

/** Supplementary info pills */
function SupplementaryPills({ supplementary }: { supplementary: Record<string, unknown> }) {
  const pills: { label: string; value: string; ok: boolean }[] = [];

  if (supplementary.lipinski_violations !== undefined) {
    const v = supplementary.lipinski_violations as number;
    pills.push({
      label: 'Lipinski',
      value: `${v} violation${v !== 1 ? 's' : ''}`,
      ok: v <= 1,
    });
  }
  if (supplementary.veber_passed !== undefined) {
    const passed = supplementary.veber_passed as boolean;
    pills.push({
      label: 'Veber',
      value: passed ? 'Pass' : 'Fail',
      ok: passed,
    });
  }

  if (pills.length === 0) return null;

  return (
    <div className="flex flex-wrap gap-2">
      {pills.map((pill) => (
        <span
          key={pill.label}
          className={cn(
            'inline-flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-xs font-medium',
            'transition-colors duration-150',
            pill.ok
              ? 'bg-emerald-500/8 text-emerald-600 dark:text-emerald-400 border border-emerald-500/20 hover:bg-emerald-500/15'
              : 'bg-red-500/8 text-red-600 dark:text-red-400 border border-red-500/20 hover:bg-red-500/15'
          )}
        >
          <div className={cn(
            'w-4 h-4 rounded flex items-center justify-center',
            pill.ok ? 'bg-emerald-500/20' : 'bg-red-500/20',
          )}>
            {pill.ok ? <Check className="w-2.5 h-2.5" strokeWidth={3} /> : <X className="w-2.5 h-2.5" strokeWidth={3} />}
          </div>
          <span className="font-semibold">{pill.label}</span>
          <span className="opacity-70">{pill.value}</span>
        </span>
      ))}
    </div>
  );
}

/**
 * ML-Readiness Score — 4-dimension scientific assessment.
 */
export function MLReadinessScore({ result, breakdownOnly = false }: MLReadinessScoreProps) {
  const { score, label, breakdown, caveats, supplementary, interpretation } = result;

  const dimensions: { key: string; dim: MLDimension }[] = [
    { key: 'structural_quality', dim: breakdown.structural_quality },
    { key: 'property_profile', dim: breakdown.property_profile },
    { key: 'complexity_feasibility', dim: breakdown.complexity_feasibility },
    { key: 'representation_quality', dim: breakdown.representation_quality },
  ];

  const calculation = `Score = Structural Quality (${fmt(breakdown.structural_quality.score)}/${fmt(breakdown.structural_quality.max_score)}) + Property Profile (${fmt(breakdown.property_profile.score)}/${fmt(breakdown.property_profile.max_score)}) + Complexity (${fmt(breakdown.complexity_feasibility.score)}/${fmt(breakdown.complexity_feasibility.max_score)}) + Representation (${fmt(breakdown.representation_quality.score)}/${fmt(breakdown.representation_quality.max_score)}) = ${fmt(score)}`;

  const overallColor = getScoreColor(score, 100);

  if (breakdownOnly) {
    return (
      <div className="space-y-4">
        {/* Interpretation banner */}
        <motion.div
          initial={{ opacity: 0, y: -10 }}
          animate={{ opacity: 1, y: 0 }}
          className={cn(
            'relative overflow-hidden rounded-2xl p-4',
            'bg-gradient-to-r from-[var(--color-primary)]/5 via-[var(--color-accent)]/5 to-[var(--color-primary)]/5',
            'border border-[var(--color-primary)]/10'
          )}
        >
          <div className="flex items-start gap-3">
            <div className="w-8 h-8 rounded-lg bg-[var(--color-primary)]/10 flex items-center justify-center flex-shrink-0">
              <Sparkles className="w-4 h-4 text-[var(--color-primary)]" />
            </div>
            <div className="flex-1 min-w-0">
              <div className="flex items-center gap-2 mb-1">
                <h4 className="text-sm font-semibold text-[var(--color-text-primary)]">
                  ML Readiness: {label}
                </h4>
                <span className={cn(
                  'text-xs font-bold px-2 py-0.5 rounded-md tabular-nums',
                  overallColor.bg, overallColor.text,
                )}>
                  {score}/100
                </span>
                <InfoTooltip
                  title="What is ML Readiness?"
                  content={
                    <div className="text-xs space-y-2">
                      <p>ML Readiness measures how suitable a molecule is for machine learning applications across 4 scientific dimensions.</p>
                      <ul className="list-disc list-inside space-y-1 text-white/70">
                        <li><strong>Structural Quality (20pts):</strong> Clean structure for ML pipelines</li>
                        <li><strong>Property Profile (35pts):</strong> Drug-like physicochemical properties</li>
                        <li><strong>Complexity & Feasibility (25pts):</strong> QED, synthesizability, 3D character</li>
                        <li><strong>Representation Quality (20pts):</strong> Descriptor/fingerprint completeness</li>
                      </ul>
                    </div>
                  }
                />
              </div>
              <p className="text-sm text-[var(--color-text-secondary)] leading-relaxed mb-2.5">
                {TIER_GUIDANCE[label] || interpretation}
              </p>
              {/* Dimension health tags */}
              <div className="flex flex-wrap gap-1.5">
                {dimensions.map(({ key, dim }) => {
                  const dimColor = getScoreColor(dim.score, dim.max_score);
                  const pct = dim.max_score > 0 ? Math.round((dim.score / dim.max_score) * 100) : 0;
                  return (
                    <span
                      key={key}
                      className="inline-flex items-center gap-1 text-[11px] px-2 py-0.5 rounded-md bg-[var(--color-surface-sunken)]"
                    >
                      <span
                        className="w-1.5 h-1.5 rounded-full flex-shrink-0"
                        style={{ backgroundColor: dimColor.stroke }}
                      />
                      <span className="text-[var(--color-text-muted)]">{DIM_SHORT[key]}</span>
                      <span className={cn('font-medium tabular-nums', dimColor.text)}>{pct}%</span>
                    </span>
                  );
                })}
              </div>
            </div>
          </div>
        </motion.div>

        {/* Caveats */}
        <CaveatsBanner caveats={caveats} />

        {/* Dimension cards */}
        <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
          {dimensions.map(({ key, dim }, i) => (
            <DimensionCard key={key} dimKey={key} dimension={dim} delay={0.1 + i * 0.05} />
          ))}
        </div>

        {/* Supplementary pills */}
        <SupplementaryPills supplementary={supplementary} />
      </div>
    );
  }

  // Full view with header + score chart
  return (
    <div className="card-chem p-6">
      {/* Header */}
      <div className="flex items-start justify-between mb-6">
        <div className="flex items-center gap-3">
          <div className="section-header-icon">
            <svg className="w-5 h-5" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <path d="M12 2v4M12 18v4M4.93 4.93l2.83 2.83M16.24 16.24l2.83 2.83M2 12h4M18 12h4M4.93 19.07l2.83-2.83M16.24 7.76l2.83-2.83" />
            </svg>
          </div>
          <div>
            <div className="flex items-center gap-1.5">
              <h3 className="text-lg font-semibold text-chem-dark">ML-Readiness Score</h3>
              <InfoTooltip
                title="What is ML Readiness?"
                content={
                  <div className="text-xs space-y-2">
                    <p>A scientifically meaningful 0-100 score across 4 dimensions assessing how suitable a molecule is for machine learning applications.</p>
                    <ul className="list-disc list-inside space-y-1 text-white/70">
                      <li><strong>Structural Quality:</strong> Clean structure checks</li>
                      <li><strong>Property Profile:</strong> Desirability-scored properties</li>
                      <li><strong>Complexity:</strong> QED, SA Score, Fsp3, stereo</li>
                      <li><strong>Representation:</strong> Descriptor/FP completeness</li>
                    </ul>
                  </div>
                }
              />
            </div>
            <p className="text-sm text-chem-dark/50">Scientific suitability for ML models</p>
          </div>
        </div>

        <ScoreChart
          score={score}
          label="ML-Readiness"
          size={120}
          calculation={calculation}
          interpretation={interpretation}
          compact
        />
      </div>

      {/* Interpretation */}
      <div className="bg-chem-primary/5 rounded-xl p-4 mb-6">
        <p className="text-sm text-chem-dark/80 mb-2.5">
          {TIER_GUIDANCE[label] || interpretation}
        </p>
        <div className="flex flex-wrap gap-1.5">
          {dimensions.map(({ key, dim }) => {
            const dimColor = getScoreColor(dim.score, dim.max_score);
            const pct = dim.max_score > 0 ? Math.round((dim.score / dim.max_score) * 100) : 0;
            return (
              <span
                key={key}
                className="inline-flex items-center gap-1 text-[11px] px-2 py-0.5 rounded-md bg-chem-primary/5 border border-chem-primary/10"
              >
                <span
                  className="w-1.5 h-1.5 rounded-full flex-shrink-0"
                  style={{ backgroundColor: dimColor.stroke }}
                />
                <span className="text-chem-dark/50">{DIM_SHORT[key]}</span>
                <span className={cn('font-medium tabular-nums', dimColor.text)}>{pct}%</span>
              </span>
            );
          })}
        </div>
      </div>

      {/* Caveats */}
      {caveats.length > 0 && (
        <div className="mb-6">
          <CaveatsBanner caveats={caveats} />
        </div>
      )}

      {/* Dimension breakdown */}
      <div className="space-y-4">
        <h4 className="text-sm font-semibold text-chem-dark/70 uppercase tracking-wide">
          Score Breakdown
        </h4>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
          {dimensions.map(({ key, dim }, i) => (
            <DimensionCard key={key} dimKey={key} dimension={dim} delay={i * 0.05} />
          ))}
        </div>

        {/* Supplementary info */}
        <SupplementaryPills supplementary={supplementary} />
      </div>
    </div>
  );
}
