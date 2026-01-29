import { useState } from 'react';
import { motion } from 'framer-motion';
import { Cpu, Fingerprint, Scale, ChevronRight, Sparkles, AlertCircle } from 'lucide-react';
import type { MLReadinessResult } from '../../types/scoring';
import { ScoreChart } from './ScoreChart';
import { InfoTooltip } from '../ui/Tooltip';
import { cn } from '../../lib/utils';

interface MLReadinessScoreProps {
  result: MLReadinessResult;
  /** When true, hides the header/score chart and only shows breakdown */
  breakdownOnly?: boolean;
}

/**
 * Get color config based on percentage
 */
function getScoreColor(percentage: number) {
  if (percentage >= 80) return {
    gradient: 'from-yellow-500 to-amber-400',
    bg: 'bg-yellow-500/10',
    text: 'text-amber-500 dark:text-yellow-400',
    border: 'border-yellow-500/20',
    glow: 'shadow-yellow-500/20',
  };
  if (percentage >= 50) return {
    gradient: 'from-orange-500 to-orange-400',
    bg: 'bg-orange-500/10',
    text: 'text-orange-500',
    border: 'border-orange-500/20',
    glow: 'shadow-orange-500/20',
  };
  return {
    gradient: 'from-red-500 to-red-400',
    bg: 'bg-red-500/10',
    text: 'text-red-500',
    border: 'border-red-500/20',
    glow: 'shadow-red-500/20',
  };
}

interface BreakdownCardProps {
  icon: React.ReactNode;
  label: string;
  score: number;
  maxScore: number;
  detail: string;
  subDetail?: string;
  delay?: number;
}

function BreakdownCard({ icon, label, score, maxScore, detail, subDetail, delay = 0 }: BreakdownCardProps) {
  const percentage = maxScore > 0 ? (score / maxScore) * 100 : 0;
  const color = getScoreColor(percentage);

  return (
    <motion.div
      initial={{ opacity: 0, y: 20 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.4, delay }}
      className={cn(
        'relative overflow-hidden rounded-2xl p-4',
        'bg-gradient-to-br from-[var(--color-surface-elevated)] to-[var(--color-surface)]',
        'border border-[var(--color-border)]',
        'hover:border-[var(--color-primary)]/30 hover:shadow-lg hover:shadow-[var(--color-primary)]/5',
        'transition-all duration-300'
      )}
    >
      {/* Background glow effect */}
      <div
        className={cn(
          'absolute -top-12 -right-12 w-32 h-32 rounded-full blur-3xl opacity-20',
          color.bg
        )}
      />

      <div className="relative">
        {/* Header */}
        <div className="flex items-center justify-between mb-3">
          <div className="flex items-center gap-3">
            <div className={cn(
              'w-10 h-10 rounded-xl flex items-center justify-center',
              color.bg, color.text
            )}>
              {icon}
            </div>
            <div>
              <h4 className="font-semibold text-[var(--color-text-primary)] text-sm">{label}</h4>
              <p className="text-xs text-[var(--color-text-muted)]">{detail}</p>
            </div>
          </div>
          <div className="text-right">
            <div className={cn('text-2xl font-bold', color.text)}>
              {score.toFixed(0)}
            </div>
            <div className="text-xs text-[var(--color-text-muted)]">/ {maxScore}</div>
          </div>
        </div>

        {/* Progress bar */}
        <div className="h-2 bg-[var(--color-surface-sunken)] rounded-full overflow-hidden">
          <motion.div
            initial={{ width: 0 }}
            animate={{ width: `${percentage}%` }}
            transition={{ duration: 0.8, delay: delay + 0.2, ease: 'easeOut' }}
            className={cn('h-full rounded-full bg-gradient-to-r', color.gradient)}
          />
        </div>

        {/* Sub detail */}
        {subDetail && (
          <p className="mt-2 text-xs text-[var(--color-text-muted)]">{subDetail}</p>
        )}
      </div>
    </motion.div>
  );
}

/**
 * Displays ML-readiness score with radial chart, breakdown bars, and informative tooltips.
 */
export function MLReadinessScore({ result, breakdownOnly = false }: MLReadinessScoreProps) {
  const [showFailedDescriptors, setShowFailedDescriptors] = useState(false);
  const { score, breakdown, interpretation, failed_descriptors } = result;

  // Calculation explanation
  const calculation = `Score = Descriptors (${breakdown.descriptors_max}pts) + Fingerprints (${breakdown.fingerprints_max}pts) + Size (${breakdown.size_max}pts)
= ${breakdown.descriptors_score.toFixed(0)} + ${breakdown.fingerprints_score.toFixed(0)} + ${breakdown.size_score.toFixed(0)} = ${score}`;

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
            <div>
              <h4 className="text-sm font-semibold text-[var(--color-text-primary)] mb-1">ML Readiness Analysis</h4>
              <p className="text-sm text-[var(--color-text-secondary)] leading-relaxed">{interpretation}</p>
            </div>
          </div>
        </motion.div>

        {/* Breakdown cards grid */}
        <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
          <BreakdownCard
            icon={<Cpu className="w-5 h-5" />}
            label="Descriptors"
            score={breakdown.descriptors_score}
            maxScore={breakdown.descriptors_max}
            detail={`${breakdown.descriptors_successful}/${breakdown.descriptors_total} calculated`}
            delay={0.1}
          />
          <BreakdownCard
            icon={<Fingerprint className="w-5 h-5" />}
            label="Fingerprints"
            score={breakdown.fingerprints_score}
            maxScore={breakdown.fingerprints_max}
            detail={breakdown.fingerprints_successful.join(', ') || 'none'}
            subDetail={breakdown.fingerprints_failed.length > 0 ? `Failed: ${breakdown.fingerprints_failed.join(', ')}` : undefined}
            delay={0.2}
          />
          <BreakdownCard
            icon={<Scale className="w-5 h-5" />}
            label="Size"
            score={breakdown.size_score}
            maxScore={breakdown.size_max}
            detail={breakdown.size_category}
            subDetail={breakdown.molecular_weight !== null ? `MW: ${breakdown.molecular_weight.toFixed(1)} Da` : undefined}
            delay={0.3}
          />
        </div>

        {/* Failed descriptors (collapsible) */}
        {failed_descriptors.length > 0 && (
          <motion.div
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            transition={{ delay: 0.4 }}
          >
            <button
              onClick={() => setShowFailedDescriptors(!showFailedDescriptors)}
              className={cn(
                'w-full flex items-center justify-between gap-2 px-4 py-3 rounded-xl',
                'bg-amber-500/5 border border-amber-500/10',
                'text-sm text-amber-600 dark:text-amber-400',
                'hover:bg-amber-500/10 transition-colors'
              )}
            >
              <div className="flex items-center gap-2">
                <AlertCircle className="w-4 h-4" />
                <span>{failed_descriptors.length} descriptors could not be calculated</span>
              </div>
              <ChevronRight className={cn(
                'w-4 h-4 transition-transform',
                showFailedDescriptors && 'rotate-90'
              )} />
            </button>
            {showFailedDescriptors && (
              <motion.div
                initial={{ opacity: 0, height: 0 }}
                animate={{ opacity: 1, height: 'auto' }}
                className="mt-2 p-3 bg-[var(--color-surface-sunken)] rounded-xl text-xs text-[var(--color-text-muted)] font-mono max-h-32 overflow-y-auto"
              >
                {failed_descriptors.join(', ')}
              </motion.div>
            )}
          </motion.div>
        )}
      </div>
    );
  }

  // Full view with header
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
            <h3 className="text-lg font-semibold text-chem-dark">ML-Readiness Score</h3>
            <p className="text-sm text-chem-dark/50">Suitability for machine learning models</p>
          </div>
        </div>

        {/* Main Score */}
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
        <p className="text-sm text-chem-dark/80">{interpretation}</p>
      </div>

      {/* Score Breakdown */}
      <div className="space-y-4">
        <h4 className="text-sm font-semibold text-chem-dark/70 uppercase tracking-wide flex items-center gap-2">
          Score Breakdown
          <InfoTooltip
            title="How ML-Readiness is Calculated"
            content={
              <div className="space-y-2 text-xs">
                <p>The ML-Readiness score measures how suitable a molecule is for machine learning models.</p>
                <ul className="list-disc list-inside space-y-1 text-white/70">
                  <li>Descriptors (40pts): % of RDKit descriptors calculable</li>
                  <li>Fingerprints (40pts): Ability to generate Morgan, MACCS, AtomPair</li>
                  <li>Size (20pts): Molecular weight and atom count constraints</li>
                </ul>
              </div>
            }
          />
        </h4>

        <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
          <BreakdownCard
            icon={<Cpu className="w-5 h-5" />}
            label="Descriptors"
            score={breakdown.descriptors_score}
            maxScore={breakdown.descriptors_max}
            detail={`${breakdown.descriptors_successful}/${breakdown.descriptors_total}`}
            delay={0}
          />
          <BreakdownCard
            icon={<Fingerprint className="w-5 h-5" />}
            label="Fingerprints"
            score={breakdown.fingerprints_score}
            maxScore={breakdown.fingerprints_max}
            detail={breakdown.fingerprints_successful.join(', ') || 'none'}
            delay={0.1}
          />
          <BreakdownCard
            icon={<Scale className="w-5 h-5" />}
            label="Size"
            score={breakdown.size_score}
            maxScore={breakdown.size_max}
            detail={breakdown.size_category}
            delay={0.2}
          />
        </div>
      </div>

      {/* Failed Descriptors */}
      {failed_descriptors.length > 0 && (
        <div className="mt-6 pt-4 border-t border-chem-dark/10">
          <button
            onClick={() => setShowFailedDescriptors(!showFailedDescriptors)}
            className="flex items-center gap-2 text-sm text-chem-dark/60 hover:text-chem-dark transition-colors"
          >
            <ChevronRight className={cn('w-4 h-4 transition-transform', showFailedDescriptors && 'rotate-90')} />
            <span>{failed_descriptors.length} descriptors could not be calculated</span>
          </button>
          {showFailedDescriptors && (
            <div className="mt-3 p-4 bg-chem-dark/5 rounded-lg text-xs text-chem-dark/60 font-mono max-h-32 overflow-y-auto">
              {failed_descriptors.join(', ')}
            </div>
          )}
        </div>
      )}
    </div>
  );
}
