import { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { CheckCircle2, XCircle, ChevronDown } from 'lucide-react';
import type { ConsensusScore } from '../../types/scoring';
import { cn } from '../../lib/utils';

interface ConsensusScoreCardProps {
  data: ConsensusScore;
}

function getScoreColor(score: number, total: number) {
  const ratio = score / total;
  if (ratio >= 0.8) return 'text-emerald-500';
  if (ratio >= 0.6) return 'text-teal-500';
  if (ratio >= 0.4) return 'text-amber-500';
  return 'text-red-500';
}

function getScoreBg(score: number, total: number) {
  const ratio = score / total;
  if (ratio >= 0.8) return 'bg-emerald-500/10 border-emerald-500/20';
  if (ratio >= 0.6) return 'bg-teal-500/10 border-teal-500/20';
  if (ratio >= 0.4) return 'bg-amber-500/10 border-amber-500/20';
  return 'bg-red-500/10 border-red-500/20';
}

export function ConsensusScoreCard({ data }: ConsensusScoreCardProps) {
  const [expandedSet, setExpandedSet] = useState<string | null>(null);

  return (
    <div className="space-y-4">
      {/* Score Display */}
      <div className={cn(
        'flex items-center justify-between p-4 rounded-xl border',
        getScoreBg(data.score, data.total)
      )}>
        <div>
          <p className="text-xs font-medium text-[var(--color-text-muted)] uppercase tracking-wider">
            Consensus Drug-Likeness
          </p>
          <p className="text-sm text-[var(--color-text-secondary)] mt-1">
            {data.interpretation}
          </p>
        </div>
        <div className="flex items-baseline gap-1">
          <span className={cn('text-4xl font-bold tracking-tight', getScoreColor(data.score, data.total))}>
            {data.score}
          </span>
          <span className="text-lg text-[var(--color-text-muted)]">/{data.total}</span>
        </div>
      </div>

      {/* Rule Sets */}
      <div className="space-y-2">
        {data.rule_sets.map((rs) => (
          <div key={rs.name} className="rounded-xl border border-[var(--color-border)] overflow-hidden">
            <button
              onClick={() => setExpandedSet(expandedSet === rs.name ? null : rs.name)}
              className="w-full flex items-center justify-between p-3 hover:bg-[var(--color-surface-sunken)]/50 transition-colors"
            >
              <div className="flex items-center gap-2">
                {rs.passed ? (
                  <CheckCircle2 className="w-4 h-4 text-emerald-500" />
                ) : (
                  <XCircle className="w-4 h-4 text-red-500" />
                )}
                <span className="text-sm font-medium text-[var(--color-text-primary)]">
                  {rs.name}
                </span>
              </div>
              <div className="flex items-center gap-2">
                <span className={cn(
                  'text-xs px-2 py-0.5 rounded-full font-medium',
                  rs.passed
                    ? 'bg-emerald-500/10 text-emerald-600 dark:text-emerald-400'
                    : 'bg-red-500/10 text-red-600 dark:text-red-400'
                )}>
                  {rs.passed ? 'PASS' : 'FAIL'}
                </span>
                <ChevronDown className={cn(
                  'w-4 h-4 text-[var(--color-text-muted)] transition-transform',
                  expandedSet === rs.name && 'rotate-180'
                )} />
              </div>
            </button>

            <AnimatePresence>
              {expandedSet === rs.name && rs.violations.length > 0 && (
                <motion.div
                  initial={{ height: 0, opacity: 0 }}
                  animate={{ height: 'auto', opacity: 1 }}
                  exit={{ height: 0, opacity: 0 }}
                  transition={{ duration: 0.2 }}
                  className="overflow-hidden"
                >
                  <div className="px-3 pb-3 space-y-1">
                    {rs.violations.map((v) => (
                      <div
                        key={v.property}
                        className="flex items-center justify-between py-1.5 px-2 rounded-lg bg-[var(--color-surface-sunken)]"
                      >
                        <span className="text-xs text-[var(--color-text-muted)]">{v.property}</span>
                        <div className="flex items-center gap-2">
                          <span className={cn(
                            'text-xs font-medium',
                            v.result === 'pass'
                              ? 'text-emerald-600 dark:text-emerald-400'
                              : 'text-red-600 dark:text-red-400'
                          )}>
                            {typeof v.value === 'number' ? (Number.isInteger(v.value) ? String(v.value) : v.value.toFixed(2)) : v.value}
                          </span>
                          <span className="text-xs text-[var(--color-text-muted)]">
                            ({v.threshold})
                          </span>
                          {v.result === 'pass' ? (
                            <CheckCircle2 className="w-3 h-3 text-emerald-500" />
                          ) : (
                            <XCircle className="w-3 h-3 text-red-500" />
                          )}
                        </div>
                      </div>
                    ))}
                  </div>
                </motion.div>
              )}
            </AnimatePresence>
          </div>
        ))}
      </div>
    </div>
  );
}
