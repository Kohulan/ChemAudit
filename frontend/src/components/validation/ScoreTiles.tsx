import { AnimatePresence, motion } from 'framer-motion';
import { Target, Zap } from 'lucide-react';

import { ScoreChart } from '../scoring/ScoreChart';
import { Badge } from '../ui/Badge';
import { getScoreLabel } from '../../lib/utils';

function IssueSeverityTags({
  issues,
  totalChecks,
}: {
  issues: { severity: string }[];
  totalChecks: number;
}) {
  const counts = { critical: 0, error: 0, warning: 0 };
  for (const issue of issues) {
    const sev = issue.severity as keyof typeof counts;
    if (sev in counts) counts[sev]++;
  }
  const passed = totalChecks - issues.length;

  const TAG_STYLES: Record<string, string> = {
    critical: 'bg-red-500/10 text-red-600 dark:text-red-400',
    error: 'bg-orange-500/10 text-orange-600 dark:text-orange-400',
    warning: 'bg-amber-500/10 text-amber-600 dark:text-amber-400',
  };

  return (
    <div className="flex flex-wrap items-center gap-1.5 mt-2">
      {totalChecks > 0 && (
        <span className="text-[10px] px-1.5 py-0.5 rounded bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)]">
          {passed}/{totalChecks} passed
        </span>
      )}
      {Object.entries(counts).map(([severity, count]) =>
        count > 0 ? (
          <span
            key={severity}
            className={`text-[10px] px-1.5 py-0.5 rounded font-medium ${TAG_STYLES[severity]}`}
          >
            {count} {severity}
          </span>
        ) : null
      )}
      {issues.length === 0 && totalChecks > 0 && (
        <span className="text-[10px] px-1.5 py-0.5 rounded bg-[rgba(251,191,36,0.18)] text-[#b45309] dark:text-[#fcd34d] font-medium">
          All clear
        </span>
      )}
    </div>
  );
}

interface ScoreTilesProps {
  show: boolean;
  qualityScore: number | null;
  mlReadyScore: number | undefined;
  mlReady: boolean | null;
  /** Validation issues for the severity tags, or null when there is no result yet. */
  issues: { severity: string }[] | null;
  totalChecks: number;
}

/** Quality + ML-readiness score tiles shown after validation/scoring. */
export function ScoreTiles({
  show,
  qualityScore,
  mlReadyScore,
  mlReady,
  issues,
  totalChecks,
}: ScoreTilesProps) {
  return (
    <AnimatePresence>
      {show && (
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          exit={{ opacity: 0, y: -20 }}
          className="grid grid-cols-2 gap-4"
        >
          {/* Quality Score */}
          <div className="card-gradient p-4 sm:p-5">
            <div className="flex items-center justify-between mb-2">
              <div className="flex items-center gap-2">
                <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-[var(--color-primary)]/10 to-[var(--color-accent)]/10 flex items-center justify-center text-[var(--color-primary)]">
                  <Target className="w-4 h-4" />
                </div>
                <span className="text-xs font-medium text-[var(--color-text-secondary)]">Quality</span>
              </div>
              {qualityScore !== null && (
                <Badge
                  variant={qualityScore >= 70 ? 'success' : qualityScore >= 40 ? 'warning' : 'error'}
                  size="sm"
                >
                  {getScoreLabel(qualityScore)}
                </Badge>
              )}
            </div>
            <div className="flex justify-center">
              {qualityScore !== null ? (
                <ScoreChart
                  score={qualityScore}
                  label="Validation Quality"
                  size={100}
                  compact
                  variant="cool"
                />
              ) : (
                <div className="w-[100px] h-[100px] flex items-center justify-center">
                  <span className="text-3xl font-bold text-[var(--color-text-muted)]">--</span>
                </div>
              )}
            </div>
            {/* Issue severity tags */}
            {issues !== null && <IssueSeverityTags issues={issues} totalChecks={totalChecks} />}
          </div>

          {/* ML Readiness */}
          <div className="card-accent p-4 sm:p-5 flex flex-col">
            <div className="flex items-center justify-between mb-2">
              <div className="flex items-center gap-2">
                <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-[var(--color-accent)]/10 to-[var(--color-primary)]/10 flex items-center justify-center text-[var(--color-accent)]">
                  <Zap className="w-4 h-4" />
                </div>
                <span className="text-xs font-medium text-[var(--color-text-secondary)]">ML Ready</span>
              </div>
              {mlReady !== null && (
                <Badge variant={mlReady ? 'success' : 'warning'} dot size="sm">
                  {mlReady ? 'Ready' : 'Review'}
                </Badge>
              )}
            </div>
            <div className="flex-1 flex items-center justify-center">
              {mlReadyScore !== undefined ? (
                <ScoreChart score={mlReadyScore} label="ML-Readiness" size={100} compact />
              ) : (
                <div className="w-[100px] h-[100px] flex items-center justify-center">
                  <span className="text-3xl font-bold text-[var(--color-text-muted)]">--</span>
                </div>
              )}
            </div>
          </div>
        </motion.div>
      )}
    </AnimatePresence>
  );
}
