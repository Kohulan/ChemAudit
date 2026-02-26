/**
 * BatchTimeline (VIZ-09)
 *
 * Batch processing timeline showing Upload -> Validation -> Analytics -> Complete.
 * Horizontal on lg screens, vertical on mobile.
 */

import React from 'react';
import { motion } from 'framer-motion';
import { Upload, Shield, BarChart3, CheckCircle } from 'lucide-react';
import type { BatchStatistics } from '../../types/batch';
import type { AnalysisStatus } from '../../types/analytics';
import { cn } from '../../lib/utils';

interface BatchTimelineProps {
  statistics: BatchStatistics;
  analyticsStatus: Record<string, AnalysisStatus> | null;
}

type PhaseStatus = 'complete' | 'computing' | 'pending' | 'failed';

interface TimelinePhase {
  label: string;
  icon: React.ReactNode;
  status: PhaseStatus;
  detail?: string;
}

function getStatusColor(status: PhaseStatus): string {
  switch (status) {
    case 'complete':
      return 'bg-emerald-500 text-white';
    case 'computing':
      return 'bg-amber-500 text-white animate-pulse';
    case 'failed':
      return 'bg-red-500 text-white';
    case 'pending':
    default:
      return 'bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)] border border-[var(--color-border)]';
  }
}

function getLineColor(status: PhaseStatus): string {
  switch (status) {
    case 'complete':
      return 'bg-emerald-500';
    case 'computing':
      return 'bg-amber-500';
    case 'failed':
      return 'bg-red-500';
    case 'pending':
    default:
      return 'bg-[var(--color-border)]';
  }
}

/**
 * Derive the analytics phase status from the response data.
 * Only considers "active" types — those that have moved beyond "pending"
 * (i.e. auto-started or user-triggered). Types that are still "pending"
 * (never started, e.g. mmp, similarity_search, rgroup) are ignored.
 */
function deriveAnalyticsPhaseStatus(
  analyticsStatus: Record<string, AnalysisStatus> | null
): PhaseStatus {
  if (!analyticsStatus) return 'pending';

  // Active types = started (non-pending)
  const activeStatuses = Object.values(analyticsStatus)
    .filter((s) => s.status !== 'pending');

  if (activeStatuses.length === 0) return 'pending';

  const anyComputing = activeStatuses.some(
    (s) => s.status === 'computing'
  );
  const anyFailed = activeStatuses.some((s) => s.status === 'failed');
  const allTerminal = activeStatuses.every(
    (s) => s.status === 'complete' || s.status === 'skipped' || s.status === 'failed'
  );

  if (allTerminal && !anyFailed) return 'complete';
  if (allTerminal && anyFailed) return 'failed';
  if (anyComputing) return 'computing';
  return 'pending';
}

/**
 * Count progress across active (non-pending) analytics types.
 * Derives the total dynamically from the response rather than a hardcoded list.
 */
function computeAnalyticsProgress(
  analyticsStatus: Record<string, AnalysisStatus> | null
): { completed: number; total: number; percent: number } {
  if (!analyticsStatus) return { completed: 0, total: 0, percent: 0 };

  let completed = 0;
  let total = 0;

  for (const s of Object.values(analyticsStatus)) {
    if (s.status === 'pending') continue;
    total++;
    if (s.status === 'complete' || s.status === 'skipped') {
      completed++;
    }
  }

  const percent = total > 0 ? Math.round((completed / total) * 100) : 0;
  return { completed, total, percent };
}

export const BatchTimeline = React.memo(function BatchTimeline({
  statistics,
  analyticsStatus,
}: BatchTimelineProps) {
  const analyticsPhase = deriveAnalyticsPhaseStatus(analyticsStatus);
  const progress = computeAnalyticsProgress(analyticsStatus);
  const allComplete = analyticsPhase === 'complete';

  const phases: TimelinePhase[] = [
    {
      label: 'Upload',
      icon: <Upload className="w-4 h-4" />,
      status: 'complete',
      detail: `${statistics.total} molecules`,
    },
    {
      label: 'Validation',
      icon: <Shield className="w-4 h-4" />,
      status: 'complete',
      detail: statistics.processing_time_seconds
        ? `${statistics.processing_time_seconds.toFixed(1)}s`
        : undefined,
    },
    {
      label: 'Analytics',
      icon: <BarChart3 className="w-4 h-4" />,
      status: analyticsPhase,
      detail:
        analyticsPhase === 'complete'
          ? 'Done'
          : analyticsPhase === 'computing'
            ? `${progress.completed}/${progress.total} complete`
            : analyticsPhase === 'failed'
              ? 'Partial failure'
              : undefined,
    },
    {
      label: 'Complete',
      icon: <CheckCircle className="w-4 h-4" />,
      status: allComplete ? 'complete' : 'pending',
    },
  ];

  // Summary stats
  const scaffoldCount = analyticsStatus?.scaffold?.status === 'complete' ? 'scaffolds analyzed' : '';
  const outlierInfo = analyticsStatus?.statistics?.status === 'complete' ? 'outliers detected' : '';

  return (
    <div className="space-y-3">
      {/* Timeline */}
      <div className="flex flex-col lg:flex-row items-start lg:items-center gap-0">
        {phases.map((phase, i) => (
          <React.Fragment key={phase.label}>
            <motion.div
              className="flex flex-row lg:flex-col items-center gap-2 lg:gap-1"
              initial={{ opacity: 0, y: 10 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ delay: i * 0.1 }}
            >
              {/* Node */}
              <div
                className={cn(
                  'w-8 h-8 rounded-full flex items-center justify-center flex-shrink-0',
                  getStatusColor(phase.status)
                )}
              >
                {phase.icon}
              </div>
              {/* Label */}
              <div className="text-center lg:min-w-[80px]">
                <p className="text-xs font-medium text-[var(--color-text-primary)]">
                  {phase.label}
                </p>
                {phase.detail && (
                  <p className="text-[10px] text-[var(--color-text-muted)]">{phase.detail}</p>
                )}
              </div>
            </motion.div>

            {/* Connecting line */}
            {i < phases.length - 1 && (
              <div
                className={cn(
                  'hidden lg:block flex-1 h-1 min-w-[40px] rounded-full mx-1',
                  getLineColor(phases[i + 1].status === 'pending' ? 'pending' : phase.status)
                )}
              />
            )}
            {i < phases.length - 1 && (
              <div
                className={cn(
                  'block lg:hidden w-1 h-4 rounded-full ml-[14px]',
                  getLineColor(phases[i + 1].status === 'pending' ? 'pending' : phase.status)
                )}
              />
            )}
          </React.Fragment>
        ))}
      </div>

      {/* Analytics progress bar */}
      {analyticsPhase === 'computing' && (
        <div className="space-y-1">
          <div className="flex items-center justify-between text-[10px] text-[var(--color-text-muted)]">
            <span>Analytics: {progress.completed}/{progress.total} complete</span>
            <span>{progress.percent}%</span>
          </div>
          <div className="h-1.5 w-full rounded-full bg-[var(--color-surface-sunken)] overflow-hidden">
            <motion.div
              className="h-full rounded-full bg-amber-500"
              initial={{ width: 0 }}
              animate={{ width: `${progress.percent}%` }}
              transition={{ duration: 0.5, ease: 'easeOut' }}
            />
          </div>
        </div>
      )}

      {/* Summary stats */}
      <div className="flex flex-wrap gap-3 text-[10px] text-[var(--color-text-muted)]">
        <span>
          {statistics.successful} of {statistics.total} validated successfully
        </span>
        {statistics.processing_time_seconds && (
          <span>in {statistics.processing_time_seconds.toFixed(1)}s</span>
        )}
        {scaffoldCount && <span>{scaffoldCount}</span>}
        {outlierInfo && <span>{outlierInfo}</span>}
      </div>
    </div>
  );
});
