import { AnimatePresence, motion } from 'framer-motion';
import { CheckCircle2 } from 'lucide-react';

import { AlertCard } from '../alerts/AlertCard';
import { Badge } from '../ui/Badge';
import type { AlertScreenResponse } from '../../types/alerts';

interface ValidationOutcomePanelsProps {
  /** Show the "All Clear" success panel (valid molecule, no issues). */
  showSuccess: boolean;
  successExecutionMs: number;
  /** Structural-alert screening result, or null if not yet screened. */
  alertResult: AlertScreenResponse | null;
  onAtomHover: (atoms: number[]) => void;
}

/** Right-column validation outcome panels: "All Clear" success + alert screening. */
export function ValidationOutcomePanels({
  showSuccess,
  successExecutionMs,
  alertResult,
  onAtomHover,
}: ValidationOutcomePanelsProps) {
  const alertIssues = alertResult?.alerts || [];

  return (
    <AnimatePresence>
      {showSuccess && (
        <motion.div
          key="validation-success"
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          exit={{ opacity: 0, y: -20 }}
          className="rounded-xl p-5 text-center bg-[rgba(251,191,36,0.18)] border border-[rgba(251,191,36,0.35)]"
        >
          <CheckCircle2
            className="w-10 h-10 mx-auto mb-2 text-[#d97706] dark:text-[#fbbf24]"
            strokeWidth={2.25}
          />
          <h3 className="text-lg font-semibold text-[#b45309] dark:text-[#fcd34d] mb-1 font-display">
            All Clear
          </h3>
          <p className="text-sm text-[var(--color-text-secondary)]">
            All validation checks passed
          </p>
          <p className="mt-3 text-xs text-[var(--color-text-muted)]">
            Completed in {successExecutionMs.toFixed(0)}ms
          </p>
        </motion.div>
      )}

      {alertResult && (
        <motion.div
          key="alert-results"
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          exit={{ opacity: 0, y: -20 }}
          className="card p-5 sm:p-6"
        >
          <div className="flex items-center justify-between mb-4">
            <div>
              <h4 className="font-semibold text-[var(--color-text-primary)] text-sm">
                Structural Alerts
              </h4>
              <p className="text-xs text-[var(--color-text-muted)] mt-0.5">
                Screened: {alertResult.screened_catalogs.join(', ')}
              </p>
            </div>
            <Badge variant={alertIssues.length === 0 ? 'success' : 'warning'}>
              {alertIssues.length} alerts
            </Badge>
          </div>

          {alertIssues.length > 0 ? (
            <div className="space-y-3 max-h-[400px] overflow-y-auto pr-2">
              {alertIssues.map((alert, index) => (
                <AlertCard
                  key={`${alert.pattern_name}-${index}`}
                  alert={alert}
                  onAtomHover={onAtomHover}
                />
              ))}
            </div>
          ) : (
            <div className="rounded-xl p-4 text-center bg-[rgba(251,191,36,0.18)] border border-[rgba(251,191,36,0.35)]">
              <CheckCircle2
                className="w-7 h-7 mx-auto mb-1 text-[#d97706] dark:text-[#fbbf24]"
                strokeWidth={2.25}
              />
              <p className="text-sm font-medium text-[#b45309] dark:text-[#fcd34d]">
                No structural alerts detected
              </p>
            </div>
          )}

          <p className="mt-4 text-xs text-[var(--color-text-muted)] text-right">
            Completed in {alertResult.execution_time_ms}ms
          </p>
        </motion.div>
      )}
    </AnimatePresence>
  );
}
