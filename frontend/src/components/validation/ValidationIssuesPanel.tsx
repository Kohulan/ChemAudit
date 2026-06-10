import { type ComponentProps } from 'react';
import { AnimatePresence, motion } from 'framer-motion';

import { Badge } from '../ui/Badge';
import { IssueCard } from './IssueCard';

type Issue = ComponentProps<typeof IssueCard>['issue'];

interface ValidationIssuesPanelProps {
  show: boolean;
  issues: Issue[];
  executionMs: number;
  highlightedAtoms: number[];
  highlightLocked: boolean;
  onAtomHover: (atoms: number[]) => void;
  onAtomLock: ComponentProps<typeof IssueCard>['onAtomLock'];
}

/** Right-column list of validation issues (validate tab) with atom highlight/lock. */
export function ValidationIssuesPanel({
  show,
  issues,
  executionMs,
  highlightedAtoms,
  highlightLocked,
  onAtomHover,
  onAtomLock,
}: ValidationIssuesPanelProps) {
  return (
    <AnimatePresence>
      {show && (
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          exit={{ opacity: 0, y: -20 }}
          className="card p-5 sm:p-6"
        >
          <div className="flex items-center justify-between mb-4">
            <h4 className="font-semibold text-[var(--color-text-primary)] text-sm">
              Validation Issues
            </h4>
            <Badge variant="warning">{issues.length} found</Badge>
          </div>
          <div className="space-y-3 max-h-[400px] overflow-y-auto pr-2">
            {issues.map((issue, index) => {
              const isThisLocked =
                highlightLocked &&
                JSON.stringify(highlightedAtoms) === JSON.stringify(issue.affected_atoms);
              return (
                <IssueCard
                  key={`${issue.check_name}-${index}`}
                  issue={issue}
                  onAtomHover={highlightLocked ? undefined : onAtomHover}
                  onAtomLock={onAtomLock}
                  isLocked={isThisLocked}
                />
              );
            })}
          </div>
          <p className="mt-4 text-xs text-[var(--color-text-muted)] text-right">
            Completed in {executionMs.toFixed(0)}ms
          </p>
        </motion.div>
      )}
    </AnimatePresence>
  );
}
