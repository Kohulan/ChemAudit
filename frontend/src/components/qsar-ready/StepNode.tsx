import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown } from 'lucide-react';
import { cn } from '../../lib/utils';
import { Badge } from '../ui/Badge';
import type { QSARStepResult } from '../../types/qsar_ready';

// =============================================================================
// Props
// =============================================================================

interface StepNodeProps {
  /** The pipeline step result to render. */
  step: QSARStepResult;
  /** Whether the collapsible detail section is open. */
  expanded: boolean;
  /** Called when the user clicks the header to toggle expansion. */
  onToggle: () => void;
}

// =============================================================================
// Status badge helper
// =============================================================================

function StatusBadge({ status }: { status: QSARStepResult['status'] }) {
  switch (status) {
    case 'applied':
      return <Badge variant="success" size="sm">Applied</Badge>;
    case 'no_change':
      return (
        <span className="inline-flex items-center px-2 py-0.5 rounded-full text-xs font-medium bg-blue-500/10 text-blue-700 dark:text-blue-400">
          No Change
        </span>
      );
    case 'skipped':
      return (
        <span className="inline-flex items-center px-2 py-0.5 rounded-full text-xs font-medium bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)]">
          Skipped
        </span>
      );
    case 'error':
      return <Badge variant="error" size="sm">Error</Badge>;
  }
}

// =============================================================================
// Timeline dot color
// =============================================================================

function dotColorClass(status: QSARStepResult['status']): string {
  switch (status) {
    case 'applied':  return 'bg-green-500 border-green-600';
    case 'no_change': return 'bg-blue-500 border-blue-600';
    case 'skipped':  return 'bg-[var(--color-border)] border-[var(--color-border)]';
    case 'error':    return 'bg-red-500 border-red-600';
  }
}

// =============================================================================
// SMILES diff highlight — bold changed tokens
// =============================================================================

function SmilesDiff({ before, after }: { before: string; after: string }) {
  if (before === after) {
    return <span className="font-mono text-xs text-[var(--color-text-muted)]">{before}</span>;
  }

  // Simple character-level diff: highlight chars that differ
  const maxLen = Math.max(before.length, after.length);
  const changedIndices = new Set<number>();
  for (let i = 0; i < maxLen; i++) {
    if (before[i] !== after[i]) changedIndices.add(i);
  }

  const renderWithHighlight = (str: string) => {
    const segments: Array<{ text: string; changed: boolean }> = [];
    let current = '';
    let lastChanged = false;

    for (let i = 0; i < str.length; i++) {
      const changed = changedIndices.has(i);
      if (changed !== lastChanged && current) {
        segments.push({ text: current, changed: lastChanged });
        current = '';
      }
      current += str[i];
      lastChanged = changed;
    }
    if (current) segments.push({ text: current, changed: lastChanged });

    return segments.map((seg, i) =>
      seg.changed ? (
        <strong key={i} className="text-[var(--color-primary)]">{seg.text}</strong>
      ) : (
        <span key={i}>{seg.text}</span>
      ),
    );
  };

  return (
    <div className="space-y-1">
      <div className="font-mono text-xs text-[var(--color-text-muted)]">
        <span className="text-[var(--color-text-muted)] mr-1.5">before:</span>
        {renderWithHighlight(before)}
      </div>
      <div className="font-mono text-xs text-[var(--color-text-muted)]">
        <span className="text-[var(--color-text-muted)] mr-1.5">after: </span>
        {renderWithHighlight(after)}
      </div>
    </div>
  );
}

// =============================================================================
// Component
// =============================================================================

/**
 * Single node in the QSAR pipeline step timeline.
 *
 * Header (clickable): step number badge + step name + status badge + chevron.
 * Collapsible body: before/after SMILES with diff highlighting + detail text.
 *
 * Expand animation: AnimatePresence + motion.div height 0→auto, 300ms ease-out.
 * Accessibility: aria-expanded on the header button.
 *
 * Per Phase 10 UI-SPEC interaction contracts.
 */
export function StepNode({ step, expanded, onToggle }: StepNodeProps) {
  const hasDetail = !!(step.before_smiles || step.after_smiles || step.detail);
  const canExpand = hasDetail;

  return (
    <div className="flex items-start gap-3">
      {/* Timeline dot */}
      <div
        className={cn(
          'flex-shrink-0 w-[18px] h-[18px] mt-[13px] rounded-full border-2 z-10',
          dotColorClass(step.status),
        )}
        aria-hidden="true"
      />

      {/* Step card */}
      <div className="flex-1 min-w-0 rounded-xl border bg-[var(--color-surface-elevated)] border-[var(--color-border)] overflow-hidden">
        {/* Clickable header */}
        <button
          type="button"
          onClick={canExpand ? onToggle : undefined}
          aria-expanded={canExpand ? expanded : undefined}
          className={cn(
            'w-full flex items-center gap-2 px-3 py-2.5 text-left',
            canExpand
              ? 'cursor-pointer hover:bg-[var(--color-surface-sunken)] transition-colors'
              : 'cursor-default',
          )}
        >
          {/* Step number badge */}
          <span className="inline-flex items-center justify-center w-5 h-5 rounded-full text-[10px] font-semibold bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)] shrink-0">
            {step.step_index}
          </span>

          {/* Step name */}
          <span className="text-xs font-semibold text-[var(--color-text-primary)] flex-1 min-w-0 truncate">
            {step.step_name}
          </span>

          {/* Status badge */}
          <StatusBadge status={step.status} />

          {/* Expand chevron */}
          {canExpand && (
            <ChevronDown
              className={cn(
                'w-3.5 h-3.5 text-[var(--color-text-muted)] transition-transform duration-200 shrink-0',
                expanded && 'rotate-180',
              )}
            />
          )}
        </button>

        {/* Collapsible detail body */}
        <AnimatePresence initial={false}>
          {expanded && canExpand && (
            <motion.div
              key="detail"
              initial={{ height: 0, opacity: 0 }}
              animate={{ height: 'auto', opacity: 1 }}
              exit={{ height: 0, opacity: 0 }}
              transition={{ duration: 0.3, ease: 'easeOut' }}
              className="overflow-hidden"
            >
              <div className="px-3 pb-3 space-y-2 border-t border-[var(--color-border)]">
                {/* SMILES diff */}
                {step.before_smiles && step.after_smiles && (
                  <div className="pt-2">
                    <SmilesDiff
                      before={step.before_smiles}
                      after={step.after_smiles}
                    />
                  </div>
                )}

                {/* Only before (no after) */}
                {step.before_smiles && !step.after_smiles && (
                  <p className="pt-2 font-mono text-xs text-[var(--color-text-muted)]">
                    <span className="mr-1.5">before:</span>
                    {step.before_smiles}
                  </p>
                )}

                {/* Detail text */}
                {step.detail && (
                  <p className="text-xs text-[var(--color-text-secondary)]">
                    {step.detail}
                  </p>
                )}
              </div>
            </motion.div>
          )}
        </AnimatePresence>
      </div>
    </div>
  );
}
