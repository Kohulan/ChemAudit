import { useState } from 'react';
import { StepNode } from './StepNode';
import type { QSARStepResult } from '../../types/qsar_ready';

// =============================================================================
// Props
// =============================================================================

interface PipelineStepTimelineProps {
  /** Array of 10 pipeline step results to render in order. */
  steps: QSARStepResult[];
}

// =============================================================================
// Component
// =============================================================================

/**
 * Vertical 10-node pipeline step timeline.
 *
 * Auto-expands any step where status === 'applied' AND before_smiles !== after_smiles.
 * All 'no_change' and 'skipped' nodes start collapsed.
 *
 * Connector line sits at left-[20px] matching ProvenanceTimeline geometry.
 * Each node has a colored dot on the connector: green=applied, blue=no_change,
 * gray=skipped, red=error.
 *
 * Per Phase 10 UI-SPEC step timeline layout.
 */
export function PipelineStepTimeline({ steps }: PipelineStepTimelineProps) {
  // Initialize expanded state: auto-expand applied steps where SMILES changed
  const [expandedSteps, setExpandedSteps] = useState<Set<number>>(() => {
    const initial = new Set<number>();
    steps.forEach((step, idx) => {
      if (
        step.status === 'applied' &&
        step.before_smiles !== null &&
        step.after_smiles !== null &&
        step.before_smiles !== step.after_smiles
      ) {
        initial.add(idx);
      }
    });
    return initial;
  });

  const toggleStep = (idx: number) => {
    setExpandedSteps((prev) => {
      const next = new Set(prev);
      if (next.has(idx)) {
        next.delete(idx);
      } else {
        next.add(idx);
      }
      return next;
    });
  };

  if (steps.length === 0) {
    return null;
  }

  return (
    <div className="relative">
      {/* Vertical connector line at left-[20px] matching ProvenanceTimeline pattern */}
      <div
        className="absolute left-[20px] top-5 bottom-5 w-0.5 bg-[var(--color-border)] z-0"
        aria-hidden="true"
      />

      {/* Step nodes */}
      <div className="space-y-2 relative z-10">
        {steps.map((step, idx) => (
          <StepNode
            key={`${step.step_name}-${step.step_index}`}
            step={step}
            expanded={expandedSteps.has(idx)}
            onToggle={() => toggleStep(idx)}
          />
        ))}
      </div>
    </div>
  );
}
