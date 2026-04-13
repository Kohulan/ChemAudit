import { useState } from 'react';
import { ChevronDown } from 'lucide-react';
import { cn } from '../../lib/utils';
import type { QSARReadyConfig } from '../../types/qsar_ready';
import { PIPELINE_STEPS } from '../../types/qsar_ready';

// =============================================================================
// Non-togglable steps (always-on)
// =============================================================================

interface FixedStep {
  index: number;
  name: string;
  description: string;
  alwaysOn: true;
}

const FIXED_STEPS: FixedStep[] = [
  { index: 1, name: 'Parse', description: 'Parse input identifier to RDKit mol', alwaysOn: true },
  { index: 9, name: 'Filter', description: 'Apply molecular property filters', alwaysOn: true },
  { index: 10, name: 'Canonical', description: 'Generate canonical SMILES output', alwaysOn: true },
];

// =============================================================================
// Props
// =============================================================================

interface StepToggleGridProps {
  /** Current pipeline configuration — drives checkbox states. */
  config: QSARReadyConfig;
  /** Called when a togglable step boolean key is changed. */
  onToggleStep: (key: string, enabled: boolean) => void;
  /** Called when a numeric/non-boolean config field changes (filter thresholds). */
  onUpdateConfig: (partial: Partial<QSARReadyConfig>) => void;
}

// =============================================================================
// Component
// =============================================================================

/**
 * 10-step toggle grid for the QSAR-ready pipeline.
 *
 * 7 togglable steps from PIPELINE_STEPS + 3 always-on fixed steps (Parse/Filter/Canonical).
 * The Filter step expands to show threshold inputs (min_heavy_atoms, max_heavy_atoms,
 * max_mw, remove_inorganics) when the chevron is clicked.
 *
 * Layout: `grid grid-cols-2 md:grid-cols-5 gap-3` per UI-SPEC.
 * Accessibility: checkbox aria-labels, 44px minimum touch target per UI-SPEC.
 */
export function StepToggleGrid({ config, onToggleStep, onUpdateConfig }: StepToggleGridProps) {
  const [filterExpanded, setFilterExpanded] = useState(false);

  // Combine all steps sorted by index for rendering
  const allSteps: Array<(typeof PIPELINE_STEPS[0]) | FixedStep> = [
    ...FIXED_STEPS,
    ...PIPELINE_STEPS,
  ].sort((a, b) => a.index - b.index);

  return (
    <div className="space-y-3">
      {/* Step toggle grid */}
      <div className="grid grid-cols-2 md:grid-cols-5 gap-3">
        {allSteps.map((step) => {
          const isFixed = 'alwaysOn' in step && step.alwaysOn;
          const isFilterStep = step.name === 'Filter';

          // For togglable steps, derive enabled state from config
          let isEnabled = true;
          if (!isFixed) {
            const toggleable = step as typeof PIPELINE_STEPS[0];
            isEnabled = config[toggleable.key] as boolean;
          }

          return (
            <div
              key={step.index}
              className={cn(
                'flex flex-col gap-1 rounded-xl p-3 border transition-opacity duration-150',
                'bg-[var(--color-surface-elevated)] border-[var(--color-border)]',
                !isEnabled && 'opacity-50',
              )}
            >
              {/* Row: checkbox + step number + name */}
              <div className="flex items-center gap-2 min-h-[44px]">
                <input
                  type="checkbox"
                  checked={isEnabled}
                  disabled={isFixed}
                  onChange={(e) => {
                    if (!isFixed) {
                      const toggleable = step as typeof PIPELINE_STEPS[0];
                      onToggleStep(toggleable.key, e.target.checked);
                    }
                  }}
                  aria-label={`${step.name} — ${isEnabled ? 'enabled' : 'disabled'}`}
                  className={cn(
                    'w-4 h-4 rounded accent-[var(--color-primary)] cursor-pointer shrink-0',
                    isFixed && 'cursor-not-allowed opacity-60',
                  )}
                />
                {/* Step number badge */}
                <span className="inline-flex items-center justify-center w-5 h-5 rounded-full text-[10px] font-semibold bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)] shrink-0">
                  {step.index}
                </span>
                {/* Step name */}
                <span className="text-xs font-semibold text-[var(--color-text-primary)] truncate">
                  {step.name}
                </span>
                {/* Filter expand chevron */}
                {isFilterStep && (
                  <button
                    type="button"
                    onClick={() => setFilterExpanded((v) => !v)}
                    className="ml-auto text-[var(--color-text-muted)] hover:text-[var(--color-text-primary)] transition-colors"
                    aria-label={filterExpanded ? 'Collapse filter options' : 'Expand filter options'}
                  >
                    <ChevronDown
                      className={cn(
                        'w-3.5 h-3.5 transition-transform duration-200',
                        filterExpanded && 'rotate-180',
                      )}
                    />
                  </button>
                )}
              </div>

              {/* Step description */}
              <p className="text-xs text-[var(--color-text-muted)] pl-[26px] leading-tight">
                {step.description}
              </p>

              {/* Filter threshold inputs (expandable, collapsed by default) */}
              {isFilterStep && filterExpanded && (
                <div className="mt-2 pt-2 border-t border-[var(--color-border)] space-y-2 pl-[26px]">
                  <label className="flex flex-col gap-1">
                    <span className="text-xs text-[var(--color-text-muted)]">Min heavy atoms</span>
                    <input
                      type="number"
                      min={1}
                      max={50}
                      value={config.min_heavy_atoms}
                      onChange={(e) => onUpdateConfig({ min_heavy_atoms: Number(e.target.value) })}
                      className="w-full px-2 py-1 text-xs rounded-lg border border-[var(--color-border)] bg-[var(--color-surface)] text-[var(--color-text-primary)] focus:outline-none focus:border-[var(--color-primary)]"
                    />
                  </label>
                  <label className="flex flex-col gap-1">
                    <span className="text-xs text-[var(--color-text-muted)]">Max heavy atoms</span>
                    <input
                      type="number"
                      min={1}
                      max={500}
                      value={config.max_heavy_atoms}
                      onChange={(e) => onUpdateConfig({ max_heavy_atoms: Number(e.target.value) })}
                      className="w-full px-2 py-1 text-xs rounded-lg border border-[var(--color-border)] bg-[var(--color-surface)] text-[var(--color-text-primary)] focus:outline-none focus:border-[var(--color-primary)]"
                    />
                  </label>
                  <label className="flex flex-col gap-1">
                    <span className="text-xs text-[var(--color-text-muted)]">Max MW (Da)</span>
                    <input
                      type="number"
                      min={1}
                      max={5000}
                      value={config.max_mw}
                      onChange={(e) => onUpdateConfig({ max_mw: Number(e.target.value) })}
                      className="w-full px-2 py-1 text-xs rounded-lg border border-[var(--color-border)] bg-[var(--color-surface)] text-[var(--color-text-primary)] focus:outline-none focus:border-[var(--color-primary)]"
                    />
                  </label>
                  <label className="flex items-center gap-2">
                    <input
                      type="checkbox"
                      checked={config.remove_inorganics}
                      onChange={(e) => onUpdateConfig({ remove_inorganics: e.target.checked })}
                      aria-label="Remove inorganics — enabled/disabled"
                      className="w-4 h-4 rounded accent-[var(--color-primary)] cursor-pointer"
                    />
                    <span className="text-xs text-[var(--color-text-muted)]">Remove inorganics</span>
                  </label>
                </div>
              )}
            </div>
          );
        })}
      </div>
    </div>
  );
}
