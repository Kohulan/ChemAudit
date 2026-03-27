import { useState, useId } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown, Zap } from 'lucide-react';
import { ClayCard } from '../ui/ClayCard';
import { ClayButton } from '../ui/ClayButton';
import type { LEResult } from '../../types/profiler';

interface LigandEfficiencyPanelProps {
  smiles: string;
  computeEfficiency: (smiles: string, activityValue: number, activityType: string) => Promise<LEResult | null>;
}

/** Activity type options per D-10. */
const ACTIVITY_TYPES = [
  { value: 'IC50_nM', label: 'IC50 (nM)' },
  { value: 'IC50_uM', label: 'IC50 (uM)' },
  { value: 'Ki_nM', label: 'Ki (nM)' },
  { value: 'pIC50', label: 'pIC50' },
  { value: 'pKd', label: 'pKd' },
] as const;

/** Metric descriptions for display in results grid. */
const METRIC_DESCRIPTIONS: Record<string, { label: string; description: string; unit: string }> = {
  LE: {
    label: 'LE',
    description: 'Ligand Efficiency',
    unit: 'kcal/mol/HA',
  },
  LLE: {
    label: 'LLE',
    description: 'Lipophilic Ligand Efficiency',
    unit: 'pIC50 - LogP',
  },
  LELP: {
    label: 'LELP',
    description: 'Ligand-Efficiency-dependent Lipophilicity',
    unit: 'LogP / LE',
  },
  BEI: {
    label: 'BEI',
    description: 'Binding Efficiency Index',
    unit: 'pIC50 × 1000 / MW',
  },
  SEI: {
    label: 'SEI',
    description: 'Surface Efficiency Index',
    unit: 'pIC50 × 100 / TPSA',
  },
};

/**
 * Ligand Efficiency Metrics panel.
 *
 * Collapsed by default (D-09). On expand: activity type selector + value input +
 * "Compute Efficiency" button. Results show LE, LLE, LELP, BEI, SEI inline.
 */
export function LigandEfficiencyPanel({ smiles, computeEfficiency }: LigandEfficiencyPanelProps) {
  const [expanded, setExpanded] = useState(false);
  const [activityType, setActivityType] = useState<string>('IC50_nM');
  const [activityValue, setActivityValue] = useState<string>('');
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<LEResult | null>(null);
  const [computeError, setComputeError] = useState<string | null>(null);
  const panelId = useId();

  const handleCompute = async () => {
    const value = parseFloat(activityValue);
    if (isNaN(value) || value <= 0) {
      setComputeError('Please enter a valid positive activity value.');
      return;
    }
    setComputeError(null);
    setLoading(true);
    try {
      const data = await computeEfficiency(smiles, value, activityType);
      setResult(data);
      if (!data) {
        setComputeError('Unable to compute ligand efficiency. Check your activity value and try again.');
      }
    } finally {
      setLoading(false);
    }
  };

  return (
    <ClayCard size="md" className="overflow-hidden">
      {/* Section header — collapsible toggle */}
      <button
        type="button"
        onClick={() => setExpanded((prev) => !prev)}
        aria-expanded={expanded}
        aria-controls={panelId}
        className="w-full text-left cursor-pointer"
      >
        <div className="flex items-center justify-between gap-3">
          <div className="flex items-center gap-2">
            <Zap className="w-4 h-4 text-text-muted flex-shrink-0" />
            <h3 className="text-lg font-semibold text-text-primary font-display">
              Ligand Efficiency Metrics
            </h3>
          </div>
          <motion.div
            animate={{ rotate: expanded ? 180 : 0 }}
            transition={{ duration: 0.2, ease: 'easeOut' }}
            className="text-text-muted flex-shrink-0"
          >
            <ChevronDown className="w-5 h-5" />
          </motion.div>
        </div>
      </button>

      {/* Expandable content */}
      <AnimatePresence initial={false}>
        {expanded && (
          <motion.div
            id={panelId}
            role="region"
            key="content"
            initial={{ height: 0, opacity: 0 }}
            animate={{ height: 'auto', opacity: 1 }}
            exit={{ height: 0, opacity: 0 }}
            transition={{ duration: 0.5, ease: 'easeOut' }}
            className="overflow-hidden"
          >
            <div className="pt-4 border-t border-[var(--color-border)] mt-4 space-y-4">
              {/* Activity input row */}
              <div className="flex flex-wrap items-end gap-3">
                {/* Activity type selector */}
                <div className="flex-1 min-w-[140px]">
                  <label className="block text-xs font-medium text-text-secondary mb-1">
                    Activity type
                  </label>
                  <select
                    value={activityType}
                    onChange={(e) => setActivityType(e.target.value)}
                    className="w-full text-sm bg-surface-elevated border border-[var(--color-border)] rounded-xl px-3 py-2 text-text-primary"
                  >
                    {ACTIVITY_TYPES.map(({ value, label }) => (
                      <option key={value} value={value}>{label}</option>
                    ))}
                  </select>
                </div>

                {/* Activity value input */}
                <div className="flex-1 min-w-[120px]">
                  <label className="block text-xs font-medium text-text-secondary mb-1">
                    Activity value
                  </label>
                  <input
                    type="number"
                    value={activityValue}
                    onChange={(e) => setActivityValue(e.target.value)}
                    placeholder="e.g. 100"
                    min={0}
                    step="any"
                    className="w-full text-sm bg-surface-elevated border border-[var(--color-border)] rounded-xl px-3 py-2 text-text-primary"
                    onKeyDown={(e) => e.key === 'Enter' && handleCompute()}
                  />
                </div>

                {/* Compute button */}
                <ClayButton
                  variant="primary"
                  size="md"
                  onClick={handleCompute}
                  loading={loading}
                  disabled={!activityValue || loading}
                >
                  Compute Efficiency
                </ClayButton>
              </div>

              {/* Input error */}
              {computeError && (
                <p className="text-xs text-status-error">{computeError}</p>
              )}

              {/* Results */}
              {result && (
                <motion.div
                  initial={{ opacity: 0, y: 4 }}
                  animate={{ opacity: 1, y: 0 }}
                  transition={{ duration: 0.3 }}
                  className="space-y-3"
                >
                  <div className="flex items-baseline gap-2">
                    <span className="text-xs font-semibold text-text-muted uppercase tracking-wide">
                      pIC50
                    </span>
                    <span className="text-lg font-semibold tabular-nums text-text-primary font-display">
                      {result.pIC50.toFixed(3)}
                    </span>
                  </div>

                  <div className="grid grid-cols-2 sm:grid-cols-3 lg:grid-cols-5 gap-3">
                    {(['LE', 'LLE', 'LELP', 'BEI', 'SEI'] as const).map((key) => {
                      const meta = METRIC_DESCRIPTIONS[key];
                      return (
                        <div key={key} className="bg-surface-sunken rounded-xl p-3">
                          <p className="text-xs font-semibold text-text-muted uppercase tracking-wide">
                            {meta.label}
                          </p>
                          <p className="text-lg font-semibold tabular-nums text-text-primary font-display mt-0.5">
                            {result[key].toFixed(3)}
                          </p>
                          <p className="text-[10px] text-text-muted mt-0.5 leading-tight">
                            {meta.description}
                          </p>
                          <p className="text-[10px] text-text-muted font-mono mt-0.5">
                            {meta.unit}
                          </p>
                        </div>
                      );
                    })}
                  </div>
                </motion.div>
              )}
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </ClayCard>
  );
}
