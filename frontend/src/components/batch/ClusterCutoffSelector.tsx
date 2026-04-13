/**
 * ClusterCutoffSelector
 *
 * Inline row with a Tanimoto cutoff selector and Cluster/Re-cluster button.
 * Follows the ActivityColumnSelector pattern for accessibility.
 */

import { Loader2 } from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';

interface ClusterCutoffSelectorProps {
  cutoff: number;
  onCutoffChange: (val: number) => void;
  onRecluster: () => void;
  isComputing: boolean;
  hasTriggered: boolean;
}

const CUTOFF_OPTIONS = [
  { value: 0.2, label: '0.2 (loose)' },
  { value: 0.3, label: '0.3' },
  { value: 0.35, label: '0.35 (default)' },
  { value: 0.4, label: '0.4' },
  { value: 0.5, label: '0.5' },
  { value: 0.6, label: '0.6 (strict)' },
];

export function ClusterCutoffSelector({
  cutoff,
  onCutoffChange,
  onRecluster,
  isComputing,
  hasTriggered,
}: ClusterCutoffSelectorProps) {
  const buttonText = hasTriggered ? 'Re-cluster' : 'Cluster';
  const buttonAriaLabel = hasTriggered
    ? 'Re-cluster with new cutoff'
    : 'Cluster molecules using Butina algorithm';

  return (
    <div className="flex items-center gap-3 flex-wrap">
      <label
        htmlFor="tanimoto-cutoff"
        className="text-sm font-medium text-[var(--color-text-secondary)]"
      >
        Tanimoto Cutoff
      </label>
      <select
        id="tanimoto-cutoff"
        value={cutoff}
        onChange={(e) => onCutoffChange(Number(e.target.value))}
        disabled={isComputing}
        className="rounded-lg border border-[var(--color-border)] bg-[var(--color-surface)] text-sm text-[var(--color-text-primary)] px-3 py-1.5 focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]/30 disabled:opacity-50"
      >
        {CUTOFF_OPTIONS.map((opt) => (
          <option key={opt.value} value={opt.value}>
            {opt.label}
          </option>
        ))}
      </select>
      <ClayButton
        variant="primary"
        size="sm"
        onClick={onRecluster}
        disabled={isComputing}
        aria-label={buttonAriaLabel}
        leftIcon={isComputing ? <Loader2 className="w-3.5 h-3.5 animate-spin" /> : undefined}
      >
        {buttonText}
      </ClayButton>
    </div>
  );
}
