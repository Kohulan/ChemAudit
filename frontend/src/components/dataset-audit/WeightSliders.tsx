import { useCallback } from 'react';
import { ClayButton } from '../ui/ClayButton';
import { InfoTooltip } from '../ui/Tooltip';
import type { DatasetWeightProfile } from '../../types/dataset_intelligence';
import { PRESET_WEIGHT_IDS } from '../../types/dataset_intelligence';

// =============================================================================
// Types
// =============================================================================

interface WeightSlidersProps {
  /** Current weight values (sum to 1.0). */
  weights: Record<string, number>;
  /** Active profile ID, or null for custom. */
  activeProfileId: string | null;
  /** All saved profiles (preset + custom). */
  savedProfiles: DatasetWeightProfile[];
  /** Update a single weight (hook handles normalization). */
  onUpdateWeight: (key: string, value: number) => void;
  /** Save current weights as a named profile. */
  onSaveProfile: (name: string) => void;
  /** Delete a custom profile by ID. */
  onDeleteProfile: (id: string) => void;
  /** Load weights from a saved profile. */
  onLoadProfile: (id: string) => void;
  /** Reset weights to literature defaults. */
  onResetToDefaults: () => void;
}

// =============================================================================
// Display names for sub-score keys
// =============================================================================

const WEIGHT_LABELS: Record<string, string> = {
  parsability: 'Parsability',
  stereo: 'Stereo Completeness',
  uniqueness: 'Uniqueness',
  alerts: 'Alert Prevalence',
  std_consistency: 'Std. Consistency',
};

/** Ordered keys for consistent rendering. */
const WEIGHT_KEYS = ['parsability', 'stereo', 'uniqueness', 'alerts', 'std_consistency'];

// =============================================================================
// Component
// =============================================================================

/**
 * 5 configurable weight sliders with profile save/load/delete.
 *
 * Per UI-SPEC:
 * - Each slider: 0 to 1.0, step 0.05
 * - Display current value as percentage
 * - Profile management: dropdown + Save/Delete/Reset buttons
 * - Save uses window.prompt per UI-SPEC interaction contract
 * - Delete uses window.confirm
 * - Preset profiles cannot be deleted (PRESET_WEIGHT_IDS guard)
 */
export function WeightSliders({
  weights,
  activeProfileId,
  savedProfiles,
  onUpdateWeight,
  onSaveProfile,
  onDeleteProfile,
  onLoadProfile,
  onResetToDefaults,
}: WeightSlidersProps) {
  const isPreset = activeProfileId ? PRESET_WEIGHT_IDS.has(activeProfileId) : false;

  const handleSave = useCallback(() => {
    const name = window.prompt('Name your custom weight profile:');
    if (name && name.trim()) {
      onSaveProfile(name.trim());
    }
  }, [onSaveProfile]);

  const handleDelete = useCallback(() => {
    if (!activeProfileId || isPreset) return;
    const profile = savedProfiles.find((p) => p.id === activeProfileId);
    const profileName = profile?.name ?? activeProfileId;
    const confirmed = window.confirm(
      `Delete '${profileName}'? This cannot be undone.`,
    );
    if (confirmed) {
      onDeleteProfile(activeProfileId);
    }
  }, [activeProfileId, isPreset, savedProfiles, onDeleteProfile]);

  const handleProfileChange = useCallback(
    (e: React.ChangeEvent<HTMLSelectElement>) => {
      onLoadProfile(e.target.value);
    },
    [onLoadProfile],
  );

  return (
    <div className="space-y-4">
      {/* Section heading */}
      <div className="flex items-center gap-1.5">
        <h3 className="text-base font-semibold font-display text-[var(--color-text-primary)]">
          Sub-Score Weights
        </h3>
        <InfoTooltip
          title="Sub-Score Weights"
          content={
            <div className="text-xs space-y-1">
              <p>Control how much each quality dimension contributes to the overall health score.</p>
              <ul className="mt-1 text-white/70 space-y-0.5">
                <li><strong>Parsability</strong> &mdash; % of SMILES that RDKit can parse</li>
                <li><strong>Stereo</strong> &mdash; % with fully defined stereochemistry</li>
                <li><strong>Uniqueness</strong> &mdash; % of non-duplicate structures</li>
                <li><strong>Alerts</strong> &mdash; % free of structural alerts</li>
                <li><strong>Std. Consistency</strong> &mdash; % where all pipelines agree</li>
              </ul>
              <p className="mt-1 text-white/60">Weights are normalized to sum to 100%. Adjust to emphasize dimensions most critical to your ML use case.</p>
            </div>
          }
          size="small"
        />
      </div>

      {/* Sliders */}
      <div className="space-y-3">
        {WEIGHT_KEYS.map((key) => {
          const value = weights[key] ?? 0;
          const label = WEIGHT_LABELS[key] ?? key;
          const percent = Math.round(value * 100);
          const sliderId = `weight-slider-${key}`;

          return (
            <div key={key}>
              <div className="flex items-center justify-between mb-1">
                <label
                  htmlFor={sliderId}
                  className="text-xs text-[var(--color-text-secondary)]"
                >
                  {label}
                </label>
                <span className="text-xs font-medium text-[var(--color-text-primary)] tabular-nums">
                  {percent}%
                </span>
              </div>
              <input
                id={sliderId}
                type="range"
                min={0}
                max={1}
                step={0.05}
                value={value}
                onChange={(e) => onUpdateWeight(key, parseFloat(e.target.value))}
                className="w-full h-2 rounded-lg appearance-none cursor-pointer accent-[var(--color-primary)]"
                aria-valuemin={0}
                aria-valuemax={1}
                aria-valuenow={value}
                aria-label={`${label} weight: ${percent}%`}
              />
            </div>
          );
        })}
      </div>

      {/* Profile management row */}
      <div className="border-t border-[var(--color-border)] pt-3 space-y-2">
        {/* Profile selector */}
        <select
          value={activeProfileId ?? ''}
          onChange={handleProfileChange}
          className="w-full text-xs px-3 py-2 rounded-lg border border-[var(--color-border)] bg-[var(--color-surface)] text-[var(--color-text-primary)] focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]"
          aria-label="Weight profile selector"
        >
          {!activeProfileId && (
            <option value="" disabled>
              Custom (unsaved)
            </option>
          )}
          {savedProfiles.map((profile) => (
            <option key={profile.id} value={profile.id}>
              {profile.name}
            </option>
          ))}
        </select>

        {/* Action buttons */}
        <div className="flex gap-2">
          <ClayButton variant="default" size="sm" onClick={handleSave} className="flex-1 text-xs">
            Save
          </ClayButton>
          {activeProfileId && !isPreset && (
            <ClayButton
              variant="danger"
              size="sm"
              onClick={handleDelete}
              className="flex-1 text-xs"
            >
              Delete
            </ClayButton>
          )}
          <ClayButton
            variant="ghost"
            size="sm"
            onClick={onResetToDefaults}
            className="flex-1 text-xs"
          >
            Reset
          </ClayButton>
        </div>
      </div>
    </div>
  );
}
