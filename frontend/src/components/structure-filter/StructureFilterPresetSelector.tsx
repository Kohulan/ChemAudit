import { useState } from 'react';
import { Trash2 } from 'lucide-react';
import { cn } from '../../lib/utils';
import { ClayButton } from '../ui/ClayButton';
import { PRESET_IDS } from '../../types/genchem';
import type { GenChemProfile } from '../../types/genchem';

// =============================================================================
// Constants
// =============================================================================

interface PresetInfo {
  id: string;
  label: string;
  subtitle: string;
}

/** 4 built-in presets (D-11, D-12 locked). */
const PRESETS: PresetInfo[] = [
  {
    id: 'drug_like',
    label: 'Drug-like',
    subtitle: 'Lipinski-compatible: MW 200–500, LogP −1–5',
  },
  {
    id: 'lead_like',
    label: 'Lead-like',
    subtitle: 'Lead optimization: MW 200–350, LogP −1–3.5',
  },
  {
    id: 'fragment_like',
    label: 'Fragment-like',
    subtitle: 'FBDD: MW 100–300, LogP −1–3, max 3 rings',
  },
  {
    id: 'permissive',
    label: 'Permissive',
    subtitle: 'Broad range: MW 100–800, alerts off',
  },
];

/** Human-readable names for preset IDs — used in "modified from" label. */
const PRESET_LABELS: Record<string, string> = {
  drug_like: 'Drug-like',
  lead_like: 'Lead-like',
  fragment_like: 'Fragment-like',
  permissive: 'Permissive',
};

// =============================================================================
// Props
// =============================================================================

interface GenChemPresetSelectorProps {
  /** The currently active preset ID, or null if a custom config is active. */
  activePreset: string | null;
  /** If non-null, shows "Custom (modified from {presetName})" label. */
  modifiedFrom: string | null;
  /** Custom profiles saved by the user. */
  profiles: GenChemProfile[];
  /** Called when a preset pill or custom profile is clicked. */
  onSelectPreset: (id: string) => void;
  /** Called when the user saves a custom configuration. */
  onSaveProfile: (name: string) => void;
  /** Called when the user deletes a custom profile. No-op for built-in presets. */
  onDeleteProfile: (id: string) => void;
}

// =============================================================================
// Component
// =============================================================================

/**
 * Preset selector for GenChem Filter configuration.
 *
 * Shows 4 built-in preset pills (Drug-like, Lead-like, Fragment-like, Permissive)
 * and custom profiles below with delete buttons.
 *
 * Active preset: ring-2 highlight + filled background.
 * Modified label: shown when a preset was edited (D-13).
 * Save: window.prompt for name (D-13 — same as PipelineConfigPanel Phase 10 pattern).
 * Delete: inline confirmation before calling onDeleteProfile.
 * Guard: delete button never rendered for built-in preset IDs (PRESET_IDS guard).
 *
 * Per Phase 11 UI-SPEC D-11, D-13 copywriting and interaction contracts.
 */
export function GenChemPresetSelector({
  activePreset,
  modifiedFrom,
  profiles,
  onSelectPreset,
  onSaveProfile,
  onDeleteProfile,
}: GenChemPresetSelectorProps) {
  const [deletingId, setDeletingId] = useState<string | null>(null);

  const handleSaveConfig = () => {
    const name = window.prompt('Name your custom configuration:');
    if (name && name.trim()) {
      onSaveProfile(name.trim());
    }
  };

  const handleDeleteClick = (id: string, name: string) => {
    setDeletingId(id);
    // Inline confirmation pattern (UI-SPEC copywriting)
    const confirmed = window.confirm(`Delete '${name}'? This cannot be undone.`);
    setDeletingId(null);
    if (confirmed) {
      // Guard: never call delete for built-in presets (PRESET_IDS no-op)
      if (!PRESET_IDS.has(id)) {
        onDeleteProfile(id);
      }
    }
  };

  return (
    <div className="space-y-3">
      {/* ── Built-in preset pills ── */}
      <div className="flex flex-wrap gap-2">
        {PRESETS.map((preset) => {
          const isActive = activePreset === preset.id;
          return (
            <button
              key={preset.id}
              type="button"
              onClick={() => onSelectPreset(preset.id)}
              className={cn(
                'px-3 py-1.5 rounded-xl border text-xs font-medium transition-all duration-150',
                isActive
                  ? 'bg-[var(--color-primary)]/10 border-[var(--color-primary)] text-[var(--color-primary)] ring-2 ring-[var(--color-primary)]/30'
                  : 'bg-[var(--color-surface)] border-[var(--color-border)] text-[var(--color-text-secondary)] hover:border-[var(--color-text-muted)] hover:text-[var(--color-text-primary)]',
              )}
              aria-pressed={isActive}
            >
              {preset.label}
            </button>
          );
        })}
      </div>

      {/* ── Preset subtitles row ── */}
      <div className="flex flex-wrap gap-2">
        {PRESETS.map((preset) => {
          const isActive = activePreset === preset.id;
          if (!isActive) return null;
          return (
            <p key={preset.id} className="text-xs text-[var(--color-text-muted)]">
              {preset.subtitle}
            </p>
          );
        })}
      </div>

      {/* ── Modified-from label (D-13) ── */}
      {modifiedFrom !== null && (
        <p className="text-xs text-[var(--color-text-muted)]">
          Custom (modified from {PRESET_LABELS[modifiedFrom] ?? modifiedFrom})
        </p>
      )}

      {/* ── Save Config button ── */}
      <ClayButton
        variant="default"
        size="sm"
        onClick={handleSaveConfig}
        className="text-xs"
      >
        Save Config
      </ClayButton>

      {/* ── Custom profiles ── */}
      {profiles.length > 0 && (
        <div className="space-y-1 pt-1">
          <p className="text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wide">
            Saved Configurations
          </p>
          <div className="flex flex-wrap gap-2">
            {profiles.map((profile) => {
              const isActiveProfile = activePreset === profile.id;
              const isDeleting = deletingId === profile.id;
              return (
                <div
                  key={profile.id}
                  className={cn(
                    'flex items-center gap-1.5 px-3 py-1.5 rounded-xl border text-xs',
                    'bg-[var(--color-surface-elevated)] border-[var(--color-border)]',
                    isActiveProfile &&
                      'border-[var(--color-primary)] text-[var(--color-primary)]',
                    isDeleting && 'opacity-50',
                  )}
                >
                  <button
                    type="button"
                    onClick={() => onSelectPreset(profile.id)}
                    className={cn(
                      'font-medium hover:text-[var(--color-primary)] transition-colors',
                      isActiveProfile
                        ? 'text-[var(--color-primary)]'
                        : 'text-[var(--color-text-primary)]',
                    )}
                  >
                    {profile.name}
                  </button>
                  {/* Guard: no delete button for built-in preset IDs */}
                  {!PRESET_IDS.has(profile.id) && (
                    <button
                      type="button"
                      onClick={() => handleDeleteClick(profile.id, profile.name)}
                      className="text-[var(--color-text-muted)] hover:text-red-500 transition-colors"
                      aria-label={`Delete configuration '${profile.name}'`}
                    >
                      <Trash2 className="w-3 h-3" />
                    </button>
                  )}
                </div>
              );
            })}
          </div>
        </div>
      )}
    </div>
  );
}
