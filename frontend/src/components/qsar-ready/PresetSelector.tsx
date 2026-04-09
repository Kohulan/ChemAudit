import { cn } from '../../lib/utils';
import { ClayButton } from '../ui/ClayButton';

// =============================================================================
// Types
// =============================================================================

interface PresetInfo {
  id: string;
  label: string;
  subtitle: string;
}

const PRESETS: PresetInfo[] = [
  {
    id: 'qsar-2d',
    label: 'QSAR-2D',
    subtitle: 'Strip stereo and isotopes for 2D descriptors',
  },
  {
    id: 'qsar-3d',
    label: 'QSAR-3D',
    subtitle: 'Preserve stereo, strip isotopes for 3D conformers',
  },
  {
    id: 'minimal',
    label: 'Minimal',
    subtitle: 'Parse, desalt, normalize, and canonicalize only',
  },
];

// Human-readable names for preset IDs (used in "modified from" label)
const PRESET_LABELS: Record<string, string> = {
  'qsar-2d': 'QSAR-2D',
  'qsar-3d': 'QSAR-3D',
  'minimal': 'Minimal',
};

// =============================================================================
// Props
// =============================================================================

interface PresetSelectorProps {
  /** The currently active preset ID, or null if using a custom config. */
  activePreset: string | null;
  /** If non-null, shows "Custom (modified from {modifiedFrom})" label. */
  modifiedFrom: string | null;
  /** Called when the user clicks a preset button. */
  onSelectPreset: (id: string) => void;
}

// =============================================================================
// Component
// =============================================================================

/**
 * Three preset buttons (QSAR-2D, QSAR-3D, Minimal) with subtitles.
 *
 * Active preset gets outline variant with primary border highlight.
 * When a preset has been modified, shows "Custom (modified from X)" below.
 *
 * Per Phase 10 UI-SPEC D-09 preset selector pattern.
 */
export function PresetSelector({ activePreset, modifiedFrom, onSelectPreset }: PresetSelectorProps) {
  return (
    <div className="space-y-3">
      {/* Preset button row */}
      <div className="flex flex-wrap gap-3">
        {PRESETS.map((preset) => {
          const isActive = activePreset === preset.id;

          return (
            <div key={preset.id} className="flex flex-col gap-1">
              <ClayButton
                variant="outline"
                size="sm"
                onClick={() => onSelectPreset(preset.id)}
                className={cn(
                  'min-w-[110px]',
                  isActive && 'border-2 border-[var(--color-primary)] text-[var(--color-primary)] font-semibold',
                )}
              >
                {preset.label}
              </ClayButton>
              <span className="text-xs text-[var(--color-text-muted)] max-w-[140px]">
                {preset.subtitle}
              </span>
            </div>
          );
        })}
      </div>

      {/* Modified-from label (shown when a preset was customized) */}
      {modifiedFrom !== null && (
        <p className="text-xs text-[var(--color-text-muted)]">
          Custom (modified from {PRESET_LABELS[modifiedFrom] ?? modifiedFrom})
        </p>
      )}
    </div>
  );
}
