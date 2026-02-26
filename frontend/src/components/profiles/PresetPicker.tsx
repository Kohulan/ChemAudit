/**
 * Card grid showing 8 preset scoring profiles with "Apply" and "Duplicate & Customize" actions.
 */

import { motion } from 'framer-motion';
import { Check, Copy, Beaker } from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';
import { Badge } from '../ui/Badge';
import { cn } from '../../lib/utils';
import type { ScoringProfile, ThresholdRange } from '../../types/workflow';

interface PresetPickerProps {
  presets: ScoringProfile[];
  activeProfileId: number | null;
  onApply: (profile: ScoringProfile) => void;
  onDuplicate: (profile: ScoringProfile) => void;
}

/** Format a threshold range for display. */
function formatRange(range: ThresholdRange | undefined, unit?: string): string {
  if (!range) return '--';
  const u = unit ?? '';
  if (range.min !== undefined && range.max !== undefined) {
    return `${range.min}${u} - ${range.max}${u}`;
  }
  if (range.max !== undefined) return `<= ${range.max}${u}`;
  if (range.min !== undefined) return `>= ${range.min}${u}`;
  return '--';
}

const KEY_LABELS: Record<string, { label: string; unit?: string }> = {
  mw: { label: 'MW', unit: '' },
  logp: { label: 'LogP' },
  hbd: { label: 'HBD' },
  hba: { label: 'HBA' },
  tpsa: { label: 'TPSA', unit: '\u00C5\u00B2' },
  rotatable_bonds: { label: 'RotB' },
  aromatic_rings: { label: 'ArRings' },
  fsp3: { label: 'Fsp3' },
};

export function PresetPicker({ presets, activeProfileId, onApply, onDuplicate }: PresetPickerProps) {
  if (presets.length === 0) {
    return (
      <div className="text-center py-8 text-[var(--color-text-muted)]">
        <Beaker className="w-8 h-8 mx-auto mb-2 opacity-50" />
        <p className="text-sm">No preset profiles available</p>
      </div>
    );
  }

  return (
    <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
      {presets.map((preset, index) => {
        const isActive = preset.id === activeProfileId;
        // Show up to 4 key thresholds
        const displayKeys = Object.keys(preset.thresholds).slice(0, 4);

        return (
          <motion.div
            key={preset.id}
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ delay: index * 0.05 }}
            className={cn(
              'rounded-xl border p-4 transition-all duration-200',
              isActive
                ? 'border-[var(--color-primary)] bg-[var(--color-primary)]/5 shadow-[0_0_16px_var(--glow-primary)]'
                : 'border-[var(--color-border)] bg-[var(--color-surface-elevated)] hover:border-[var(--color-border-strong)]'
            )}
          >
            {/* Header */}
            <div className="flex items-start justify-between mb-2">
              <div className="flex-1 min-w-0">
                <div className="flex items-center gap-2">
                  <h4 className="font-semibold text-sm text-[var(--color-text-primary)] truncate">
                    {preset.name}
                  </h4>
                  {isActive && (
                    <Badge variant="success" size="sm" icon={<Check className="w-3 h-3" />}>
                      Active
                    </Badge>
                  )}
                </div>
                {preset.description && (
                  <p className="text-xs text-[var(--color-text-muted)] mt-0.5 line-clamp-2">
                    {preset.description}
                  </p>
                )}
              </div>
            </div>

            {/* Key thresholds preview */}
            <div className="grid grid-cols-2 gap-1.5 mt-3 mb-3">
              {displayKeys.map((key) => {
                const meta = KEY_LABELS[key];
                return (
                  <div key={key} className="text-[10px] bg-[var(--color-surface-sunken)] rounded-md px-2 py-1">
                    <span className="text-[var(--color-text-muted)]">{meta?.label ?? key}: </span>
                    <span className="text-[var(--color-text-primary)] font-medium">
                      {formatRange(preset.thresholds[key], meta?.unit)}
                    </span>
                  </div>
                );
              })}
            </div>

            {/* Actions */}
            <div className="flex items-center gap-2">
              <ClayButton
                size="sm"
                variant={isActive ? 'default' : 'primary'}
                onClick={() => onApply(preset)}
                className="flex-1"
              >
                {isActive ? 'Applied' : 'Apply'}
              </ClayButton>
              <ClayButton
                size="sm"
                onClick={() => onDuplicate(preset)}
                leftIcon={<Copy className="w-3.5 h-3.5" />}
              >
                Customize
              </ClayButton>
            </div>
          </motion.div>
        );
      })}
    </div>
  );
}
