/**
 * Visual scoring profile builder with sliders for thresholds and number
 * inputs for weights. Supports create, edit, export, and import modes.
 */

import { useState, useEffect, useCallback } from 'react';
import { motion } from 'framer-motion';
import { Save, X, Upload, Download, RotateCcw } from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';
import { cn } from '../../lib/utils';
import type {
  ScoringProfile,
  ScoringProfileCreate,
  ScoringProfileExport,
  ThresholdRange,
} from '../../types/workflow';

interface ThresholdConfig {
  key: string;
  label: string;
  min: number;
  max: number;
  step: number;
  unit: string;
}

const THRESHOLD_CONFIGS: ThresholdConfig[] = [
  { key: 'mw', label: 'Molecular Weight', min: 0, max: 1000, step: 10, unit: 'Da' },
  { key: 'logp', label: 'LogP', min: -5, max: 10, step: 0.5, unit: '' },
  { key: 'tpsa', label: 'TPSA', min: 0, max: 300, step: 5, unit: 'A\u00B2' },
  { key: 'hbd', label: 'H-Bond Donors', min: 0, max: 10, step: 1, unit: '' },
  { key: 'hba', label: 'H-Bond Acceptors', min: 0, max: 20, step: 1, unit: '' },
  { key: 'rotatable_bonds', label: 'Rotatable Bonds', min: 0, max: 20, step: 1, unit: '' },
  { key: 'aromatic_rings', label: 'Aromatic Rings', min: 0, max: 10, step: 1, unit: '' },
  { key: 'fsp3', label: 'Fsp3', min: 0, max: 1, step: 0.05, unit: '' },
];

interface ProfileBuilderProps {
  /** Existing profile to edit (null for create mode). */
  profile?: ScoringProfile | null;
  /** Called when the user saves the profile. */
  onSave: (data: ScoringProfileCreate) => Promise<void>;
  /** Called when the user cancels. */
  onCancel: () => void;
  /** Called to export a profile as JSON. */
  onExport?: (id: number) => Promise<ScoringProfileExport>;
  /** Called to import a profile from JSON. */
  onImport?: (data: ScoringProfileExport) => Promise<void>;
  /** Is saving in progress? */
  isSaving?: boolean;
}

export function ProfileBuilder({
  profile,
  onSave,
  onCancel,
  onExport,
  onImport,
  isSaving = false,
}: ProfileBuilderProps) {
  const [name, setName] = useState(profile?.name ?? '');
  const [description, setDescription] = useState(profile?.description ?? '');
  const [thresholds, setThresholds] = useState<Record<string, ThresholdRange>>(() => {
    if (profile?.thresholds) return { ...profile.thresholds };
    // Default thresholds (Lipinski-like)
    return {
      mw: { min: 150, max: 500 },
      logp: { min: -2, max: 5 },
      tpsa: { min: 0, max: 140 },
      hbd: { min: 0, max: 5 },
      hba: { min: 0, max: 10 },
      rotatable_bonds: { min: 0, max: 10 },
      aromatic_rings: { min: 0, max: 4 },
      fsp3: { min: 0, max: 1 },
    };
  });
  const [weights, setWeights] = useState<Record<string, number>>(() => {
    if (profile?.weights) return { ...profile.weights };
    const defaults: Record<string, number> = {};
    THRESHOLD_CONFIGS.forEach((c) => {
      defaults[c.key] = 1.0;
    });
    return defaults;
  });

  // Reset form when profile prop changes
  useEffect(() => {
    if (profile) {
      setName(profile.name);
      setDescription(profile.description ?? '');
      if (profile.thresholds) setThresholds({ ...profile.thresholds });
      if (profile.weights) setWeights({ ...profile.weights });
    }
  }, [profile]);

  const updateThreshold = useCallback(
    (key: string, field: 'min' | 'max', value: number) => {
      setThresholds((prev) => ({
        ...prev,
        [key]: { ...prev[key], [field]: value },
      }));
    },
    []
  );

  const updateWeight = useCallback((key: string, value: number) => {
    setWeights((prev) => ({ ...prev, [key]: value }));
  }, []);

  const handleSave = async () => {
    if (!name.trim()) return;
    await onSave({
      name: name.trim(),
      description: description.trim() || undefined,
      thresholds,
      weights,
    });
  };

  const handleReset = () => {
    setName(profile?.name ?? '');
    setDescription(profile?.description ?? '');
    if (profile?.thresholds) setThresholds({ ...profile.thresholds });
    if (profile?.weights) {
      setWeights({ ...profile.weights });
    }
  };

  const handleExportClick = async () => {
    if (!onExport || !profile?.id) return;
    const data = await onExport(profile.id);
    const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `profile_${name.replace(/\s+/g, '_').toLowerCase()}.json`;
    a.click();
    URL.revokeObjectURL(url);
  };

  const handleImportClick = () => {
    const input = document.createElement('input');
    input.type = 'file';
    input.accept = '.json';
    input.onchange = async (e) => {
      const file = (e.target as HTMLInputElement).files?.[0];
      if (!file || !onImport) return;
      const text = await file.text();
      try {
        const data = JSON.parse(text) as ScoringProfileExport;
        await onImport(data);
      } catch {
        // Invalid JSON
      }
    };
    input.click();
  };

  return (
    <motion.div
      initial={{ opacity: 0, y: 10 }}
      animate={{ opacity: 1, y: 0 }}
      className="card p-6 space-y-6"
    >
      {/* Header */}
      <div className="flex items-center justify-between">
        <h3 className="text-lg font-semibold text-[var(--color-text-primary)] font-display">
          {profile ? 'Edit Profile' : 'Create Profile'}
        </h3>
        <div className="flex items-center gap-2">
          {onExport && profile?.id && (
            <ClayButton size="sm" onClick={handleExportClick} leftIcon={<Download className="w-3.5 h-3.5" />}>
              Export
            </ClayButton>
          )}
          {onImport && (
            <ClayButton size="sm" onClick={handleImportClick} leftIcon={<Upload className="w-3.5 h-3.5" />}>
              Import
            </ClayButton>
          )}
        </div>
      </div>

      {/* Name and description */}
      <div className="grid grid-cols-1 gap-4">
        <div>
          <label className="block text-sm font-medium text-[var(--color-text-secondary)] mb-1.5">
            Profile Name
          </label>
          <input
            type="text"
            value={name}
            onChange={(e) => setName(e.target.value)}
            placeholder="e.g., My Custom Filter"
            className={cn(
              'w-full px-3 py-2 rounded-lg text-sm',
              'bg-[var(--color-surface-sunken)] border border-[var(--color-border)]',
              'text-[var(--color-text-primary)] placeholder:text-[var(--color-text-muted)]',
              'focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]/30 focus:border-[var(--color-primary)]'
            )}
          />
        </div>
        <div>
          <label className="block text-sm font-medium text-[var(--color-text-secondary)] mb-1.5">
            Description
          </label>
          <textarea
            value={description}
            onChange={(e) => setDescription(e.target.value)}
            placeholder="Optional description of this scoring profile..."
            rows={2}
            className={cn(
              'w-full px-3 py-2 rounded-lg text-sm resize-none',
              'bg-[var(--color-surface-sunken)] border border-[var(--color-border)]',
              'text-[var(--color-text-primary)] placeholder:text-[var(--color-text-muted)]',
              'focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]/30 focus:border-[var(--color-primary)]'
            )}
          />
        </div>
      </div>

      {/* Threshold sliders */}
      <div>
        <h4 className="text-sm font-medium text-[var(--color-text-primary)] mb-3">Thresholds</h4>
        <div className="space-y-4">
          {THRESHOLD_CONFIGS.map((config) => {
            const range = thresholds[config.key] ?? { min: config.min, max: config.max };
            const weight = weights[config.key] ?? 1.0;
            return (
              <div key={config.key} className="bg-[var(--color-surface-sunken)] rounded-xl p-3">
                <div className="flex items-center justify-between mb-2">
                  <span className="text-sm font-medium text-[var(--color-text-primary)]">{config.label}</span>
                  <div className="flex items-center gap-3">
                    <span className="text-xs text-[var(--color-text-muted)]">
                      {range.min ?? config.min}{config.unit} - {range.max ?? config.max}{config.unit}
                    </span>
                    <div className="flex items-center gap-1">
                      <span className="text-[10px] text-[var(--color-text-muted)] uppercase">W:</span>
                      <input
                        type="number"
                        min={0}
                        max={2}
                        step={0.1}
                        value={weight}
                        onChange={(e) => updateWeight(config.key, parseFloat(e.target.value) || 0)}
                        className={cn(
                          'w-14 px-1.5 py-0.5 rounded text-xs text-center',
                          'bg-[var(--color-surface-elevated)] border border-[var(--color-border)]',
                          'text-[var(--color-text-primary)]',
                          'focus:outline-none focus:ring-1 focus:ring-[var(--color-primary)]/30'
                        )}
                      />
                    </div>
                  </div>
                </div>
                <div className="flex items-center gap-3">
                  <span className="text-[10px] text-[var(--color-text-muted)] w-8 text-right">
                    {config.min}
                  </span>
                  <div className="flex-1 relative">
                    {/* Min slider */}
                    <input
                      type="range"
                      min={config.min}
                      max={config.max}
                      step={config.step}
                      value={range.min ?? config.min}
                      onChange={(e) => updateThreshold(config.key, 'min', parseFloat(e.target.value))}
                      className="w-full h-1.5 appearance-none bg-[var(--color-border)] rounded-full cursor-pointer accent-[var(--color-primary)]"
                    />
                    {/* Max slider */}
                    <input
                      type="range"
                      min={config.min}
                      max={config.max}
                      step={config.step}
                      value={range.max ?? config.max}
                      onChange={(e) => updateThreshold(config.key, 'max', parseFloat(e.target.value))}
                      className="w-full h-1.5 appearance-none bg-transparent rounded-full cursor-pointer accent-[var(--color-accent)] -mt-1.5"
                    />
                  </div>
                  <span className="text-[10px] text-[var(--color-text-muted)] w-8">
                    {config.max}
                  </span>
                </div>
              </div>
            );
          })}
        </div>
      </div>

      {/* Actions */}
      <div className="flex items-center justify-between pt-2 border-t border-[var(--color-border)]">
        <ClayButton size="sm" onClick={handleReset} leftIcon={<RotateCcw className="w-3.5 h-3.5" />}>
          Reset
        </ClayButton>
        <div className="flex items-center gap-3">
          <ClayButton onClick={onCancel} leftIcon={<X className="w-4 h-4" />}>
            Cancel
          </ClayButton>
          <ClayButton
            variant="primary"
            onClick={handleSave}
            disabled={!name.trim() || isSaving}
            loading={isSaving}
            leftIcon={<Save className="w-4 h-4" />}
          >
            {profile ? 'Update Profile' : 'Create Profile'}
          </ClayButton>
        </div>
      </div>
    </motion.div>
  );
}
