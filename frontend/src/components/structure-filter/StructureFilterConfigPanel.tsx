import { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown } from 'lucide-react';
import { cn } from '../../lib/utils';
import { ClayCard } from '../ui/ClayCard';
import { StructureFilterPresetSelector } from './StructureFilterPresetSelector';
import type { FilterConfig, StructureFilterProfile } from '../../types/structure_filter';

// =============================================================================
// Props
// =============================================================================

interface StructureFilterConfigPanelProps {
  /** Current filter configuration. */
  config: FilterConfig;
  /** ID of the active preset (null = custom). */
  activePreset: string | null;
  /** Preset this config was modified from (null = not modified). */
  modifiedFrom: string | null;
  /** Custom profiles saved to localStorage. */
  profiles: StructureFilterProfile[];
  /** Called on any config field change. */
  onUpdateConfig: (partial: Partial<FilterConfig>) => void;
  /** Called when a preset is selected. */
  onSelectPreset: (id: string) => void;
  /** Called when a custom profile is saved. */
  onSaveProfile: (name: string) => void;
  /** Called when a custom profile is deleted. */
  onDeleteProfile: (id: string) => void;
}

// =============================================================================
// Helpers
// =============================================================================

interface LabeledInputProps {
  id: string;
  label: string;
  value: number;
  onChange: (value: number) => void;
  step?: number;
  min?: number;
  max?: number;
  unit?: string;
}

/** Labeled number input with associated label element (UI-SPEC accessibility). */
function LabeledInput({
  id,
  label,
  value,
  onChange,
  step = 1,
  min,
  max,
  unit,
}: LabeledInputProps) {
  return (
    <div className="flex flex-col gap-1">
      <label
        htmlFor={id}
        className="text-xs font-medium text-[var(--color-text-secondary)]"
      >
        {label}
        {unit && <span className="text-[var(--color-text-muted)] ml-1">({unit})</span>}
      </label>
      <input
        id={id}
        type="number"
        value={value}
        step={step}
        min={min}
        max={max}
        onChange={(e) => onChange(parseFloat(e.target.value))}
        className="w-full border border-[var(--color-border)] bg-[var(--color-surface-elevated)] text-[var(--color-text-primary)] rounded-lg px-3 py-1.5 text-sm focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)] focus:border-transparent"
      />
    </div>
  );
}

interface CheckboxRowProps {
  id: string;
  label: string;
  checked: boolean;
  onChange: (checked: boolean) => void;
  description?: string;
}

/** Checkbox toggle row with associated label element. */
function CheckboxRow({ id, label, checked, onChange, description }: CheckboxRowProps) {
  return (
    <div className="flex items-start gap-2">
      <input
        id={id}
        type="checkbox"
        checked={checked}
        onChange={(e) => onChange(e.target.checked)}
        className="mt-0.5 w-4 h-4 accent-[var(--color-primary)] rounded cursor-pointer"
      />
      <div>
        <label
          htmlFor={id}
          className="text-xs font-medium text-[var(--color-text-primary)] cursor-pointer"
        >
          {label}
        </label>
        {description && (
          <p className="text-xs text-[var(--color-text-muted)] mt-0.5">{description}</p>
        )}
      </div>
    </div>
  );
}

// =============================================================================
// Component
// =============================================================================

/**
 * Configuration panel for the Structure Filter page.
 *
 * Shows preset selector, property threshold inputs, alert catalog toggles,
 * and novelty stage controls in a collapsible ClayCard.
 *
 * All inputs have associated label elements (htmlFor) for accessibility.
 * Each onChange triggers onUpdateConfig(partial), which updates the hook state
 * and triggers modified-from detection in useStructureFilterConfig.
 *
 * Per Phase 11 UI-SPEC D-12 config panel spec.
 */
export function StructureFilterConfigPanel({
  config,
  activePreset,
  modifiedFrom,
  profiles,
  onUpdateConfig,
  onSelectPreset,
  onSaveProfile,
  onDeleteProfile,
}: StructureFilterConfigPanelProps) {
  const [collapsed, setCollapsed] = useState(false);
  const [ringLimitEnabled, setRingLimitEnabled] = useState(config.max_rings !== null);

  const handleRingLimitToggle = (enabled: boolean) => {
    setRingLimitEnabled(enabled);
    onUpdateConfig({ max_rings: enabled ? (config.max_rings ?? 6) : null });
  };

  return (
    <ClayCard variant="default" size="md" className="p-0">
      {/* ── Header row ── */}
      <div className="flex items-center justify-between px-6 py-4">
        <div className="flex-1 min-w-0">
          <h2 className="text-base font-semibold font-display text-[var(--color-text-primary)]">
            Configuration
          </h2>
          <p className="text-sm text-[var(--color-text-secondary)] mt-0.5">
            Select a preset or customize thresholds and alerts
          </p>
        </div>

        {/* Collapse toggle */}
        <button
          type="button"
          onClick={() => setCollapsed((v) => !v)}
          className="text-[var(--color-text-muted)] hover:text-[var(--color-text-primary)] transition-colors p-1 ml-4"
          aria-label={collapsed ? 'Expand configuration panel' : 'Collapse configuration panel'}
        >
          <ChevronDown
            className={cn(
              'w-5 h-5 transition-transform duration-250',
              collapsed && 'rotate-180',
            )}
          />
        </button>
      </div>

      {/* ── Collapsible body ── */}
      <AnimatePresence initial={false}>
        {!collapsed && (
          <motion.div
            key="config-panel-body"
            initial={{ height: 0, opacity: 0 }}
            animate={{ height: 'auto', opacity: 1 }}
            exit={{ height: 0, opacity: 0 }}
            transition={{ duration: 0.25, ease: 'easeIn' }}
            className="overflow-hidden"
          >
            <div className="px-6 pb-6 space-y-5">
              <div className="border-t border-[var(--color-border)]" />

              {/* ── Preset selector ── */}
              <StructureFilterPresetSelector
                activePreset={activePreset}
                modifiedFrom={modifiedFrom}
                profiles={profiles}
                onSelectPreset={onSelectPreset}
                onSaveProfile={onSaveProfile}
                onDeleteProfile={onDeleteProfile}
              />

              {/* ── Property Thresholds ── */}
              <div>
                <p className="text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wide mb-3">
                  Property Thresholds
                </p>
                <div className="grid grid-cols-2 gap-3">
                  {/* MW: min + max side by side */}
                  <LabeledInput
                    id="config-min-mw"
                    label="MW min"
                    value={config.min_mw}
                    onChange={(v) => onUpdateConfig({ min_mw: v })}
                    step={10}
                    min={0}
                    unit="Da"
                  />
                  <LabeledInput
                    id="config-max-mw"
                    label="MW max"
                    value={config.max_mw}
                    onChange={(v) => onUpdateConfig({ max_mw: v })}
                    step={10}
                    min={0}
                    unit="Da"
                  />

                  {/* LogP: min + max side by side */}
                  <LabeledInput
                    id="config-min-logp"
                    label="LogP min"
                    value={config.min_logp}
                    onChange={(v) => onUpdateConfig({ min_logp: v })}
                    step={0.5}
                  />
                  <LabeledInput
                    id="config-max-logp"
                    label="LogP max"
                    value={config.max_logp}
                    onChange={(v) => onUpdateConfig({ max_logp: v })}
                    step={0.5}
                  />

                  {/* TPSA: single full-width */}
                  <div className="col-span-2">
                    <LabeledInput
                      id="config-max-tpsa"
                      label="TPSA max"
                      value={config.max_tpsa}
                      onChange={(v) => onUpdateConfig({ max_tpsa: v })}
                      step={10}
                      min={0}
                      unit="Å²"
                    />
                  </div>

                  {/* Rotatable bonds */}
                  <LabeledInput
                    id="config-max-rot-bonds"
                    label="Rotatable bonds max"
                    value={config.max_rot_bonds}
                    onChange={(v) => onUpdateConfig({ max_rot_bonds: Math.round(v) })}
                    step={1}
                    min={0}
                  />

                  {/* SA Score */}
                  <LabeledInput
                    id="config-max-sa-score"
                    label="SA Score max"
                    value={config.max_sa_score}
                    onChange={(v) => onUpdateConfig({ max_sa_score: v })}
                    step={0.1}
                    min={1.0}
                    max={10.0}
                  />
                </div>

                {/* Ring limit (optional — checkbox to enable) */}
                <div className="mt-3 space-y-2">
                  <CheckboxRow
                    id="config-ring-limit-enabled"
                    label="Enable ring count limit"
                    checked={ringLimitEnabled}
                    onChange={handleRingLimitToggle}
                  />
                  {ringLimitEnabled && (
                    <div className="pl-6">
                      <LabeledInput
                        id="config-max-rings"
                        label="Rings max"
                        value={config.max_rings ?? 6}
                        onChange={(v) => onUpdateConfig({ max_rings: Math.round(v) })}
                        step={1}
                        min={1}
                      />
                    </div>
                  )}
                </div>
              </div>

              {/* ── Alert Catalogs ── */}
              <div>
                <p className="text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wide mb-3">
                  Alert Catalogs
                </p>
                <div className="space-y-2">
                  <CheckboxRow
                    id="config-use-pains"
                    label="PAINS"
                    checked={config.use_pains}
                    onChange={(checked) => onUpdateConfig({ use_pains: checked })}
                    description="Pan-assay interference compounds (480 patterns)"
                  />
                  <CheckboxRow
                    id="config-use-brenk"
                    label="Brenk"
                    checked={config.use_brenk}
                    onChange={(checked) => onUpdateConfig({ use_brenk: checked })}
                    description="Undesirable functional groups (105 patterns)"
                  />
                  <CheckboxRow
                    id="config-use-kazius"
                    label="Kazius"
                    checked={config.use_kazius}
                    onChange={(checked) => onUpdateConfig({ use_kazius: checked })}
                    description="Mutagenicity toxicophores (29 patterns)"
                  />
                  <CheckboxRow
                    id="config-use-nibr"
                    label="NIBR"
                    checked={config.use_nibr}
                    onChange={(checked) => onUpdateConfig({ use_nibr: checked })}
                    description="Novartis screening deck filters (444 patterns)"
                  />
                </div>
              </div>

              {/* ── Novelty Stage ── */}
              <div>
                <p className="text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wide mb-3">
                  Novelty Stage
                </p>
                <div className="space-y-3">
                  <CheckboxRow
                    id="config-enable-novelty"
                    label="Enable novelty check"
                    checked={config.enable_novelty}
                    onChange={(checked) => onUpdateConfig({ enable_novelty: checked })}
                  />

                  {config.enable_novelty ? (
                    <div className="pl-6 space-y-3">
                      <LabeledInput
                        id="config-novelty-threshold"
                        label="Tanimoto threshold (reject if >)"
                        value={config.novelty_threshold}
                        onChange={(v) => onUpdateConfig({ novelty_threshold: v })}
                        step={0.05}
                        min={0.5}
                        max={1.0}
                      />

                      {/* Reference set selector (informational) */}
                      <div className="space-y-1.5">
                        <p className="text-xs font-medium text-[var(--color-text-secondary)]">
                          Reference set
                        </p>
                        <label className="flex items-center gap-2 cursor-pointer">
                          <input
                            type="radio"
                            name="novelty-reference"
                            defaultChecked
                            className="accent-[var(--color-primary)]"
                          />
                          <span className="text-xs text-[var(--color-text-primary)]">
                            ChEMBL approved drugs (built-in)
                          </span>
                        </label>
                        <label className="flex items-center gap-2 cursor-pointer">
                          <input
                            type="radio"
                            name="novelty-reference"
                            className="accent-[var(--color-primary)]"
                          />
                          <span className="text-xs text-[var(--color-text-primary)]">
                            Upload reference SMILES (.txt)
                          </span>
                        </label>
                      </div>
                    </div>
                  ) : (
                    <p className="pl-6 text-xs text-[var(--color-text-muted)]">
                      Novelty check is disabled. Enable it in the config panel to compare against
                      ChEMBL approved drugs.
                    </p>
                  )}
                </div>
              </div>
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </ClayCard>
  );
}
