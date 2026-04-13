import { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown, Trash2, Save } from 'lucide-react';
import { cn } from '../../lib/utils';
import { ClayCard } from '../ui/ClayCard';
import { ClayButton } from '../ui/ClayButton';
import { PresetSelector } from './PresetSelector';
import { StepToggleGrid } from './StepToggleGrid';
import type { UseQSARPipelineConfigReturn } from '../../hooks/useQSARPipelineConfig';
import type { QSARReadyConfig } from '../../types/qsar_ready';

// =============================================================================
// Props
// =============================================================================

/**
 * All props from useQSARPipelineConfig return value — spread directly from the hook.
 */
type PipelineConfigPanelProps = UseQSARPipelineConfigReturn;

// =============================================================================
// Component
// =============================================================================

/**
 * Pipeline Configuration Panel — wraps PresetSelector + StepToggleGrid.
 *
 * Shows 3 built-in preset buttons and a 10-step toggle grid in a collapsible
 * ClayCard. Custom profiles can be saved (window.prompt for name), loaded, and
 * deleted (window.confirm). Built-in presets have no save/delete actions.
 *
 * Collapse animation: AnimatePresence + motion.div height auto→0, 250ms ease-in.
 * Per Phase 10 UI-SPEC copywriting and interaction contracts.
 */
export function PipelineConfigPanel({
  config,
  activePreset,
  modifiedFrom,
  profiles,
  selectPreset,
  updateConfig,
  saveProfile,
  deleteProfile,
  loadProfile,
  isPreset,
}: PipelineConfigPanelProps) {
  const [collapsed, setCollapsed] = useState(false);

  // Handle step boolean toggles
  const handleToggleStep = (key: string, enabled: boolean) => {
    updateConfig({ [key]: enabled } as Partial<QSARReadyConfig>);
  };

  // Handle filter numeric/boolean threshold updates
  const handleUpdateConfig = (partial: Partial<QSARReadyConfig>) => {
    updateConfig(partial);
  };

  // Save current config as a named custom profile
  const handleSaveConfig = () => {
    const name = window.prompt('Enter a name for this configuration:');
    if (name && name.trim()) {
      saveProfile(name.trim());
    }
  };

  // Delete a custom preset with confirmation
  const handleDeleteProfile = (id: string, name: string) => {
    const confirmed = window.confirm(`Delete '${name}'? This cannot be undone.`);
    if (confirmed) {
      deleteProfile(id);
    }
  };

  return (
    <ClayCard variant="default" size="md" className="p-0">
      {/* ── Header row ── */}
      <div className="flex items-center justify-between px-6 py-4">
        <div className="flex-1 min-w-0">
          <h2 className="text-lg font-semibold font-display text-[var(--color-text-primary)]">
            Pipeline Configuration
          </h2>
          <p className="text-sm text-[var(--color-text-secondary)] mt-0.5">
            Select a preset or customize individual steps
          </p>
        </div>

        {/* Right side: Save Config + collapse toggle */}
        <div className="flex items-center gap-2 shrink-0 ml-4">
          <ClayButton
            variant="default"
            size="sm"
            leftIcon={<Save className="w-3.5 h-3.5" />}
            onClick={handleSaveConfig}
          >
            Save Config
          </ClayButton>

          <button
            type="button"
            onClick={() => setCollapsed((v) => !v)}
            className="text-[var(--color-text-muted)] hover:text-[var(--color-text-primary)] transition-colors p-1"
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
      </div>

      {/* ── Collapsible body ── */}
      <AnimatePresence initial={false}>
        {!collapsed && (
          <motion.div
            key="panel-body"
            initial={{ height: 0, opacity: 0 }}
            animate={{ height: 'auto', opacity: 1 }}
            exit={{ height: 0, opacity: 0 }}
            transition={{ duration: 0.25, ease: 'easeIn' }}
            className="overflow-hidden"
          >
            <div className="px-6 pb-6 space-y-5">
              {/* Divider */}
              <div className="border-t border-[var(--color-border)]" />

              {/* Preset selector */}
              <PresetSelector
                activePreset={activePreset}
                modifiedFrom={modifiedFrom}
                onSelectPreset={selectPreset}
              />

              {/* Step toggle grid */}
              <StepToggleGrid
                config={config}
                onToggleStep={handleToggleStep}
                onUpdateConfig={handleUpdateConfig}
              />

              {/* Custom saved profiles list */}
              {profiles.length > 0 && (
                <div className="space-y-2">
                  <p className="text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wide">
                    Saved Configurations
                  </p>
                  <div className="flex flex-wrap gap-2">
                    {profiles.map((profile) => {
                      const isActiveProfile = activePreset === profile.id;
                      return (
                        <div
                          key={profile.id}
                          className={cn(
                            'flex items-center gap-1.5 px-3 py-1.5 rounded-xl border text-xs',
                            'bg-[var(--color-surface-elevated)] border-[var(--color-border)]',
                            isActiveProfile && 'border-[var(--color-primary)] text-[var(--color-primary)]',
                          )}
                        >
                          <button
                            type="button"
                            onClick={() => loadProfile(profile.id)}
                            className={cn(
                              'font-medium hover:text-[var(--color-primary)] transition-colors',
                              isActiveProfile
                                ? 'text-[var(--color-primary)]'
                                : 'text-[var(--color-text-primary)]',
                            )}
                          >
                            {profile.name}
                          </button>
                          {!isPreset(profile.id) && (
                            <button
                              type="button"
                              onClick={() => handleDeleteProfile(profile.id, profile.name)}
                              className="text-[var(--color-text-muted)] hover:text-red-500 transition-colors"
                              aria-label={`Delete preset '${profile.name}'`}
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
          </motion.div>
        )}
      </AnimatePresence>
    </ClayCard>
  );
}
