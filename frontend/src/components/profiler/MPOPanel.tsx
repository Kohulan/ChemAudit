import { useState, useCallback } from 'react';
import { Settings2 } from 'lucide-react';
import { ClayCard } from '../ui/ClayCard';
import { ClayButton } from '../ui/ClayButton';
import { Badge } from '../ui/Badge';
import { DesirabilityEditorModal } from './DesirabilityEditorModal';
import { useMPOProfiles, type MPOProfile } from '../../hooks/useMPOProfiles';
import type { CNSMPOResult, CustomMPOResult, MPOProperty } from '../../types/profiler';

interface MPOPanelProps {
  smiles: string;
  cnsMPO: CNSMPOResult;
  computeCustomMPO: (smiles: string, profile: MPOProperty[]) => Promise<CustomMPOResult | null>;
}

/**
 * Custom Multi-Parameter Optimization (MPO) panel.
 *
 * Per D-11: preset-first approach. Defaults to CNS MPO.
 * Allows selecting from available presets or user profiles.
 * "Configure MPO" button (ClayButton accent) opens DesirabilityEditorModal for editing.
 * On selecting a custom profile: calls computeCustomMPO and displays result.
 */
export function MPOPanel({ smiles, cnsMPO, computeCustomMPO }: MPOPanelProps) {
  const { profiles, saveProfile, isPreset } = useMPOProfiles();

  // Default to CNS MPO preset
  const [selectedProfileId, setSelectedProfileId] = useState<string>('cns-mpo');
  const [customResult, setCustomResult] = useState<CustomMPOResult | null>(null);
  const [isComputing, setIsComputing] = useState(false);
  const [isEditorOpen, setIsEditorOpen] = useState(false);

  const selectedProfile = profiles.find((p) => p.id === selectedProfileId) ?? profiles[0];

  // Determine displayed score
  const isCNSMPO = selectedProfileId === 'cns-mpo';
  const displayScore = isCNSMPO ? cnsMPO.score : (customResult?.score ?? null);
  const displayMax = isCNSMPO ? cnsMPO.max_score : (customResult?.max_score ?? null);
  const displayNormalized = isCNSMPO
    ? cnsMPO.max_score > 0 ? cnsMPO.score / cnsMPO.max_score : 0
    : (customResult?.normalized ?? null);

  // CNS MPO components or custom components
  const components = isCNSMPO
    ? Object.entries(cnsMPO.components).map(([property, desirability]) => ({
        property,
        raw_value: 0, // not returned by CNS endpoint
        desirability: desirability as number,
        weight: 1,
      }))
    : (customResult?.components ?? []);

  const handleProfileChange = useCallback(
    async (profileId: string) => {
      setSelectedProfileId(profileId);
      if (profileId === 'cns-mpo') {
        setCustomResult(null);
        return;
      }
      const profile = profiles.find((p) => p.id === profileId);
      if (!profile) return;

      setIsComputing(true);
      try {
        const result = await computeCustomMPO(smiles, profile.properties);
        setCustomResult(result);
      } finally {
        setIsComputing(false);
      }
    },
    [profiles, smiles, computeCustomMPO]
  );

  const handleSaveProfile = useCallback(
    async (saved: MPOProfile) => {
      saveProfile(saved);
      // Auto-switch to the saved profile and compute
      setSelectedProfileId(saved.id);
      setIsComputing(true);
      try {
        const result = await computeCustomMPO(smiles, saved.properties);
        setCustomResult(result);
      } finally {
        setIsComputing(false);
      }
    },
    [saveProfile, smiles, computeCustomMPO]
  );

  // Score classification for CNS MPO (max 4)
  const getScoreVariant = (normalized: number): 'success' | 'warning' | 'error' => {
    if (normalized >= 0.75) return 'success';
    if (normalized >= 0.5) return 'warning';
    return 'error';
  };

  const newProfileTemplate: MPOProfile = {
    id: String(Date.now()),
    name: 'My MPO Profile',
    properties: selectedProfile?.properties ?? [],
    createdAt: Date.now(),
  };

  return (
    <>
      <ClayCard size="md">
        <div className="flex items-start justify-between gap-3 mb-4">
          <div>
            <h3 className="text-lg font-semibold text-text-primary font-display">
              Custom Multi-Parameter Optimization (MPO)
            </h3>
            <p className="text-xs text-text-muted mt-0.5">
              Composite desirability score across multiple molecular properties.
            </p>
          </div>
          <ClayButton
            variant="accent"
            size="sm"
            onClick={() => setIsEditorOpen(true)}
            leftIcon={<Settings2 className="w-3.5 h-3.5" />}
          >
            Configure MPO
          </ClayButton>
        </div>

        {/* Preset selector */}
        <div className="mb-4">
          <label className="block text-xs font-medium text-text-secondary mb-1">
            Active profile
          </label>
          <select
            value={selectedProfileId}
            onChange={(e) => handleProfileChange(e.target.value)}
            className="text-sm bg-surface-elevated border border-[var(--color-border)] rounded-xl px-3 py-2 text-text-primary w-full max-w-xs"
          >
            {profiles.map((p) => (
              <option key={p.id} value={p.id}>
                {p.name}{isPreset(p.id) ? ' (preset)' : ''}
              </option>
            ))}
          </select>
        </div>

        {/* Score display */}
        {isComputing ? (
          <div className="flex items-center gap-2 text-sm text-text-muted py-2">
            <span className="inline-block w-4 h-4 border-2 border-[var(--color-primary)] border-t-transparent rounded-full animate-spin" />
            Computing MPO score...
          </div>
        ) : displayScore !== null && displayMax !== null && displayNormalized !== null ? (
          <div className="space-y-3">
            {/* Main score */}
            <div className="flex items-baseline gap-3">
              <span className="text-3xl font-semibold tabular-nums text-text-primary font-display">
                {displayScore.toFixed(2)}
              </span>
              <span className="text-sm text-text-muted">/ {displayMax.toFixed(2)}</span>
              <Badge
                variant={getScoreVariant(displayNormalized)}
                size="sm"
              >
                {(displayNormalized * 100).toFixed(0)}%
              </Badge>
            </div>

            {/* Component breakdown */}
            {components.length > 0 && (
              <div className="space-y-1.5">
                <p className="text-xs font-semibold text-text-muted uppercase tracking-wide">
                  Component breakdown
                </p>
                {components.map((c) => (
                  <div key={c.property} className="flex items-center gap-2">
                    <span className="text-xs text-text-secondary w-16 shrink-0">{c.property}</span>
                    {/* Desirability bar */}
                    <div className="flex-1 h-1.5 bg-surface-sunken rounded-full overflow-hidden">
                      <div
                        className="h-full bg-[var(--color-primary)] rounded-full transition-all duration-500"
                        style={{ width: `${Math.max(0, Math.min(1, c.desirability)) * 100}%` }}
                      />
                    </div>
                    <span className="text-xs tabular-nums text-text-muted w-8 text-right">
                      {c.desirability.toFixed(2)}
                    </span>
                  </div>
                ))}
              </div>
            )}
          </div>
        ) : (
          <p className="text-sm text-text-muted">Select a profile to compute MPO score.</p>
        )}
      </ClayCard>

      {/* Desirability editor modal */}
      <DesirabilityEditorModal
        isOpen={isEditorOpen}
        onClose={() => setIsEditorOpen(false)}
        profile={
          isPreset(selectedProfileId)
            ? (selectedProfile ?? newProfileTemplate)
            : (selectedProfile ?? newProfileTemplate)
        }
        onSave={handleSaveProfile}
      />
    </>
  );
}
