import { useState, useEffect } from 'react';
import { ChevronDown, ChevronRight, Beaker, Check } from 'lucide-react';
import { profilesApi } from '../../services/api';
import { ClayButton } from '../ui/ClayButton';
import { cn } from '../../lib/utils';
import type { ScoringProfile, ThresholdRange } from '../../types/workflow';

interface ProfileSidebarProps {
  selectedProfileId: number | null;
  onProfileChange: (profileId: number | null) => void;
  disabled?: boolean;
}

function formatRange(range: ThresholdRange | undefined): string {
  if (!range) return '--';
  if (range.min !== undefined && range.max !== undefined) return `${range.min}-${range.max}`;
  if (range.max !== undefined) return `\u2264${range.max}`;
  if (range.min !== undefined) return `\u2265${range.min}`;
  return '--';
}

export function ProfileSidebar({ selectedProfileId, onProfileChange, disabled }: ProfileSidebarProps) {
  const [isExpanded, setIsExpanded] = useState(selectedProfileId != null);
  const [profiles, setProfiles] = useState<ScoringProfile[]>([]);
  const [isLoading, setIsLoading] = useState(false);

  useEffect(() => {
    if (isExpanded && profiles.length === 0) {
      setIsLoading(true);
      profilesApi.getProfiles()
        .then(data => setProfiles(data))
        .catch(() => {})
        .finally(() => setIsLoading(false));
    }
  }, [isExpanded, profiles.length]);

  const handleToggle = () => {
    const next = !isExpanded;
    setIsExpanded(next);
    if (!next) onProfileChange(null);
  };

  const selected = profiles.find(p => p.id === selectedProfileId);

  return (
    <div className="bg-[var(--color-surface-elevated)] rounded-xl border border-[var(--color-border)] overflow-hidden">
      <button
        onClick={handleToggle}
        disabled={disabled}
        className="w-full p-4 flex items-center justify-between bg-[var(--color-surface-sunken)]/50 hover:bg-[var(--color-surface-sunken)] transition-colors"
      >
        <div className="flex items-center gap-2">
          <Beaker className="w-4 h-4 text-[var(--color-primary)]" />
          <span className="text-sm font-medium text-[var(--color-text-primary)]">Scoring Profile</span>
          {selected && (
            <span className="text-xs px-2 py-0.5 rounded-full bg-[var(--color-primary)]/10 text-[var(--color-primary)] font-medium">
              {selected.name}
            </span>
          )}
        </div>
        {isExpanded ? (
          <ChevronDown className="w-4 h-4 text-[var(--color-text-muted)]" />
        ) : (
          <ChevronRight className="w-4 h-4 text-[var(--color-text-muted)]" />
        )}
      </button>

      {isExpanded && (
        <div className="p-4 space-y-3">
          <p className="text-xs text-[var(--color-text-muted)]">
            Select a profile to score each molecule against property thresholds.
          </p>

          {isLoading ? (
            <div className="flex justify-center py-4">
              <div className="animate-spin w-6 h-6 border-2 border-[var(--color-primary)] border-t-transparent rounded-full" />
            </div>
          ) : (
            <div className="grid grid-cols-1 gap-2">
              {profiles.map(profile => {
                const isSelected = profile.id === selectedProfileId;
                const keys = Object.keys(profile.thresholds).slice(0, 4);
                return (
                  <button
                    key={profile.id}
                    onClick={() => onProfileChange(isSelected ? null : profile.id)}
                    disabled={disabled}
                    className={cn(
                      'w-full text-left p-3 rounded-lg border transition-all',
                      isSelected
                        ? 'border-[var(--color-primary)] bg-[var(--color-primary)]/5 ring-1 ring-[var(--color-primary)]'
                        : 'border-[var(--color-border)] hover:border-[var(--color-text-muted)] bg-[var(--color-surface-elevated)]'
                    )}
                  >
                    <div className="flex items-center justify-between">
                      <span className="text-sm font-medium text-[var(--color-text-primary)]">{profile.name}</span>
                      {isSelected && <Check className="w-4 h-4 text-[var(--color-primary)]" />}
                    </div>
                    <div className="flex flex-wrap gap-1.5 mt-1.5">
                      {keys.map(key => (
                        <span key={key} className="text-[10px] px-1.5 py-0.5 rounded bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)]">
                          {key}: {formatRange(profile.thresholds[key])}
                        </span>
                      ))}
                    </div>
                  </button>
                );
              })}
            </div>
          )}

          {selectedProfileId != null && (
            <ClayButton
              variant="ghost"
              size="sm"
              onClick={() => onProfileChange(null)}
              className="w-full"
            >
              Clear selection
            </ClayButton>
          )}
        </div>
      )}
    </div>
  );
}
