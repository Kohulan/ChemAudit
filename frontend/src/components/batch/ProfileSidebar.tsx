import { useState, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown, Beaker, Check, FlaskConical, Dna, Pill, Microscope, Atom, TestTube, X } from 'lucide-react';
import { profilesApi } from '../../services/api';
import { ClayButton } from '../ui/ClayButton';
import { cn } from '../../lib/utils';
import type { ScoringProfile, ThresholdRange } from '../../types/workflow';

interface ProfileSidebarProps {
  selectedProfileId: number | null;
  onProfileChange: (profileId: number | null) => void;
  disabled?: boolean;
}

const KEY_META: Record<string, { label: string; unit?: string }> = {
  mw: { label: 'MW' },
  logp: { label: 'LogP' },
  hbd: { label: 'HBD' },
  hba: { label: 'HBA' },
  tpsa: { label: 'TPSA', unit: '\u00C5\u00B2' },
  rotatable_bonds: { label: 'RotB' },
  aromatic_rings: { label: 'ArRings' },
  fsp3: { label: 'Fsp3' },
};

function formatRange(range: ThresholdRange | undefined, unit?: string): string {
  if (!range) return '--';
  const u = unit ?? '';
  if (range.min !== undefined && range.max !== undefined) return `${range.min}\u2013${range.max}${u}`;
  if (range.max !== undefined) return `\u2264${range.max}${u}`;
  if (range.min !== undefined) return `\u2265${range.min}${u}`;
  return '--';
}

const PROFILE_ICONS = [FlaskConical, Dna, Pill, Microscope, Atom, TestTube, Beaker, FlaskConical];

const ACCENT = [
  { bg: 'bg-rose-500/12', text: 'text-rose-500' },
  { bg: 'bg-violet-500/12', text: 'text-violet-500' },
  { bg: 'bg-emerald-500/12', text: 'text-emerald-500' },
  { bg: 'bg-sky-500/12', text: 'text-sky-500' },
  { bg: 'bg-amber-500/12', text: 'text-amber-600' },
  { bg: 'bg-pink-500/12', text: 'text-pink-500' },
  { bg: 'bg-cyan-500/12', text: 'text-cyan-500' },
  { bg: 'bg-lime-500/12', text: 'text-lime-600' },
];

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
    <div className="mb-2">
      {/* Toggle row — lightweight, no box */}
      <button
        onClick={handleToggle}
        disabled={disabled}
        className={cn(
          'w-full group cursor-pointer flex items-center gap-2.5 py-2 transition-colors',
          disabled && 'opacity-50 pointer-events-none'
        )}
      >
        <div className={cn(
          'w-7 h-7 rounded-lg flex items-center justify-center transition-all duration-200',
          isExpanded
            ? 'bg-gradient-to-br from-[var(--color-primary)]/20 to-[var(--color-accent)]/10'
            : 'bg-[var(--color-surface-sunken)] group-hover:bg-[var(--color-primary)]/10'
        )}>
          <Beaker className={cn(
            'w-3.5 h-3.5 transition-colors',
            isExpanded ? 'text-[var(--color-primary)]' : 'text-[var(--color-text-muted)] group-hover:text-[var(--color-primary)]'
          )} />
        </div>
        <span className="text-sm font-semibold text-[var(--color-text-primary)] font-display">
          Scoring Profile
        </span>
        <AnimatePresence>
          {selected && (
            <motion.span
              initial={{ opacity: 0, scale: 0.85 }}
              animate={{ opacity: 1, scale: 1 }}
              exit={{ opacity: 0, scale: 0.85 }}
              className="inline-flex items-center gap-1 text-xs px-2 py-0.5 rounded-full bg-[var(--color-primary)]/10 text-[var(--color-primary)] font-medium"
            >
              <Check className="w-3 h-3" />
              {selected.name}
            </motion.span>
          )}
        </AnimatePresence>
        <motion.div
          animate={{ rotate: isExpanded ? 180 : 0 }}
          transition={{ duration: 0.25 }}
          className="ml-auto text-[var(--color-text-muted)]"
        >
          <ChevronDown className="w-4 h-4" />
        </motion.div>
      </button>

      {/* Cards */}
      <AnimatePresence>
        {isExpanded && (
          <motion.div
            initial={{ height: 0, opacity: 0 }}
            animate={{ height: 'auto', opacity: 1 }}
            exit={{ height: 0, opacity: 0 }}
            transition={{ duration: 0.25, ease: [0.4, 0, 0.2, 1] }}
            className="overflow-hidden"
          >
            {isLoading ? (
              <div className="flex items-center gap-2 py-4 justify-center">
                <div className="w-5 h-5 rounded-full border-2 border-transparent border-t-[var(--color-primary)] animate-spin" />
                <span className="text-sm text-[var(--color-text-muted)]">Loading profiles...</span>
              </div>
            ) : (
              <motion.div
                className="grid grid-cols-1 sm:grid-cols-2 gap-2"
                initial="initial"
                animate="animate"
                variants={{
                  initial: {},
                  animate: { transition: { staggerChildren: 0.04 } },
                }}
              >
                {profiles.map((profile, index) => {
                  const isSelected = profile.id === selectedProfileId;
                  const keys = Object.keys(profile.thresholds).slice(0, 4);
                  const Icon = PROFILE_ICONS[index % PROFILE_ICONS.length];
                  const accent = ACCENT[index % ACCENT.length];

                  return (
                    <motion.button
                      key={profile.id}
                      variants={{
                        initial: { opacity: 0, y: 8 },
                        animate: { opacity: 1, y: 0, transition: { type: 'spring', stiffness: 400, damping: 28 } },
                      }}
                      onClick={() => onProfileChange(isSelected ? null : profile.id)}
                      disabled={disabled}
                      whileHover={{ y: -1, transition: { duration: 0.15 } }}
                      whileTap={{ scale: 0.98 }}
                      className={cn(
                        'relative w-full text-left rounded-xl overflow-hidden cursor-pointer',
                        'border transition-all duration-200',
                        isSelected
                          ? 'border-[var(--color-primary)] bg-[var(--color-primary)]/[0.03] shadow-[0_0_14px_var(--glow-primary)]'
                          : 'border-[var(--color-border)] bg-[var(--color-surface-elevated)] hover:border-[var(--color-border-strong)] hover:shadow-[var(--shadow-sm)]'
                      )}
                    >
                      {/* Active top bar */}
                      {isSelected && (
                        <motion.div
                          layoutId="profile-bar"
                          className="absolute top-0 inset-x-0 h-[2px] bg-gradient-to-r from-[var(--color-primary)] to-[var(--color-accent)]"
                          transition={{ type: 'spring', stiffness: 500, damping: 30 }}
                        />
                      )}

                      <div className="relative px-3 py-2.5">
                        {/* Row 1: icon + name + description + check */}
                        <div className="flex items-center gap-2 mb-2">
                          <div className={cn(
                            'w-7 h-7 rounded-lg flex-shrink-0 flex items-center justify-center',
                            isSelected ? 'bg-[var(--color-primary)]/15' : accent.bg
                          )}>
                            <Icon className={cn(
                              'w-3.5 h-3.5',
                              isSelected ? 'text-[var(--color-primary)]' : accent.text
                            )} />
                          </div>
                          <div className="flex-1 min-w-0">
                            <span className="text-sm font-semibold text-[var(--color-text-primary)] font-display truncate block leading-tight">
                              {profile.name}
                            </span>
                            {profile.description && (
                              <span className="text-xs text-[var(--color-text-muted)] truncate block leading-tight">
                                {profile.description}
                              </span>
                            )}
                          </div>
                          {isSelected && (
                            <motion.span
                              initial={{ scale: 0 }}
                              animate={{ scale: 1 }}
                              className="flex-shrink-0 w-5 h-5 rounded-full bg-[var(--color-primary)] flex items-center justify-center shadow-[0_0_6px_var(--glow-primary)]"
                            >
                              <Check className="w-3 h-3 text-white" />
                            </motion.span>
                          )}
                        </div>

                        {/* Row 2: threshold pills in 2x2 grid */}
                        <div className="grid grid-cols-2 gap-x-2 gap-y-1">
                          {keys.map((key) => {
                            const meta = KEY_META[key];
                            return (
                              <div
                                key={key}
                                className={cn(
                                  'flex items-center justify-between rounded-md px-2 py-0.5',
                                  isSelected
                                    ? 'bg-[var(--color-primary)]/[0.06]'
                                    : 'bg-[var(--color-surface-sunken)]/60'
                                )}
                              >
                                <span className="text-xs text-[var(--color-text-muted)]">
                                  {meta?.label ?? key}
                                </span>
                                <span className={cn(
                                  'text-xs font-semibold font-mono',
                                  isSelected ? 'text-[var(--color-primary)]' : 'text-[var(--color-text-secondary)]'
                                )}>
                                  {formatRange(profile.thresholds[key], meta?.unit)}
                                </span>
                              </div>
                            );
                          })}
                        </div>
                      </div>
                    </motion.button>
                  );
                })}
              </motion.div>
            )}

            {/* Clear selection */}
            <AnimatePresence>
              {selectedProfileId != null && (
                <motion.div
                  initial={{ opacity: 0 }}
                  animate={{ opacity: 1 }}
                  exit={{ opacity: 0 }}
                  className="mt-2"
                >
                  <ClayButton
                    variant="ghost"
                    size="sm"
                    onClick={() => onProfileChange(null)}
                    className="w-full"
                    leftIcon={<X className="w-3.5 h-3.5" />}
                  >
                    Clear selection
                  </ClayButton>
                </motion.div>
              )}
            </AnimatePresence>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}
