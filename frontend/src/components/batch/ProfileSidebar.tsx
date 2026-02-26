import { useState, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown, Beaker, Check, FlaskConical, Dna, Pill, Microscope, Atom, TestTube, X, SlidersHorizontal } from 'lucide-react';
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
  { bg: 'bg-rose-500/12', text: 'text-rose-500', border: 'border-rose-500/20' },
  { bg: 'bg-violet-500/12', text: 'text-violet-500', border: 'border-violet-500/20' },
  { bg: 'bg-emerald-500/12', text: 'text-emerald-500', border: 'border-emerald-500/20' },
  { bg: 'bg-sky-500/12', text: 'text-sky-500', border: 'border-sky-500/20' },
  { bg: 'bg-amber-500/12', text: 'text-amber-600', border: 'border-amber-500/20' },
  { bg: 'bg-pink-500/12', text: 'text-pink-500', border: 'border-pink-500/20' },
  { bg: 'bg-cyan-500/12', text: 'text-cyan-500', border: 'border-cyan-500/20' },
  { bg: 'bg-lime-500/12', text: 'text-lime-600', border: 'border-lime-500/20' },
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
    <div className="my-3">
      {/* Toggle — prominent card with gradient background */}
      <button
        onClick={handleToggle}
        disabled={disabled}
        className={cn(
          'w-full group cursor-pointer relative rounded-xl transition-all duration-300 overflow-hidden',
          'border px-4 py-3.5',
          disabled && 'opacity-50 pointer-events-none',
          isExpanded
            ? 'border-[var(--color-primary)]/30 shadow-[0_0_20px_var(--glow-soft)]'
            : selected
              ? 'border-[var(--color-primary)]/30 shadow-[var(--shadow-sm)]'
              : 'border-[var(--color-accent)]/25 hover:border-[var(--color-primary)]/40 hover:shadow-[0_0_20px_var(--glow-soft)]'
        )}
        style={{
          background: isExpanded
            ? 'var(--gradient-primary-subtle)'
            : undefined,
        }}
      >
        {/* Subtle gradient background when collapsed */}
        {!isExpanded && (
          <div
            className="absolute inset-0 opacity-60 group-hover:opacity-100 transition-opacity duration-300"
            style={{ background: 'var(--gradient-primary-subtle)' }}
          />
        )}

        {/* Gradient accent line at top when expanded */}
        <AnimatePresence>
          {isExpanded && (
            <motion.div
              initial={{ scaleX: 0, opacity: 0 }}
              animate={{ scaleX: 1, opacity: 1 }}
              exit={{ scaleX: 0, opacity: 0 }}
              transition={{ duration: 0.3, ease: [0.4, 0, 0.2, 1] }}
              className="absolute top-0 left-0 right-0 h-[2px] origin-left"
              style={{ background: 'var(--gradient-primary)' }}
            />
          )}
        </AnimatePresence>

        <div className="relative flex items-center gap-3">
          {/* Icon with gradient fill */}
          <div className="relative flex-shrink-0">
            <div className={cn(
              'w-10 h-10 rounded-xl flex items-center justify-center transition-all duration-300',
              isExpanded || selected
                ? 'shadow-[0_0_12px_var(--glow-primary)]'
                : 'group-hover:shadow-[0_0_10px_var(--glow-accent)]'
            )}
              style={{
                background: isExpanded || selected
                  ? 'var(--gradient-primary)'
                  : 'linear-gradient(135deg, rgba(var(--color-accent-rgb), 0.15) 0%, rgba(var(--color-primary-rgb), 0.1) 100%)',
              }}
            >
              <SlidersHorizontal className={cn(
                'w-4.5 h-4.5 transition-colors duration-200',
                isExpanded || selected
                  ? 'text-white'
                  : 'text-[var(--color-accent)] group-hover:text-[var(--color-primary)]'
              )} />
            </div>
            {/* Subtle pulse ring when collapsed and no selection */}
            {!isExpanded && !selected && (
              <span
                className="absolute -inset-1 rounded-xl animate-[ping_3s_ease-in-out_infinite] pointer-events-none"
                style={{ border: '1px solid rgba(var(--color-accent-rgb), 0.25)' }}
              />
            )}
          </div>

          {/* Title + hint */}
          <div className="flex-1 text-left min-w-0">
            <div className="flex items-center gap-2">
              <span className="text-sm font-semibold text-[var(--color-text-primary)] font-display">
                Scoring Profile
              </span>
              {!isExpanded && !selected && (
                <span
                  className="text-[10px] uppercase tracking-wider font-bold px-2 py-0.5 rounded-full"
                  style={{
                    background: 'rgba(var(--color-accent-rgb), 0.12)',
                    color: 'var(--color-accent)',
                  }}
                >
                  Optional
                </span>
              )}
            </div>
            <AnimatePresence mode="wait">
              {selected ? (
                <motion.span
                  key="selected"
                  initial={{ opacity: 0, y: 4 }}
                  animate={{ opacity: 1, y: 0 }}
                  exit={{ opacity: 0, y: -4 }}
                  className="inline-flex items-center gap-1.5 text-xs font-medium mt-0.5"
                  style={{ color: 'var(--color-primary)' }}
                >
                  <span className="w-4 h-4 rounded-full flex items-center justify-center" style={{ background: 'var(--gradient-primary)' }}>
                    <Check className="w-2.5 h-2.5 text-white" />
                  </span>
                  {selected.name}
                </motion.span>
              ) : (
                <motion.span
                  key="hint"
                  initial={{ opacity: 0, y: 4 }}
                  animate={{ opacity: 1, y: 0 }}
                  exit={{ opacity: 0, y: -4 }}
                  className="text-xs mt-0.5 block"
                  style={{ color: 'var(--color-text-muted)' }}
                >
                  Tailor scoring thresholds for your target domain
                </motion.span>
              )}
            </AnimatePresence>
          </div>

          {/* Chevron with accent background */}
          <motion.div
            animate={{ rotate: isExpanded ? 180 : 0 }}
            transition={{ duration: 0.3, ease: [0.4, 0, 0.2, 1] }}
            className={cn(
              'flex-shrink-0 w-8 h-8 rounded-lg flex items-center justify-center transition-all duration-200',
              isExpanded
                ? 'text-white shadow-[0_0_8px_var(--glow-primary)]'
                : 'text-[var(--color-accent)] group-hover:text-white group-hover:shadow-[0_0_8px_var(--glow-primary)]'
            )}
            style={{
              background: isExpanded
                ? 'var(--gradient-primary)'
                : undefined,
            }}
          >
            {/* Hover-only gradient fill for chevron */}
            {!isExpanded && (
              <div
                className="absolute inset-0 rounded-lg opacity-0 group-hover:opacity-100 transition-opacity duration-200"
                style={{ background: 'var(--gradient-primary)' }}
              />
            )}
            <ChevronDown className="w-4 h-4 relative z-10" />
          </motion.div>
        </div>
      </button>

      {/* Expanded profile cards */}
      <AnimatePresence>
        {isExpanded && (
          <motion.div
            initial={{ height: 0, opacity: 0 }}
            animate={{ height: 'auto', opacity: 1 }}
            exit={{ height: 0, opacity: 0 }}
            transition={{ duration: 0.3, ease: [0.4, 0, 0.2, 1] }}
            className="overflow-hidden"
          >
            <div className="pt-3">
              {isLoading ? (
                <div className="flex items-center gap-2 py-4 justify-center">
                  <div
                    className="w-5 h-5 rounded-full border-2 border-transparent animate-spin"
                    style={{ borderTopColor: 'var(--color-primary)' }}
                  />
                  <span className="text-sm" style={{ color: 'var(--color-text-muted)' }}>Loading profiles...</span>
                </div>
              ) : (
                <motion.div
                  className="grid grid-cols-1 sm:grid-cols-2 gap-2.5"
                  initial="initial"
                  animate="animate"
                  variants={{
                    initial: {},
                    animate: { transition: { staggerChildren: 0.05 } },
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
                        whileHover={{ y: -2, transition: { duration: 0.15 } }}
                        whileTap={{ scale: 0.98 }}
                        className={cn(
                          'relative w-full text-left rounded-xl overflow-hidden cursor-pointer',
                          'border transition-all duration-200',
                          isSelected
                            ? 'border-[var(--color-primary)] shadow-[0_0_18px_var(--glow-primary)]'
                            : cn('bg-[var(--color-surface-elevated)] hover:shadow-[var(--shadow-md)]', accent.border, 'hover:border-[var(--color-border-strong)]')
                        )}
                      >
                        {/* Background gradient when selected */}
                        {isSelected && (
                          <div
                            className="absolute inset-0"
                            style={{ background: 'var(--gradient-primary-subtle)' }}
                          />
                        )}

                        {/* Active top bar */}
                        {isSelected && (
                          <motion.div
                            layoutId="profile-bar"
                            className="absolute top-0 inset-x-0 h-[2px]"
                            style={{ background: 'var(--gradient-primary)' }}
                            transition={{ type: 'spring', stiffness: 500, damping: 30 }}
                          />
                        )}

                        <div className="relative px-3 py-2.5">
                          {/* Row 1: icon + name + description + check */}
                          <div className="flex items-center gap-2 mb-2">
                            <div className={cn(
                              'w-8 h-8 rounded-lg flex-shrink-0 flex items-center justify-center transition-all duration-200',
                              isSelected ? 'shadow-[0_0_8px_var(--glow-primary)]' : '',
                              !isSelected && accent.bg
                            )}
                              style={isSelected ? { background: 'var(--gradient-primary)' } : undefined}
                            >
                              <Icon className={cn(
                                'w-4 h-4',
                                isSelected ? 'text-white' : accent.text
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
                                className="flex-shrink-0 w-5 h-5 rounded-full flex items-center justify-center text-white shadow-[0_0_8px_var(--glow-primary)]"
                                style={{ background: 'var(--gradient-primary)' }}
                              >
                                <Check className="w-3 h-3" />
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
                                      ? 'bg-[var(--color-primary)]/[0.08]'
                                      : 'bg-[var(--color-surface-sunken)]/60'
                                  )}
                                >
                                  <span className="text-xs text-[var(--color-text-muted)]">
                                    {meta?.label ?? key}
                                  </span>
                                  <span className={cn(
                                    'text-xs font-semibold font-mono',
                                    isSelected ? 'text-[var(--color-primary)]' : accent.text
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
                    className="mt-2.5"
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
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}
