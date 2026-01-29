import { motion, AnimatePresence } from 'framer-motion';
import { Sun, Moon, Monitor } from 'lucide-react';
import { useThemeContext } from '../../contexts/ThemeContext';
import type { Theme } from '../../hooks/useTheme';
import { cn } from '../../lib/utils';

interface ThemeToggleProps {
  className?: string;
  showLabel?: boolean;
  variant?: 'button' | 'dropdown';
}

const themes: { value: Theme; label: string; icon: typeof Sun }[] = [
  { value: 'light', label: 'Light', icon: Sun },
  { value: 'dark', label: 'Dark', icon: Moon },
  { value: 'system', label: 'System', icon: Monitor },
];

/**
 * Theme toggle button with animated icon transition
 */
export function ThemeToggle({ className, showLabel = false, variant = 'button' }: ThemeToggleProps) {
  const { theme, setTheme, resolvedTheme } = useThemeContext();

  if (variant === 'dropdown') {
    return <ThemeDropdown className={className} />;
  }

  const cycleTheme = () => {
    const currentIndex = themes.findIndex(t => t.value === theme);
    const nextIndex = (currentIndex + 1) % themes.length;
    setTheme(themes[nextIndex].value);
  };

  const Icon = resolvedTheme === 'dark' ? Moon : Sun;

  return (
    <motion.button
      onClick={cycleTheme}
      className={cn(
        'relative p-2.5 rounded-xl',
        'bg-[var(--color-surface-sunken)] hover:bg-[var(--color-surface-elevated)]',
        'border border-[var(--color-border)]',
        'transition-all duration-200',
        className
      )}
      whileHover={{ scale: 1.05 }}
      whileTap={{ scale: 0.95 }}
      aria-label={`Current theme: ${theme}. Click to change.`}
    >
      <AnimatePresence mode="wait" initial={false}>
        <motion.div
          key={resolvedTheme}
          initial={{ rotate: -90, opacity: 0 }}
          animate={{ rotate: 0, opacity: 1 }}
          exit={{ rotate: 90, opacity: 0 }}
          transition={{ duration: 0.2 }}
        >
          <Icon className="w-5 h-5 text-[var(--color-text-primary)]" />
        </motion.div>
      </AnimatePresence>

      {showLabel && (
        <span className="ml-2 text-sm font-medium text-[var(--color-text-primary)] capitalize">
          {theme}
        </span>
      )}
    </motion.button>
  );
}

/**
 * Theme dropdown selector
 */
function ThemeDropdown({ className }: { className?: string }) {
  const { theme, setTheme } = useThemeContext();

  return (
    <div className={cn('relative', className)}>
      <div className={cn(
        'flex items-center gap-1 p-1 rounded-xl',
        'bg-[var(--color-surface-sunken)]',
        'border border-[var(--color-border)]'
      )}>
        {themes.map(({ value, icon: Icon }) => (
          <motion.button
            key={value}
            onClick={() => setTheme(value)}
            className={cn(
              'relative p-2 rounded-lg transition-colors duration-200',
              theme === value
                ? 'text-[var(--color-text-primary)]'
                : 'text-[var(--color-text-muted)] hover:text-[var(--color-text-secondary)]'
            )}
            whileHover={{ scale: 1.05 }}
            whileTap={{ scale: 0.95 }}
            aria-label={`Switch to ${value} theme`}
          >
            {theme === value && (
              <motion.div
                layoutId="theme-indicator"
                className="absolute inset-0 bg-[var(--color-primary)]/10 rounded-lg"
                transition={{ type: 'spring', stiffness: 500, damping: 30 }}
              />
            )}
            <Icon className="w-4 h-4 relative z-10" />
          </motion.button>
        ))}
      </div>
    </div>
  );
}

/**
 * Full theme selector with labels (for settings pages)
 */
export function ThemeSelector({ className }: { className?: string }) {
  const { theme, setTheme } = useThemeContext();

  return (
    <div className={cn('space-y-2', className)}>
      <label className="text-sm font-medium text-[var(--color-text-primary)]">Theme</label>
      <div className="flex flex-wrap gap-2">
        {themes.map(({ value, label, icon: Icon }) => (
          <motion.button
            key={value}
            onClick={() => setTheme(value)}
            className={cn(
              'flex items-center gap-2 px-4 py-2.5 rounded-xl',
              'border-2 transition-all duration-200',
              theme === value
                ? 'border-[var(--color-primary)] bg-[var(--color-primary)]/10 text-[var(--color-primary)]'
                : 'border-transparent bg-[var(--color-surface-elevated)] text-[var(--color-text-secondary)] hover:bg-[var(--color-primary)]/5'
            )}
            whileHover={{ scale: 1.02 }}
            whileTap={{ scale: 0.98 }}
          >
            <Icon className="w-4 h-4" />
            <span className="text-sm font-medium">{label}</span>
          </motion.button>
        ))}
      </div>
    </div>
  );
}
