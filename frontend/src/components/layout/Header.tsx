import { motion } from 'framer-motion';
import { BookOpen, ExternalLink } from 'lucide-react';
import { ThemeToggle } from '../ui/ThemeToggle';
import { cn } from '../../lib/utils';

/**
 * Premium header with glassmorphism and gradient accent
 */
export function Header() {
  return (
    <header className="sticky top-0 z-50">
      {/* Glass background */}
      <div
        className={cn(
          'absolute inset-0',
          'bg-[var(--color-surface-elevated)]/80 dark:bg-[var(--color-surface)]/80',
          'backdrop-blur-xl',
          'border-b border-[var(--color-border)]'
        )}
      />

      {/* Gradient accent line */}
      <div className="absolute bottom-0 left-0 right-0 h-[1px] bg-gradient-to-r from-transparent via-[var(--color-primary)]/40 to-transparent" />

      <div className="relative max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
        <div className="flex items-center justify-between h-16">
          {/* Logo and branding */}
          <motion.a
            href="/"
            className="flex items-center gap-3 group"
            whileHover={{ scale: 1.02 }}
            whileTap={{ scale: 0.98 }}
          >
            {/* Molecule icon with gradient */}
            <div
              className={cn(
                'w-10 h-10 rounded-xl flex items-center justify-center',
                'bg-gradient-to-br from-[var(--color-primary)] to-[var(--color-accent)]',
                'shadow-[0_2px_8px_var(--glow-primary)]',
                'group-hover:shadow-[0_4px_16px_var(--glow-primary)]',
                'transition-shadow duration-300'
              )}
            >
              <svg
                className="w-5 h-5 text-white"
                viewBox="0 0 24 24"
                fill="currentColor"
              >
                {/* Benzene-like hexagon with atoms */}
                <circle cx="12" cy="12" r="3" />
                <circle cx="12" cy="4" r="2" />
                <circle cx="18.9" cy="8" r="2" />
                <circle cx="18.9" cy="16" r="2" />
                <circle cx="12" cy="20" r="2" />
                <circle cx="5.1" cy="16" r="2" />
                <circle cx="5.1" cy="8" r="2" />
                {/* Bonds */}
                <line x1="12" y1="9" x2="12" y2="6" stroke="currentColor" strokeWidth="1.5" opacity="0.7" />
                <line x1="14.6" y1="10.5" x2="17" y2="9" stroke="currentColor" strokeWidth="1.5" opacity="0.7" />
                <line x1="14.6" y1="13.5" x2="17" y2="15" stroke="currentColor" strokeWidth="1.5" opacity="0.7" />
                <line x1="12" y1="15" x2="12" y2="18" stroke="currentColor" strokeWidth="1.5" opacity="0.7" />
                <line x1="9.4" y1="13.5" x2="7" y2="15" stroke="currentColor" strokeWidth="1.5" opacity="0.7" />
                <line x1="9.4" y1="10.5" x2="7" y2="9" stroke="currentColor" strokeWidth="1.5" opacity="0.7" />
              </svg>
            </div>
            <div>
              <h1 className="text-lg font-bold text-[var(--color-text-primary)] tracking-tight">
                ChemVault
              </h1>
              <p className="text-[10px] text-[var(--color-text-muted)] -mt-0.5 tracking-wide uppercase">
                Structure Validation
              </p>
            </div>
          </motion.a>

          {/* Navigation */}
          <nav className="flex items-center gap-1">
            <NavLink href="/" active>
              Validate
            </NavLink>
            <NavLink
              href="http://localhost:8000/docs"
              external
              icon={<BookOpen className="w-4 h-4" />}
            >
              API Docs
            </NavLink>

            {/* Divider */}
            <div className="w-px h-6 bg-[var(--color-border)] mx-2" />

            {/* Theme Toggle */}
            <ThemeToggle />
          </nav>
        </div>
      </div>
    </header>
  );
}

interface NavLinkProps {
  href: string;
  children: React.ReactNode;
  active?: boolean;
  external?: boolean;
  icon?: React.ReactNode;
}

function NavLink({ href, children, active, external, icon }: NavLinkProps) {
  return (
    <motion.a
      href={href}
      target={external ? '_blank' : undefined}
      rel={external ? 'noopener noreferrer' : undefined}
      className={cn(
        'px-3 py-2 text-sm font-medium rounded-lg transition-all duration-200',
        'flex items-center gap-2',
        active
          ? 'bg-[var(--color-primary)]/10 text-[var(--color-primary)]'
          : 'text-[var(--color-text-secondary)] hover:text-[var(--color-text-primary)] hover:bg-[var(--color-surface-sunken)]'
      )}
      whileHover={{ scale: 1.02 }}
      whileTap={{ scale: 0.98 }}
    >
      {icon}
      {children}
      {external && <ExternalLink className="w-3 h-3 opacity-50" />}
    </motion.a>
  );
}
