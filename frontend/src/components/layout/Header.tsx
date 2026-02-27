import { useState, useEffect } from 'react';
import { NavLink, useLocation } from 'react-router-dom';
import { motion, AnimatePresence } from 'framer-motion';
import {
  Atom, Grid3X3, Info, BookOpen, FileText,
  ExternalLink, Star, Clock, Menu, X,
} from 'lucide-react';

import { cn } from '../../lib/utils';
import { API_DOCS_URL } from '../../services/api';
import { ThemeToggle } from '../ui/ThemeToggle';

const navItems = [
  { to: '/', label: 'Single Validation', icon: Atom },
  { to: '/batch', label: 'Batch Validation', icon: Grid3X3 },
  { to: '/bookmarks', label: 'Bookmarks', icon: Star },
  { to: '/history', label: 'History', icon: Clock },
  { to: '/about', label: 'About', icon: Info },
] as const;

const externalLinks = [
  { href: 'https://kohulan.github.io/ChemAudit/', label: 'Docs', icon: FileText },
  { href: API_DOCS_URL, label: 'API', icon: BookOpen },
] as const;

/**
 * Floating glass-pill header with sliding focus highlight and spring animations.
 */
export function Header() {
  const location = useLocation();
  const [hoveredNav, setHoveredNav] = useState<string | null>(null);
  const [mobileOpen, setMobileOpen] = useState(false);
  const [scrolled, setScrolled] = useState(false);

  // Strengthen glass effect on scroll
  useEffect(() => {
    const onScroll = () => setScrolled(window.scrollY > 16);
    window.addEventListener('scroll', onScroll, { passive: true });
    return () => window.removeEventListener('scroll', onScroll);
  }, []);

  // Close mobile menu on route change
  useEffect(() => { setMobileOpen(false); }, [location.pathname]);

  function isActive(path: string): boolean {
    if (path === '/') return location.pathname === '/';
    return location.pathname.startsWith(path);
  }

  // The focus pill follows hover; when nothing is hovered it snaps back to active route
  const activeNav = navItems.find(item => isActive(item.to))?.to;
  const focusTarget = hoveredNav ?? activeNav;

  function navTextColor(active: boolean, focused: boolean): string {
    if (active) return 'text-[var(--color-primary)]';
    if (focused) return 'text-[var(--color-text-primary)]';
    return 'text-[var(--color-text-secondary)]';
  }

  return (
      <header className="sticky top-0 z-50 flex justify-center px-3 sm:px-5 pt-3 sm:pt-4 pb-1">
        {/* ── Floating pill container ── */}
        <motion.div
          initial={{ y: -28, opacity: 0, scale: 0.97 }}
          animate={{ y: 0, opacity: 1, scale: 1 }}
          transition={{ duration: 0.8, ease: [0.22, 1, 0.36, 1] }}
          className={cn(
            'relative w-full max-w-[1200px] rounded-2xl',
            'transition-shadow duration-500 ease-out',
            scrolled
              ? 'shadow-[0_8px_40px_rgba(0,0,0,0.08),0_2px_6px_rgba(0,0,0,0.04)] dark:shadow-[0_8px_40px_rgba(0,0,0,0.45)]'
              : 'shadow-[0_4px_20px_rgba(0,0,0,0.03)] dark:shadow-[0_4px_20px_rgba(0,0,0,0.15)]',
          )}
        >
          {/* Glass background — top-to-bottom opacity gradient, fades on scroll */}
          <motion.div
            className="absolute inset-0 rounded-2xl backdrop-blur-2xl backdrop-saturate-[1.8]"
            animate={{ opacity: scrolled ? 0.6 : 1 }}
            transition={{ duration: 0.5, ease: 'easeOut' }}
            style={{
              background: 'linear-gradient(to bottom, var(--color-surface-overlay) 0%, transparent 100%)',
            }}
          />

          {/* Top-edge white shimmer */}
          <div className="absolute inset-x-6 top-0 h-px overflow-hidden rounded-t-2xl">
            <div className="h-full bg-gradient-to-r from-transparent via-white/50 dark:via-white/[0.08] to-transparent" />
          </div>

          {/* Bottom accent gradient — stronger on scroll */}
          <motion.div
            className="absolute inset-x-10 bottom-0 h-px overflow-hidden"
            animate={{ opacity: scrolled ? 0.7 : 0.2 }}
            transition={{ duration: 0.5 }}
          >
            <div className="h-full bg-gradient-to-r from-transparent via-[var(--color-primary)] to-transparent" />
          </motion.div>

          {/* Border ring */}
          <div className="absolute inset-0 rounded-2xl border border-[var(--color-border)] pointer-events-none" />

          {/* ── Content row ── */}
          <div className="relative flex items-center justify-between h-14 sm:h-[60px] px-4 sm:px-6">

            {/* Logo */}
            <NavLink to="/" className="flex items-center gap-2.5 group shrink-0">
              <motion.div
                className={cn(
                  'w-8 h-8 sm:w-9 sm:h-9 rounded-xl flex items-center justify-center overflow-hidden',
                  'bg-gradient-to-br from-[var(--color-surface-elevated)] to-[var(--color-surface-sunken)]',
                  'border border-[var(--color-border)] shadow-sm',
                  'group-hover:shadow-[0_4px_20px_var(--glow-soft)]',
                  'group-hover:border-[rgba(var(--color-primary-rgb),0.2)]',
                  'transition-all duration-300',
                )}
                whileHover={{ scale: 1.08, rotate: -3 }}
                whileTap={{ scale: 0.93 }}
              >
                <img src="/logo.png" alt="ChemAudit" className="w-full h-full object-contain" />
              </motion.div>
              <motion.div
                className="hidden sm:block"
                whileHover={{ x: 2 }}
                transition={{ type: 'spring', stiffness: 300, damping: 20 }}
              >
                <h1 className="text-lg font-bold text-[var(--color-text-primary)] tracking-tight leading-none font-display">
                  <span className="font-extrabold text-[var(--color-primary)]">Chem</span>
                  <span className="font-semibold">Audit</span>
                </h1>
              </motion.div>
            </NavLink>

            {/* ── Desktop right group: nav + external + theme ── */}
            <div className="hidden lg:flex items-center gap-0.5">
              {/* Primary navigation */}
              <nav
                className="flex items-center gap-0.5"
                onMouseLeave={() => setHoveredNav(null)}
              >
                {navItems.map(item => {
                  const active = isActive(item.to);
                  const focused = focusTarget === item.to;

                  return (
                    <NavLink
                      key={item.to}
                      to={item.to}
                      onMouseEnter={() => setHoveredNav(item.to)}
                      className="relative"
                    >
                      <motion.div
                        className={cn(
                          'relative px-3 py-2 text-[13px] font-medium rounded-xl cursor-pointer',
                          'flex items-center gap-1.5 transition-colors duration-200',
                          navTextColor(active, focused),
                        )}
                        whileTap={{ scale: 0.96 }}
                      >
                        {/* Sliding focus pill — follows hover, returns to active */}
                        {focused && (
                          <motion.div
                            layoutId="navFocus"
                            className={cn(
                              'absolute inset-0 rounded-xl transition-colors duration-200',
                              active
                                ? 'bg-[rgba(var(--color-primary-rgb),0.08)] shadow-[0_0_16px_var(--glow-soft)]'
                                : 'bg-[var(--color-surface-sunken)]',
                            )}
                            initial={false}
                            transition={{ type: 'spring', stiffness: 380, damping: 28, mass: 0.8 }}
                          />
                        )}

                        {/* Active dot */}
                        {active && (
                          <motion.span
                            className="relative z-10 w-[5px] h-[5px] rounded-full bg-[var(--color-primary)] shrink-0"
                            initial={{ scale: 0, opacity: 0 }}
                            animate={{ scale: 1, opacity: 1 }}
                            transition={{ type: 'spring', stiffness: 500, damping: 25 }}
                          />
                        )}

                        <item.icon className="relative z-10 w-3.5 h-3.5 shrink-0" />
                        <span className="relative z-10 font-display whitespace-nowrap">{item.label}</span>
                      </motion.div>
                    </NavLink>
                  );
                })}
              </nav>

              {/* Divider */}
              <div className="w-px h-5 bg-[var(--color-border-strong)] mx-2.5" />

              {/* External links */}
              {externalLinks.map(item => (
                <motion.a
                  key={item.href}
                  href={item.href}
                  target="_blank"
                  rel="noopener noreferrer"
                  className={cn(
                    'relative px-2.5 py-1.5 text-[13px] font-medium rounded-lg cursor-pointer',
                    'flex items-center gap-1.5',
                    'text-[var(--color-text-muted)] hover:text-[var(--color-text-primary)]',
                    'transition-colors duration-200',
                  )}
                  whileHover={{ y: -1 }}
                  whileTap={{ scale: 0.96 }}
                >
                  <item.icon className="w-3.5 h-3.5" />
                  <span className="font-display">{item.label}</span>
                  <ExternalLink className="w-2.5 h-2.5 opacity-30" />
                </motion.a>
              ))}

              {/* Divider */}
              <div className="w-px h-5 bg-[var(--color-border-strong)] mx-2.5" />

              <ThemeToggle />
            </div>

            {/* ── Mobile controls ── */}
            <div className="flex lg:hidden items-center gap-1.5">
              <ThemeToggle />
              <motion.button
                onClick={() => setMobileOpen(v => !v)}
                className={cn(
                  'p-2 rounded-xl cursor-pointer',
                  'text-[var(--color-text-secondary)]',
                  'hover:bg-[var(--color-surface-sunken)]',
                  'transition-colors duration-150',
                )}
                whileTap={{ scale: 0.9 }}
                aria-label="Toggle navigation menu"
              >
                <AnimatePresence mode="wait" initial={false}>
                  <motion.div
                    key={mobileOpen ? 'close' : 'open'}
                    initial={{ rotate: mobileOpen ? -90 : 90, opacity: 0 }}
                    animate={{ rotate: 0, opacity: 1 }}
                    exit={{ rotate: mobileOpen ? 90 : -90, opacity: 0 }}
                    transition={{ duration: 0.15 }}
                  >
                    {mobileOpen ? <X className="w-5 h-5" /> : <Menu className="w-5 h-5" />}
                  </motion.div>
                </AnimatePresence>
              </motion.button>
            </div>
          </div>

          {/* ── Mobile dropdown ── */}
          <AnimatePresence>
            {mobileOpen && (
              <motion.div
                initial={{ height: 0, opacity: 0 }}
                animate={{ height: 'auto', opacity: 1 }}
                exit={{ height: 0, opacity: 0 }}
                transition={{ duration: 0.3, ease: [0.25, 0.46, 0.45, 0.94] }}
                className="overflow-hidden lg:hidden"
              >
                <div className="border-t border-[var(--color-border)] mx-4" />
                <div className="px-3 py-3 flex flex-col gap-0.5 bg-[var(--color-surface-elevated)]/95 backdrop-blur-lg rounded-b-2xl">
                  {navItems.map((item, i) => {
                    const active = isActive(item.to);
                    return (
                      <motion.div
                        key={item.to}
                        initial={{ x: -16, opacity: 0 }}
                        animate={{ x: 0, opacity: 1 }}
                        transition={{ delay: i * 0.04, duration: 0.2, ease: 'easeOut' }}
                      >
                        <NavLink
                          to={item.to}
                          className={cn(
                            'flex items-center gap-3 px-4 py-2.5 rounded-xl text-sm font-medium font-display cursor-pointer',
                            'transition-colors duration-150',
                            active
                              ? 'text-[var(--color-primary)] bg-[rgba(var(--color-primary-rgb),0.08)]'
                              : 'text-[var(--color-text-secondary)] hover:bg-[var(--color-surface-sunken)]',
                          )}
                        >
                          {active && (
                            <span className="w-1.5 h-1.5 rounded-full bg-[var(--color-primary)] shrink-0" />
                          )}
                          <item.icon className="w-4 h-4 shrink-0" />
                          {item.label}
                        </NavLink>
                      </motion.div>
                    );
                  })}

                  <div className="h-px bg-[var(--color-border)] mx-3 my-1.5" />

                  {externalLinks.map((item, i) => (
                    <motion.a
                      key={item.href}
                      href={item.href}
                      target="_blank"
                      rel="noopener noreferrer"
                      initial={{ x: -16, opacity: 0 }}
                      animate={{ x: 0, opacity: 1 }}
                      transition={{ delay: (navItems.length + i) * 0.04, duration: 0.2, ease: 'easeOut' }}
                      className={cn(
                        'flex items-center gap-3 px-4 py-2.5 rounded-xl text-sm font-medium font-display cursor-pointer',
                        'text-[var(--color-text-muted)] hover:bg-[var(--color-surface-sunken)]',
                        'transition-colors duration-150',
                      )}
                    >
                      <item.icon className="w-4 h-4 shrink-0" />
                      {item.label}
                      <ExternalLink className="w-3 h-3 ml-auto opacity-30" />
                    </motion.a>
                  ))}
                </div>
              </motion.div>
            )}
          </AnimatePresence>
        </motion.div>
      </header>
  );
}
