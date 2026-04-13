import { useState, useEffect, useLayoutEffect, useRef, useCallback } from 'react';
import { NavLink, useLocation, useNavigate } from 'react-router-dom';
import { motion, AnimatePresence } from 'framer-motion';
import {
  Atom,
  Layers,
  Info,
  BookOpen,
  FileText,
  ExternalLink,
  Bookmark,
  Clock,
  Menu,
  X,
  Beaker,
  Filter,
  Database,
  ChevronDown,
  FlaskConical,
  Stethoscope,
  Library,
} from 'lucide-react';

import { cn } from '../../lib/utils';
import { API_DOCS_URL } from '../../services/api';
import { ThemeToggle } from '../ui/ThemeToggle';

/* ─── Nav item types ─── */

interface NavItem {
  to: string;
  label: string;
  icon: React.FC<React.SVGProps<SVGSVGElement> & { size?: string | number }>;
  description: string;
}

interface NavGroup {
  label: string;
  icon: React.FC<React.SVGProps<SVGSVGElement> & { size?: string | number }>;
  description: string;
  items: NavItem[];
}

type NavEntry = NavItem | NavGroup;

function isGroup(entry: NavEntry): entry is NavGroup {
  return 'items' in entry;
}

/* ─── Navigation config ─── */

const navEntries: NavEntry[] = [
  {
    to: '/',
    label: 'Single Validation',
    icon: Atom,
    description: 'Validate, score, and profile individual chemical structures',
  },
  {
    to: '/batch',
    label: 'Batch Validation',
    icon: Layers,
    description: 'Process thousands of molecules with analytics and export',
  },
  {
    label: 'Data Preparation',
    icon: FlaskConical,
    description: 'Curate datasets for ML and generative workflows',
    items: [
      {
        to: '/qsar-ready',
        label: 'QSAR-Ready',
        icon: Beaker,
        description: '10-step structure curation pipeline for ML-ready SMILES',
      },
      {
        to: '/structure-filter',
        label: 'Structure Filter',
        icon: Filter,
        description: '6-stage validation funnel for generative model output',
      },
      {
        to: '/dataset-audit',
        label: 'Dataset Audit',
        icon: Database,
        description: 'Health scoring, contradiction detection, and curation reports',
      },
      {
        to: '/diagnostics',
        label: 'Diagnostics',
        icon: Stethoscope,
        description: 'SMILES error detection, InChI diff, round-trip checks, file validation',
      },
    ],
  },
  {
    label: 'Library',
    icon: Library,
    description: 'Your saved molecules and validation history',
    items: [
      {
        to: '/bookmarks',
        label: 'Bookmarks',
        icon: Bookmark,
        description: 'Saved molecules and validation snapshots',
      },
      {
        to: '/history',
        label: 'History',
        icon: Clock,
        description: 'Recently validated molecules and past results',
      },
    ],
  },
  {
    to: '/about',
    label: 'About',
    icon: Info,
    description: 'About ChemAudit, credits, and methodology',
  },
];

const externalLinks = [
  { href: 'https://kohulan.github.io/ChemAudit/', label: 'Docs', icon: FileText },
  { href: API_DOCS_URL, label: 'API', icon: BookOpen },
] as const;

/* ─── Tooltip component ─── */

function NavTooltip({ text, visible }: { text: string; visible: boolean }) {
  return (
    <AnimatePresence>
      {visible && (
        <motion.div
          initial={{ opacity: 0, y: 6, scale: 0.96 }}
          animate={{ opacity: 1, y: 0, scale: 1 }}
          exit={{ opacity: 0, y: 4, scale: 0.98 }}
          transition={{ duration: 0.15, ease: 'easeOut' }}
          className={cn(
            'absolute top-full left-1/2 -translate-x-1/2 mt-2 z-50',
            'px-3 py-1.5 rounded-lg',
            'bg-[var(--color-surface-elevated)] border border-[var(--color-border)]',
            'shadow-lg shadow-black/8 dark:shadow-black/30',
            'text-[11px] text-[var(--color-text-secondary)] font-medium leading-snug',
            'whitespace-nowrap pointer-events-none',
          )}
        >
          {text}
          {/* Arrow */}
          <div className="absolute -top-1 left-1/2 -translate-x-1/2 w-2 h-2 rotate-45 bg-[var(--color-surface-elevated)] border-l border-t border-[var(--color-border)]" />
        </motion.div>
      )}
    </AnimatePresence>
  );
}

/* ─── Dropdown menu for grouped items ─── */

function NavDropdown({
  group,
  isActive,
  hoveredNav,
  setHoveredNav,
  navTextColor,
  navigate,
  registerRef,
}: {
  group: NavGroup;
  isActive: (path: string) => boolean;
  hoveredNav: string | null;
  setHoveredNav: React.Dispatch<React.SetStateAction<string | null>>;
  navTextColor: (active: boolean, focused: boolean) => string;
  navigate: ReturnType<typeof useNavigate>;
  registerRef: (key: string) => (el: HTMLElement | null) => void;
}) {
  const [open, setOpen] = useState(false);
  const [tooltipVisible, setTooltipVisible] = useState(false);
  const timeoutRef = useRef<ReturnType<typeof setTimeout>>();
  const tooltipTimeoutRef = useRef<ReturnType<typeof setTimeout>>();
  const containerRef = useRef<HTMLDivElement>(null);

  const groupActive = group.items.some(i => isActive(i.to));
  const groupFocused = hoveredNav === `group:${group.label}`;

  const handleEnter = () => {
    clearTimeout(timeoutRef.current);
    setOpen(true);
    setHoveredNav(`group:${group.label}`);
    setTooltipVisible(false);
  };

  const handleLeave = () => {
    tooltipTimeoutRef.current && clearTimeout(tooltipTimeoutRef.current);
    setTooltipVisible(false);
    const groupKey = `group:${group.label}`;
    timeoutRef.current = setTimeout(() => {
      setOpen(false);
      // Only clear if still pointing at this group — prevents overwriting
      // a valid hover state set by another item's onMouseEnter
      setHoveredNav(prev => prev === groupKey ? null : prev);
    }, 150);
  };

  return (
    <div
      ref={containerRef}
      className="relative"
      onMouseEnter={handleEnter}
      onMouseLeave={handleLeave}
    >
      <motion.button
        ref={registerRef(`group:${group.label}`)}
        type="button"
        className={cn(
          'relative px-3 py-2 text-[13px] font-medium rounded-xl cursor-pointer',
          'flex items-center gap-1.5 transition-colors duration-200',
          navTextColor(groupActive, groupFocused),
        )}
        whileTap={{ scale: 0.96 }}
        onClick={() => {
          navigate(group.items[0].to);
          setOpen(false);
          setHoveredNav(null);
        }}
        onMouseEnter={() => {
          tooltipTimeoutRef.current = setTimeout(() => {
            if (!open) setTooltipVisible(true);
          }, 400);
        }}
      >
        {groupActive && (
          <motion.span
            className="relative z-10 w-[5px] h-[5px] rounded-full bg-[var(--color-primary)] shrink-0"
            initial={{ scale: 0 }}
            animate={{ scale: 1 }}
            transition={{ type: 'spring', stiffness: 500, damping: 25 }}
          />
        )}
        <group.icon className="relative z-10 w-3.5 h-3.5 shrink-0" />
        <span className="relative z-10 font-display whitespace-nowrap">{group.label}</span>
        <motion.div
          className="relative z-10"
          animate={{ rotate: open ? 180 : 0 }}
          transition={{ duration: 0.2 }}
        >
          <ChevronDown className="w-3 h-3 opacity-50" />
        </motion.div>
      </motion.button>

      {/* Tooltip (only when not open) */}
      {!open && <NavTooltip text={group.description} visible={tooltipVisible} />}

      {/* Dropdown panel */}
      <AnimatePresence>
        {open && (
          <motion.div
            initial={{ opacity: 0, y: -4, scale: 0.97 }}
            animate={{ opacity: 1, y: 0, scale: 1 }}
            exit={{ opacity: 0, y: -4, scale: 0.97 }}
            transition={{ duration: 0.18, ease: [0.25, 0.46, 0.45, 0.94] }}
            className={cn(
              'absolute top-full left-0 mt-2 z-50',
              'min-w-[280px] rounded-xl overflow-hidden',
              'bg-[var(--color-surface-elevated)] border border-[var(--color-border)]',
              'shadow-xl shadow-black/10 dark:shadow-black/40',
              'backdrop-blur-xl',
            )}
          >
            {/* Dropdown inner glow */}
            <div className="absolute inset-0 rounded-xl bg-gradient-to-b from-white/40 dark:from-white/[0.03] to-transparent pointer-events-none" />

            <div className="relative py-1.5">
              {group.items.map((item) => {
                const active = isActive(item.to);
                return (
                  <NavLink
                    key={item.to}
                    to={item.to}
                    className={cn(
                      'flex items-start gap-3 px-4 py-3 mx-1.5 rounded-lg',
                      'transition-all duration-150 group/item',
                      active
                        ? 'bg-[rgba(var(--color-primary-rgb),0.06)]'
                        : 'hover:bg-[var(--color-surface-sunken)]',
                    )}
                  >
                    <div
                      className={cn(
                        'mt-0.5 w-8 h-8 rounded-lg flex items-center justify-center shrink-0',
                        'transition-all duration-200',
                        active
                          ? 'bg-[rgba(var(--color-primary-rgb),0.1)] text-[var(--color-primary)]'
                          : 'bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)] group-hover/item:text-[var(--color-text-primary)] group-hover/item:bg-[rgba(var(--color-primary-rgb),0.06)]',
                      )}
                    >
                      <item.icon className="w-4 h-4" />
                    </div>
                    <div className="min-w-0">
                      <div className={cn(
                        'text-[13px] font-semibold font-display leading-snug',
                        active ? 'text-[var(--color-primary)]' : 'text-[var(--color-text-primary)]',
                      )}>
                        {item.label}
                      </div>
                      <div className="text-[11px] text-[var(--color-text-muted)] leading-snug mt-0.5">
                        {item.description}
                      </div>
                    </div>
                    {active && (
                      <span className="mt-2 w-1.5 h-1.5 rounded-full bg-[var(--color-primary)] shrink-0 ml-auto" />
                    )}
                  </NavLink>
                );
              })}
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}

/* ─── Header ─── */

export function Header() {
  const location = useLocation();
  const navigate = useNavigate();
  const [hoveredNav, setHoveredNav] = useState<string | null>(null);
  const [hoveredTooltip, setHoveredTooltip] = useState<string | null>(null);
  const [mobileOpen, setMobileOpen] = useState(false);
  const [scrolled, setScrolled] = useState(false);
  const tooltipTimerRef = useRef<ReturnType<typeof setTimeout>>();
  const navRef = useRef<HTMLDivElement>(null);
  const itemRefs = useRef<Map<string, HTMLElement>>(new Map());
  const registerItemRef = useCallback((key: string) => (el: HTMLElement | null) => {
    if (el) itemRefs.current.set(key, el);
    else itemRefs.current.delete(key);
  }, []);

  useEffect(() => {
    const onScroll = () => setScrolled(window.scrollY > 16);
    window.addEventListener('scroll', onScroll, { passive: true });
    return () => window.removeEventListener('scroll', onScroll);
  }, []);

  useEffect(() => { setMobileOpen(false); }, [location.pathname]);

  function isActive(path: string): boolean {
    if (path === '/') return location.pathname === '/';
    return location.pathname.startsWith(path);
  }

  // Compute active nav target — maps group children to group key for pill placement
  const allDirectItems = navEntries.filter((e): e is NavItem => !isGroup(e));
  const activeNav: string | null = allDirectItems.find(item => isActive(item.to))?.to
    ?? (() => {
      const g = navEntries.filter(isGroup).find(g => g.items.some(i => isActive(i.to)));
      return g ? `group:${g.label}` : null;
    })();
  // Unified pill target: hovered item takes priority, falls back to active route
  const pillTarget = hoveredNav ?? activeNav;
  const pillIsActive = pillTarget !== null && pillTarget === activeNav;

  // Measure pill position from the target element's ref
  const [pillRect, setPillRect] = useState<{ left: number; top: number; width: number; height: number } | null>(null);

  useLayoutEffect(() => {
    if (!pillTarget || !navRef.current) {
      setPillRect(null);
      return;
    }
    const el = itemRefs.current.get(pillTarget);
    if (!el) {
      setPillRect(null);
      return;
    }
    const navRect = navRef.current.getBoundingClientRect();
    const elRect = el.getBoundingClientRect();
    setPillRect({
      left: elRect.left - navRect.left,
      top: elRect.top - navRect.top,
      width: elRect.width,
      height: elRect.height,
    });
  }, [pillTarget]);

  function navTextColor(active: boolean, focused: boolean): string {
    if (active) return 'text-[var(--color-primary)]';
    if (focused) return 'text-[var(--color-text-primary)]';
    return 'text-[var(--color-text-secondary)]';
  }

  // All items flat for mobile menu
  const allMobileItems: (NavItem & { section?: string })[] = [];
  for (const entry of navEntries) {
    if (isGroup(entry)) {
      for (const item of entry.items) {
        allMobileItems.push({ ...item, section: entry.label });
      }
    } else {
      allMobileItems.push(entry);
    }
  }

  return (
    <header className="sticky top-0 z-50 flex justify-center px-3 sm:px-5 pt-3 sm:pt-4 pb-1">
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
        {/* Glass background */}
        <motion.div
          className="absolute inset-0 rounded-2xl backdrop-blur-2xl backdrop-saturate-[1.8]"
          animate={{ opacity: scrolled ? 0.6 : 1 }}
          transition={{ duration: 0.5, ease: 'easeOut' }}
          style={{
            background: 'linear-gradient(to bottom, var(--color-surface-overlay) 0%, transparent 100%)',
          }}
        />

        {/* Top-edge shimmer */}
        <div className="absolute inset-x-6 top-0 h-px overflow-hidden rounded-t-2xl">
          <div className="h-full bg-gradient-to-r from-transparent via-white/50 dark:via-white/[0.08] to-transparent" />
        </div>

        {/* Bottom accent gradient */}
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

          {/* ── Desktop nav ── */}
          <div className="hidden lg:flex items-center gap-0.5">
            <nav
              ref={navRef}
              className="relative flex items-center gap-0.5"
              onMouseLeave={() => {
                setHoveredNav(null);
                setHoveredTooltip(null);
                tooltipTimerRef.current && clearTimeout(tooltipTimerRef.current);
              }}
            >
              {/* Single sliding focus pill — never unmounts/remounts per item */}
              {pillRect && (
                <motion.div
                  className={cn(
                    'absolute top-0 left-0 rounded-xl pointer-events-none z-0',
                    'transition-[background-color,box-shadow] duration-200',
                    pillIsActive
                      ? 'bg-[rgba(var(--color-primary-rgb),0.08)] shadow-[0_0_16px_var(--glow-soft)]'
                      : 'bg-[var(--color-surface-sunken)]',
                  )}
                  initial={false}
                  animate={{
                    x: pillRect.left,
                    y: pillRect.top,
                    width: pillRect.width,
                    height: pillRect.height,
                  }}
                  transition={{ type: 'spring', bounce: 0, duration: 0.25 }}
                />
              )}

              {navEntries.map((entry) => {
                if (isGroup(entry)) {
                  return (
                    <NavDropdown
                      key={entry.label}
                      group={entry}
                      isActive={isActive}
                      hoveredNav={hoveredNav}
                      setHoveredNav={setHoveredNav}
                      navTextColor={navTextColor}
                      navigate={navigate}
                      registerRef={registerItemRef}
                    />
                  );
                }

                const item = entry;
                const active = isActive(item.to);
                const focused = pillTarget === item.to;

                return (
                  <NavLink
                    key={item.to}
                    to={item.to}
                    onMouseEnter={() => {
                      setHoveredNav(item.to);
                      tooltipTimerRef.current && clearTimeout(tooltipTimerRef.current);
                      tooltipTimerRef.current = setTimeout(() => setHoveredTooltip(item.to), 500);
                    }}
                    onMouseLeave={() => {
                      tooltipTimerRef.current && clearTimeout(tooltipTimerRef.current);
                      setHoveredTooltip(null);
                    }}
                    className="relative"
                  >
                    <motion.div
                      ref={registerItemRef(item.to)}
                      className={cn(
                        'relative px-3 py-2 text-[13px] font-medium rounded-xl cursor-pointer',
                        'flex items-center gap-1.5 transition-colors duration-200',
                        navTextColor(active, focused),
                      )}
                      whileTap={{ scale: 0.96 }}
                    >
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
                    <NavTooltip text={item.description} visible={hoveredTooltip === item.to} />
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
                {(() => {
                  let lastSection: string | undefined;
                  return allMobileItems.map((item, i) => {
                    const active = isActive(item.to);
                    const showSection = item.section && item.section !== lastSection;
                    lastSection = item.section;
                    return (
                      <div key={item.to}>
                        {showSection && (
                          <motion.div
                            initial={{ opacity: 0 }}
                            animate={{ opacity: 1 }}
                            transition={{ delay: i * 0.04 }}
                            className="px-4 pt-3 pb-1 text-[10px] font-bold uppercase tracking-widest text-[var(--color-text-muted)]"
                          >
                            {item.section}
                          </motion.div>
                        )}
                        <motion.div
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
                            <div>
                              <div>{item.label}</div>
                              <div className="text-[10px] text-[var(--color-text-muted)] font-normal">{item.description}</div>
                            </div>
                          </NavLink>
                        </motion.div>
                      </div>
                    );
                  });
                })()}

                <div className="h-px bg-[var(--color-border)] mx-3 my-1.5" />

                {externalLinks.map((item, i) => (
                  <motion.a
                    key={item.href}
                    href={item.href}
                    target="_blank"
                    rel="noopener noreferrer"
                    initial={{ x: -16, opacity: 0 }}
                    animate={{ x: 0, opacity: 1 }}
                    transition={{ delay: (allMobileItems.length + i) * 0.04, duration: 0.2, ease: 'easeOut' }}
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
