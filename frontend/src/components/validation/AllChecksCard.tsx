import { useEffect, useMemo, useState } from 'react';
import { motion } from 'framer-motion';
import {
  Atom,
  AlertTriangle,
  CheckCircle2,
  ChevronDown,
  FlaskConical,
  GitMerge,
  Hexagon,
  Layers,
  Lock,
  type LucideIcon,
  Share2,
} from 'lucide-react';
import { cn } from '../../lib/utils';
import type { CheckResult } from '../../types/validation';

const CHECK_DESCRIPTIONS: Record<string, string> = {
  // Basic checks
  parsability: 'Verifies the input string can be parsed into a valid molecular structure by RDKit.',
  sanitization: 'Checks if RDKit can sanitize the molecule (assign aromaticity, add implicit hydrogens, validate bonds).',
  valence: 'Validates that all atoms have chemically valid valence states (e.g., carbon with 4 bonds, nitrogen with 3).',
  aromaticity: 'Confirms aromatic ring systems are properly defined and assigned by RDKit\'s aromaticity model.',
  connectivity: 'Checks molecular connectivity - ensures the structure is a single connected component without fragments.',
  // Stereo checks
  undefined_stereocenters: 'Identifies chiral centers (sp3 carbons with 4 different substituents) that lack R/S stereochemistry assignment.',
  undefined_doublebond_stereo: 'Finds double bonds that could have E/Z isomerism but lack defined geometry.',
  conflicting_stereo: 'Detects contradictory stereochemistry assignments that cannot exist in a real molecule.',
  // Representation checks
  smiles_roundtrip: 'Tests if converting SMILES → molecule → SMILES preserves the structure identity.',
  inchi_generation: 'Verifies that a valid InChI identifier can be generated for the molecule.',
  inchi_roundtrip: 'Tests if converting to InChI and back preserves the molecular structure.',
  // Deep: Stereo & Tautomers
  stereoisomer_enumeration: 'Finds undefined stereocenters and enumerates possible stereoisomers (up to 128). Helps identify ambiguous chirality.',
  tautomer_detection: 'Detects tautomeric forms and identifies the canonical tautomer. Reports whether the input matches the canonical form.',
  aromatic_system_validation: 'Checks for unusual aromatic ring sizes (not 5 or 6 membered) and charged aromatic atoms that may indicate issues.',
  coordinate_dimension: 'Reports whether the molecule has 2D coordinates, 3D coordinates, or no coordinate information.',
  // Deep: Chemical Composition
  mixture_detection: 'Detects multi-fragment inputs (dot-separated SMILES) and classifies each fragment as drug, salt, solvent, or unknown.',
  solvent_contamination: 'Screens for common lab solvents (water, DMSO, DMF, methanol, etc.) that may contaminate the input structure.',
  inorganic_filter: 'Flags molecules lacking carbon (inorganic) or containing metal atoms (organometallic) that may not suit standard validation.',
  radical_detection: 'Identifies atoms with unpaired radical electrons that may indicate unstable or reactive species.',
  isotope_label_detection: 'Detects isotope-labeled atoms (deuterium, carbon-13, etc.) often used in pharmacokinetic studies.',
  trivial_molecule: 'Flags molecules with 3 or fewer heavy atoms as too small for meaningful chemical validation.',
  // Deep: Structural Complexity
  hypervalent_atoms: 'Detects atoms exceeding their normal valence limits, which may indicate unusual bonding or input errors.',
  polymer_detection: 'Identifies possible polymers via SGroup markers, molecular weight above 1500 Da, or dummy atom attachment points.',
  ring_strain: 'Flags 3-membered (cyclopropane) and 4-membered (cyclobutane) rings that have significant ring strain.',
  macrocycle_detection: 'Identifies macrocyclic rings with more than 12 atoms, common in natural products and cyclic peptides.',
  charged_species: 'Reports formal charges, identifies zwitterions (net charge zero with both positive and negative atoms).',
  explicit_hydrogen_audit: 'Reports atoms with explicit hydrogen specifications and detects H atom objects from AddHs() processing.',
};

// Per DESIGN.md "Warm-Status Rule": pass / warning / error stay in the warm
// spectrum (amber-gold → flame → fire). Differentiation comes from icon +
// label first, with hue and saturation reinforcing the gradient.
const CHECK_SEVERITY_STYLES: Record<string, string> = {
  pass: 'bg-[rgba(251,191,36,0.18)] border border-[rgba(251,191,36,0.35)] text-[#b45309] dark:text-[#fcd34d]',
  critical: 'bg-red-500/10 text-red-600 dark:text-red-400',
  error: 'bg-orange-500/10 text-orange-600 dark:text-orange-400',
  warning: 'bg-amber-500/10 text-amber-600 dark:text-amber-400',
  info: 'bg-sky-500/10 text-sky-600 dark:text-sky-400',
};

// Group the 27 individual checks into 7 chemistry-meaningful families so
// the expanded "All Checks" panel reads as scannable sections, not as a
// 27-card grid. Each category carries its own warm-tinted accent and
// lucide icon — restrained color (still well under the "10% of any
// screen" budget) that gives every section a visual key. Order follows
// reading flow: parsing first, identifiers last.
type CategoryAccent = {
  icon: LucideIcon;
  // Section header icon chip
  textClass: string;
  bgClass: string;
  // Section box (card-inset container)
  sectionBgClass: string;
  sectionBorderClass: string;
  sectionHoverBorderClass: string;
  // Pass-count chip when allPassed
  countBgClass: string;
  countTextClass: string;
};

const CHECK_CATEGORIES: ReadonlyArray<{
  key: string;
  label: string;
  members: ReadonlyArray<string>;
  accent: CategoryAccent;
}> = [
  {
    key: 'parsing',
    label: 'Parsing & sanitization',
    members: ['parsability', 'sanitization', 'valence', 'aromaticity', 'connectivity'],
    accent: {
      icon: FlaskConical,
      textClass: 'text-[#c41e3a] dark:text-[#f87171]',
      bgClass: 'bg-[rgba(196,30,58,0.12)]',
      sectionBgClass: 'bg-[rgba(196,30,58,0.03)] dark:bg-[rgba(248,113,113,0.04)]',
      sectionBorderClass: 'border-[rgba(196,30,58,0.18)] dark:border-[rgba(248,113,113,0.20)]',
      sectionHoverBorderClass: 'hover:border-[rgba(196,30,58,0.32)] dark:hover:border-[rgba(248,113,113,0.36)]',
      countBgClass: 'bg-[rgba(196,30,58,0.14)] dark:bg-[rgba(248,113,113,0.18)]',
      countTextClass: 'text-[#c41e3a] dark:text-[#f87171]',
    },
  },
  {
    key: 'atomic',
    label: 'Atomic composition',
    members: [
      'hypervalent_atoms',
      'charged_species',
      'explicit_hydrogen_audit',
      'radical_detection',
      'isotope_label_detection',
    ],
    accent: {
      icon: Atom,
      textClass: 'text-[#d97706] dark:text-[#fbbf24]',
      bgClass: 'bg-[rgba(217,119,6,0.14)]',
      sectionBgClass: 'bg-[rgba(217,119,6,0.04)] dark:bg-[rgba(251,191,36,0.05)]',
      sectionBorderClass: 'border-[rgba(217,119,6,0.20)] dark:border-[rgba(251,191,36,0.22)]',
      sectionHoverBorderClass: 'hover:border-[rgba(217,119,6,0.34)] dark:hover:border-[rgba(251,191,36,0.38)]',
      countBgClass: 'bg-[rgba(217,119,6,0.16)] dark:bg-[rgba(251,191,36,0.20)]',
      countTextClass: 'text-[#b45309] dark:text-[#fbbf24]',
    },
  },
  {
    key: 'topology',
    label: 'Topology & rings',
    members: ['polymer_detection', 'ring_strain', 'macrocycle_detection', 'aromatic_system_validation'],
    accent: {
      icon: Hexagon,
      textClass: 'text-[#b45309] dark:text-[#fcd34d]',
      bgClass: 'bg-[rgba(180,83,9,0.12)]',
      sectionBgClass: 'bg-[rgba(180,83,9,0.03)] dark:bg-[rgba(252,211,77,0.05)]',
      sectionBorderClass: 'border-[rgba(180,83,9,0.20)] dark:border-[rgba(252,211,77,0.20)]',
      sectionHoverBorderClass: 'hover:border-[rgba(180,83,9,0.34)] dark:hover:border-[rgba(252,211,77,0.36)]',
      countBgClass: 'bg-[rgba(180,83,9,0.16)] dark:bg-[rgba(252,211,77,0.20)]',
      countTextClass: 'text-[#b45309] dark:text-[#fcd34d]',
    },
  },
  {
    key: 'composition',
    label: 'Composition & contaminants',
    members: ['mixture_detection', 'solvent_contamination', 'inorganic_filter', 'trivial_molecule'],
    accent: {
      icon: Layers,
      textClass: 'text-[#ea580c] dark:text-[#fdba74]',
      bgClass: 'bg-[rgba(234,88,12,0.12)]',
      sectionBgClass: 'bg-[rgba(234,88,12,0.03)] dark:bg-[rgba(253,186,116,0.05)]',
      sectionBorderClass: 'border-[rgba(234,88,12,0.20)] dark:border-[rgba(253,186,116,0.20)]',
      sectionHoverBorderClass: 'hover:border-[rgba(234,88,12,0.34)] dark:hover:border-[rgba(253,186,116,0.36)]',
      countBgClass: 'bg-[rgba(234,88,12,0.16)] dark:bg-[rgba(253,186,116,0.20)]',
      countTextClass: 'text-[#ea580c] dark:text-[#fdba74]',
    },
  },
  {
    key: 'stereo',
    label: 'Stereochemistry',
    members: [
      'stereoisomer_enumeration',
      'undefined_stereocenters',
      'undefined_doublebond_stereo',
      'conflicting_stereo',
    ],
    accent: {
      icon: Share2,
      textClass: 'text-[#e11d48] dark:text-[#fb7185]',
      bgClass: 'bg-[rgba(225,29,72,0.12)]',
      sectionBgClass: 'bg-[rgba(225,29,72,0.03)] dark:bg-[rgba(251,113,133,0.05)]',
      sectionBorderClass: 'border-[rgba(225,29,72,0.18)] dark:border-[rgba(251,113,133,0.20)]',
      sectionHoverBorderClass: 'hover:border-[rgba(225,29,72,0.32)] dark:hover:border-[rgba(251,113,133,0.36)]',
      countBgClass: 'bg-[rgba(225,29,72,0.14)] dark:bg-[rgba(251,113,133,0.20)]',
      countTextClass: 'text-[#e11d48] dark:text-[#fb7185]',
    },
  },
  {
    key: 'tautomer',
    label: 'Tautomers & coordinates',
    members: ['tautomer_detection', 'coordinate_dimension'],
    accent: {
      icon: GitMerge,
      textClass: 'text-[#9d1830] dark:text-[#f87171]',
      bgClass: 'bg-[rgba(157,24,48,0.12)]',
      sectionBgClass: 'bg-[rgba(157,24,48,0.03)] dark:bg-[rgba(248,113,113,0.05)]',
      sectionBorderClass: 'border-[rgba(157,24,48,0.18)] dark:border-[rgba(248,113,113,0.20)]',
      sectionHoverBorderClass: 'hover:border-[rgba(157,24,48,0.32)] dark:hover:border-[rgba(248,113,113,0.36)]',
      countBgClass: 'bg-[rgba(157,24,48,0.14)] dark:bg-[rgba(248,113,113,0.18)]',
      countTextClass: 'text-[#9d1830] dark:text-[#f87171]',
    },
  },
  {
    key: 'identifier',
    label: 'Identifiers & roundtrips',
    members: ['smiles_roundtrip', 'inchi_generation', 'inchi_roundtrip'],
    accent: {
      icon: Lock,
      textClass: 'text-[#92400e] dark:text-[#fbbf24]',
      bgClass: 'bg-[rgba(146,64,14,0.12)]',
      sectionBgClass: 'bg-[rgba(146,64,14,0.03)] dark:bg-[rgba(251,191,36,0.05)]',
      sectionBorderClass: 'border-[rgba(146,64,14,0.20)] dark:border-[rgba(251,191,36,0.20)]',
      sectionHoverBorderClass: 'hover:border-[rgba(146,64,14,0.34)] dark:hover:border-[rgba(251,191,36,0.36)]',
      countBgClass: 'bg-[rgba(146,64,14,0.16)] dark:bg-[rgba(251,191,36,0.20)]',
      countTextClass: 'text-[#92400e] dark:text-[#fbbf24]',
    },
  },
];

const FALLBACK_CATEGORY_ACCENT: CategoryAccent = {
  icon: Layers,
  textClass: 'text-[var(--color-text-secondary)]',
  bgClass: 'bg-[var(--color-surface-sunken)]',
  sectionBgClass: 'bg-[var(--color-surface-sunken)]',
  sectionBorderClass: 'border-[var(--color-border)]/50',
  sectionHoverBorderClass: 'hover:border-[var(--color-border)]',
  countBgClass: 'bg-[var(--color-surface-sunken)]',
  countTextClass: 'text-[var(--color-text-secondary)]',
};

type AnyCheck = { check_name: string; passed: boolean; severity: string; message?: string | null };

function groupChecksByCategory<C extends AnyCheck>(checks: ReadonlyArray<C>) {
  const byName = new Map(checks.map((c) => [c.check_name, c] as const));
  const seen = new Set<string>();
  const groups: { label: string; items: C[]; accent: CategoryAccent }[] = [];
  for (const cat of CHECK_CATEGORIES) {
    const items: C[] = [];
    for (const name of cat.members) {
      const c = byName.get(name);
      if (c) {
        items.push(c);
        seen.add(name);
      }
    }
    if (items.length > 0) groups.push({ label: cat.label, items, accent: cat.accent });
  }
  const leftovers = checks.filter((c) => !seen.has(c.check_name));
  if (leftovers.length > 0)
    groups.push({ label: 'Other', items: leftovers, accent: FALLBACK_CATEGORY_ACCENT });
  return groups;
}

// Maximum columns per section at the current viewport. Keeping the cap
// in sync with the responsive sm/md/lg/xl breakpoints means a 5-item
// section can fill a single row at xl (5 cols allowed), while smaller
// viewports degrade gracefully without crushing card width.
function useChecksMaxCols(): number {
  const compute = () => {
    if (typeof window === 'undefined') return 5;
    const w = window.innerWidth;
    if (w >= 1280) return 5;
    if (w >= 1024) return 4;
    if (w >= 768) return 3;
    if (w >= 640) return 2;
    return 1;
  };
  const [cols, setCols] = useState(compute);
  useEffect(() => {
    if (typeof window === 'undefined') return;
    const handler = () => setCols(compute());
    handler();
    window.addEventListener('resize', handler);
    return () => window.removeEventListener('resize', handler);
  }, []);
  return cols;
}

/** Compact passed/flagged summary tags shown in the collapsed header. */
function renderAllChecksSummary(checks: CheckResult[]) {
  const passed = checks.filter((c) => c.passed).length;
  const flagged = checks.length - passed;
  return (
    <div className="flex items-center gap-1.5 ml-1">
      <span className="text-[10px] px-1.5 py-0.5 rounded bg-[rgba(251,191,36,0.18)] text-[#b45309] dark:text-[#fcd34d] font-medium">
        {passed} passed
      </span>
      {flagged > 0 && (
        <span className="text-[10px] px-1.5 py-0.5 rounded bg-amber-500/10 text-amber-600 dark:text-amber-400 font-medium">
          {flagged} flagged
        </span>
      )}
    </div>
  );
}

export type FloatRect = { top: number; left: number; width: number; height: number };
export type AllChecksBounds = { collapsed: FloatRect; expanded: FloatRect };
export type AllChecksPhase = 'collapsed' | 'expanding' | 'expanded' | 'collapsing';

interface AllChecksCardProps {
  /** Current animation phase of the floating card. */
  phase: AllChecksPhase;
  /** Measured collapsed/expanded anchor rectangles (relative to the page container). */
  bounds: AllChecksBounds;
  /** Natural expanded height (measured content height + header strip). */
  expandedHeight: number;
  /** The full set of validation checks to render in the categorized grid. */
  checks: CheckResult[];
  /** Callback ref attached to the inner content div so the parent's ResizeObserver can measure it. */
  contentRef: (node: HTMLDivElement | null) => void;
  /** Toggle handler — expands when collapsed, collapses when expanded. */
  onToggle: () => void;
}

/**
 * The single floating "All Checks" card. Rendered once as a position:absolute
 * motion.div that animates top/left/width/height between two layout anchors
 * (a collapsed anchor in the right column and an expanded anchor below the
 * grid) owned by the parent page. Same DOM element throughout — no
 * mount/unmount, no fade — sequenced so it MOVES first then GROWS on expand,
 * and SHRINKS first then MOVES on collapse.
 */
export function AllChecksCard({
  phase,
  bounds,
  expandedHeight,
  checks,
  contentRef,
  onToggle,
}: AllChecksCardProps) {
  const groups = useMemo(() => groupChecksByCategory(checks), [checks]);
  const maxCols = useChecksMaxCols();
  const expanded = phase === 'expanded' || phase === 'expanding';

  return (
    <motion.div
      className="absolute card overflow-hidden z-10"
      style={{
        marginTop: 0,
        boxShadow: expanded ? 'var(--shadow-lg)' : 'var(--shadow-sm)',
      }}
      initial={false}
      animate={
        expanded
          ? {
              top: bounds.expanded.top,
              left: bounds.expanded.left,
              width: bounds.expanded.width,
              height: bounds.expanded.top !== bounds.collapsed.top
                ? expandedHeight
                : bounds.collapsed.height,
            }
          : {
              top: bounds.collapsed.top,
              left: bounds.collapsed.left,
              width: bounds.collapsed.width,
              height: bounds.collapsed.height || 62,
            }
      }
      transition={
        phase === 'expanding'
          ? {
              // Phase 1 — move STRAIGHT DOWN: only `top` animates (350ms).
              // Phase 2 — grow LEFTWARD: `left` and `width` animate
              // together so the right edge stays anchored and the card
              // expands to the left. `height` reveals the body grid.
              top: { duration: 0.35, ease: [0.4, 0, 0.2, 1] },
              left: { duration: 0.4, ease: [0.4, 0, 0.2, 1], delay: 0.3 },
              width: { duration: 0.4, ease: [0.4, 0, 0.2, 1], delay: 0.3 },
              height: { duration: 0.4, ease: [0.4, 0, 0.2, 1], delay: 0.3 },
            }
          : phase === 'collapsing'
            ? {
                // Phase 1 — shrink RIGHTWARD: `left` slides right and
                // `width` shrinks together (right edge anchored).
                // Phase 2 — move STRAIGHT UP: `top` animates last.
                left: { duration: 0.4, ease: [0.4, 0, 0.2, 1] },
                width: { duration: 0.4, ease: [0.4, 0, 0.2, 1] },
                height: { duration: 0.4, ease: [0.4, 0, 0.2, 1] },
                top: { duration: 0.35, ease: [0.4, 0, 0.2, 1], delay: 0.35 },
              }
            : { duration: 0.2 }
      }
    >
      <button
        onClick={onToggle}
        aria-expanded={expanded}
        aria-controls="all-checks-grid"
        className="w-full flex items-center justify-between text-left px-5 py-4 sm:px-6 sm:py-5 hover:bg-[var(--color-surface-sunken)]/40 transition-colors"
      >
        <div className="flex items-center gap-3">
          <h4 className="font-semibold text-[var(--color-text-primary)] text-sm font-display">
            All Checks
          </h4>
          <span className="text-xs text-[var(--color-text-muted)]">
            {checks.length} checks
          </span>
          {renderAllChecksSummary(checks)}
        </div>
        <motion.div
          animate={{ rotate: expanded ? 180 : 0 }}
          transition={{ duration: 0.3 }}
        >
          <ChevronDown className="w-5 h-5 text-[var(--color-text-muted)]" />
        </motion.div>
      </button>

      <motion.div
        id="all-checks-grid"
        animate={{
          opacity: phase === 'expanded' ? 1 : 0,
        }}
        initial={false}
        transition={{
          duration: 0.25,
          delay: phase === 'expanded' ? 0.6 : 0,
        }}
        className="border-t border-[var(--color-border)]/40"
      >
        {/* Categorized card grid. Replaces the flat 27-card grid
             (the "Identical card grids" anti-pattern in DESIGN.md)
             by wrapping cards in 7 chemistry-meaningful sections,
             each with its own warm-tinted icon. Per-section column
             count = min(items.length, breakpointMax) so 5-item
             sections fill a single row at xl (no orphan cells), and
             small sections are capped at CARD_MAX width so cards
             don't balloon. ResizeObserver on this div drives the
             floating-card height — exact, no estimation. */}
        <div ref={contentRef} className="px-5 sm:px-6 pt-4 pb-5 sm:pb-6 space-y-4">
          {groups.map((group) => {
            const passedCount = group.items.filter((c) => c.passed).length;
            const allPassed = passedCount === group.items.length;
            const cols = Math.max(1, Math.min(group.items.length, maxCols));
            const CARD_MAX = 320;
            const GAP = 10;
            const sectionMaxWidth = cols * CARD_MAX + (cols - 1) * GAP;
            const Icon = group.accent.icon;
            return (
              <section
                key={group.label}
                className={cn(
                  'rounded-2xl p-4 sm:p-5 border transition-colors duration-300',
                  group.accent.sectionBgClass,
                  group.accent.sectionBorderClass,
                  group.accent.sectionHoverBorderClass,
                )}
                style={{ boxShadow: 'inset 0 1px 2px rgba(26, 24, 21, 0.03)' }}
              >
                <header className="flex items-center gap-3 mb-3.5">
                  <span
                    className={cn(
                      'flex items-center justify-center w-8 h-8 rounded-lg flex-shrink-0',
                      group.accent.bgClass,
                    )}
                  >
                    <Icon
                      className={cn('w-4 h-4', group.accent.textClass)}
                      strokeWidth={2.25}
                    />
                  </span>
                  <h5 className="text-[0.9375rem] font-semibold text-[var(--color-text-primary)] font-display tracking-tight whitespace-nowrap flex-1 min-w-0 truncate">
                    {group.label}
                  </h5>
                  <span
                    className={cn(
                      'text-xs font-semibold tabular-nums px-2.5 py-1 rounded-full whitespace-nowrap',
                      allPassed
                        ? cn(group.accent.countBgClass, group.accent.countTextClass)
                        : 'bg-[rgba(220,38,38,0.14)] text-[#dc2626] dark:bg-[rgba(248,113,113,0.20)] dark:text-[#fb7185]',
                    )}
                  >
                    {passedCount} / {group.items.length}
                  </span>
                </header>
                <div
                  className="grid auto-rows-fr"
                  style={{
                    gridTemplateColumns: `repeat(${cols}, minmax(0, 1fr))`,
                    gap: `${GAP}px`,
                    maxWidth: `${sectionMaxWidth}px`,
                  }}
                >
                  {group.items.map((check, index) => {
                    const severity = check.passed ? 'pass' : check.severity;
                    const severityClass =
                      CHECK_SEVERITY_STYLES[severity] ?? CHECK_SEVERITY_STYLES.info;
                    const description = CHECK_DESCRIPTIONS[check.check_name];
                    const prettyName = check.check_name.replace(/_/g, ' ');
                    return (
                      <div
                        key={`${check.check_name}-${index}`}
                        className={cn(
                          'rounded-lg p-3 border flex flex-col gap-1.5',
                          'transition-all duration-200',
                          'hover:-translate-y-0.5 hover:shadow-[0_4px_12px_rgba(26,24,21,0.08)]',
                          check.passed
                            ? 'bg-[var(--color-surface-elevated)] border-[var(--color-border)]/50'
                            : 'bg-[rgba(220,38,38,0.06)] dark:bg-[rgba(248,113,113,0.10)] border-[rgba(220,38,38,0.30)] dark:border-[rgba(248,113,113,0.35)]',
                        )}
                      >
                        <div className="flex items-center gap-1.5 min-w-0">
                          {check.passed ? (
                            <CheckCircle2
                              className="w-3.5 h-3.5 flex-shrink-0 text-[#d97706] dark:text-[#fbbf24]"
                              strokeWidth={2.25}
                            />
                          ) : (
                            <AlertTriangle
                              className="w-3.5 h-3.5 flex-shrink-0 text-[#dc2626] dark:text-[#fb7185]"
                              strokeWidth={2.5}
                            />
                          )}
                          <span
                            className="text-xs font-semibold text-[var(--color-text-primary)] truncate flex-1 min-w-0"
                            title={prettyName}
                          >
                            {prettyName}
                          </span>
                          <span
                            className={cn(
                              'text-[10px] px-1.5 py-0.5 rounded font-semibold tracking-wide flex-shrink-0',
                              severityClass,
                            )}
                          >
                            {check.passed ? 'PASS' : check.severity.toUpperCase()}
                          </span>
                        </div>
                        {description && (
                          <p className="text-[11px] text-[var(--color-text-secondary)] leading-snug line-clamp-3">
                            {description}
                          </p>
                        )}
                        {check.message && (
                          <p className="text-[11px] text-[var(--color-text-muted)] leading-snug line-clamp-2 font-medium mt-auto">
                            {check.message}
                          </p>
                        )}
                      </div>
                    );
                  })}
                </div>
              </section>
            );
          })}
        </div>
      </motion.div>
    </motion.div>
  );
}
