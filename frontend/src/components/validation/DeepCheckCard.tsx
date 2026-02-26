/**
 * DeepCheckCard Component
 *
 * Renders a single deep validation check result with:
 * - Severity badge (reflects user-overridden effective severity)
 * - Human-readable check name and message
 * - Clickable atom indices for molecule viewer cross-linking
 * - Domain-specific structured detail rendering (stereoisomers, fragments, etc.)
 * - Generic key-value fallback for other checks
 * - Framer Motion collapse/expand animation
 */
import { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown, CheckCircle2, XCircle } from 'lucide-react';
import { Badge } from '../ui/Badge';
import { InfoTooltip } from '../ui/Tooltip';
import { StereoisomerList } from './StereoisomerList';
import { FragmentClassificationTable } from './FragmentClassificationTable';
import { cn } from '../../lib/utils';

/** Tooltip descriptions for all 16 deep validation checks */
const DEEP_CHECK_DESCRIPTIONS: Record<string, string> = {
  // Stereo & Tautomers
  stereoisomer_enumeration: 'Finds undefined stereocenters and enumerates possible stereoisomers (up to 128). Helps identify ambiguous chirality.',
  tautomer_detection: 'Detects tautomeric forms and identifies the canonical tautomer. Reports whether the input matches the canonical form.',
  aromatic_system_validation: 'Checks for unusual aromatic ring sizes (not 5 or 6 membered) and charged aromatic atoms that may indicate issues.',
  coordinate_dimension: 'Reports whether the molecule has 2D coordinates, 3D coordinates, or no coordinate information.',
  // Chemical Composition
  mixture_detection: 'Detects multi-fragment inputs (dot-separated SMILES) and classifies each fragment as drug, salt, solvent, or unknown.',
  solvent_contamination: 'Screens for common lab solvents (water, DMSO, DMF, methanol, etc.) that may contaminate the input structure.',
  inorganic_filter: 'Flags molecules lacking carbon (inorganic) or containing metal atoms (organometallic) that may not suit standard validation.',
  radical_detection: 'Identifies atoms with unpaired radical electrons that may indicate unstable or reactive species.',
  isotope_label_detection: 'Detects isotope-labeled atoms (deuterium, carbon-13, etc.) often used in pharmacokinetic studies.',
  trivial_molecule: 'Flags molecules with 3 or fewer heavy atoms as too small for meaningful chemical validation.',
  // Structural Complexity
  hypervalent_atoms: 'Detects atoms exceeding their normal valence limits, which may indicate unusual bonding or input errors.',
  polymer_detection: 'Identifies possible polymers via SGroup markers, molecular weight above 1500 Da, or dummy atom attachment points.',
  ring_strain: 'Flags 3-membered (cyclopropane) and 4-membered (cyclobutane) rings that have significant ring strain.',
  macrocycle_detection: 'Identifies macrocyclic rings with more than 12 atoms, common in natural products and cyclic peptides.',
  charged_species: 'Reports formal charges, identifies zwitterions (net charge zero with both positive and negative atoms).',
  explicit_hydrogen_audit: 'Reports atoms with explicit hydrogen specifications and detects H atom objects from AddHs() processing.',
};
import type {
  CheckResult,
  StereoisomerDetail,
  TautomerDetail,
  MixtureDetail,
  ChargedSpeciesDetail,
  RingStrainDetail,
  MacrocycleDetail,
  AromaticSystemDetail,
  SolventDetail,
  InorganicDetail,
  RadicalDetail,
  IsotopeDetail,
  TrivialMoleculeDetail,
  HypervalentDetail,
  PolymerDetail,
  CoordinateDimensionDetail,
  ExplicitHydrogenDetail,
} from '../../types/validation';

interface DeepCheckCardProps {
  check: CheckResult;
  effectiveSeverity: string;
  onHighlightAtoms: (atoms: number[]) => void;
}

/** Maps severity string to Badge variant */
function getSeverityVariant(
  severity: string
): 'error' | 'warning' | 'info' | 'success' | 'default' {
  switch (severity) {
    case 'critical':
    case 'error':
      return 'error';
    case 'warning':
      return 'warning';
    case 'info':
      return 'info';
    case 'pass':
      return 'success';
    default:
      return 'default';
  }
}

/** Maps severity string to icon color */
function getSeverityBg(severity: string): string {
  switch (severity) {
    case 'critical':
    case 'error':
      return 'bg-red-500/10 border-red-500/20';
    case 'warning':
      return 'bg-amber-500/10 border-amber-500/20';
    case 'info':
      return 'bg-sky-500/10 border-sky-500/20';
    default:
      return 'bg-[var(--color-surface-elevated)] border-[var(--color-border)]';
  }
}

/** Converts snake_case to Title Case */
function formatCheckName(name: string): string {
  return name
    .split('_')
    .map((word) => word.charAt(0).toUpperCase() + word.slice(1))
    .join(' ');
}

/** Clickable atom index badge */
function AtomBadge({
  atomIdx,
  onClick,
}: {
  atomIdx: number;
  onClick: () => void;
}) {
  return (
    <button
      onClick={onClick}
      className={cn(
        'inline-flex items-center px-1.5 py-0.5 rounded text-[10px] font-mono font-medium',
        'bg-[var(--color-primary)]/10 text-[var(--color-primary)]',
        'border border-[var(--color-primary)]/20',
        'hover:bg-[var(--color-primary)]/20 transition-colors cursor-pointer',
        'mr-1 mb-1'
      )}
    >
      #{atomIdx}
    </button>
  );
}

/** Generic key-value detail display */
function GenericDetails({ details }: { details: Record<string, unknown> }) {
  const entries = Object.entries(details);
  if (entries.length === 0) return null;

  return (
    <div className="mt-2 space-y-1">
      {entries.map(([key, value]) => {
        let displayValue: string;
        if (value === null || value === undefined) {
          displayValue = '—';
        } else if (typeof value === 'boolean') {
          displayValue = value ? 'Yes' : 'No';
        } else if (Array.isArray(value)) {
          displayValue = value.length > 0 ? `[${value.length} items]` : '[]';
        } else if (typeof value === 'object') {
          displayValue = JSON.stringify(value);
        } else {
          displayValue = String(value);
        }
        return (
          <div key={key} className="flex items-baseline gap-2 text-xs">
            <span className="text-[var(--color-text-muted)] w-40 flex-shrink-0 font-medium">
              {formatCheckName(key)}:
            </span>
            <span className="text-[var(--color-text-secondary)] break-all">{displayValue}</span>
          </div>
        );
      })}
    </div>
  );
}

/** Render check-specific structured detail section */
function CheckDetails({
  check,
  onHighlightAtoms,
}: {
  check: CheckResult;
  onHighlightAtoms: (atoms: number[]) => void;
}) {
  const details = check.details as Record<string, unknown>;

  switch (check.check_name) {
    case 'stereoisomer_enumeration': {
      const d = details as unknown as StereoisomerDetail;
      return (
        <div className="mt-2 text-xs text-[var(--color-text-secondary)]">
          <div className="flex gap-4 flex-wrap mb-1">
            <span>
              <span className="text-[var(--color-text-muted)]">Undefined centers: </span>
              <strong>{d.undefined_count}</strong>
            </span>
            <span>
              <span className="text-[var(--color-text-muted)]">Total centers: </span>
              <strong>{d.total_centers}</strong>
            </span>
          </div>
          <StereoisomerList
            smiles={d.stereoisomer_smiles ?? []}
            cap={d.enumeration_cap ?? 128}
            capExceeded={d.cap_exceeded ?? false}
          />
        </div>
      );
    }

    case 'mixture_detection': {
      const d = details as unknown as MixtureDetail;
      return (
        <div className="mt-2 text-xs text-[var(--color-text-secondary)]">
          <p>
            <span className="text-[var(--color-text-muted)]">Fragments found: </span>
            <strong>{d.num_fragments}</strong>
          </p>
          {d.fragments && d.fragments.length > 0 && (
            <FragmentClassificationTable fragments={d.fragments} />
          )}
        </div>
      );
    }

    case 'tautomer_detection': {
      const d = details as unknown as TautomerDetail;
      return (
        <div className="mt-2 text-xs text-[var(--color-text-secondary)] space-y-1">
          <div className="flex gap-4 flex-wrap">
            <span>
              <span className="text-[var(--color-text-muted)]">Tautomers: </span>
              <strong>{d.tautomer_count}</strong>
            </span>
            <span>
              <span className="text-[var(--color-text-muted)]">Is canonical: </span>
              <strong>{d.is_canonical_form ? 'Yes' : 'No'}</strong>
            </span>
          </div>
          {d.canonical_smiles && (
            <div>
              <span className="text-[var(--color-text-muted)]">Canonical SMILES: </span>
              <code className="font-mono text-[var(--color-text-secondary)] break-all">
                {d.canonical_smiles}
              </code>
            </div>
          )}
        </div>
      );
    }

    case 'charged_species': {
      const d = details as unknown as ChargedSpeciesDetail;
      return (
        <div className="mt-2 text-xs text-[var(--color-text-secondary)] space-y-1">
          <div className="flex gap-4 flex-wrap">
            <span>
              <span className="text-[var(--color-text-muted)]">Net charge: </span>
              <strong>{d.net_charge >= 0 ? `+${d.net_charge}` : d.net_charge}</strong>
            </span>
            <span>
              <span className="text-[var(--color-text-muted)]">Zwitterion: </span>
              <strong>{d.is_zwitterion ? 'Yes' : 'No'}</strong>
            </span>
            <span>
              <span className="text-[var(--color-text-muted)]">Charged atoms: </span>
              <strong>{d.total_charged_atoms}</strong>
            </span>
          </div>
          {d.positive_atoms && d.positive_atoms.length > 0 && (
            <div>
              <span className="text-[var(--color-text-muted)]">Positive: </span>
              {d.positive_atoms.map((a) => (
                <AtomBadge
                  key={a.atom_idx}
                  atomIdx={a.atom_idx}
                  onClick={() => onHighlightAtoms([a.atom_idx])}
                />
              ))}
            </div>
          )}
          {d.negative_atoms && d.negative_atoms.length > 0 && (
            <div>
              <span className="text-[var(--color-text-muted)]">Negative: </span>
              {d.negative_atoms.map((a) => (
                <AtomBadge
                  key={a.atom_idx}
                  atomIdx={a.atom_idx}
                  onClick={() => onHighlightAtoms([a.atom_idx])}
                />
              ))}
            </div>
          )}
        </div>
      );
    }

    case 'ring_strain': {
      const d = details as unknown as RingStrainDetail;
      return (
        <div className="mt-2 text-xs text-[var(--color-text-secondary)] space-y-1">
          <span>
            <span className="text-[var(--color-text-muted)]">Strained rings: </span>
            <strong>{d.total_strained_rings}</strong>
          </span>
          {d.strained_rings &&
            d.strained_rings.map((ring, idx) => (
              <div key={idx} className="flex items-center gap-2 flex-wrap">
                <span className="text-[var(--color-text-muted)]">
                  {ring.ring_size}-membered ring:
                </span>
                {ring.atom_indices.map((atomIdx) => (
                  <AtomBadge
                    key={atomIdx}
                    atomIdx={atomIdx}
                    onClick={() => onHighlightAtoms([atomIdx])}
                  />
                ))}
                <button
                  className="text-[10px] text-[var(--color-primary)] hover:underline ml-1"
                  onClick={() => onHighlightAtoms(ring.atom_indices)}
                >
                  highlight all
                </button>
              </div>
            ))}
        </div>
      );
    }

    case 'macrocycle_detection': {
      const d = details as unknown as MacrocycleDetail;
      return (
        <div className="mt-2 text-xs text-[var(--color-text-secondary)] space-y-1">
          <div className="flex gap-4 flex-wrap">
            <span>
              <span className="text-[var(--color-text-muted)]">Macrocycles: </span>
              <strong>{d.total_macrocycles}</strong>
            </span>
          </div>
          {d.macrocycles &&
            d.macrocycles.map((mc, idx) => (
              <div key={idx} className="flex items-center gap-2 flex-wrap">
                <span className="text-[var(--color-text-muted)]">
                  {mc.ring_size}-membered macrocycle:
                </span>
                {mc.atom_indices.slice(0, 8).map((atomIdx) => (
                  <AtomBadge
                    key={atomIdx}
                    atomIdx={atomIdx}
                    onClick={() => onHighlightAtoms([atomIdx])}
                  />
                ))}
                {mc.atom_indices.length > 8 && (
                  <span className="text-[var(--color-text-muted)]">+{mc.atom_indices.length - 8} more</span>
                )}
                <button
                  className="text-[10px] text-[var(--color-primary)] hover:underline ml-1"
                  onClick={() => onHighlightAtoms(mc.atom_indices)}
                >
                  highlight all
                </button>
              </div>
            ))}
          {d.sssr_note && (
            <p className="text-[var(--color-text-muted)] text-[10px] italic">{d.sssr_note}</p>
          )}
        </div>
      );
    }

    case 'aromatic_system_validation': {
      const d = details as unknown as AromaticSystemDetail;
      return (
        <div className="mt-2 text-xs text-[var(--color-text-secondary)] space-y-1">
          {d.unusual_ring_sizes && d.unusual_ring_sizes.length > 0 && (
            <div>
              <span className="text-[var(--color-text-muted)]">Unusual ring sizes: </span>
              {d.unusual_ring_sizes.map((ring, idx) => (
                <span key={idx}>
                  <span className="font-medium">{ring.ring_size}-membered</span>
                  {' '}
                  <button
                    className="text-[10px] text-[var(--color-primary)] hover:underline"
                    onClick={() => onHighlightAtoms(ring.atom_indices)}
                  >
                    highlight
                  </button>{' '}
                </span>
              ))}
            </div>
          )}
          {d.charged_aromatics && d.charged_aromatics.length > 0 && (
            <div className="flex items-center gap-1 flex-wrap">
              <span className="text-[var(--color-text-muted)]">Charged aromatic atoms: </span>
              {d.charged_aromatics.map((a) => (
                <AtomBadge
                  key={a.atom_idx}
                  atomIdx={a.atom_idx}
                  onClick={() => onHighlightAtoms([a.atom_idx])}
                />
              ))}
            </div>
          )}
        </div>
      );
    }

    case 'solvent_contamination': {
      const d = details as unknown as SolventDetail;
      return (
        <div className="mt-2 text-xs text-[var(--color-text-secondary)] space-y-1">
          <span>
            <span className="text-[var(--color-text-muted)]">Pure solvent: </span>
            <strong>{d.is_pure_solvent ? 'Yes' : 'No'}</strong>
          </span>
          {d.solvents_found && d.solvents_found.length > 0 && (
            <ul className="mt-1 space-y-0.5">
              {d.solvents_found.map((s, idx) => (
                <li key={idx} className="flex gap-2">
                  <span className="font-medium">{s.name}</span>
                  <span className="text-[var(--color-text-muted)]">
                    ({s.molecular_weight.toFixed(1)} g/mol)
                  </span>
                </li>
              ))}
            </ul>
          )}
        </div>
      );
    }

    case 'inorganic_filter': {
      const d = details as unknown as InorganicDetail;
      return (
        <div className="mt-2 text-xs text-[var(--color-text-secondary)] space-y-1">
          <div className="flex gap-4 flex-wrap">
            <span>
              <span className="text-[var(--color-text-muted)]">Has carbon: </span>
              <strong>{d.has_carbon ? 'Yes' : 'No'}</strong>
            </span>
            <span>
              <span className="text-[var(--color-text-muted)]">Inorganic: </span>
              <strong>{d.is_inorganic ? 'Yes' : 'No'}</strong>
            </span>
            <span>
              <span className="text-[var(--color-text-muted)]">Organometallic: </span>
              <strong>{d.is_organometallic ? 'Yes' : 'No'}</strong>
            </span>
          </div>
          {d.metal_atoms && d.metal_atoms.length > 0 && (
            <div className="flex items-center gap-1 flex-wrap">
              <span className="text-[var(--color-text-muted)]">Metal atoms: </span>
              {d.metal_atoms.map((a) => (
                <AtomBadge
                  key={a.atom_idx}
                  atomIdx={a.atom_idx}
                  onClick={() => onHighlightAtoms([a.atom_idx])}
                />
              ))}
            </div>
          )}
        </div>
      );
    }

    case 'radical_detection': {
      const d = details as unknown as RadicalDetail;
      return (
        <div className="mt-2 text-xs text-[var(--color-text-secondary)] space-y-1">
          <span>
            <span className="text-[var(--color-text-muted)]">Total radical electrons: </span>
            <strong>{d.total_radical_electrons}</strong>
          </span>
          {d.radical_atoms && d.radical_atoms.length > 0 && (
            <div className="flex items-center gap-1 flex-wrap">
              <span className="text-[var(--color-text-muted)]">Radical atoms: </span>
              {d.radical_atoms.map((a) => (
                <AtomBadge
                  key={a.atom_idx}
                  atomIdx={a.atom_idx}
                  onClick={() => onHighlightAtoms([a.atom_idx])}
                />
              ))}
            </div>
          )}
        </div>
      );
    }

    case 'isotope_label_detection': {
      const d = details as unknown as IsotopeDetail;
      return (
        <div className="mt-2 text-xs text-[var(--color-text-secondary)] space-y-1">
          <span>
            <span className="text-[var(--color-text-muted)]">Labeled atoms: </span>
            <strong>{d.total_labeled}</strong>
          </span>
          {d.labeled_atoms && d.labeled_atoms.length > 0 && (
            <div className="flex items-center gap-1 flex-wrap">
              {d.labeled_atoms.map((a) => (
                <span key={a.atom_idx} className="flex items-center gap-1">
                  <AtomBadge
                    atomIdx={a.atom_idx}
                    onClick={() => onHighlightAtoms([a.atom_idx])}
                  />
                  {a.common_name && (
                    <span className="text-[var(--color-text-muted)] text-[10px]">
                      ({a.common_name})
                    </span>
                  )}
                </span>
              ))}
            </div>
          )}
        </div>
      );
    }

    case 'trivial_molecule': {
      const d = details as unknown as TrivialMoleculeDetail;
      return (
        <div className="mt-2 text-xs text-[var(--color-text-secondary)] space-y-1">
          <div className="flex gap-4 flex-wrap">
            <span>
              <span className="text-[var(--color-text-muted)]">Heavy atoms: </span>
              <strong>{d.heavy_atom_count}</strong>
            </span>
            <span>
              <span className="text-[var(--color-text-muted)]">Threshold: </span>
              <strong>&le; {d.threshold}</strong>
            </span>
            <span>
              <span className="text-[var(--color-text-muted)]">Single atom: </span>
              <strong>{d.is_single_atom ? 'Yes' : 'No'}</strong>
            </span>
          </div>
        </div>
      );
    }

    case 'hypervalent_atoms': {
      const d = details as unknown as HypervalentDetail;
      return (
        <div className="mt-2 text-xs text-[var(--color-text-secondary)] space-y-1">
          {d.hypervalent_atoms && d.hypervalent_atoms.length > 0 ? (
            <div className="space-y-1">
              {d.hypervalent_atoms.map((a) => (
                <div key={a.atom_idx} className="flex items-center gap-2 flex-wrap">
                  <AtomBadge
                    atomIdx={a.atom_idx}
                    onClick={() => onHighlightAtoms([a.atom_idx])}
                  />
                  <span className="font-medium">{a.symbol}</span>
                  <span className="text-[var(--color-text-muted)]">
                    valence {a.actual_valence} (allowed: {a.allowed_valences.join(', ')})
                  </span>
                </div>
              ))}
            </div>
          ) : (
            <span className="text-[var(--color-text-muted)]">No hypervalent atoms detected.</span>
          )}
        </div>
      );
    }

    case 'polymer_detection': {
      const d = details as unknown as PolymerDetail;
      return (
        <div className="mt-2 text-xs text-[var(--color-text-secondary)] space-y-1">
          <div className="flex gap-4 flex-wrap">
            <span>
              <span className="text-[var(--color-text-muted)]">SGroup markers: </span>
              <strong>{d.has_sgroup_markers ? 'Yes' : 'No'}</strong>
            </span>
            <span>
              <span className="text-[var(--color-text-muted)]">MW: </span>
              <strong>{d.molecular_weight.toFixed(1)}</strong>
            </span>
            <span>
              <span className="text-[var(--color-text-muted)]">Exceeds MW: </span>
              <strong>{d.exceeds_mw_threshold ? 'Yes' : 'No'}</strong>
            </span>
            <span>
              <span className="text-[var(--color-text-muted)]">Dummy atoms: </span>
              <strong>{d.dummy_atom_count}</strong>
            </span>
          </div>
          {d.sgroup_types && d.sgroup_types.length > 0 && (
            <div>
              <span className="text-[var(--color-text-muted)]">SGroup types: </span>
              <span>{d.sgroup_types.join(', ')}</span>
            </div>
          )}
        </div>
      );
    }

    case 'coordinate_dimension': {
      const d = details as unknown as CoordinateDimensionDetail;
      return (
        <div className="mt-2 text-xs text-[var(--color-text-secondary)]">
          <div className="flex gap-4 flex-wrap">
            <span>
              <span className="text-[var(--color-text-muted)]">Dimension: </span>
              <strong>{d.dimension}</strong>
            </span>
            <span>
              <span className="text-[var(--color-text-muted)]">Conformers: </span>
              <strong>{d.num_conformers}</strong>
            </span>
          </div>
        </div>
      );
    }

    case 'explicit_hydrogen_audit': {
      const d = details as unknown as ExplicitHydrogenDetail;
      return (
        <div className="mt-2 text-xs text-[var(--color-text-secondary)] space-y-1">
          <div className="flex gap-4 flex-wrap">
            <span>
              <span className="text-[var(--color-text-muted)]">Explicit H count: </span>
              <strong>{d.total_explicit_h}</strong>
            </span>
            <span>
              <span className="text-[var(--color-text-muted)]">H atom objects: </span>
              <strong>{d.h_atom_object_count}</strong>
            </span>
          </div>
          {d.atoms_with_explicit_h && d.atoms_with_explicit_h.length > 0 && (
            <div className="flex items-center gap-1 flex-wrap">
              <span className="text-[var(--color-text-muted)]">Atoms with explicit H: </span>
              {d.atoms_with_explicit_h.map((a) => (
                <AtomBadge
                  key={a.atom_idx}
                  atomIdx={a.atom_idx}
                  onClick={() => onHighlightAtoms([a.atom_idx])}
                />
              ))}
            </div>
          )}
        </div>
      );
    }

    default:
      return <GenericDetails details={details} />;
  }
}

/**
 * Renders a single deep validation check result card.
 */
export function DeepCheckCard({ check, effectiveSeverity, onHighlightAtoms }: DeepCheckCardProps) {
  const [detailsOpen, setDetailsOpen] = useState(false);

  const isOverridden = effectiveSeverity !== check.severity;
  const hasDetails = Object.keys(check.details).length > 0;
  const hasAffectedAtoms = check.affected_atoms.length > 0;

  // Determine card background based on effective severity (not original)
  const cardBg = check.passed
    ? 'bg-[var(--color-surface-elevated)] border-[var(--color-border)]'
    : getSeverityBg(effectiveSeverity);

  return (
    <div
      className={cn(
        'rounded-xl border p-4 transition-all',
        cardBg
      )}
    >
      {/* Header */}
      <div className="flex items-start gap-3">
        {/* Pass/Fail icon */}
        <div className="flex-shrink-0 mt-0.5">
          {check.passed ? (
            <CheckCircle2 className="w-4 h-4 text-emerald-500" />
          ) : (
            <XCircle className="w-4 h-4 text-red-500" />
          )}
        </div>

        {/* Content */}
        <div className="flex-1 min-w-0">
          <div className="flex flex-wrap items-center gap-2 mb-1">
            <span className="text-sm font-semibold text-[var(--color-text-primary)]">
              {formatCheckName(check.check_name)}
            </span>

            {/* Check description tooltip */}
            {DEEP_CHECK_DESCRIPTIONS[check.check_name] && (
              <InfoTooltip
                content={DEEP_CHECK_DESCRIPTIONS[check.check_name]}
                position="right"
              />
            )}

            {/* Effective severity badge */}
            <Badge variant={getSeverityVariant(effectiveSeverity)} size="sm">
              {effectiveSeverity.toUpperCase()}
            </Badge>

            {/* Show original if overridden */}
            {isOverridden && (
              <span className="text-[10px] text-[var(--color-text-muted)] italic">
                (orig: {check.severity})
              </span>
            )}

            {/* PASS badge */}
            {check.passed && (
              <Badge variant="success" size="sm">
                PASS
              </Badge>
            )}
          </div>

          {/* Message */}
          <p className="text-sm text-[var(--color-text-secondary)]">{check.message}</p>

          {/* Affected atoms */}
          {hasAffectedAtoms && (
            <div className="mt-2">
              <div className="flex items-center gap-1 flex-wrap">
                <span className="text-[10px] text-[var(--color-text-muted)] mr-1">Affected atoms:</span>
                {check.affected_atoms.map((atomIdx) => (
                  <AtomBadge
                    key={atomIdx}
                    atomIdx={atomIdx}
                    onClick={() => onHighlightAtoms([atomIdx])}
                  />
                ))}
                <button
                  onClick={() => onHighlightAtoms(check.affected_atoms)}
                  className="text-[10px] text-[var(--color-primary)] hover:underline ml-1"
                >
                  highlight all
                </button>
              </div>
            </div>
          )}

          {/* Details section toggle */}
          {hasDetails && (
            <div className="mt-2">
              <button
                onClick={() => setDetailsOpen((prev) => !prev)}
                className={cn(
                  'flex items-center gap-1 text-[10px] font-medium',
                  'text-[var(--color-text-muted)] hover:text-[var(--color-text-secondary)]',
                  'transition-colors'
                )}
              >
                <ChevronDown
                  className={cn(
                    'w-3 h-3 transition-transform duration-200',
                    detailsOpen && 'rotate-180'
                  )}
                />
                {detailsOpen ? 'Hide details' : 'Show details'}
              </button>

              <AnimatePresence>
                {detailsOpen && (
                  <motion.div
                    initial={{ opacity: 0, height: 0 }}
                    animate={{ opacity: 1, height: 'auto' }}
                    exit={{ opacity: 0, height: 0 }}
                    transition={{ duration: 0.2 }}
                    className="overflow-hidden"
                  >
                    <div className="mt-2 p-3 bg-[var(--color-surface-sunken)] rounded-lg border border-[var(--color-border)]">
                      <CheckDetails check={check} onHighlightAtoms={onHighlightAtoms} />
                    </div>
                  </motion.div>
                )}
              </AnimatePresence>
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
