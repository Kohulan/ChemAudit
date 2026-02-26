import { useState } from 'react';
import { ChevronDown } from 'lucide-react';
import type { ValidationResponse } from '../../types/validation';
import { ScoreGauge } from './ScoreGauge';
import { IssueCard } from './IssueCard';
import { CopyButton } from '../ui/CopyButton';
import { Tooltip } from '../ui/Tooltip';
import { cn } from '../../lib/utils';

/** Tooltip descriptions for all validation checks */
const CHECK_DESCRIPTIONS: Record<string, string> = {
  // Basic checks
  parsability: 'Verifies the input string can be parsed into a valid molecular graph by RDKit.',
  sanitization: 'Runs RDKit sanitization to fix implicit valences, set aromaticity, and normalize the structure.',
  valence: 'Checks that every atom has a valid number of bonds according to its element type.',
  aromaticity: 'Validates aromatic ring systems can be properly kekulized (assigned alternating single/double bonds).',
  connectivity: 'Detects disconnected fragments and verifies the molecule is a single connected component.',
  // Stereo checks
  undefined_stereocenters: 'Identifies chiral centers (sp3 atoms with 4 different substituents) lacking R/S specification.',
  undefined_doublebond_stereo: 'Finds double bonds with two different substituents on each end that lack E/Z specification.',
  conflicting_stereo: 'Detects contradictory stereochemistry assignments that cannot coexist in a valid 3D structure.',
  // Representation checks
  smiles_roundtrip: 'Converts SMILES → molecule → SMILES and checks if the round-trip produces an identical canonical string.',
  inchi_generation: 'Tests whether a valid InChI identifier can be generated from the molecular structure.',
  inchi_roundtrip: 'Converts to InChI and back to a molecule, then compares with the original to detect information loss.',
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

interface ValidationResultsProps {
  result: ValidationResponse;
  onHighlightAtoms?: (atoms: number[]) => void;
  className?: string;
}

export function ValidationResults({
  result,
  onHighlightAtoms,
  className = '',
}: ValidationResultsProps) {
  const [showAllChecks, setShowAllChecks] = useState(false);

  const { molecule_info, overall_score, issues, all_checks, execution_time_ms } = result;

  return (
    <div className={cn('space-y-6', className)}>
      {/* Score Summary */}
      <div className="card p-6">
        <h3 className="text-lg font-semibold text-[var(--color-text-primary)] mb-4 text-center">
          Validation Score
        </h3>
        <ScoreGauge score={overall_score} size={140} className="mx-auto" />
        <div className="mt-4 text-center text-sm text-[var(--color-text-muted)]">
          Completed in {execution_time_ms.toFixed(0)}ms
        </div>
      </div>

      {/* Molecule Information */}
      <div className="card p-6">
        <h3 className="text-lg font-semibold text-[var(--color-text-primary)] mb-4">
          Molecule Information
        </h3>
        <div className="space-y-2 text-sm">
          {molecule_info.canonical_smiles && (
            <div className="flex items-start">
              <span className="font-medium text-[var(--color-text-secondary)] w-32 shrink-0">SMILES:</span>
              <code className="flex-1 text-[var(--color-text-secondary)] font-mono text-xs break-all">
                {molecule_info.canonical_smiles}
              </code>
              <CopyButton text={molecule_info.canonical_smiles} className="ml-2 shrink-0" />
            </div>
          )}
          {molecule_info.molecular_formula && (
            <div className="flex items-center">
              <span className="font-medium text-[var(--color-text-secondary)] w-32 shrink-0">Formula:</span>
              <span className="text-[var(--color-text-secondary)]">{molecule_info.molecular_formula}</span>
            </div>
          )}
          {molecule_info.molecular_weight && (
            <div className="flex items-center">
              <span className="font-medium text-[var(--color-text-secondary)] w-32 shrink-0">Mol. Weight:</span>
              <span className="text-[var(--color-text-secondary)]">
                {molecule_info.molecular_weight.toFixed(2)} g/mol
              </span>
            </div>
          )}
          {molecule_info.inchi && (
            <div className="flex items-start">
              <span className="font-medium text-[var(--color-text-secondary)] w-32 shrink-0">InChI:</span>
              <code className="flex-1 text-[var(--color-text-secondary)] font-mono text-xs break-all">
                {molecule_info.inchi}
              </code>
              <CopyButton text={molecule_info.inchi} className="ml-2 shrink-0" />
            </div>
          )}
          {molecule_info.inchikey && (
            <div className="flex items-start">
              <span className="font-medium text-[var(--color-text-secondary)] w-32 shrink-0">InChIKey:</span>
              <code className="flex-1 text-[var(--color-text-secondary)] font-mono text-xs break-all">
                {molecule_info.inchikey}
              </code>
              <CopyButton text={molecule_info.inchikey} className="ml-2 shrink-0" />
            </div>
          )}
          {molecule_info.num_atoms !== null && (
            <div className="flex items-center">
              <span className="font-medium text-[var(--color-text-secondary)] w-32 shrink-0">Atoms:</span>
              <span className="text-[var(--color-text-secondary)]">{molecule_info.num_atoms}</span>
            </div>
          )}
        </div>
      </div>

      {/* Issues */}
      {issues.length > 0 ? (
        <div className="card p-6">
          <h3 className="text-lg font-semibold text-[var(--color-text-primary)] mb-4">
            Issues Found ({issues.length})
          </h3>
          <div className="space-y-3">
            {issues.map((issue, index) => (
              <IssueCard
                key={`${issue.check_name}-${index}`}
                issue={issue}
                onAtomHover={onHighlightAtoms}
              />
            ))}
          </div>
        </div>
      ) : (
        <div className="rounded-xl p-6 text-center bg-yellow-500/10 border border-yellow-500/20">
          <div className="text-4xl mb-2">✓</div>
          <h3 className="text-lg font-semibold text-amber-600 dark:text-yellow-400 mb-1">
            No Issues Found
          </h3>
          <p className="text-sm text-amber-600/80 dark:text-yellow-400/80">
            All validation checks passed successfully
          </p>
        </div>
      )}

      {/* All Checks (collapsible) */}
      <div className="card p-6">
        <button
          onClick={() => setShowAllChecks(!showAllChecks)}
          className="w-full flex items-center justify-between text-left"
        >
          <h3 className="text-lg font-semibold text-[var(--color-text-primary)]">
            All Checks ({all_checks.length})
          </h3>
          <ChevronDown
            className={cn(
              'w-5 h-5 text-[var(--color-text-muted)] transition-transform',
              showAllChecks && 'rotate-180'
            )}
          />
        </button>

        {showAllChecks && (
          <div className="mt-4 space-y-2">
            {all_checks.map((check, index) => (
              <div
                key={`${check.check_name}-${index}`}
                className="flex items-center justify-between py-2 px-3 bg-[var(--color-surface-sunken)] rounded-lg"
              >
                <div className="flex items-center gap-2">
                  <span className={check.passed ? 'text-amber-500 dark:text-yellow-400' : 'text-red-500'}>
                    {check.passed ? '✓' : '✗'}
                  </span>
                  {CHECK_DESCRIPTIONS[check.check_name] ? (
                    <Tooltip
                      content={CHECK_DESCRIPTIONS[check.check_name]}
                      position="top"
                      maxWidth={300}
                    >
                      <span className="text-sm font-medium text-[var(--color-text-primary)] cursor-help border-b border-dashed border-[var(--color-text-muted)]/40 hover:border-[var(--color-primary)] transition-colors">
                        {check.check_name.replace(/_/g, ' ')}
                      </span>
                    </Tooltip>
                  ) : (
                    <span className="text-sm font-medium text-[var(--color-text-primary)]">
                      {check.check_name.replace(/_/g, ' ')}
                    </span>
                  )}
                </div>
                <span
                  className={cn(
                    'text-xs px-2 py-1 rounded-md font-medium',
                    check.passed
                      ? 'bg-yellow-500/10 text-amber-600 dark:text-yellow-400'
                      : check.severity === 'critical'
                      ? 'bg-red-500/10 text-red-600 dark:text-red-400'
                      : check.severity === 'error'
                      ? 'bg-orange-500/10 text-orange-600 dark:text-orange-400'
                      : check.severity === 'warning'
                      ? 'bg-amber-500/10 text-amber-600 dark:text-amber-400'
                      : 'bg-sky-500/10 text-sky-600 dark:text-sky-400'
                  )}
                >
                  {check.passed ? 'PASS' : check.severity.toUpperCase()}
                </span>
              </div>
            ))}
          </div>
        )}
      </div>
    </div>
  );
}
