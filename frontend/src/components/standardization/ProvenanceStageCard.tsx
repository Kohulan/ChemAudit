/**
 * ProvenanceStageCard Component
 *
 * Expandable card for a single provenance pipeline stage.
 * Shows before/after SMILES, status icon, change count badge, and
 * expandable detail sections for charge, bond, ring, radical changes,
 * fragment removals, and DVAL cross-references.
 *
 * On-demand structure rendering via "Show structure" button (not eagerly loaded).
 * Color-coded atom highlighting: blue for charge changes, green for bond/ring changes.
 */
import { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { CheckCircle2, MinusCircle, XCircle, ChevronDown, ChevronRight } from 'lucide-react';
import { MoleculeViewer } from '../molecules/MoleculeViewer';
import type { ProvStageRecord } from '../../types/standardization';

interface ProvenanceStageCardProps {
  stage: ProvStageRecord;
  isExpanded: boolean;
  onToggle: () => void;
}

/** Format pipeline stage name for display. */
function formatStageName(name: string): string {
  switch (name) {
    case 'checker': return 'Checker';
    case 'standardizer': return 'Standardizer';
    case 'get_parent': return 'Parent Extraction';
    case 'tautomer_canonicalization': return 'Tautomer Canonicalization';
    default:
      return name.split('_').map((w) => w.charAt(0).toUpperCase() + w.slice(1)).join(' ');
  }
}

/** Count total notable changes in a stage record. */
function countChanges(stage: ProvStageRecord): number {
  return (
    stage.charge_changes.length +
    stage.bond_changes.length +
    stage.ring_changes.length +
    stage.radical_changes.length +
    stage.fragment_removals.length +
    stage.dval_cross_refs.length
  );
}

/** Role badge for fragment removals. */
function RoleBadge({ role }: { role: string }) {
  const colors: Record<string, string> = {
    salt: 'bg-blue-100 text-blue-700',
    solvent: 'bg-yellow-100 text-yellow-700',
    counterion: 'bg-purple-100 text-purple-700',
    unknown: 'bg-gray-100 text-gray-600',
  };
  return (
    <span className={`inline-flex px-1.5 py-0.5 rounded text-xs font-medium ${colors[role] ?? colors.unknown}`}>
      {role}
    </span>
  );
}

export function ProvenanceStageCard({ stage, isExpanded, onToggle }: ProvenanceStageCardProps) {
  const [showStructure, setShowStructure] = useState(false);

  const changeCount = countChanges(stage);
  const hasChanges = changeCount > 0;
  const smilesDiffer = stage.input_smiles !== stage.output_smiles;

  // Build highlight atom lists for "Show structure"
  const chargeAtomIndices = stage.charge_changes.map((c) => c.atom_idx);
  const bondAtomIndices = stage.bond_changes.flatMap((b) => [b.atom1_idx, b.atom2_idx]);
  const ringAtomIndices = stage.ring_changes.flatMap((r) => r.ring_atoms);
  // Combine all highlight atoms (blue for charge, green for bond/ring — MoleculeViewer
  // uses a single list, so we pass all affected atom indices)
  const highlightAtoms = [
    ...new Set([...chargeAtomIndices, ...bondAtomIndices, ...ringAtomIndices]),
  ];

  // Status icon
  let StatusIcon: typeof CheckCircle2 | typeof MinusCircle | typeof XCircle;
  let iconColor: string;
  if (!stage.applied) {
    StatusIcon = XCircle;
    iconColor = 'text-gray-400';
  } else if (hasChanges) {
    StatusIcon = CheckCircle2;
    iconColor = 'text-green-500';
  } else {
    StatusIcon = MinusCircle;
    iconColor = 'text-gray-400';
  }

  return (
    <div className="border border-gray-200 rounded-lg overflow-hidden bg-white">
      {/* Header */}
      <button
        onClick={onToggle}
        className="w-full flex items-center gap-3 px-4 py-3 hover:bg-gray-50 transition-colors text-left"
        aria-expanded={isExpanded}
      >
        <StatusIcon className={`w-5 h-5 flex-shrink-0 ${iconColor}`} />

        <span className="flex-1 text-sm font-semibold text-gray-900">
          {formatStageName(stage.stage_name)}
        </span>

        {/* Change count badge */}
        {hasChanges && (
          <span className="flex-shrink-0 text-xs px-2 py-0.5 rounded-full bg-amber-100 text-amber-700 font-medium">
            {changeCount} change{changeCount !== 1 ? 's' : ''}
          </span>
        )}

        {/* Skipped badge */}
        {!stage.applied && (
          <span className="flex-shrink-0 text-xs px-2 py-0.5 rounded-full bg-gray-100 text-gray-500">
            skipped
          </span>
        )}

        {/* Expand/collapse chevron */}
        {isExpanded ? (
          <ChevronDown className="w-4 h-4 text-gray-400 flex-shrink-0" />
        ) : (
          <ChevronRight className="w-4 h-4 text-gray-400 flex-shrink-0" />
        )}
      </button>

      {/* Before → After SMILES (collapsed one-liner, shown when SMILES differ) */}
      {smilesDiffer && !isExpanded && (
        <div className="px-4 pb-3">
          <div className="flex items-start gap-2 text-xs font-mono bg-gray-50 rounded p-2 border border-gray-100">
            <span className="text-gray-500 flex-shrink-0 pt-0.5">before:</span>
            <span className="text-gray-700 break-all">{stage.input_smiles}</span>
          </div>
          <div className="flex items-start gap-2 text-xs font-mono bg-green-50 rounded p-2 border border-green-100 mt-1">
            <span className="text-green-600 flex-shrink-0 pt-0.5">after:</span>
            <span className="text-green-800 break-all">{stage.output_smiles}</span>
          </div>
        </div>
      )}

      {/* Expanded content */}
      <AnimatePresence initial={false}>
        {isExpanded && (
          <motion.div
            initial={{ opacity: 0, height: 0 }}
            animate={{ opacity: 1, height: 'auto' }}
            exit={{ opacity: 0, height: 0 }}
            transition={{ duration: 0.2 }}
            className="overflow-hidden"
          >
            <div className="px-4 pb-4 space-y-4 border-t border-gray-100">

              {/* Before / After SMILES */}
              <div className="mt-3 space-y-1">
                <div className="flex items-start gap-2 text-xs font-mono bg-gray-50 rounded p-2 border border-gray-100">
                  <span className="text-gray-500 flex-shrink-0 pt-0.5">before:</span>
                  <span className="text-gray-700 break-all">{stage.input_smiles || '—'}</span>
                </div>
                <div className="flex items-start gap-2 text-xs font-mono bg-green-50 rounded p-2 border border-green-100">
                  <span className="text-green-600 flex-shrink-0 pt-0.5">after:</span>
                  <span className="text-green-800 break-all">{stage.output_smiles || '—'}</span>
                </div>
              </div>

              {/* Charge changes */}
              {stage.charge_changes.length > 0 && (
                <div>
                  <h5 className="text-xs font-semibold text-gray-600 uppercase tracking-wide mb-2">
                    Charge Changes ({stage.charge_changes.length})
                  </h5>
                  <table className="w-full text-xs border-collapse">
                    <thead>
                      <tr className="bg-gray-50">
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Atom</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Element</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Before</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">After</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Rule</th>
                      </tr>
                    </thead>
                    <tbody>
                      {stage.charge_changes.map((c, i) => (
                        <tr key={i} className="border-b border-gray-100">
                          <td className="p-1.5">
                            <span className="inline-flex items-center justify-center w-6 h-6 rounded bg-blue-100 text-blue-700 font-mono font-medium text-[11px]">
                              {c.atom_idx}
                            </span>
                          </td>
                          <td className="p-1.5 font-medium text-gray-800">{c.element}</td>
                          <td className="p-1.5 font-mono text-gray-600">{c.before_charge > 0 ? `+${c.before_charge}` : c.before_charge}</td>
                          <td className="p-1.5 font-mono text-gray-600">{c.after_charge > 0 ? `+${c.after_charge}` : c.after_charge}</td>
                          <td className="p-1.5 text-gray-500">{c.rule_name}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}

              {/* Bond changes */}
              {stage.bond_changes.length > 0 && (
                <div>
                  <h5 className="text-xs font-semibold text-gray-600 uppercase tracking-wide mb-2">
                    Bond Changes ({stage.bond_changes.length})
                  </h5>
                  <table className="w-full text-xs border-collapse">
                    <thead>
                      <tr className="bg-gray-50">
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Bond</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Atoms</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Before</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">After</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Rule</th>
                      </tr>
                    </thead>
                    <tbody>
                      {stage.bond_changes.map((b, i) => (
                        <tr key={i} className="border-b border-gray-100">
                          <td className="p-1.5 font-mono text-gray-600">{b.bond_idx}</td>
                          <td className="p-1.5 text-gray-600">{b.atom1_idx}–{b.atom2_idx}</td>
                          <td className="p-1.5 font-mono text-gray-600">{b.before_type}</td>
                          <td className="p-1.5 font-mono text-gray-600">{b.after_type}</td>
                          <td className="p-1.5 text-gray-500">{b.rule_name}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}

              {/* Ring changes */}
              {stage.ring_changes.length > 0 && (
                <div>
                  <h5 className="text-xs font-semibold text-gray-600 uppercase tracking-wide mb-2">
                    Ring Aromaticity Changes ({stage.ring_changes.length})
                  </h5>
                  <table className="w-full text-xs border-collapse">
                    <thead>
                      <tr className="bg-gray-50">
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Ring Size</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Atoms</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Before Type</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">After Type</th>
                      </tr>
                    </thead>
                    <tbody>
                      {stage.ring_changes.map((r, i) => (
                        <tr key={i} className="border-b border-gray-100">
                          <td className="p-1.5 text-gray-800">{r.ring_size}-membered</td>
                          <td className="p-1.5 font-mono text-gray-600 text-[11px]">
                            [{r.ring_atoms.join(', ')}]
                          </td>
                          <td className="p-1.5">
                            <span className={`px-1.5 py-0.5 rounded text-[11px] font-medium ${
                              r.before_type === 'aromatic' ? 'bg-purple-100 text-purple-700' : 'bg-gray-100 text-gray-600'
                            }`}>
                              {r.before_type}
                            </span>
                          </td>
                          <td className="p-1.5">
                            <span className={`px-1.5 py-0.5 rounded text-[11px] font-medium ${
                              r.after_type === 'aromatic' ? 'bg-purple-100 text-purple-700' : 'bg-gray-100 text-gray-600'
                            }`}>
                              {r.after_type}
                            </span>
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}

              {/* Radical changes */}
              {stage.radical_changes.length > 0 && (
                <div>
                  <h5 className="text-xs font-semibold text-gray-600 uppercase tracking-wide mb-2">
                    Radical Changes ({stage.radical_changes.length})
                  </h5>
                  <table className="w-full text-xs border-collapse">
                    <thead>
                      <tr className="bg-gray-50">
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Atom</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Element</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Before</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">After</th>
                      </tr>
                    </thead>
                    <tbody>
                      {stage.radical_changes.map((r, i) => (
                        <tr key={i} className="border-b border-gray-100">
                          <td className="p-1.5 font-mono text-gray-600">{r.atom_idx}</td>
                          <td className="p-1.5 font-medium text-gray-800">{r.element}</td>
                          <td className="p-1.5 text-gray-600">{r.before_radicals} e⁻</td>
                          <td className="p-1.5 text-gray-600">{r.after_radicals} e⁻</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}

              {/* Fragment removals */}
              {stage.fragment_removals.length > 0 && (
                <div>
                  <h5 className="text-xs font-semibold text-gray-600 uppercase tracking-wide mb-2">
                    Removed Fragments ({stage.fragment_removals.length})
                  </h5>
                  <table className="w-full text-xs border-collapse">
                    <thead>
                      <tr className="bg-gray-50">
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Fragment</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Name</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">Role</th>
                        <th className="text-left p-1.5 text-gray-500 font-medium border-b border-gray-200">MW</th>
                      </tr>
                    </thead>
                    <tbody>
                      {stage.fragment_removals.map((f, i) => (
                        <tr key={i} className="border-b border-gray-100">
                          <td className="p-1.5 font-mono text-gray-700 text-[11px] break-all max-w-[120px]">
                            {f.smiles}
                          </td>
                          <td className="p-1.5 text-gray-600">{f.name ?? '—'}</td>
                          <td className="p-1.5">
                            <RoleBadge role={f.role} />
                          </td>
                          <td className="p-1.5 text-gray-600">{f.mw.toFixed(2)}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}

              {/* DVAL cross-references */}
              {stage.dval_cross_refs.length > 0 && (
                <div>
                  <h5 className="text-xs font-semibold text-gray-600 uppercase tracking-wide mb-2">
                    Deep Validation Cross-References
                  </h5>
                  <ul className="space-y-1">
                    {stage.dval_cross_refs.map((ref, i) => (
                      <li key={i} className="text-xs text-gray-600 bg-blue-50 rounded px-2 py-1 border border-blue-100">
                        {ref}
                      </li>
                    ))}
                  </ul>
                </div>
              )}

              {/* On-demand structure rendering */}
              {stage.output_smiles && (
                <div>
                  {!showStructure ? (
                    <button
                      onClick={() => setShowStructure(true)}
                      className="text-xs text-blue-600 hover:text-blue-800 underline underline-offset-2 transition-colors"
                    >
                      Show structure
                    </button>
                  ) : (
                    <div className="space-y-2">
                      <div className="flex items-center justify-between">
                        <span className="text-xs font-medium text-gray-600">Output structure</span>
                        <button
                          onClick={() => setShowStructure(false)}
                          className="text-xs text-gray-400 hover:text-gray-600 transition-colors"
                        >
                          Hide
                        </button>
                      </div>
                      <MoleculeViewer
                        smiles={stage.output_smiles}
                        highlightAtoms={highlightAtoms}
                        width={280}
                        height={180}
                        className="border border-gray-200 rounded-lg"
                      />
                      {highlightAtoms.length > 0 && (
                        <p className="text-[11px] text-gray-400">
                          Highlighted atoms: {highlightAtoms.join(', ')} (affected by changes)
                        </p>
                      )}
                    </div>
                  )}
                </div>
              )}
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}
