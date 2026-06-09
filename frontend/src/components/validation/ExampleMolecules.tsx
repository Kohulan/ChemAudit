import { type ComponentProps } from 'react';
import { motion } from 'framer-motion';

import { cn } from '../../lib/utils';
import { RecentMolecules } from '../molecules/RecentMolecules';

const EXAMPLE_MOLECULES = [
  { name: 'Aspirin', smiles: 'CC(=O)Oc1ccccc1C(=O)O' },
  { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' },
  { name: 'L-Alanine (chiral)', smiles: 'C[C@H](N)C(=O)O' },
  { name: 'E-Stilbene (E/Z)', smiles: 'C(/c1ccccc1)=C/c1ccccc1' },
  { name: 'Morphine', smiles: 'CN1CCC23C4=C5C=CC(O)=C4OC2C(O)C=CC3C1C5' },
  { name: 'Rhodanine (PAINS)', smiles: 'O=C1NC(=S)SC1' },
  { name: 'Amine HCl (salt)', smiles: 'CCN.Cl' },
];

const IDENTIFIER_EXAMPLES = [
  { name: 'ibuprofen', label: 'ibuprofen' },
  { name: 'CHEMBL25', label: 'CHEMBL25' },
  { name: '50-78-2', label: 'CAS 50-78-2' },
  { name: 'DB00945', label: 'DrugBank' },
];

type RecentProps = ComponentProps<typeof RecentMolecules>;

interface ExampleMoleculesProps {
  /** Called with a SMILES or identifier string when an example is clicked. */
  onExampleClick: (value: string) => void;
  recent: RecentProps['recent'];
  onSelectRecent: RecentProps['onSelect'];
  onRemoveRecent: RecentProps['onRemove'];
  onClearRecent: RecentProps['onClear'];
}

/**
 * Quick-start example molecules + identifier examples + the recent-molecules
 * dropdown. Extracted from SingleValidation to keep that page focused.
 */
export function ExampleMolecules({
  onExampleClick,
  recent,
  onSelectRecent,
  onRemoveRecent,
  onClearRecent,
}: ExampleMoleculesProps) {
  return (
    <motion.div
      className="flex flex-wrap gap-2 justify-center items-center"
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ delay: 0.2, duration: 0.5 }}
    >
      <span className="text-sm text-[var(--color-text-muted)] self-center mr-2">Try:</span>
      {EXAMPLE_MOLECULES.map((example, i) => (
        <motion.button
          key={example.name}
          onClick={() => onExampleClick(example.smiles)}
          className={cn(
            'px-3 py-1.5 text-sm rounded-full transition-all',
            'bg-[var(--color-surface-elevated)] border border-[var(--color-border)]',
            'text-[var(--color-text-secondary)]',
            'hover:border-[var(--color-primary)] hover:text-[var(--color-primary)]',
            'hover:shadow-[0_0_12px_var(--glow-primary)]'
          )}
          initial={{ opacity: 0, scale: 0.9 }}
          animate={{ opacity: 1, scale: 1 }}
          transition={{ delay: 0.3 + i * 0.05 }}
          whileHover={{ scale: 1.05 }}
          whileTap={{ scale: 0.95 }}
        >
          {example.name}
        </motion.button>
      ))}

      {/* Identifier examples separator */}
      <div className="hidden sm:block w-px h-6 bg-[var(--color-border)] mx-1" />
      {IDENTIFIER_EXAMPLES.map((example, i) => (
        <motion.button
          key={example.name}
          onClick={() => onExampleClick(example.name)}
          className={cn(
            'px-3 py-1.5 text-sm rounded-full transition-all',
            'bg-emerald-500/5 border border-emerald-500/20',
            'text-emerald-600 dark:text-emerald-400',
            'hover:border-emerald-500/40 hover:bg-emerald-500/10',
            'hover:shadow-[0_0_12px_rgba(16,185,129,0.15)]'
          )}
          initial={{ opacity: 0, scale: 0.9 }}
          animate={{ opacity: 1, scale: 1 }}
          transition={{ delay: 0.3 + (EXAMPLE_MOLECULES.length + i) * 0.05 }}
          whileHover={{ scale: 1.05 }}
          whileTap={{ scale: 0.95 }}
        >
          {example.label}
        </motion.button>
      ))}

      {/* Separator */}
      {recent.length > 0 && (
        <div className="hidden sm:block w-px h-6 bg-[var(--color-border)] mx-2" />
      )}

      {/* Recent molecules dropdown */}
      <RecentMolecules
        recent={recent}
        onSelect={onSelectRecent}
        onRemove={onRemoveRecent}
        onClear={onClearRecent}
      />
    </motion.div>
  );
}
