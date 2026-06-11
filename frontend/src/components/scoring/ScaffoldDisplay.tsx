import { useState } from 'react';
import type { ReactElement } from 'react';
import { MoleculeViewer } from '../molecules/MoleculeViewer';
import { CopyButton } from '../ui/CopyButton';
import { InfoTooltip, DoiLink } from '../ui/Tooltip';
import type { ScaffoldResult } from '../../types/scoring';

interface ScaffoldDisplayProps {
  result: ScaffoldResult;
}

type ScaffoldView = 'murcko' | 'generic';

/**
 * Displays Murcko scaffold with toggle between standard and generic forms.
 */
export function ScaffoldDisplay({ result }: ScaffoldDisplayProps): ReactElement {
  const [activeView, setActiveView] = useState<ScaffoldView>('murcko');

  const { scaffold_smiles, generic_scaffold_smiles, has_scaffold, message, details } = result;

  // Get the current scaffold to display
  const currentSmiles = activeView === 'murcko' ? scaffold_smiles : generic_scaffold_smiles;

  return (
    <div className="bg-[var(--color-surface-elevated)] rounded-lg border border-[var(--color-border)] p-6">
      <div className="flex items-center justify-between mb-4">
        <div className="flex items-center gap-2">
          <h3 className="text-lg font-semibold text-[var(--color-text-primary)]">Murcko Scaffold</h3>
          <InfoTooltip
            title="Murcko Scaffold Analysis"
            content={
              <div className="text-xs space-y-2">
                <p>The Murcko scaffold is the core ring structure of a molecule with linkers, removing all side chains.</p>
                <ul className="list-disc list-inside space-y-1 text-white/80">
                  <li><strong>Standard:</strong> Preserves atom types and bond orders</li>
                  <li><strong>Generic:</strong> Converts all atoms to carbon, all bonds to single (framework)</li>
                </ul>
                <p className="text-white/60">Useful for analyzing structural similarity and compound series.</p>
                <div className="mt-2 pt-2 border-t border-white/20 text-white/60 space-y-1">
                  <p>📖 Bemis & Murcko. J Med Chem (1996)</p>
                  <DoiLink doi="10.1021/jm9602928" />
                </div>
              </div>
            }
          />
        </div>
        {has_scaffold && (
          <span className="text-sm text-[var(--color-text-secondary)]">
            {Number(details.scaffold_rings) || 0} ring(s)
          </span>
        )}
      </div>

      {has_scaffold ? (
        <>
          {/* Toggle buttons */}
          <div className="flex gap-2 mb-4">
            <button
              onClick={() => setActiveView('murcko')}
              className={`px-3 py-1.5 text-sm rounded-md transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-[var(--color-primary)] focus-visible:ring-offset-2 ${activeView === 'murcko'
                  ? 'bg-chem-primary-600 text-white'
                  : 'bg-[var(--color-surface-sunken)] text-[var(--color-text-secondary)] hover:bg-[var(--color-border-strong)]'
                }`}
            >
              Standard
            </button>
            <button
              onClick={() => setActiveView('generic')}
              className={`px-3 py-1.5 text-sm rounded-md transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-[var(--color-primary)] focus-visible:ring-offset-2 ${activeView === 'generic'
                  ? 'bg-chem-primary-600 text-white'
                  : 'bg-[var(--color-surface-sunken)] text-[var(--color-text-secondary)] hover:bg-[var(--color-border-strong)]'
                }`}
            >
              Generic (Framework)
            </button>
          </div>

          {/* Scaffold viewer */}
          <div className="flex justify-center w-full min-h-[300px]">
            <MoleculeViewer
              smiles={currentSmiles}
              width={500}
              height={350}
            />
          </div>

          {/* SMILES display */}
          <div className="mt-4 p-2 bg-[var(--color-surface-sunken)] rounded text-xs font-mono text-[var(--color-text-secondary)] break-all flex items-center gap-2">
            <span className="flex-1">{currentSmiles}</span>
            <CopyButton text={currentSmiles || ''} size={14} />
          </div>

          {/* Message */}
          <p className="mt-3 text-sm text-[var(--color-text-secondary)]">{message}</p>
        </>
      ) : (
        /* No scaffold message */
        <div className="text-center py-8">
          <div className="text-4xl mb-3">&#128279;</div>
          <p className="text-[var(--color-text-secondary)]">{message}</p>
          <p className="text-sm text-[var(--color-text-muted)] mt-2">
            Acyclic molecules do not have ring scaffolds.
          </p>
        </div>
      )}
    </div>
  );
}
