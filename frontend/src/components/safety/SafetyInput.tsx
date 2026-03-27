import { useState, useEffect } from 'react';
import { useSearchParams } from 'react-router-dom';
import { motion } from 'framer-motion';
import { ClayButton } from '../ui/ClayButton';
import { integrationsApi } from '../../services/api';
import { cn } from '../../lib/utils';

interface SafetyInputProps {
  onScreen: (smiles: string) => void;
  isLoading: boolean;
}

/**
 * Detect whether the input looks like a SMILES string.
 * Returns true for strings containing typical SMILES characters.
 */
function looksLikeSMILES(input: string): boolean {
  const trimmed = input.trim();
  if (trimmed.startsWith('InChI=') || trimmed.startsWith('InChI ')) return false;
  const identifierPrefixes = [
    /^\d{2,7}-\d{2}-\d$/,         // CAS number (e.g. 50-78-2)
    /^CHEMBL\d+$/i,                // ChEMBL ID
    /^CID\s*\d+$/i,                // PubChem CID with prefix
    /^\d{4,}$/,                    // Numeric-only PubChem CID
    /^DB\d{5}$/i,                  // DrugBank ID
    /^Q\d+$/,                      // Wikidata QID
  ];
  if (identifierPrefixes.some(p => p.test(trimmed))) return false;
  const smilesChars = /[CNOSPFIBcnosp()[\]=@#%+\-1-9]/;
  return smilesChars.test(trimmed) && !trimmed.includes(' ');
}

/**
 * Safety page input.
 *
 * Accepts SMILES, InChI, CAS, ChEMBL ID, PubChem CID, or DrugBank ID.
 * Non-SMILES identifiers are resolved via /api/v1/resolve before screening.
 * Auto-populates and auto-submits from `?smiles=` query parameter (D-04).
 */
export function SafetyInput({ onScreen, isLoading }: SafetyInputProps) {
  const [searchParams] = useSearchParams();
  const [input, setInput] = useState(() => searchParams.get('smiles') || '');
  const [resolving, setResolving] = useState(false);
  const [resolveError, setResolveError] = useState<string | null>(null);

  // Auto-submit if query param present on mount (D-04)
  useEffect(() => {
    const smilesParam = searchParams.get('smiles');
    if (smilesParam) {
      onScreen(smilesParam);
    }
  }, []); // eslint-disable-line react-hooks/exhaustive-deps

  const isWorking = isLoading || resolving;

  async function handleSubmit() {
    const trimmed = input.trim();
    if (!trimmed) return;

    setResolveError(null);

    if (looksLikeSMILES(trimmed)) {
      onScreen(trimmed);
      return;
    }

    // Resolve identifier via API
    setResolving(true);
    try {
      const result = await integrationsApi.resolveIdentifier({ identifier: trimmed });
      if (result.resolved && result.canonical_smiles) {
        onScreen(result.canonical_smiles);
      } else {
        setResolveError(
          `Could not resolve "${trimmed}" to a structure. Check the identifier and try again.`
        );
      }
    } catch {
      setResolveError('Identifier resolution failed. Try entering the SMILES directly.');
    } finally {
      setResolving(false);
    }
  }

  function handleKeyDown(e: React.KeyboardEvent<HTMLInputElement>) {
    if (e.key === 'Enter' && !isWorking) {
      handleSubmit();
    }
  }

  return (
    <div className="w-full">
      <div className="flex gap-3 items-stretch">
        {/* Input field */}
        <div className="flex-1 relative">
          <input
            type="text"
            value={input}
            onChange={e => setInput(e.target.value)}
            onKeyDown={handleKeyDown}
            disabled={isWorking}
            placeholder="Enter SMILES, InChI, CAS number, ChEMBL ID, PubChem CID, or DrugBank ID..."
            className={cn(
              'w-full h-full px-4 py-3 rounded-2xl font-mono text-sm',
              'bg-[var(--color-surface-elevated)] dark:bg-[var(--color-surface-sunken)]',
              'border border-[var(--color-border-strong)]',
              'text-[var(--color-text-primary)]',
              'placeholder:text-[var(--color-text-muted)] placeholder:font-sans',
              'focus:outline-none focus:border-[var(--color-primary)] focus:ring-2 focus:ring-[var(--glow-primary)]',
              'disabled:opacity-50 disabled:cursor-not-allowed',
              'transition-all duration-200',
            )}
            spellCheck={false}
            autoComplete="off"
            autoCorrect="off"
          />
        </div>

        {/* Screen Molecule button */}
        <ClayButton
          variant="primary"
          size="lg"
          loading={isWorking}
          disabled={!input.trim() || isWorking}
          onClick={handleSubmit}
          className="shrink-0"
        >
          {resolving ? 'Resolving...' : 'Screen Molecule'}
        </ClayButton>
      </div>

      {/* Resolve error */}
      {resolveError && (
        <motion.p
          initial={{ opacity: 0, y: -4 }}
          animate={{ opacity: 1, y: 0 }}
          className="mt-2 text-sm text-status-error"
        >
          {resolveError}
        </motion.p>
      )}
    </div>
  );
}
