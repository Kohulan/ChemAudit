import { useState, useRef, useCallback, useEffect } from 'react';
import { useSearchParams } from 'react-router-dom';
import { motion } from 'framer-motion';
import { Upload } from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';
import { MoleculeLoader } from '../ui/MoleculeLoader';
import { integrationsApi } from '../../services/api';
import { cn } from '../../lib/utils';

interface DiagnosticsInputProps {
  onAnalyze: (smiles: string) => void;
  onFileUpload: (file: File) => void;
  isLoading: boolean;
  fileLoading: boolean;
  initialSmiles?: string;
}

const ACCEPTED_EXTENSIONS = ['.sdf', '.sd', '.csv'];

/**
 * Detect whether the input looks like a SMILES string (not an identifier).
 */
function looksLikeSMILES(input: string): boolean {
  const trimmed = input.trim();
  if (trimmed.startsWith('InChI=') || trimmed.startsWith('InChI ')) return false;
  const identifierPrefixes = [
    /^\d{2,7}-\d{2}-\d$/,   // CAS number
    /^CHEMBL\d+$/i,           // ChEMBL ID
    /^CID\s*\d+$/i,           // PubChem CID with prefix
    /^\d{4,}$/,               // Numeric-only PubChem CID
    /^DB\d{5}$/i,             // DrugBank ID
    /^Q\d+$/,                 // Wikidata QID
  ];
  if (identifierPrefixes.some(p => p.test(trimmed))) return false;
  const smilesChars = /[CNOSPFIBcnosp()[\]=@#%+\-1-9]/;
  return smilesChars.test(trimmed) && !trimmed.includes(' ');
}

/**
 * DiagnosticsInput — unified input area for the Diagnostics page.
 *
 * Contains two sections (always visible, not tabbed):
 * 1. Molecule identifier input with identifier resolution and "Analyze" button.
 * 2. SDF / CSV file drag-and-drop upload zone for file pre-validation.
 *
 * Supports ?smiles= URL query param for cross-link from SingleValidation (D-02).
 */
export function DiagnosticsInput({
  onAnalyze,
  onFileUpload,
  isLoading,
  fileLoading,
  initialSmiles,
}: DiagnosticsInputProps) {
  const [searchParams] = useSearchParams();
  const [input, setInput] = useState(() => initialSmiles ?? searchParams.get('smiles') ?? '');
  const [resolving, setResolving] = useState(false);
  const [resolveError, setResolveError] = useState<string | null>(null);

  // File drop zone state
  const [isDraggingOver, setIsDraggingOver] = useState(false);
  const [dragLabel, setDragLabel] = useState<string | null>(null);
  const [fileError, setFileError] = useState<string | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  const isWorking = isLoading || resolving;

  // Auto-submit if ?smiles= query param present on mount (D-02 cross-link support)
  useEffect(() => {
    const smilesParam = searchParams.get('smiles') ?? initialSmiles;
    if (smilesParam) {
      setInput(smilesParam);
      onAnalyze(smilesParam);
    }
  }, []); // eslint-disable-line react-hooks/exhaustive-deps

  async function handleSubmit() {
    const trimmed = input.trim();
    if (!trimmed) return;
    setResolveError(null);

    if (looksLikeSMILES(trimmed)) {
      onAnalyze(trimmed);
      return;
    }

    // Resolve identifier via the existing integrations API
    setResolving(true);
    try {
      const result = await integrationsApi.resolveIdentifier({ identifier: trimmed });
      if (result.resolved && result.canonical_smiles) {
        onAnalyze(result.canonical_smiles);
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

  // ── File drop zone handlers ──

  function validateExtension(filename: string): boolean {
    const lower = filename.toLowerCase();
    return ACCEPTED_EXTENSIONS.some(ext => lower.endsWith(ext));
  }

  const handleFile = useCallback(
    (file: File) => {
      setFileError(null);
      if (!validateExtension(file.name)) {
        setFileError(`Unsupported file type "${file.name}". Please upload an SDF or CSV file.`);
        return;
      }
      onFileUpload(file);
    },
    [onFileUpload],
  );

  function handleDragOver(e: React.DragEvent<HTMLDivElement>) {
    e.preventDefault();
    setIsDraggingOver(true);
  }

  function handleDragEnter(e: React.DragEvent<HTMLDivElement>) {
    e.preventDefault();
    setIsDraggingOver(true);
    setDragLabel('Release to upload');
  }

  function handleDragLeave(e: React.DragEvent<HTMLDivElement>) {
    e.preventDefault();
    setIsDraggingOver(false);
    setDragLabel(null);
  }

  function handleDrop(e: React.DragEvent<HTMLDivElement>) {
    e.preventDefault();
    setIsDraggingOver(false);
    setDragLabel(null);
    const file = e.dataTransfer.files[0];
    if (file) handleFile(file);
  }

  function handleFileInputChange(e: React.ChangeEvent<HTMLInputElement>) {
    const file = e.target.files?.[0];
    if (file) handleFile(file);
    // Reset input so the same file can be re-uploaded
    if (fileInputRef.current) fileInputRef.current.value = '';
  }

  function handleZoneKeyDown(e: React.KeyboardEvent<HTMLDivElement>) {
    if (e.key === 'Enter' || e.key === ' ') {
      e.preventDefault();
      fileInputRef.current?.click();
    }
  }

  function handleZoneClick() {
    fileInputRef.current?.click();
  }

  return (
    <div className="w-full space-y-0">
      {/* ── Section 1: Molecule identifier input ── */}
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

          {/* Analyze button */}
          <ClayButton
            variant="primary"
            size="lg"
            loading={isWorking}
            disabled={!input.trim() || isWorking}
            onClick={handleSubmit}
            className="shrink-0"
          >
            {resolving ? 'Resolving...' : 'Analyze'}
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

      {/* ── Divider ── */}
      <div className="border-b border-[var(--color-border)] pt-6 mb-6" />

      {/* ── Section 2: File upload zone ── */}
      <div>
        <p className="text-xs font-semibold text-[var(--color-text-secondary)] uppercase tracking-wide mb-3">
          File Pre-Validation
        </p>

        {/* Hidden file input */}
        <input
          ref={fileInputRef}
          type="file"
          accept=".sdf,.sd,.csv"
          className="sr-only"
          onChange={handleFileInputChange}
          tabIndex={-1}
          aria-hidden="true"
        />

        {/* Drag-and-drop zone */}
        <div
          role="button"
          tabIndex={0}
          aria-label="Drop SDF or CSV file here to upload, or press Enter to open file browser"
          onClick={handleZoneClick}
          onKeyDown={handleZoneKeyDown}
          onDragOver={handleDragOver}
          onDragEnter={handleDragEnter}
          onDragLeave={handleDragLeave}
          onDrop={handleDrop}
          className={cn(
            'relative flex flex-col items-center justify-center',
            'min-h-[120px] rounded-2xl cursor-pointer',
            'border-2 border-dashed transition-all duration-200',
            'focus:outline-none focus:ring-2 focus:ring-[var(--glow-primary)]',
            isDraggingOver
              ? 'border-chem-primary-600/40 bg-[rgba(var(--color-primary-rgb),0.04)]'
              : 'border-[var(--color-border)] hover:border-[var(--color-border-strong)]',
            fileLoading ? 'pointer-events-none' : '',
          )}
        >
          {fileLoading ? (
            <MoleculeLoader size="sm" text="Validating file..." />
          ) : (
            <div className="flex flex-col items-center gap-2 text-center px-6">
              <Upload
                className={cn(
                  'w-6 h-6 transition-colors duration-200',
                  isDraggingOver
                    ? 'text-[var(--color-primary)]'
                    : 'text-[var(--color-text-muted)]',
                )}
              />
              <p className="text-sm text-[var(--color-text-secondary)]">
                {dragLabel ?? 'Drop an SDF or CSV file here, or click to browse'}
              </p>
              <p className="text-xs text-[var(--color-text-muted)]">
                Supported formats: .sdf, .sd, .csv
              </p>
            </div>
          )}
        </div>

        {/* File extension error */}
        {fileError && (
          <motion.p
            initial={{ opacity: 0, y: -4 }}
            animate={{ opacity: 1, y: 0 }}
            className="mt-2 text-sm text-status-error"
          >
            {fileError}
          </motion.p>
        )}
      </div>
    </div>
  );
}
