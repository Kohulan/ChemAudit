import { useState, useCallback, useRef } from 'react';
import { Upload } from 'lucide-react';
import { ClayCard } from '../ui/ClayCard';
import { ClayButton } from '../ui/ClayButton';

// =============================================================================
// Props
// =============================================================================

interface GenChemInputProps {
  /** Current textarea SMILES value. */
  smilesInput: string;
  /** Called when the textarea value changes. */
  onSmilesChange: (value: string) => void;
  /** Currently uploaded file (or null). */
  uploadedFile: File | null;
  /** Called when the file selection changes (null = cleared). */
  onFileChange: (file: File | null) => void;
  /** Called when the Run Filter button is clicked. */
  onRunFilter: () => void;
  /** Disables the Run Filter button and file drop zone while filtering. */
  isLoading: boolean;
}

// =============================================================================
// Tooltip descriptions for each funnel stage (used in FunnelChart — exported)
// =============================================================================

export const STAGE_CRITERIA: Record<string, string> = {
  parse: 'rejects molecules that cannot be parsed by RDKit',
  valence: 'rejects molecules with valence errors',
  alerts: 'rejects molecules matching PAINS, Brenk, Kazius, or NIBR alert patterns',
  property: 'rejects molecules outside MW, LogP, TPSA, rotatable bonds, or ring count thresholds',
  sa_score: 'rejects molecules with SA score above the configured maximum',
  dedup: 'rejects duplicate molecules by InChIKey',
  novelty: 'rejects molecules too similar to known drugs (Tanimoto similarity above threshold)',
};

// =============================================================================
// Component
// =============================================================================

/**
 * GenChem bulk SMILES input area.
 *
 * Provides two input modes:
 * 1. Textarea for pasting SMILES (one per line, monospace font per D-03)
 * 2. Drag-and-drop / click-to-browse file upload zone (.txt and .csv only)
 *
 * The two modes are separated by an OR divider. Only one is used per run:
 * if a file is selected, the file takes precedence over textarea content.
 *
 * Accessibility: drop zone has role="button", aria-label, tabIndex=0.
 * Per UI-SPEC D-02, D-03 copywriting and layout contracts.
 */
export function GenChemInput({
  smilesInput,
  onSmilesChange,
  uploadedFile,
  onFileChange,
  onRunFilter,
  isLoading,
}: GenChemInputProps) {
  const [isDragging, setIsDragging] = useState(false);
  const [fileError, setFileError] = useState<string | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  // ── Drag handlers ──

  const handleDragEnter = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(true);
  }, []);

  const handleDragLeave = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(false);
  }, []);

  const handleDragOver = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
  }, []);

  // ── File validation helper ──

  const processFile = useCallback(
    (selectedFile: File) => {
      const name = selectedFile.name.toLowerCase();
      if (!name.endsWith('.txt') && !name.endsWith('.csv')) {
        setFileError('Invalid file type. Please upload a .txt or .csv file.');
        return;
      }
      setFileError(null);
      onFileChange(selectedFile);
    },
    [onFileChange],
  );

  const handleDrop = useCallback(
    (e: React.DragEvent) => {
      e.preventDefault();
      e.stopPropagation();
      setIsDragging(false);

      if (isLoading) return;

      const files = e.dataTransfer.files;
      if (files.length > 0) {
        processFile(files[0]);
      }
    },
    [isLoading, processFile],
  );

  const handleFileSelect = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const files = e.target.files;
      if (files && files.length > 0) {
        processFile(files[0]);
      }
      // Reset input so re-selecting same file triggers onChange
      if (fileInputRef.current) fileInputRef.current.value = '';
    },
    [processFile],
  );

  const handleDropZoneClick = useCallback(() => {
    if (!isLoading) {
      fileInputRef.current?.click();
    }
  }, [isLoading]);

  const handleDropZoneKeyDown = useCallback(
    (e: React.KeyboardEvent) => {
      if (e.key === 'Enter' || e.key === ' ') {
        e.preventDefault();
        handleDropZoneClick();
      }
    },
    [handleDropZoneClick],
  );

  const canRun = !isLoading && (smilesInput.trim().length > 0 || uploadedFile !== null);

  return (
    <ClayCard variant="default" size="md" className="p-4">
      {/* ── SMILES textarea ── */}
      <div>
        <label
          htmlFor="genchem-smiles-input"
          className="block text-sm font-medium text-[var(--color-text-primary)] mb-2"
        >
          Paste SMILES
        </label>
        <textarea
          id="genchem-smiles-input"
          rows={6}
          value={smilesInput}
          onChange={(e) => onSmilesChange(e.target.value)}
          disabled={isLoading}
          placeholder="Paste SMILES here, one per line, or drop a .txt or .csv file to filter generative model output."
          className="w-full border border-[var(--color-border)] bg-[var(--color-surface-elevated)] text-[var(--color-text-primary)] rounded-lg px-4 py-3 text-sm font-mono resize-none focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)] focus:border-transparent placeholder:text-[var(--color-text-muted)] placeholder:font-sans disabled:opacity-50"
        />
        <p className="text-xs text-[var(--color-text-muted)] mt-1">
          One SMILES per line
        </p>
      </div>

      {/* ── OR divider ── */}
      <div className="flex items-center gap-2 my-4">
        <div className="flex-1 h-px bg-[var(--color-border)]" />
        <span className="text-xs text-[var(--color-text-muted)] font-medium uppercase tracking-wide">
          OR
        </span>
        <div className="flex-1 h-px bg-[var(--color-border)]" />
      </div>

      {/* ── File drop zone ── */}
      <div
        role="button"
        aria-label="Upload SMILES file"
        tabIndex={0}
        onDragEnter={handleDragEnter}
        onDragLeave={handleDragLeave}
        onDragOver={handleDragOver}
        onDrop={handleDrop}
        onClick={handleDropZoneClick}
        onKeyDown={handleDropZoneKeyDown}
        className={[
          'border-2 border-dashed rounded-lg p-8 text-center cursor-pointer transition-all duration-200',
          isDragging
            ? 'border-[var(--color-primary)] bg-[var(--color-primary)]/5'
            : 'border-[var(--color-border-strong)] hover:border-[var(--color-text-muted)] hover:bg-[var(--color-surface-sunken)]/50',
          isLoading ? 'opacity-50 cursor-not-allowed' : '',
        ]
          .filter(Boolean)
          .join(' ')}
      >
        <input
          ref={fileInputRef}
          type="file"
          accept=".txt,.csv"
          onChange={handleFileSelect}
          className="hidden"
          disabled={isLoading}
          aria-hidden="true"
        />

        <div className="flex flex-col items-center gap-3 text-[var(--color-text-secondary)]">
          <div className="w-12 h-12 rounded-xl bg-[var(--color-surface-sunken)] flex items-center justify-center">
            <Upload className="w-6 h-6 text-[var(--color-text-muted)]" />
          </div>
          {uploadedFile ? (
            <>
              <p className="text-sm font-medium text-[var(--color-text-primary)]">
                {uploadedFile.name}
              </p>
              <p className="text-xs text-[var(--color-text-muted)]">
                Click to change file
              </p>
            </>
          ) : (
            <>
              <p className="text-sm font-medium text-[var(--color-text-primary)]">
                {isDragging
                  ? 'Drop file here'
                  : 'Drop a .txt or .csv file here, or click to browse'}
              </p>
              <p className="text-xs text-[var(--color-text-muted)]">
                Accepts <span className="font-medium">.txt</span> and{' '}
                <span className="font-medium">.csv</span> files
              </p>
            </>
          )}
        </div>
      </div>

      {/* ── File error message ── */}
      {fileError && (
        <div className="mt-3 bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-800 rounded-lg px-4 py-3 text-sm text-red-700 dark:text-red-400">
          {fileError}
        </div>
      )}

      {/* ── Run Filter CTA — right-aligned ── */}
      <div className="flex justify-end mt-4">
        <ClayButton
          variant="primary"
          size="md"
          onClick={onRunFilter}
          disabled={!canRun}
          loading={isLoading}
        >
          Run Filter
        </ClayButton>
      </div>
    </ClayCard>
  );
}
