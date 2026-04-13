import { useState, useCallback, useRef } from 'react';
import { Upload, ChevronDown, CheckCircle2, AlertCircle } from 'lucide-react';
import { ClayCard } from '../ui/ClayCard';
import { ClayButton } from '../ui/ClayButton';
import { batchApi } from '../../services/api';
import type { CSVColumnsResponse } from '../../types/batch';

// =============================================================================
// Props
// =============================================================================

interface StructureFilterInputProps {
  /** Current textarea SMILES value. */
  smilesInput: string;
  /** Called when the textarea value changes. */
  onSmilesChange: (value: string) => void;
  /** Called when the Run Filter button is clicked with parsed SMILES list. */
  onRunFilter: (smilesList: string[]) => void;
  /** Called when an SDF file needs to be uploaded directly to the backend. */
  onRunFilterFile: (file: File) => void;
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
  sa: 'rejects molecules with SA score above the configured maximum',
  sa_score: 'rejects molecules with SA score above the configured maximum',
  dedup: 'rejects duplicate molecules by InChIKey',
  novelty: 'rejects molecules too similar to known drugs (Tanimoto similarity above threshold)',
};

// Supported file extensions
const DELIMITED_FORMATS = ['.csv', '.tsv', '.txt'];
const ALL_FORMATS = [...DELIMITED_FORMATS, '.sdf'];

function isDelimitedFormat(filename: string): boolean {
  const lower = filename.toLowerCase();
  return DELIMITED_FORMATS.some(ext => lower.endsWith(ext));
}

function isSdfFormat(filename: string): boolean {
  return filename.toLowerCase().endsWith('.sdf');
}

// =============================================================================
// Component
// =============================================================================

export function StructureFilterInput({
  smilesInput,
  onSmilesChange,
  onRunFilter,
  onRunFilterFile,
  isLoading,
}: StructureFilterInputProps) {
  const [isDragging, setIsDragging] = useState(false);
  const [fileError, setFileError] = useState<string | null>(null);
  const [selectedFile, setSelectedFile] = useState<File | null>(null);
  const [isAnalyzing, setIsAnalyzing] = useState(false);

  // CSV column detection state
  const [csvColumns, setCsvColumns] = useState<CSVColumnsResponse | null>(null);
  const [selectedSmilesColumn, setSelectedSmilesColumn] = useState('');

  // SDF state
  const [sdfMoleculeCount, setSdfMoleculeCount] = useState<number | null>(null);

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

  // ── File processing ──

  const processFile = useCallback(
    async (file: File) => {
      const name = file.name.toLowerCase();
      if (!ALL_FORMATS.some(ext => name.endsWith(ext))) {
        setFileError('Invalid file type. Please upload a CSV, TSV, TXT, or SDF file.');
        return;
      }

      setFileError(null);
      setSelectedFile(file);
      setCsvColumns(null);
      setSelectedSmilesColumn('');
      setSdfMoleculeCount(null);

      if (isDelimitedFormat(file.name)) {
        // For CSV/TSV: detect columns and suggest SMILES column
        setIsAnalyzing(true);
        try {
          const columns = await batchApi.detectColumnsLocal(file);
          setCsvColumns(columns);
          setSelectedSmilesColumn(columns.suggested_smiles || columns.columns[0] || '');
        } catch (e: any) {
          // If column detection fails, treat as plain SMILES list (one per line)
          setCsvColumns(null);
        } finally {
          setIsAnalyzing(false);
        }
      } else if (isSdfFormat(file.name)) {
        // For SDF: estimate molecule count from $$$$
        setIsAnalyzing(true);
        try {
          const text = await file.text();
          const count = (text.match(/\$\$\$\$/g) || []).length;
          setSdfMoleculeCount(count);
        } catch {
          setSdfMoleculeCount(null);
        } finally {
          setIsAnalyzing(false);
        }
      }
    },
    [],
  );

  const handleDrop = useCallback(
    async (e: React.DragEvent) => {
      e.preventDefault();
      e.stopPropagation();
      setIsDragging(false);
      if (isLoading || isAnalyzing) return;
      const files = e.dataTransfer.files;
      if (files.length > 0) await processFile(files[0]);
    },
    [isLoading, isAnalyzing, processFile],
  );

  const handleFileSelect = useCallback(
    async (e: React.ChangeEvent<HTMLInputElement>) => {
      const files = e.target.files;
      if (files && files.length > 0) await processFile(files[0]);
      if (fileInputRef.current) fileInputRef.current.value = '';
    },
    [processFile],
  );

  const handleDropZoneClick = useCallback(() => {
    if (!isLoading && !isAnalyzing) fileInputRef.current?.click();
  }, [isLoading, isAnalyzing]);

  const handleDropZoneKeyDown = useCallback(
    (e: React.KeyboardEvent) => {
      if (e.key === 'Enter' || e.key === ' ') {
        e.preventDefault();
        handleDropZoneClick();
      }
    },
    [handleDropZoneClick],
  );

  // ── Run Filter: extract SMILES and call parent ──

  const handleRunFilter = useCallback(() => {
    if (selectedFile) {
      const reader = new FileReader();
      reader.onload = (e) => {
        const text = e.target?.result as string;

        if (isSdfFormat(selectedFile.name)) {
          // SDF requires RDKit parsing on the backend — send the raw file
          onRunFilterFile(selectedFile);
          return;
        }

        if (csvColumns && selectedSmilesColumn) {
          // CSV with column detection: extract SMILES from the selected column
          const lines = text.split('\n').filter(line => line.trim());
          if (lines.length < 2) return;

          // Detect delimiter
          const firstLine = lines[0];
          const delimiter = firstLine.includes('\t') ? '\t' : ',';
          const header = parseDelimitedLine(firstLine, delimiter);
          const colIndex = header.findIndex(
            col => col.toLowerCase() === selectedSmilesColumn.toLowerCase()
          );

          if (colIndex < 0) return;

          const smilesList: string[] = [];
          for (let i = 1; i < lines.length; i++) {
            const values = parseDelimitedLine(lines[i], delimiter);
            const smi = values[colIndex]?.trim();
            if (smi) smilesList.push(smi);
          }
          onRunFilter(smilesList);
        } else {
          // Plain text / .txt: one SMILES per line
          const smilesList = text.split('\n').map(s => s.trim()).filter(Boolean);
          onRunFilter(smilesList);
        }
      };
      reader.readAsText(selectedFile);
    } else {
      // Textarea mode
      const smilesList = smilesInput.split('\n').map(s => s.trim()).filter(Boolean);
      if (smilesList.length > 0) onRunFilter(smilesList);
    }
  }, [selectedFile, csvColumns, selectedSmilesColumn, smilesInput, onRunFilter]);

  const canRun =
    !isLoading &&
    !isAnalyzing &&
    (smilesInput.trim().length > 0 ||
      (selectedFile !== null && (!csvColumns || selectedSmilesColumn)));

  // Whether to show the column selector (CSV with detected columns)
  const showColumnSelector = selectedFile && csvColumns && !isAnalyzing;
  const isSdf = selectedFile && isSdfFormat(selectedFile.name);

  return (
    <ClayCard variant="default" size="md" className="p-4">
      {/* ── SMILES textarea ── */}
      <div>
        <label
          htmlFor="structure-filter-smiles-input"
          className="block text-sm font-medium text-[var(--color-text-primary)] mb-2"
        >
          Paste SMILES
        </label>
        <textarea
          id="structure-filter-smiles-input"
          rows={6}
          value={smilesInput}
          onChange={(e) => onSmilesChange(e.target.value)}
          disabled={isLoading}
          placeholder="Paste SMILES here, one per line, or upload a CSV/SDF file below."
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
        aria-label="Upload file"
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
          isLoading || isAnalyzing ? 'opacity-50 cursor-not-allowed' : '',
        ]
          .filter(Boolean)
          .join(' ')}
      >
        <input
          ref={fileInputRef}
          type="file"
          accept=".txt,.csv,.tsv,.sdf"
          onChange={handleFileSelect}
          className="hidden"
          disabled={isLoading || isAnalyzing}
          aria-hidden="true"
        />

        <div className="flex flex-col items-center gap-3 text-[var(--color-text-secondary)]">
          <div className="w-12 h-12 rounded-xl bg-[var(--color-surface-sunken)] flex items-center justify-center">
            <Upload className="w-6 h-6 text-[var(--color-text-muted)]" />
          </div>
          {isAnalyzing ? (
            <>
              <p className="text-sm font-medium text-[var(--color-text-primary)]">
                Analyzing file...
              </p>
              <div className="w-5 h-5 border-2 border-[var(--color-primary)] border-t-transparent rounded-full animate-spin" />
            </>
          ) : selectedFile ? (
            <>
              <p className="text-sm font-medium text-[var(--color-text-primary)]">
                {selectedFile.name}
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
                  : 'Drop a CSV, TXT, or SDF file here, or click to browse'}
              </p>
              <p className="text-xs text-[var(--color-text-muted)]">
                Accepts <span className="font-medium">.csv</span>,{' '}
                <span className="font-medium">.tsv</span>,{' '}
                <span className="font-medium">.txt</span>, and{' '}
                <span className="font-medium">.sdf</span> files
              </p>
            </>
          )}
        </div>
      </div>

      {/* ── CSV column selector ── */}
      {showColumnSelector && (
        <div className="mt-4 space-y-3">
          <div className="flex items-start gap-2 p-3 rounded-lg bg-[var(--color-surface-sunken)]">
            <AlertCircle className="w-4 h-4 text-[var(--color-primary)] mt-0.5 flex-shrink-0" />
            <div className="text-sm text-[var(--color-text-secondary)]">
              <p>
                Detected <strong>{csvColumns.columns.length}</strong> columns and ~<strong>{csvColumns.row_count_estimate.toLocaleString()}</strong> rows.
                Select the column containing SMILES:
              </p>
            </div>
          </div>

          <div>
            <label className="block text-sm font-medium text-[var(--color-text-primary)] mb-2">
              SMILES Column <span className="text-red-500">*</span>
            </label>
            <div className="relative">
              <select
                value={selectedSmilesColumn}
                onChange={(e) => setSelectedSmilesColumn(e.target.value)}
                disabled={isLoading}
                className="w-full appearance-none border border-[var(--color-border)] bg-[var(--color-surface-elevated)] text-[var(--color-text-primary)] rounded-lg px-4 py-2.5 pr-10 text-sm focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)] focus:border-transparent"
              >
                <option value="">Select SMILES column...</option>
                {csvColumns.columns.map((col) => (
                  <option key={col} value={col}>
                    {col}
                    {col === csvColumns.suggested_smiles ? ' (suggested)' : ''}
                    {csvColumns.column_samples[col] ? ` — e.g., "${csvColumns.column_samples[col]}"` : ''}
                  </option>
                ))}
              </select>
              <ChevronDown className="absolute right-3 top-1/2 -translate-y-1/2 w-4 h-4 text-[var(--color-text-muted)] pointer-events-none" />
            </div>
            {selectedSmilesColumn && csvColumns.column_samples[selectedSmilesColumn] && (
              <p className="mt-1 text-xs text-[var(--color-text-muted)]">
                Sample: <code className="bg-[var(--color-surface-sunken)] px-1 py-0.5 rounded">{csvColumns.column_samples[selectedSmilesColumn]}</code>
              </p>
            )}
          </div>

          {selectedSmilesColumn && (
            <div className="flex items-start gap-2 p-3 rounded-lg bg-green-500/10 border border-green-500/20">
              <CheckCircle2 className="w-4 h-4 text-green-500 mt-0.5 flex-shrink-0" />
              <p className="text-sm text-green-600 dark:text-green-400">
                Ready to filter ~{csvColumns.row_count_estimate.toLocaleString()} molecules from column &ldquo;{selectedSmilesColumn}&rdquo;
              </p>
            </div>
          )}
        </div>
      )}

      {/* ── SDF info ── */}
      {isSdf && !isAnalyzing && (
        <div className="mt-4">
          <div className="flex items-start gap-2 p-3 rounded-lg bg-green-500/10 border border-green-500/20">
            <CheckCircle2 className="w-4 h-4 text-green-500 mt-0.5 flex-shrink-0" />
            <p className="text-sm text-green-600 dark:text-green-400">
              SDF file ready to process
              {sdfMoleculeCount !== null && ` (~${sdfMoleculeCount.toLocaleString()} molecules)`}
            </p>
          </div>
        </div>
      )}

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
          onClick={handleRunFilter}
          disabled={!canRun}
          loading={isLoading}
        >
          Run Filter
        </ClayButton>
      </div>
    </ClayCard>
  );
}

// =============================================================================
// Helpers (CSV parsing — matches the logic in api.ts detectColumnsLocal)
// =============================================================================

function parseDelimitedLine(line: string, delimiter: string): string[] {
  const result: string[] = [];
  let current = '';
  let inQuotes = false;

  for (let i = 0; i < line.length; i++) {
    const char = line[i];
    if (char === '"') {
      if (inQuotes && line[i + 1] === '"') {
        current += '"';
        i++;
      } else {
        inQuotes = !inQuotes;
      }
    } else if (char === delimiter && !inQuotes) {
      result.push(current.trim());
      current = '';
    } else {
      current += char;
    }
  }
  result.push(current.trim());
  return result;
}
