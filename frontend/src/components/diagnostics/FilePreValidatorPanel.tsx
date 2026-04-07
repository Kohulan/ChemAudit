import { useRef, useCallback, useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { Upload, AlertTriangle } from 'lucide-react';
import { ClayCard } from '../ui/ClayCard';
import { ClayButton } from '../ui/ClayButton';
import { Badge } from '../ui/Badge';
import { MoleculeLoader } from '../ui/MoleculeLoader';
import { FileIssueList } from './FileIssueList';
import type { FilePreValidationResponse } from '../../types/diagnostics';
import { cn } from '../../lib/utils';

interface FilePreValidatorPanelProps {
  result: FilePreValidationResponse | null;
  isLoading: boolean;
  error: string | null;
  onFileUpload: (file: File) => void;
  onClear: () => void;
}

const ACCEPTED_EXTENSIONS = ['.sdf', '.sd', '.csv'];
// Max file size in MB — aligns with existing env config (default 50MB)
const MAX_FILE_SIZE_MB = 50;

function validateExtension(filename: string): boolean {
  const lower = filename.toLowerCase();
  return ACCEPTED_EXTENSIONS.some((ext) => lower.endsWith(ext));
}

/**
 * Determine verdict badge variant and text based on file validation result.
 *
 * Per UI-SPEC copywriting (D-10):
 * - No issues: "Ready to parse" (success)
 * - Issues without error severity: "{M} issue(s) require attention" (warning)
 * - Has error-severity issues: "File has critical issues — do not parse" (error)
 */
function getVerdictInfo(result: FilePreValidationResponse): {
  variant: 'success' | 'warning' | 'error';
  text: string;
} {
  const hasErrorSeverity = result.issues.some(
    (i) => i.severity.toLowerCase() === 'error'
  );

  if (result.valid && result.issue_count === 0) {
    return { variant: 'success', text: 'Ready to parse' };
  }
  if (hasErrorSeverity) {
    return { variant: 'error', text: 'File has critical issues — do not parse' };
  }
  return {
    variant: 'warning',
    text: `${result.issue_count} issue${result.issue_count !== 1 ? 's' : ''} require attention`,
  };
}

/**
 * Get the total record count from a file validation result.
 * SDF uses total_blocks, CSV uses total_rows.
 */
function getTotalRecords(result: FilePreValidationResponse): number | null {
  return result.file_type === 'CSV' ? result.total_rows : result.total_blocks;
}

/**
 * File Pre-Validator panel (DIAG-05, D-10, D-11, D-12).
 *
 * States:
 * 1. No file: drag-and-drop zone
 * 2. Loading: MoleculeLoader spinner centered in zone
 * 3. Error: error card with Retry button
 * 4. Result: summary card with verdict badge + FileIssueList
 *
 * This panel manages its own drag-and-drop zone (separate from DiagnosticsInput).
 * Per D-11: standalone — no integration with Batch upload page.
 * Per D-12: file size limit inherits from existing env config.
 */
export function FilePreValidatorPanel({
  result,
  isLoading,
  error,
  onFileUpload,
  onClear,
}: FilePreValidatorPanelProps) {
  const fileInputRef = useRef<HTMLInputElement>(null);
  const isDraggingRef = useRef(false);
  const [isDraggingOver, setIsDraggingOver] = useState(false);
  const [localError, setLocalError] = useState<string | null>(null);

  const handleFile = useCallback(
    (file: File) => {
      setLocalError(null);

      if (!validateExtension(file.name)) {
        setLocalError(
          `Unsupported file type "${file.name}". Please upload an SDF or CSV file.`
        );
        return;
      }

      const sizeMB = file.size / (1024 * 1024);
      if (sizeMB > MAX_FILE_SIZE_MB) {
        setLocalError(
          `File exceeds ${MAX_FILE_SIZE_MB}MB limit. Use the batch API for large files.`
        );
        return;
      }

      onFileUpload(file);
    },
    [onFileUpload],
  );

  function handleDragOver(e: React.DragEvent<HTMLDivElement>) {
    e.preventDefault();
    if (!isDraggingRef.current) {
      isDraggingRef.current = true;
      setIsDraggingOver(true);
    }
  }

  function handleDragEnter(e: React.DragEvent<HTMLDivElement>) {
    e.preventDefault();
    setIsDraggingOver(true);
  }

  function handleDragLeave(e: React.DragEvent<HTMLDivElement>) {
    e.preventDefault();
    isDraggingRef.current = false;
    setIsDraggingOver(false);
  }

  function handleDrop(e: React.DragEvent<HTMLDivElement>) {
    e.preventDefault();
    isDraggingRef.current = false;
    setIsDraggingOver(false);
    const file = e.dataTransfer.files[0];
    if (file) handleFile(file);
  }

  function handleFileInputChange(e: React.ChangeEvent<HTMLInputElement>) {
    const file = e.target.files?.[0];
    if (file) handleFile(file);
    if (fileInputRef.current) fileInputRef.current.value = '';
  }

  function handleZoneClick() {
    fileInputRef.current?.click();
  }

  function handleZoneKeyDown(e: React.KeyboardEvent<HTMLDivElement>) {
    if (e.key === 'Enter' || e.key === ' ') {
      e.preventDefault();
      fileInputRef.current?.click();
    }
  }

  // State: result available
  if (result) {
    const totalRecords = getTotalRecords(result);
    const verdict = getVerdictInfo(result);
    const summaryText =
      totalRecords !== null
        ? `${totalRecords} record${totalRecords !== 1 ? 's' : ''} — ${
            result.issue_count === 0
              ? 'no issues found'
              : `${result.issue_count} issue${result.issue_count !== 1 ? 's' : ''} found`
          }`
        : result.issue_count === 0
          ? 'No issues found'
          : `${result.issue_count} issue${result.issue_count !== 1 ? 's' : ''} found`;

    return (
      <AnimatePresence>
        <motion.div
          key="result"
          initial={{ opacity: 0, y: 16 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, ease: 'easeOut' }}
          className="space-y-3"
        >
          <ClayCard variant="default" size="sm">
            {/* Header: file type badge + summary + clear button */}
            <div className="flex items-start justify-between gap-3 mb-3">
              <div className="flex items-center gap-2 flex-wrap">
                <Badge variant="default" size="md">
                  {result.file_type}
                </Badge>
                <span className="text-sm text-[var(--color-text-secondary)]">
                  {summaryText}
                </span>
              </div>
              <button
                type="button"
                onClick={onClear}
                className="text-xs text-[var(--color-text-muted)] hover:text-status-error transition-colors shrink-0"
              >
                Remove file
              </button>
            </div>

            {/* Verdict badge */}
            <div className="flex items-center gap-2">
              <Badge variant={verdict.variant}>{verdict.text}</Badge>
              {result.encoding && (
                <span className="text-xs text-[var(--color-text-muted)]">
                  Encoding: {result.encoding}
                </span>
              )}
            </div>
          </ClayCard>

          {/* Issue list */}
          {result.issues.length > 0 && (
            <FileIssueList issues={result.issues} />
          )}
        </motion.div>
      </AnimatePresence>
    );
  }

  // State: error
  if (error) {
    return (
      <motion.div
        initial={{ opacity: 0, y: 8 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.25 }}
      >
        <ClayCard variant="flat" size="sm" className="border border-[var(--color-border)]">
          <div className="flex items-start gap-3">
            <AlertTriangle className="w-4 h-4 text-status-error mt-0.5 shrink-0" />
            <div className="flex-1 min-w-0">
              <p className="text-sm font-semibold text-[var(--color-text-primary)] font-display">
                File pre-validation failed
              </p>
              <p className="text-sm text-[var(--color-text-secondary)] mt-0.5">
                File pre-validation failed. Check the file format and retry.
              </p>
            </div>
            <ClayButton variant="ghost" size="sm" onClick={onClear} className="shrink-0">
              Retry
            </ClayButton>
          </div>
        </ClayCard>
      </motion.div>
    );
  }

  // State: loading
  if (isLoading) {
    return (
      <div
        className={cn(
          'flex flex-col items-center justify-center',
          'min-h-[120px] rounded-2xl',
          'border-2 border-dashed border-[var(--color-border)]',
        )}
      >
        <MoleculeLoader size="sm" text="Validating file..." />
      </div>
    );
  }

  // State: no file — drag-and-drop zone
  return (
    <div className="space-y-2">
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
          'flex flex-col items-center justify-center text-center',
          'min-h-[120px] rounded-2xl cursor-pointer px-6',
          'border-2 border-dashed transition-all duration-200',
          'focus:outline-none focus:ring-2 focus:ring-[var(--glow-primary)]',
          isDraggingOver
            ? 'border-chem-primary-600/40 bg-[rgba(var(--color-primary-rgb),0.04)]'
            : 'border-[var(--color-border)] hover:border-[var(--color-border-strong)]',
        )}
      >
        <Upload
          className={cn(
            'w-6 h-6 mb-2 transition-colors duration-200',
            isDraggingOver ? 'text-[var(--color-primary)]' : 'text-[var(--color-text-muted)]',
          )}
        />
        <p className="text-sm text-[var(--color-text-secondary)]">
          {isDraggingOver
            ? 'Release to upload'
            : 'Drop an SDF or CSV file here, or click to browse'}
        </p>
        <p className="text-xs text-[var(--color-text-muted)] mt-1">
          Supported formats: .sdf, .sd, .csv
        </p>
      </div>

      {/* Local validation errors */}
      {localError && (
        <motion.p
          initial={{ opacity: 0, y: -4 }}
          animate={{ opacity: 1, y: 0 }}
          className="text-sm text-status-error"
        >
          {localError}
        </motion.p>
      )}
    </div>
  );
}

