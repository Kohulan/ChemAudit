import { useState, useCallback, useRef } from 'react';
import { Upload, AlertCircle } from 'lucide-react';
import type { DatasetDiffResults } from '../../types/dataset_intelligence';
import { DiffSummary } from './DiffSummary';
import { DiffMoleculeTable } from './DiffMoleculeTable';

// =============================================================================
// Types
// =============================================================================

interface DatasetDiffTabProps {
  /** Diff results, null before comparison file uploaded. */
  diffResults: DatasetDiffResults | null;
  /** Whether a diff upload is in progress. */
  diffLoading: boolean;
  /** Error message for diff upload. */
  diffError: string | null;
  /** Callback when a comparison file is selected. */
  onUploadDiffFile: (file: File) => void;
  /** Primary job ID (null if no primary dataset loaded). */
  primaryJobId: string | null;
}

// =============================================================================
// Helpers
// =============================================================================

const VALID_EXTENSIONS = ['.csv', '.sdf'];

function isValidExtension(filename: string): boolean {
  const lower = filename.toLowerCase();
  return VALID_EXTENSIONS.some((ext) => lower.endsWith(ext));
}

// =============================================================================
// Component
// =============================================================================

/**
 * Dataset Diff tab content for the Dataset Audit page.
 *
 * Per UI-SPEC D-12, D-13:
 * - Second upload zone for comparison file (reuses DatasetUploadZone visual pattern)
 * - Before upload: prompt message
 * - During processing: loading spinner
 * - After results: DiffSummary + DiffMoleculeTable
 * - Disabled if no primary dataset loaded
 */
export function DatasetDiffTab({
  diffResults,
  diffLoading,
  diffError,
  onUploadDiffFile,
  primaryJobId,
}: DatasetDiffTabProps) {
  const [isDragOver, setIsDragOver] = useState(false);
  const fileInputRef = useRef<HTMLInputElement>(null);
  const isDisabled = !primaryJobId;

  const handleFileChange = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const file = e.target.files?.[0];
      if (file && isValidExtension(file.name)) {
        onUploadDiffFile(file);
      }
      // Reset input so the same file can be re-selected
      e.target.value = '';
    },
    [onUploadDiffFile],
  );

  const handleDrop = useCallback(
    (e: React.DragEvent) => {
      e.preventDefault();
      setIsDragOver(false);
      if (isDisabled || diffLoading) return;
      const file = e.dataTransfer.files[0];
      if (file && isValidExtension(file.name)) {
        onUploadDiffFile(file);
      }
    },
    [isDisabled, diffLoading, onUploadDiffFile],
  );

  const handleDragOver = useCallback(
    (e: React.DragEvent) => {
      e.preventDefault();
      if (!isDisabled && !diffLoading) {
        setIsDragOver(true);
      }
    },
    [isDisabled, diffLoading],
  );

  const handleDragLeave = useCallback(() => {
    setIsDragOver(false);
  }, []);

  return (
    <div className="min-h-[200px] space-y-6">
      {/* Second upload zone for comparison file */}
      <div
        className={[
          'relative border-2 border-dashed rounded-xl p-8 text-center transition-colors duration-200',
          isDisabled
            ? 'border-[var(--color-border)]/50 opacity-50 cursor-not-allowed'
            : isDragOver
              ? 'border-[var(--color-primary)] bg-[var(--color-primary)]/5'
              : 'border-[var(--color-border)] hover:border-[var(--color-primary)]/50 cursor-pointer',
        ].join(' ')}
        onDrop={handleDrop}
        onDragOver={handleDragOver}
        onDragLeave={handleDragLeave}
        onClick={() => {
          if (!isDisabled && !diffLoading) {
            fileInputRef.current?.click();
          }
        }}
        role="button"
        tabIndex={isDisabled ? -1 : 0}
        aria-label="Upload comparison dataset"
        aria-disabled={isDisabled}
      >
        <input
          ref={fileInputRef}
          type="file"
          accept=".csv,.sdf"
          onChange={handleFileChange}
          className="hidden"
          disabled={isDisabled || diffLoading}
        />
        <div className="flex flex-col items-center gap-3">
          <div className="w-12 h-12 rounded-xl bg-[var(--color-surface-sunken)] flex items-center justify-center">
            <Upload className="w-6 h-6 text-[var(--color-text-muted)]" />
          </div>
          <div>
            <p className="text-sm font-medium text-[var(--color-text-primary)]">
              Upload comparison dataset
            </p>
            <p className="text-xs text-[var(--color-text-muted)] mt-1">
              Drop a second CSV or SDF file here to compare against the primary
              dataset
            </p>
          </div>
        </div>
      </div>

      {/* Error display */}
      {diffError && (
        <div className="flex items-start gap-3 p-4 rounded-xl bg-red-500/10 border border-red-500/20">
          <AlertCircle className="w-5 h-5 text-red-500 shrink-0 mt-0.5" />
          <p className="text-sm text-red-600 dark:text-red-400">{diffError}</p>
        </div>
      )}

      {/* Loading state */}
      {diffLoading && (
        <div className="flex items-center justify-center py-12">
          <div className="flex items-center gap-3">
            <div className="w-5 h-5 border-2 border-[var(--color-primary)] border-t-transparent rounded-full animate-spin" />
            <span className="text-sm text-[var(--color-text-secondary)]">
              Comparing datasets...
            </span>
          </div>
        </div>
      )}

      {/* Pre-comparison prompt */}
      {!diffResults && !diffLoading && !diffError && (
        <div className="flex items-center justify-center py-12">
          <p className="text-sm text-[var(--color-text-muted)]">
            Upload a comparison dataset above to see what changed between
            versions.
          </p>
        </div>
      )}

      {/* Diff results */}
      {diffResults && !diffLoading && (
        <>
          {/* Summary badges */}
          <DiffSummary
            addedCount={diffResults.added_count}
            removedCount={diffResults.removed_count}
            modifiedCount={diffResults.modified_count}
            unchangedCount={diffResults.unchanged_count}
          />

          {/* Molecule table with category filters */}
          <DiffMoleculeTable
            molecules={[
              ...diffResults.added,
              ...diffResults.removed,
              ...diffResults.modified,
            ]}
            category="all"
            diffResults={diffResults}
          />
        </>
      )}
    </div>
  );
}
