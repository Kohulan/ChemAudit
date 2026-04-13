import { useState, useCallback, useRef } from 'react';
import { Upload, FileSpreadsheet, Database } from 'lucide-react';
import type { DatasetAuditStatus } from '../../types/dataset_intelligence';
import { ClayButton } from '../ui/ClayButton';

// =============================================================================
// Types
// =============================================================================

interface DatasetUploadZoneProps {
  /** Callback when a valid file is selected. */
  onFileSelect: (file: File) => void;
  /** Whether the zone is disabled (e.g. during processing). */
  disabled: boolean;
  /** Currently selected file, or null. */
  file: File | null;
  /** Current audit lifecycle status. */
  status: DatasetAuditStatus;
  /** Processing progress 0-100 (shown when status='processing'). */
  progress?: number;
}

// =============================================================================
// Helpers
// =============================================================================

const VALID_EXTENSIONS = ['.csv', '.sdf'];

function isValidExtension(filename: string): boolean {
  const lower = filename.toLowerCase();
  return VALID_EXTENSIONS.some((ext) => lower.endsWith(ext));
}

function formatFileSize(bytes: number): string {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
  return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
}

// =============================================================================
// Component
// =============================================================================

/**
 * Drag-and-drop file upload zone for dataset audit.
 *
 * Follows the BatchUpload.tsx pattern from Phase 1:
 * - Drop zone with dashed border
 * - States: idle, dragging, file-selected, processing
 * - Accepts CSV and SDF files only
 * - Accessible: role="button", aria-label, tabIndex, Enter/Space triggers
 * - "Upload & Analyze" ClayButton CTA below the zone
 *
 * Per UI-SPEC: p-8 internal padding, dashed border, primary accent on hover/focus.
 */
export function DatasetUploadZone({
  onFileSelect,
  disabled,
  file,
  status,
  progress = 0,
}: DatasetUploadZoneProps) {
  const [isDragging, setIsDragging] = useState(false);
  const [fileError, setFileError] = useState<string | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  const isProcessing = status === 'uploading' || status === 'processing';
  const isDisabled = disabled || isProcessing;

  // Drag handlers
  const handleDragEnter = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    if (!isDisabled) setIsDragging(true);
  }, [isDisabled]);

  const handleDragLeave = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(false);
  }, []);

  const handleDragOver = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
  }, []);

  const processFile = useCallback(
    (f: File) => {
      setFileError(null);
      if (!isValidExtension(f.name)) {
        setFileError('Invalid file type. Please upload a CSV or SDF file.');
        return;
      }
      onFileSelect(f);
    },
    [onFileSelect],
  );

  const handleDrop = useCallback(
    (e: React.DragEvent) => {
      e.preventDefault();
      e.stopPropagation();
      setIsDragging(false);
      if (isDisabled) return;
      const files = e.dataTransfer.files;
      if (files.length > 0) {
        processFile(files[0]);
      }
    },
    [isDisabled, processFile],
  );

  const handleFileChange = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const files = e.target.files;
      if (files && files.length > 0) {
        processFile(files[0]);
      }
      // Reset input so re-selecting the same file triggers onChange
      if (fileInputRef.current) fileInputRef.current.value = '';
    },
    [processFile],
  );

  const handleClick = useCallback(() => {
    if (!isDisabled) fileInputRef.current?.click();
  }, [isDisabled]);

  const handleKeyDown = useCallback(
    (e: React.KeyboardEvent) => {
      if ((e.key === 'Enter' || e.key === ' ') && !isDisabled) {
        e.preventDefault();
        fileInputRef.current?.click();
      }
    },
    [isDisabled],
  );

  const isSdf = file?.name.toLowerCase().endsWith('.sdf');

  return (
    <div className="space-y-4">
      {/* Drop zone */}
      <div
        role="button"
        aria-label="Upload dataset file"
        tabIndex={isDisabled ? -1 : 0}
        onDragEnter={handleDragEnter}
        onDragLeave={handleDragLeave}
        onDragOver={handleDragOver}
        onDrop={handleDrop}
        onClick={handleClick}
        onKeyDown={handleKeyDown}
        className={[
          'border-2 border-dashed rounded-xl p-8 text-center transition-all duration-200',
          isDragging
            ? 'border-[var(--color-primary)] bg-[var(--color-primary)]/10 scale-[1.02]'
            : 'border-[var(--color-border)] hover:border-[var(--color-primary)]',
          isDisabled ? 'opacity-50 cursor-not-allowed' : 'cursor-pointer',
          'focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-[var(--color-primary)] focus-visible:ring-offset-2',
        ].join(' ')}
      >
        <input
          ref={fileInputRef}
          type="file"
          accept=".csv,.sdf"
          onChange={handleFileChange}
          className="hidden"
          disabled={isDisabled}
        />

        {/* Idle state */}
        {!file && !isProcessing && (
          <div className="text-[var(--color-text-secondary)]">
            <div className="mx-auto w-16 h-16 rounded-2xl bg-gradient-to-br from-[var(--color-primary)]/10 to-[var(--color-accent)]/10 flex items-center justify-center mb-4">
              <Upload className="w-8 h-8 text-[var(--color-primary)]" />
            </div>
            <p className="text-lg font-medium text-[var(--color-text-primary)]">
              {isDragging
                ? 'Drop file here'
                : 'Drop a CSV or SDF file here, or click to browse'}
            </p>
            <p className="text-sm text-[var(--color-text-muted)] mt-2">
              Supports <span className="font-medium">CSV</span> and{' '}
              <span className="font-medium">SDF</span> files
            </p>
          </div>
        )}

        {/* File selected state */}
        {file && !isProcessing && (
          <div className="flex items-center justify-center gap-3">
            <div className="w-10 h-10 rounded-lg bg-gradient-to-br from-[var(--color-primary)]/10 to-[var(--color-accent)]/10 flex items-center justify-center">
              {isSdf ? (
                <Database className="w-5 h-5 text-[var(--color-primary)]" />
              ) : (
                <FileSpreadsheet className="w-5 h-5 text-[var(--color-primary)]" />
              )}
            </div>
            <div className="text-left">
              <p className="font-medium text-[var(--color-text-primary)]">
                {file.name}{' '}
                <span className="text-[var(--color-text-muted)] font-normal">
                  ({formatFileSize(file.size)})
                </span>
              </p>
              <p className="text-xs text-[var(--color-text-muted)]">
                Click to change file
              </p>
            </div>
          </div>
        )}

        {/* Processing state */}
        {file && isProcessing && (
          <div className="text-[var(--color-text-secondary)]">
            <div className="animate-spin w-8 h-8 border-2 border-[var(--color-primary)] border-t-transparent rounded-full mx-auto mb-3" />
            <p className="text-sm font-medium text-[var(--color-text-primary)]">
              Analyzing {file.name}... {Math.round(progress)}%
            </p>
          </div>
        )}
      </div>

      {/* File validation error */}
      {fileError && (
        <p className="text-sm text-[var(--status-error)]">{fileError}</p>
      )}

      {/* Upload & Analyze button */}
      {file && status !== 'complete' && (
        <div className="flex justify-end">
          <ClayButton
            variant="primary"
            onClick={(e) => {
              e.stopPropagation();
              onFileSelect(file);
            }}
            disabled={!file || isProcessing}
            loading={isProcessing}
          >
            Upload &amp; Analyze
          </ClayButton>
        </div>
      )}
    </div>
  );
}
