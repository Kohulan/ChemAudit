import { useState, useCallback, useRef } from 'react';
import { Upload } from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';
import { diagnosticsApi } from '../../services/api';
import { FilePreValidationWarning } from './FilePreValidationWarning';
import type { QSARReadyConfig } from '../../types/qsar_ready';
import type { FilePreValidationResponse } from '../../types/diagnostics';

interface QSARBatchInputProps {
  /** Called when the user submits the batch input (file or SMILES text). */
  onSubmit: (file: File | null, smilesText: string | null, config: QSARReadyConfig) => void;
  config: QSARReadyConfig;
  loading: boolean;
}

/**
 * Batch input component for the QSAR-Ready Pipeline.
 *
 * Provides two input modes:
 * 1. Textarea for pasting SMILES (one per line)
 * 2. Drag-and-drop / click-to-browse file upload zone (.csv, .sdf)
 *
 * On file selection, auto-runs Phase 9 pre-validator (D-13):
 * - critical issues (error severity) → show FilePreValidationWarning with Proceed/Cancel
 * - warnings only → show amber info banner briefly, proceed automatically
 * - no issues → submit immediately
 */
export function QSARBatchInput({ onSubmit, config, loading }: QSARBatchInputProps) {
  const [smilesText, setSmilesText] = useState('');
  const [file, setFile] = useState<File | null>(null);
  const [isDragging, setIsDragging] = useState(false);
  const [showPrevalidation, setShowPrevalidation] = useState(false);
  const [prevalidationResult, setPrevalidationResult] = useState<FilePreValidationResponse | null>(null);
  const [pendingSubmit, setPendingSubmit] = useState(false);
  const [prevalidating, setPrevalidating] = useState(false);
  const [warningMessage, setWarningMessage] = useState<string | null>(null);

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

  // ── File processing: pre-validate then decide ──

  const processFile = useCallback(
    async (selectedFile: File) => {
      const name = selectedFile.name.toLowerCase();
      if (!name.endsWith('.csv') && !name.endsWith('.sdf')) {
        setWarningMessage('Invalid file type. Please upload a CSV or SDF file.');
        return;
      }

      setFile(selectedFile);
      setShowPrevalidation(false);
      setPrevalidationResult(null);
      setPendingSubmit(false);
      setWarningMessage(null);
      setPrevalidating(true);

      try {
        const result = await diagnosticsApi.filePrevalidate(selectedFile);
        setPrevalidationResult(result);

        const hasCritical = result.issues.some((i) => i.severity.toLowerCase() === 'error');

        if (hasCritical) {
          // Stop — show warning banner with Proceed/Cancel
          setShowPrevalidation(true);
          setPendingSubmit(true);
        } else if (result.issue_count > 0) {
          // Warnings only — show brief info banner but continue
          setWarningMessage(
            `Note: ${result.issue_count} warning(s) found in file. Proceeding with upload.`,
          );
          onSubmit(selectedFile, null, config);
        } else {
          // No issues — submit immediately
          onSubmit(selectedFile, null, config);
        }
      } catch {
        // Pre-validation failed — proceed anyway (non-blocking)
        onSubmit(selectedFile, null, config);
      } finally {
        setPrevalidating(false);
      }
    },
    [config, onSubmit],
  );

  const handleDrop = useCallback(
    async (e: React.DragEvent) => {
      e.preventDefault();
      e.stopPropagation();
      setIsDragging(false);

      if (loading || prevalidating) return;

      const files = e.dataTransfer.files;
      if (files.length > 0) {
        await processFile(files[0]);
      }
    },
    [loading, prevalidating, processFile],
  );

  const handleFileSelect = useCallback(
    async (e: React.ChangeEvent<HTMLInputElement>) => {
      const files = e.target.files;
      if (files && files.length > 0) {
        await processFile(files[0]);
      }
      // Reset input so re-selecting same file triggers onChange
      if (fileInputRef.current) fileInputRef.current.value = '';
    },
    [processFile],
  );

  // ── Submit handlers ──

  const handleSubmit = () => {
    if (loading || prevalidating) return;
    if (file) {
      onSubmit(file, null, config);
    } else if (smilesText.trim()) {
      onSubmit(null, smilesText.trim(), config);
    }
  };

  const handleProceedAnyway = () => {
    setShowPrevalidation(false);
    setPendingSubmit(false);
    if (file) {
      onSubmit(file, null, config);
    }
  };

  const handleCancelUpload = () => {
    setShowPrevalidation(false);
    setPendingSubmit(false);
    setFile(null);
    setPrevalidationResult(null);
    setWarningMessage(null);
  };

  const canSubmit =
    !loading &&
    !prevalidating &&
    !pendingSubmit &&
    (!!file || smilesText.trim().length > 0);

  const criticalIssues = prevalidationResult?.issues.filter(
    (i) => i.severity.toLowerCase() === 'error',
  ) ?? [];

  return (
    <div className="space-y-4">
      {/* ── SMILES textarea ── */}
      <div>
        <label
          htmlFor="qsar-smiles-input"
          className="block text-sm font-medium text-[var(--color-text-primary)] mb-2"
        >
          Paste SMILES
        </label>
        <textarea
          id="qsar-smiles-input"
          rows={6}
          value={smilesText}
          onChange={(e) => setSmilesText(e.target.value)}
          disabled={loading}
          placeholder="Paste SMILES below or drop a CSV or SDF file to begin batch curation."
          className="w-full border border-[var(--color-border)] bg-[var(--color-surface-elevated)] text-[var(--color-text-primary)] rounded-lg px-4 py-3 text-sm font-mono resize-none focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)] focus:border-transparent placeholder:text-[var(--color-text-muted)] placeholder:font-sans disabled:opacity-50"
        />
        <p className="text-xs text-[var(--color-text-muted)] mt-1">
          One SMILES per line
        </p>
      </div>

      {/* ── OR divider ── */}
      <div className="flex items-center gap-3">
        <div className="flex-1 h-px bg-[var(--color-border)]" />
        <span className="text-xs text-[var(--color-text-muted)] font-medium uppercase tracking-wide">
          or
        </span>
        <div className="flex-1 h-px bg-[var(--color-border)]" />
      </div>

      {/* ── File drop zone ── */}
      <div
        onDragEnter={handleDragEnter}
        onDragLeave={handleDragLeave}
        onDragOver={handleDragOver}
        onDrop={handleDrop}
        onClick={() => !loading && !prevalidating && fileInputRef.current?.click()}
        className={[
          'border-2 border-dashed rounded-lg p-8 text-center cursor-pointer transition-all duration-200',
          isDragging
            ? 'border-[var(--color-primary)] bg-[var(--color-primary)]/5'
            : 'border-[var(--color-border-strong)] hover:border-[var(--color-text-muted)] hover:bg-[var(--color-surface-sunken)]/50',
          loading || prevalidating ? 'opacity-50 cursor-not-allowed' : '',
        ]
          .filter(Boolean)
          .join(' ')}
      >
        <input
          ref={fileInputRef}
          type="file"
          accept=".csv,.sdf"
          onChange={handleFileSelect}
          className="hidden"
          disabled={loading || prevalidating}
        />

        <div className="flex flex-col items-center gap-3 text-[var(--color-text-secondary)]">
          <div className="w-12 h-12 rounded-xl bg-[var(--color-surface-sunken)] flex items-center justify-center">
            <Upload className="w-6 h-6 text-[var(--color-text-muted)]" />
          </div>
          {prevalidating ? (
            <>
              <p className="text-sm font-medium text-[var(--color-text-primary)]">
                Validating file…
              </p>
              <div className="w-5 h-5 border-2 border-[var(--color-primary)] border-t-transparent rounded-full animate-spin" />
            </>
          ) : file ? (
            <>
              <p className="text-sm font-medium text-[var(--color-text-primary)]">
                {file.name}
              </p>
              <p className="text-xs text-[var(--color-text-muted)]">
                Click to change file
              </p>
            </>
          ) : (
            <>
              <p className="text-sm font-medium text-[var(--color-text-primary)]">
                {isDragging ? 'Drop file here' : 'Drop a CSV or SDF file here, or click to browse'}
              </p>
              <p className="text-xs text-[var(--color-text-muted)]">
                Accepts <span className="font-medium">.csv</span> and{' '}
                <span className="font-medium">.sdf</span> files
              </p>
            </>
          )}
        </div>
      </div>

      {/* ── Warnings-only info banner ── */}
      {warningMessage && (
        <div className="bg-amber-50 border border-amber-200 rounded-lg px-4 py-3 text-sm text-amber-700">
          {warningMessage}
        </div>
      )}

      {/* ── Critical issues warning ── */}
      {showPrevalidation && prevalidationResult && (
        <FilePreValidationWarning
          issues={criticalIssues}
          issueCount={criticalIssues.length}
          onProceed={handleProceedAnyway}
          onCancel={handleCancelUpload}
        />
      )}

      {/* ── Submit button ── */}
      <ClayButton
        variant="primary"
        size="lg"
        onClick={handleSubmit}
        disabled={!canSubmit}
        loading={loading}
        className="w-full"
      >
        Upload and Curate
      </ClayButton>
    </div>
  );
}
