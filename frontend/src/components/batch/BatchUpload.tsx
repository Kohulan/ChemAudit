import { useState, useCallback, useRef } from 'react';
import { Upload, X } from 'lucide-react';
import { batchApi } from '../../services/api';
import { useLimits } from '../../context/ConfigContext';
import { ClayButton } from '../ui/ClayButton';
import type { CSVColumnsResponse } from '../../types/batch';

interface BatchUploadProps {
  onUploadSuccess: (jobId: string, totalMolecules: number) => void;
  onUploadError: (error: string) => void;
  disabled?: boolean;
}

/**
 * File upload component with drag-and-drop support.
 * Accepts SDF and CSV files, auto-detects CSV columns.
 */
export function BatchUpload({
  onUploadSuccess,
  onUploadError,
  disabled = false,
}: BatchUploadProps) {
  const limits = useLimits();
  const [isDragging, setIsDragging] = useState(false);
  const [selectedFile, setSelectedFile] = useState<File | null>(null);
  const [isUploading, setIsUploading] = useState(false);
  const [csvColumns, setCsvColumns] = useState<CSVColumnsResponse | null>(null);
  const [selectedSmilesColumn, setSelectedSmilesColumn] = useState<string>('');

  const fileInputRef = useRef<HTMLInputElement>(null);

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

  const validateFile = (file: File): string | null => {
    const name = file.name.toLowerCase();
    if (!name.endsWith('.sdf') && !name.endsWith('.csv')) {
      return 'Invalid file type. Please upload an SDF or CSV file.';
    }
    if (file.size > limits.max_file_size_bytes) {
      return `File too large. Maximum size is ${limits.max_file_size_mb}MB.`;
    }
    return null;
  };

  const processFile = async (file: File) => {
    const error = validateFile(file);
    if (error) {
      onUploadError(error);
      return;
    }

    setSelectedFile(file);

    // If CSV, detect columns for SMILES selection
    if (file.name.toLowerCase().endsWith('.csv')) {
      try {
        const columns = await batchApi.detectColumns(file);
        setCsvColumns(columns);
        setSelectedSmilesColumn(columns.suggested_smiles || columns.columns[0] || '');
      } catch (e) {
        onUploadError('Failed to read CSV columns');
        setSelectedFile(null);
      }
    } else {
      setCsvColumns(null);
    }
  };

  const handleDrop = useCallback(
    async (e: React.DragEvent) => {
      e.preventDefault();
      e.stopPropagation();
      setIsDragging(false);

      if (disabled || isUploading) return;

      const files = e.dataTransfer.files;
      if (files.length > 0) {
        await processFile(files[0]);
      }
    },
    [disabled, isUploading]
  );

  const handleFileSelect = useCallback(
    async (e: React.ChangeEvent<HTMLInputElement>) => {
      const files = e.target.files;
      if (files && files.length > 0) {
        await processFile(files[0]);
      }
    },
    []
  );

  const handleUpload = async () => {
    if (!selectedFile || isUploading) return;

    setIsUploading(true);

    try {
      const response = await batchApi.uploadBatch(
        selectedFile,
        csvColumns ? selectedSmilesColumn : undefined
      );
      onUploadSuccess(response.job_id, response.total_molecules);
    } catch (e: any) {
      const errorMessage =
        e.response?.data?.detail || e.message || 'Upload failed';
      onUploadError(errorMessage);
    } finally {
      setIsUploading(false);
    }
  };

  const handleReset = () => {
    setSelectedFile(null);
    setCsvColumns(null);
    setSelectedSmilesColumn('');
    if (fileInputRef.current) {
      fileInputRef.current.value = '';
    }
  };

  const formatFileSize = (bytes: number): string => {
    if (bytes < 1024) return `${bytes} B`;
    if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
    return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
  };

  return (
    <div className="space-y-4">
      {/* Drop zone */}
      <div
        onDragEnter={handleDragEnter}
        onDragLeave={handleDragLeave}
        onDragOver={handleDragOver}
        onDrop={handleDrop}
        onClick={() => !disabled && !isUploading && fileInputRef.current?.click()}
        className={`
          border-2 border-dashed rounded-lg p-8 text-center cursor-pointer
          transition-colors duration-200
          ${isDragging
            ? 'border-[var(--color-primary)] bg-[var(--color-primary)]/10'
            : 'border-[var(--color-border-strong)] hover:border-[var(--color-text-muted)]'}
          ${disabled || isUploading ? 'opacity-50 cursor-not-allowed' : ''}
        `}
      >
        <input
          ref={fileInputRef}
          type="file"
          accept=".sdf,.csv"
          onChange={handleFileSelect}
          className="hidden"
          disabled={disabled || isUploading}
        />

        <div className="text-[var(--color-text-secondary)]">
          <svg
            className="mx-auto h-12 w-12 text-[var(--color-text-muted)] mb-4"
            stroke="currentColor"
            fill="none"
            viewBox="0 0 48 48"
          >
            <path
              strokeLinecap="round"
              strokeLinejoin="round"
              strokeWidth={2}
              d="M8 14v20c0 4.418 7.163 8 16 8s16-3.582 16-8V14m-32 0c0 4.418 7.163 8 16 8s16-3.582 16-8m-32 0c0-4.418 7.163-8 16-8s16 3.582 16 8"
            />
          </svg>
          <p className="text-lg font-medium text-[var(--color-text-primary)]">
            {isDragging ? 'Drop file here' : 'Drop file here or click to browse'}
          </p>
          <p className="text-sm text-[var(--color-text-muted)] mt-1">
            Supports SDF and CSV files (up to {limits.max_batch_size.toLocaleString()} molecules)
          </p>
        </div>
      </div>

      {/* Selected file info */}
      {selectedFile && (
        <div className="bg-[var(--color-surface-sunken)] rounded-lg p-4">
          <div className="flex items-center justify-between">
            <div>
              <p className="font-medium text-[var(--color-text-primary)]">{selectedFile.name}</p>
              <p className="text-sm text-[var(--color-text-muted)]">
                {formatFileSize(selectedFile.size)}
                {csvColumns && ` | ~${csvColumns.row_count_estimate} molecules`}
              </p>
            </div>
            <ClayButton
              variant="ghost"
              size="icon"
              onClick={handleReset}
              disabled={isUploading}
            >
              <X className="h-5 w-5" />
            </ClayButton>
          </div>

          {/* CSV column selector */}
          {csvColumns && (
            <div className="mt-4">
              <label className="block text-sm font-medium text-[var(--color-text-secondary)] mb-1">
                Select SMILES column:
              </label>
              <select
                value={selectedSmilesColumn}
                onChange={(e) => setSelectedSmilesColumn(e.target.value)}
                disabled={isUploading}
                className="w-full border border-[var(--color-border)] bg-[var(--color-surface-elevated)] text-[var(--color-text-primary)] rounded-md px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]"
              >
                {csvColumns.columns.map((col) => (
                  <option key={col} value={col}>
                    {col}
                    {col === csvColumns.suggested_smiles ? ' (suggested)' : ''}
                  </option>
                ))}
              </select>
            </div>
          )}
        </div>
      )}

      {/* Upload button */}
      <ClayButton
        variant="primary"
        size="lg"
        onClick={handleUpload}
        disabled={!selectedFile || isUploading || disabled}
        loading={isUploading}
        leftIcon={!isUploading ? <Upload className="w-5 h-5" /> : undefined}
        className="w-full"
      >
        {isUploading ? 'Uploading...' : 'Upload and Process'}
      </ClayButton>
    </div>
  );
}
