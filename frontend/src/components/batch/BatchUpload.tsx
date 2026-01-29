import { useState, useCallback, useRef } from 'react';
import { batchApi } from '../../services/api';
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
    // Max 1GB
    if (file.size > 1024 * 1024 * 1024) {
      return 'File too large. Maximum size is 1GB.';
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
          ${isDragging ? 'border-blue-500 bg-blue-50' : 'border-gray-300 hover:border-gray-400'}
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

        <div className="text-gray-600">
          <svg
            className="mx-auto h-12 w-12 text-gray-400 mb-4"
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
          <p className="text-lg font-medium">
            {isDragging ? 'Drop file here' : 'Drop file here or click to browse'}
          </p>
          <p className="text-sm text-gray-500 mt-1">
            Supports SDF and CSV files (up to 1,000,000 molecules)
          </p>
        </div>
      </div>

      {/* Selected file info */}
      {selectedFile && (
        <div className="bg-gray-50 rounded-lg p-4">
          <div className="flex items-center justify-between">
            <div>
              <p className="font-medium text-gray-900">{selectedFile.name}</p>
              <p className="text-sm text-gray-500">
                {formatFileSize(selectedFile.size)}
                {csvColumns && ` | ~${csvColumns.row_count_estimate} molecules`}
              </p>
            </div>
            <button
              onClick={handleReset}
              disabled={isUploading}
              className="text-gray-400 hover:text-gray-600 disabled:opacity-50"
            >
              <svg
                className="h-5 w-5"
                fill="none"
                strokeLinecap="round"
                strokeLinejoin="round"
                strokeWidth={2}
                viewBox="0 0 24 24"
                stroke="currentColor"
              >
                <path d="M6 18L18 6M6 6l12 12" />
              </svg>
            </button>
          </div>

          {/* CSV column selector */}
          {csvColumns && (
            <div className="mt-4">
              <label className="block text-sm font-medium text-gray-700 mb-1">
                Select SMILES column:
              </label>
              <select
                value={selectedSmilesColumn}
                onChange={(e) => setSelectedSmilesColumn(e.target.value)}
                disabled={isUploading}
                className="w-full border border-gray-300 rounded-md px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-blue-500"
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
      <button
        onClick={handleUpload}
        disabled={!selectedFile || isUploading || disabled}
        className={`
          w-full py-3 px-4 rounded-lg font-medium text-white
          transition-colors duration-200
          ${
            selectedFile && !isUploading && !disabled
              ? 'bg-blue-600 hover:bg-blue-700'
              : 'bg-gray-300 cursor-not-allowed'
          }
        `}
      >
        {isUploading ? (
          <span className="flex items-center justify-center">
            <svg
              className="animate-spin -ml-1 mr-3 h-5 w-5 text-white"
              fill="none"
              viewBox="0 0 24 24"
            >
              <circle
                className="opacity-25"
                cx="12"
                cy="12"
                r="10"
                stroke="currentColor"
                strokeWidth="4"
              />
              <path
                className="opacity-75"
                fill="currentColor"
                d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
              />
            </svg>
            Uploading...
          </span>
        ) : (
          'Upload and Process'
        )}
      </button>
    </div>
  );
}
