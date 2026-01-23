import { useState } from 'react';

/**
 * Export format type matching backend ExportFormat enum
 */
export type ExportFormat = 'csv' | 'excel' | 'sdf' | 'json' | 'pdf';

interface ExportDialogProps {
  jobId: string;
  isOpen: boolean;
  onClose: () => void;
}

/**
 * Dialog for selecting export format and downloading batch results.
 *
 * Supports five formats:
 * - CSV: Plain text data
 * - Excel: Formatted spreadsheet with conditional coloring
 * - SDF: Chemical structure file with properties
 * - JSON: Programmatic access with full metadata
 * - PDF: Professional report with charts and images
 */
export function ExportDialog({ jobId, isOpen, onClose }: ExportDialogProps) {
  const [selectedFormat, setSelectedFormat] = useState<ExportFormat>('csv');
  const [isExporting, setIsExporting] = useState(false);
  const [error, setError] = useState<string | null>(null);

  if (!isOpen) return null;

  const handleExport = async () => {
    setIsExporting(true);
    setError(null);

    try {
      // Build export URL
      const url = `/api/v1/batch/${jobId}/export?format=${selectedFormat}`;

      // Fetch file
      const response = await fetch(url);

      if (!response.ok) {
        const data = await response.json();
        throw new Error(data.detail || 'Export failed');
      }

      // Get filename from Content-Disposition header or generate default
      const contentDisposition = response.headers.get('Content-Disposition');
      let filename = `batch_${jobId}.${selectedFormat}`;
      if (contentDisposition) {
        const match = contentDisposition.match(/filename="?([^"]+)"?/);
        if (match) filename = match[1];
      }

      // Download file
      const blob = await response.blob();
      const downloadUrl = window.URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = downloadUrl;
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      window.URL.revokeObjectURL(downloadUrl);

      // Close dialog on success
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Export failed');
    } finally {
      setIsExporting(false);
    }
  };

  const formats: Array<{
    value: ExportFormat;
    label: string;
    description: string;
    icon: string;
  }> = [
    {
      value: 'csv',
      label: 'CSV',
      description: 'Plain text data - Excel, Google Sheets compatible',
      icon: '📄',
    },
    {
      value: 'excel',
      label: 'Excel',
      description: 'Formatted spreadsheet with colors and summary sheet',
      icon: '📊',
    },
    {
      value: 'sdf',
      label: 'SDF',
      description: 'Chemical structure file with attached properties',
      icon: '🧪',
    },
    {
      value: 'json',
      label: 'JSON',
      description: 'Programmatic access with full metadata',
      icon: '{ }',
    },
    {
      value: 'pdf',
      label: 'PDF',
      description: 'Professional report with charts and molecule images',
      icon: '📑',
    },
  ];

  return (
    <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50">
      <div className="bg-white rounded-lg shadow-xl max-w-2xl w-full mx-4 max-h-[90vh] overflow-y-auto">
        {/* Header */}
        <div className="border-b border-gray-200 px-6 py-4">
          <h2 className="text-xl font-semibold text-gray-900">Export Results</h2>
          <p className="text-sm text-gray-500 mt-1">
            Choose a format to download batch validation results
          </p>
        </div>

        {/* Format selection */}
        <div className="px-6 py-4 space-y-3">
          {formats.map((format) => (
            <label
              key={format.value}
              className={`
                flex items-start p-4 border-2 rounded-lg cursor-pointer
                transition-colors hover:bg-gray-50
                ${
                  selectedFormat === format.value
                    ? 'border-blue-500 bg-blue-50'
                    : 'border-gray-200'
                }
              `}
            >
              <input
                type="radio"
                name="format"
                value={format.value}
                checked={selectedFormat === format.value}
                onChange={(e) => setSelectedFormat(e.target.value as ExportFormat)}
                className="mt-1 text-blue-600 focus:ring-blue-500"
              />
              <div className="ml-3 flex-1">
                <div className="flex items-center">
                  <span className="text-2xl mr-2">{format.icon}</span>
                  <span className="font-medium text-gray-900">{format.label}</span>
                </div>
                <p className="text-sm text-gray-600 mt-1">{format.description}</p>
              </div>
            </label>
          ))}
        </div>

        {/* Error message */}
        {error && (
          <div className="mx-6 mb-4 p-3 bg-red-50 border border-red-200 rounded-lg">
            <p className="text-sm text-red-800">{error}</p>
          </div>
        )}

        {/* Actions */}
        <div className="border-t border-gray-200 px-6 py-4 flex justify-end space-x-3">
          <button
            onClick={onClose}
            disabled={isExporting}
            className="px-4 py-2 text-sm font-medium text-gray-700 bg-white border border-gray-300 rounded-lg hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed"
          >
            Cancel
          </button>
          <button
            onClick={handleExport}
            disabled={isExporting}
            className="px-4 py-2 text-sm font-medium text-white bg-blue-600 rounded-lg hover:bg-blue-700 disabled:opacity-50 disabled:cursor-not-allowed flex items-center"
          >
            {isExporting ? (
              <>
                <svg
                  className="animate-spin -ml-1 mr-2 h-4 w-4 text-white"
                  xmlns="http://www.w3.org/2000/svg"
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
                Exporting...
              </>
            ) : (
              <>Download {selectedFormat.toUpperCase()}</>
            )}
          </button>
        </div>
      </div>
    </div>
  );
}
