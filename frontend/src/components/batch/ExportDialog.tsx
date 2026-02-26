import React, { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { Download, X, FileText, Table2, FlaskConical, Braces, FileBarChart, Fingerprint, Copy, Network, LayoutGrid, Image } from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';
import { Badge } from '../ui/Badge';
import { EXPORT_FORMATS, PDF_SECTION_OPTIONS } from '../../types/workflow';
import type { ExportFormat } from '../../types/workflow';

const FORMAT_ICONS: Record<ExportFormat, React.ReactNode> = {
  csv: <Table2 className="w-5 h-5" />,
  excel: <FileBarChart className="w-5 h-5" />,
  sdf: <FlaskConical className="w-5 h-5" />,
  json: <Braces className="w-5 h-5" />,
  pdf: <FileText className="w-5 h-5" />,
  fingerprint: <Fingerprint className="w-5 h-5" />,
  dedup: <Copy className="w-5 h-5" />,
  scaffold: <Network className="w-5 h-5" />,
  property_matrix: <LayoutGrid className="w-5 h-5" />,
};

interface ExportDialogProps {
  jobId: string;
  isOpen: boolean;
  onClose: () => void;
  selectedIndices?: Set<number>;
}

/**
 * Enhanced export dialog with all 9 formats and PDF section selection.
 *
 * Formats: CSV, Excel, SDF, JSON, PDF, Fingerprints, Deduplicated,
 * Scaffold-Grouped, Property Matrix.
 */
export function ExportDialog({ jobId, isOpen, onClose, selectedIndices }: ExportDialogProps) {
  const [selectedFormat, setSelectedFormat] = useState<ExportFormat>('csv');
  const [isExporting, setIsExporting] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [includeImages, setIncludeImages] = useState(false);
  const [pdfSections, setPdfSections] = useState<Set<string>>(
    new Set(PDF_SECTION_OPTIONS.filter((s) => s.defaultChecked).map((s) => s.id))
  );

  if (!isOpen) return null;

  const isSelectedMode = selectedIndices && selectedIndices.size > 0;
  const selectedFormatInfo = EXPORT_FORMATS.find((f) => f.value === selectedFormat);
  const dateStr = new Date().toISOString().slice(0, 10).replace(/-/g, '');
  const filePreview = `chemaudit_${selectedFormat}_${dateStr}_${jobId.slice(0, 8)}.${selectedFormatInfo?.extension ?? 'dat'}`;

  const togglePdfSection = (id: string) => {
    setPdfSections((prev) => {
      const next = new Set(prev);
      if (next.has(id)) {
        next.delete(id);
      } else {
        next.add(id);
      }
      return next;
    });
  };

  const handleExport = async () => {
    setIsExporting(true);
    setError(null);

    try {
      const params = new URLSearchParams({ format: selectedFormat });

      if (selectedFormat === 'pdf' && pdfSections.size > 0) {
        params.set('sections', Array.from(pdfSections).join(','));
      }
      if (selectedFormat === 'excel' && includeImages) {
        params.set('include_images', 'true');
      }

      const indicesArray = isSelectedMode
        ? Array.from(selectedIndices).sort((a, b) => a - b)
        : null;

      let fetchOptions: RequestInit | undefined;
      if (indicesArray && indicesArray.length > 200) {
        // Large selection: send indices in POST body to avoid URL length limits
        fetchOptions = {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ indices: indicesArray }),
        };
      } else if (indicesArray) {
        params.set('indices', indicesArray.join(','));
      }

      const exportUrl = `/api/v1/batch/${jobId}/export?${params.toString()}`;
      const response = await fetch(exportUrl, fetchOptions);

      if (!response.ok) {
        const data = await response.json();
        throw new Error(data.detail || 'Export failed');
      }

      const contentDisposition = response.headers.get('Content-Disposition');
      let filename = filePreview;
      if (contentDisposition) {
        const match = contentDisposition.match(/filename="?([^"]+)"?/);
        if (match) filename = match[1];
      }

      const blob = await response.blob();
      const downloadUrl = window.URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = downloadUrl;
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      window.URL.revokeObjectURL(downloadUrl);

      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Export failed');
    } finally {
      setIsExporting(false);
    }
  };

  return (
    <div className="fixed inset-0 bg-black/50 flex items-center justify-center z-50" onClick={onClose}>
      <motion.div
        initial={{ opacity: 0, scale: 0.95, y: 10 }}
        animate={{ opacity: 1, scale: 1, y: 0 }}
        exit={{ opacity: 0, scale: 0.95, y: 10 }}
        transition={{ duration: 0.2 }}
        className="bg-[var(--color-surface-elevated)] rounded-2xl shadow-xl max-w-2xl w-full mx-4 max-h-[90vh] overflow-hidden border border-[var(--color-border)] flex flex-col"
        onClick={(e) => e.stopPropagation()}
      >
        {/* Header */}
        <div className="border-b border-[var(--color-border)] px-6 py-4 flex items-center justify-between flex-shrink-0">
          <div>
            <h2 className="text-xl font-semibold text-[var(--color-text-primary)] font-display">Export Results</h2>
            <p className="text-sm text-[var(--color-text-muted)] mt-1">
              {isSelectedMode ? (
                <>Export {selectedIndices.size} selected molecules</>
              ) : (
                <>Choose a format to download batch results</>
              )}
            </p>
          </div>
          <button
            onClick={onClose}
            className="p-2 rounded-lg hover:bg-[var(--color-surface-sunken)] transition-colors"
          >
            <X className="w-5 h-5 text-[var(--color-text-muted)]" />
          </button>
        </div>

        {/* Scrollable content */}
        <div className="relative overflow-y-auto flex-1 px-6 py-4">
          {/* Format selection grid */}
          <div className="grid grid-cols-1 gap-2">
            {EXPORT_FORMATS.map((format) => (
              <React.Fragment key={format.value}>
                <label
                  className={`
                    flex items-center p-3.5 border-2 rounded-xl cursor-pointer
                    transition-all duration-200
                    ${
                      selectedFormat === format.value
                        ? 'border-[var(--color-primary)] bg-[var(--color-primary)]/5 shadow-[0_0_12px_var(--glow-primary)]'
                        : 'border-[var(--color-border)] hover:border-[var(--color-border-strong)] hover:bg-[var(--color-surface-sunken)]'
                    }
                  `}
                >
                  <input
                    type="radio"
                    name="format"
                    value={format.value}
                    checked={selectedFormat === format.value}
                    onChange={(e) => setSelectedFormat(e.target.value as ExportFormat)}
                    className="sr-only"
                  />
                  <div className={`
                    w-10 h-10 rounded-xl flex items-center justify-center flex-shrink-0 mr-3
                    ${selectedFormat === format.value
                      ? 'bg-[var(--color-primary)]/15 text-[var(--color-primary)]'
                      : 'bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)]'
                    }
                  `}>
                    {FORMAT_ICONS[format.value]}
                  </div>
                  <div className="flex-1 min-w-0">
                    <div className="flex items-center gap-2">
                      <span className="font-medium text-[var(--color-text-primary)] text-sm">{format.label}</span>
                      <span className="text-[10px] text-[var(--color-text-muted)] font-mono uppercase">.{format.extension}</span>
                      {format.isNew && <Badge variant="info" size="sm">New</Badge>}
                    </div>
                    <p className="text-xs text-[var(--color-text-secondary)] mt-0.5 truncate">{format.description}</p>
                  </div>
                  {selectedFormat === format.value && (
                    <div className="w-5 h-5 rounded-full bg-[var(--color-primary)] flex items-center justify-center flex-shrink-0 ml-2">
                      <svg className="w-3 h-3 text-white" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="3">
                        <polyline points="20 6 9 17 4 12" />
                      </svg>
                    </div>
                  )}
                </label>

                {/* Excel structure images option — inline, after the Excel format card */}
                {format.value === 'excel' && selectedFormat === 'excel' && (
                  <AnimatePresence>
                    <motion.div
                      key="excel-images"
                      initial={{ opacity: 0, height: 0 }}
                      animate={{ opacity: 1, height: 'auto' }}
                      exit={{ opacity: 0, height: 0 }}
                      transition={{ duration: 0.2 }}
                    >
                      <div className="border border-[var(--color-border)] rounded-xl p-4 bg-[var(--color-surface-sunken)]/50">
                        <label className="flex items-start gap-3 cursor-pointer">
                          <input
                            type="checkbox"
                            checked={includeImages}
                            onChange={(e) => setIncludeImages(e.target.checked)}
                            className="mt-0.5 rounded border-[var(--color-border-strong)] text-[var(--color-primary)] focus:ring-[var(--color-primary)]/30"
                          />
                          <div className="flex-1">
                            <div className="flex items-center gap-2">
                              <Image className="w-4 h-4 text-[var(--color-text-secondary)]" />
                              <span className="text-sm font-medium text-[var(--color-text-primary)]">Include 2D structure images</span>
                            </div>
                            <p className="text-[10px] text-[var(--color-text-muted)] leading-tight mt-1 ml-6">
                              Embeds RDKit-rendered molecule images in each row. Increases file size and export time for large batches.
                            </p>
                          </div>
                        </label>
                      </div>
                    </motion.div>
                  </AnimatePresence>
                )}

                {/* PDF Section Selection — inline, immediately after the PDF format card */}
                {format.value === 'pdf' && selectedFormat === 'pdf' && (
                  <AnimatePresence>
                    <motion.div
                      key="pdf-sections"
                      initial={{ opacity: 0, height: 0 }}
                      animate={{ opacity: 1, height: 'auto' }}
                      exit={{ opacity: 0, height: 0 }}
                      transition={{ duration: 0.2 }}
                    >
                      <div className="border border-[var(--color-border)] rounded-xl p-4 bg-[var(--color-surface-sunken)]/50">
                        <div className="flex items-center justify-between mb-3">
                          <h3 className="text-sm font-medium text-[var(--color-text-primary)]">Report Sections</h3>
                          <button
                            onClick={() => {
                              if (pdfSections.size === PDF_SECTION_OPTIONS.length) {
                                setPdfSections(new Set());
                              } else {
                                setPdfSections(new Set(PDF_SECTION_OPTIONS.map((s) => s.id)));
                              }
                            }}
                            className="text-xs text-[var(--color-primary)] hover:underline"
                          >
                            {pdfSections.size === PDF_SECTION_OPTIONS.length ? 'Deselect All' : 'Select All'}
                          </button>
                        </div>
                        <div className="grid grid-cols-1 sm:grid-cols-2 gap-2">
                          {PDF_SECTION_OPTIONS.map((section) => (
                            <label
                              key={section.id}
                              className="flex items-start gap-2 p-2 rounded-lg hover:bg-[var(--color-surface-elevated)] transition-colors cursor-pointer"
                            >
                              <input
                                type="checkbox"
                                checked={pdfSections.has(section.id)}
                                onChange={() => togglePdfSection(section.id)}
                                className="mt-0.5 rounded border-[var(--color-border-strong)] text-[var(--color-primary)] focus:ring-[var(--color-primary)]/30"
                              />
                              <div>
                                <span className="text-sm font-medium text-[var(--color-text-primary)]">{section.label}</span>
                                <p className="text-[10px] text-[var(--color-text-muted)] leading-tight mt-0.5">{section.description}</p>
                              </div>
                            </label>
                          ))}
                        </div>
                      </div>
                    </motion.div>
                  </AnimatePresence>
                )}
              </React.Fragment>
            ))}
          </div>

          {/* File name preview */}
          <div className="mt-4 flex items-center gap-2 px-3 py-2 bg-[var(--color-surface-sunken)] rounded-lg">
            <Download className="w-4 h-4 text-[var(--color-text-muted)] flex-shrink-0" />
            <code className="text-xs text-[var(--color-text-secondary)] font-mono truncate">{filePreview}</code>
          </div>

          {/* Error message */}
          {error && (
            <div className="mt-4 p-3 bg-red-500/10 border border-red-500/30 rounded-lg">
              <p className="text-sm text-red-600 dark:text-red-400">{error}</p>
            </div>
          )}
        </div>

        {/* Actions */}
        <div className="border-t border-[var(--color-border)] px-6 py-4 flex justify-end gap-3 flex-shrink-0">
          <ClayButton variant="default" onClick={onClose} disabled={isExporting}>
            Cancel
          </ClayButton>
          <ClayButton
            variant="primary"
            onClick={handleExport}
            disabled={isExporting || (selectedFormat === 'pdf' && pdfSections.size === 0)}
            loading={isExporting}
            leftIcon={!isExporting ? <Download className="w-4 h-4" /> : undefined}
          >
            {isExporting ? 'Exporting...' : `Download ${selectedFormatInfo?.label ?? selectedFormat.toUpperCase()}`}
          </ClayButton>
        </div>
      </motion.div>
    </div>
  );
}
