import { Download } from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';

interface QSARDownloadPanelProps {
  /** Called with the selected download format. */
  onDownload: (format: 'csv' | 'sdf' | 'json') => void;
  /** When true, all download buttons are disabled (batch still processing). */
  disabled: boolean;
}

/**
 * Download panel for QSAR-Ready batch results.
 *
 * Per UI-SPEC D-12 and copywriting contract:
 * - Three ClayButton components aligned right
 * - Labels: "Download CSV", "Download SDF", "Download JSON (full provenance)"
 * - Left icon: Download from lucide-react
 * - aria-label on each button for accessibility
 * - All disabled while batch is still processing
 */
export function QSARDownloadPanel({ onDownload, disabled }: QSARDownloadPanelProps) {
  return (
    <div className="flex items-center justify-end gap-2 flex-wrap">
      <ClayButton
        variant="default"
        size="md"
        leftIcon={<Download className="w-4 h-4" />}
        onClick={() => onDownload('csv')}
        disabled={disabled}
        aria-label="Download results as csv"
      >
        Download CSV
      </ClayButton>

      <ClayButton
        variant="default"
        size="md"
        leftIcon={<Download className="w-4 h-4" />}
        onClick={() => onDownload('sdf')}
        disabled={disabled}
        aria-label="Download results as sdf"
      >
        Download SDF
      </ClayButton>

      <ClayButton
        variant="default"
        size="md"
        leftIcon={<Download className="w-4 h-4" />}
        onClick={() => onDownload('json')}
        disabled={disabled}
        aria-label="Download results as json"
      >
        Download JSON (full provenance)
      </ClayButton>
    </div>
  );
}
