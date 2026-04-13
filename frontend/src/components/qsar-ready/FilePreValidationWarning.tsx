import { ClayButton } from '../ui/ClayButton';
import type { FileIssue } from '../../types/diagnostics';

interface FilePreValidationWarningProps {
  /** List of critical (error-severity) issues to display. */
  issues: FileIssue[];
  /** Total count of critical issues for the heading. */
  issueCount: number;
  /** Called when the user chooses to proceed despite the issues. */
  onProceed: () => void;
  /** Called when the user cancels the upload. */
  onCancel: () => void;
}

/**
 * Warning banner displayed when the Phase 9 file pre-validator detects
 * critical issues (error-severity) in the uploaded file.
 *
 * Per UI-SPEC D-13:
 * - Amber left-border style matches ProvenanceTimeline stereo warning pattern
 * - role="alert" for immediate screen reader announcement
 * - Two actions: "Proceed Anyway" (default) and "Cancel Upload" (ghost)
 */
export function FilePreValidationWarning({
  issues,
  issueCount,
  onProceed,
  onCancel,
}: FilePreValidationWarningProps) {
  return (
    <div
      role="alert"
      className="bg-amber-50 border-l-4 border-amber-400 p-3 rounded-r-lg"
    >
      <div className="flex items-start gap-2 mb-3">
        <span className="text-amber-500 font-bold flex-shrink-0 mt-0.5">!</span>
        <div className="flex-1 min-w-0">
          <p className="text-sm font-semibold text-amber-800">
            Pre-validation found {issueCount} critical issue{issueCount !== 1 ? 's' : ''} in this
            file. Review below before proceeding.
          </p>

          {issues.length > 0 && (
            <ul className="mt-2 space-y-1">
              {issues.map((issue, idx) => (
                <li
                  key={idx}
                  className="text-xs text-amber-700 flex items-start gap-1.5"
                >
                  <span className="mt-0.5 flex-shrink-0">•</span>
                  <span>
                    {issue.block !== null && (
                      <span className="font-medium">Block {issue.block}: </span>
                    )}
                    {issue.description}
                    {issue.severity && (
                      <span className="ml-1 text-amber-600/80">({issue.severity})</span>
                    )}
                  </span>
                </li>
              ))}
            </ul>
          )}
        </div>
      </div>

      <div className="flex items-center gap-2 justify-end">
        <ClayButton variant="ghost" size="sm" onClick={onCancel}>
          Cancel Upload
        </ClayButton>
        <ClayButton variant="default" size="sm" onClick={onProceed}>
          Proceed Anyway
        </ClayButton>
      </div>
    </div>
  );
}
