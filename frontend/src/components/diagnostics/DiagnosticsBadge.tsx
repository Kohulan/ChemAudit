import { Stethoscope } from 'lucide-react';

interface DiagnosticsBadgeProps {
  smiles: string;
}

/**
 * Cross-link badge from SingleValidation to /diagnostics (D-02).
 *
 * Per UI-SPEC Cross-Link contract:
 * - Compact inline element with stethoscope icon, text, and link
 * - Navigates to /diagnostics?smiles={encodedSmiles}
 * - Renders only when smiles is a non-empty string
 */
export function DiagnosticsBadge({ smiles }: DiagnosticsBadgeProps) {
  if (!smiles) return null;

  const href = `/diagnostics?smiles=${encodeURIComponent(smiles)}`;

  return (
    <div className="flex items-center gap-2 p-3 rounded-lg border border-[var(--color-border)] bg-[var(--color-surface-sunken)]">
      {/* Left: stethoscope icon */}
      <Stethoscope className="w-4 h-4 text-chem-primary-600 shrink-0" />

      {/* Center: label */}
      <span className="text-sm text-[var(--color-text-secondary)] flex-1 min-w-0">
        Structure diagnostics available
      </span>

      {/* Right: Diagnose link */}
      <a
        href={href}
        className="text-sm text-chem-primary-600 hover:underline shrink-0 whitespace-nowrap"
      >
        Diagnose &rarr;
      </a>
    </div>
  );
}
