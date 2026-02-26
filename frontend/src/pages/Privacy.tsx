import { useState } from 'react';
import { motion } from 'framer-motion';
import { Shield, Database, Cookie, Eye, Lock, Server, Trash2, AlertTriangle, Clock } from 'lucide-react';
import { isAxiosError } from 'axios';
import { sessionApi } from '../services/api';
import { clearAllSnapshots } from '../lib/bookmarkStore';

function extractErrorMessage(err: unknown): string {
  if (isAxiosError(err)) {
    if (err.response?.data?.detail) {
      return err.response.data.detail;
    }
    if (err.response?.status === 404) {
      return 'No active session found. Nothing to delete.';
    }
  }
  if (err instanceof Error) {
    return err.message;
  }
  return 'Failed to purge data';
}

/**
 * Privacy Policy Page
 * Transparent disclosure of data practices including session-scoped storage.
 */
export function PrivacyPage() {
  const [purging, setPurging] = useState(false);
  const [purgeResult, setPurgeResult] = useState<{
    bookmarks: number;
    history: number;
  } | null>(null);
  const [purgeError, setPurgeError] = useState<string | null>(null);
  const [showConfirm, setShowConfirm] = useState(false);

  const handlePurge = async () => {
    setPurging(true);
    setPurgeError(null);
    setPurgeResult(null);
    try {
      const result = await sessionApi.purgeMyData();
      await clearAllSnapshots();
      setPurgeResult(result.deleted);
      setShowConfirm(false);
    } catch (err: unknown) {
      const msg = extractErrorMessage(err);
      setPurgeError(msg);
      setShowConfirm(false);
    } finally {
      setPurging(false);
    }
  };

  return (
    <div className="max-w-4xl mx-auto px-4 sm:px-6 py-8">
      {/* Header */}
      <motion.div
        className="text-center mb-12"
        initial={{ opacity: 0, y: -20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.5 }}
      >
        <div className="inline-flex items-center justify-center w-16 h-16 rounded-2xl bg-gradient-to-br from-[var(--color-primary)]/20 to-[var(--color-accent)]/20 mb-4">
          <Shield className="w-8 h-8 text-[var(--color-primary)]" />
        </div>
        <h1 className="text-3xl sm:text-4xl font-bold text-[var(--color-text-primary)] tracking-tight font-display mb-3">
          Privacy Policy
        </h1>
        <p className="text-[var(--color-text-secondary)] max-w-2xl mx-auto">
          ChemAudit is designed with privacy in mind. We believe in transparency about how your data is handled.
        </p>
      </motion.div>

      {/* Content */}
      <motion.div
        className="space-y-8"
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.5, delay: 0.1 }}
      >
        {/* TL;DR Section */}
        <div className="card-glass p-6 border-l-4 border-l-[var(--color-primary)]">
          <h2 className="text-lg font-semibold text-[var(--color-text-primary)] mb-3 flex items-center gap-2">
            <Eye className="w-5 h-5 text-[var(--color-primary)]" />
            TL;DR - The Short Version
          </h2>
          <ul className="space-y-2 text-[var(--color-text-secondary)]">
            <li className="flex items-start gap-2">
              <span className="text-green-500 mt-1">&#10003;</span>
              <span><strong>No tracking</strong> - No analytics, no advertising cookies, no third-party monitoring</span>
            </li>
            <li className="flex items-start gap-2">
              <span className="text-green-500 mt-1">&#10003;</span>
              <span><strong>No accounts</strong> - No registration or personal data collection required</span>
            </li>
            <li className="flex items-start gap-2">
              <span className="text-green-500 mt-1">&#10003;</span>
              <span><strong>Session-scoped storage</strong> - Bookmarks and validation history are tied to an anonymous session cookie and automatically purged after 30 days</span>
            </li>
            <li className="flex items-start gap-2">
              <span className="text-green-500 mt-1">&#10003;</span>
              <span><strong>You control your data</strong> - Delete all your data at any time using the button below</span>
            </li>
          </ul>
        </div>

        {/* What We Store */}
        <div className="card p-6">
          <h2 className="text-lg font-semibold text-[var(--color-text-primary)] mb-4 flex items-center gap-2">
            <Database className="w-5 h-5 text-[var(--color-primary)]" />
            What We Store
          </h2>
          <div className="space-y-4">
            {/* Session Cookie */}
            <div className="bg-[var(--color-surface-sunken)] rounded-xl p-4">
              <div className="flex items-center gap-3 mb-3">
                <div className="w-10 h-10 rounded-lg bg-[var(--color-primary)]/10 flex items-center justify-center">
                  <Cookie className="w-5 h-5 text-[var(--color-primary)]" />
                </div>
                <div>
                  <h3 className="font-medium text-[var(--color-text-primary)]">Session Cookie</h3>
                  <p className="text-sm text-[var(--color-text-muted)]">Cookie: <code className="bg-[var(--color-surface-elevated)] px-1.5 py-0.5 rounded text-xs">chemaudit_sid</code></p>
                </div>
              </div>
              <p className="text-sm text-[var(--color-text-secondary)]">
                A random identifier (UUID) stored as an HttpOnly cookie. This is <strong>not</strong> used for tracking
                &mdash; it links your bookmarks and history to your browser session so only you can see them.
                Expires after 30 days of inactivity.
              </p>
            </div>

            {/* Theme Preference */}
            <div className="bg-[var(--color-surface-sunken)] rounded-xl p-4">
              <div className="flex items-center gap-3 mb-3">
                <div className="w-10 h-10 rounded-lg bg-[var(--color-primary)]/10 flex items-center justify-center">
                  <Eye className="w-5 h-5 text-[var(--color-primary)]" />
                </div>
                <div>
                  <h3 className="font-medium text-[var(--color-text-primary)]">Theme Preference</h3>
                  <p className="text-sm text-[var(--color-text-muted)]">localStorage key: <code className="bg-[var(--color-surface-elevated)] px-1.5 py-0.5 rounded text-xs">chemaudit-theme</code></p>
                </div>
              </div>
              <p className="text-sm text-[var(--color-text-secondary)]">
                Stores your display preference (light, dark, or system). Purely functional, contains no personal information.
              </p>
            </div>

            {/* Bookmarks */}
            <div className="bg-[var(--color-surface-sunken)] rounded-xl p-4">
              <div className="flex items-center gap-3 mb-3">
                <div className="w-10 h-10 rounded-lg bg-amber-500/10 flex items-center justify-center">
                  <Database className="w-5 h-5 text-amber-500" />
                </div>
                <div>
                  <h3 className="font-medium text-[var(--color-text-primary)]">Bookmarks</h3>
                  <p className="text-sm text-[var(--color-text-muted)]">Server-side, scoped to your session</p>
                </div>
              </div>
              <p className="text-sm text-[var(--color-text-secondary)]">
                When you bookmark a molecule, we store the SMILES string, molecule name, tags, and notes on the server.
                Result snapshots are stored locally in your browser (IndexedDB).
                All bookmark data is scoped to your session and invisible to other users.
              </p>
            </div>

            {/* Validation History */}
            <div className="bg-[var(--color-surface-sunken)] rounded-xl p-4">
              <div className="flex items-center gap-3 mb-3">
                <div className="w-10 h-10 rounded-lg bg-blue-500/10 flex items-center justify-center">
                  <Clock className="w-5 h-5 text-blue-500" />
                </div>
                <div>
                  <h3 className="font-medium text-[var(--color-text-primary)]">Validation History</h3>
                  <p className="text-sm text-[var(--color-text-muted)]">Server-side, auto-purged after 30 days</p>
                </div>
              </div>
              <p className="text-sm text-[var(--color-text-secondary)]">
                Each validation you run is logged with the SMILES string, validation score, and outcome (pass/fail).
                This powers the History page and is scoped to your session.
                Entries are automatically deleted after 30 days.
              </p>
            </div>
          </div>
        </div>

        {/* What We Don't Do */}
        <div className="card p-6">
          <h2 className="text-lg font-semibold text-[var(--color-text-primary)] mb-4 flex items-center gap-2">
            <Lock className="w-5 h-5 text-[var(--color-primary)]" />
            What We Don't Do
          </h2>
          <div className="grid gap-3">
            {[
              'No cookies for tracking or advertising',
              'No Google Analytics or similar services',
              'No third-party scripts that collect data',
              'No user accounts or personal information collection',
              'No sharing of any data with third parties',
            ].map((text, i) => (
              <div key={i} className="flex items-center gap-3 p-3 rounded-lg bg-[var(--color-surface-sunken)]">
                <span className="text-lg">&#128683;</span>
                <span className="text-[var(--color-text-secondary)]">{text}</span>
              </div>
            ))}
          </div>
        </div>

        {/* Data Processing */}
        <div className="card p-6">
          <h2 className="text-lg font-semibold text-[var(--color-text-primary)] mb-4 flex items-center gap-2">
            <Server className="w-5 h-5 text-[var(--color-primary)]" />
            Chemical Structure Processing
          </h2>
          <p className="text-[var(--color-text-secondary)] mb-4">
            When you submit a chemical structure for validation:
          </p>
          <ol className="space-y-3 text-[var(--color-text-secondary)]">
            <li className="flex items-start gap-3">
              <span className="w-6 h-6 rounded-full bg-[var(--color-primary)]/10 text-[var(--color-primary)] text-sm font-medium flex items-center justify-center flex-shrink-0">1</span>
              <span>Your structure is sent to our server for processing</span>
            </li>
            <li className="flex items-start gap-3">
              <span className="w-6 h-6 rounded-full bg-[var(--color-primary)]/10 text-[var(--color-primary)] text-sm font-medium flex items-center justify-center flex-shrink-0">2</span>
              <span>Validation, scoring, and standardization are performed in memory</span>
            </li>
            <li className="flex items-start gap-3">
              <span className="w-6 h-6 rounded-full bg-[var(--color-primary)]/10 text-[var(--color-primary)] text-sm font-medium flex items-center justify-center flex-shrink-0">3</span>
              <span>Results are returned to your browser immediately</span>
            </li>
            <li className="flex items-start gap-3">
              <span className="w-6 h-6 rounded-full bg-[var(--color-primary)]/10 text-[var(--color-primary)] text-sm font-medium flex items-center justify-center flex-shrink-0">4</span>
              <span>A summary (SMILES, score, outcome) is logged in your validation history, scoped to your session</span>
            </li>
          </ol>
        </div>

        {/* GDPR Compliance */}
        <div className="card p-6">
          <h2 className="text-lg font-semibold text-[var(--color-text-primary)] mb-4">
            GDPR Compliance
          </h2>
          <p className="text-[var(--color-text-secondary)] mb-3">
            ChemAudit does not collect personal data such as names, emails, or IP addresses.
            The session cookie is a &ldquo;strictly necessary&rdquo; functional cookie under GDPR Article 6(1)(f)
            and does not require consent.
          </p>
          <p className="text-[var(--color-text-secondary)] mb-3">
            Bookmarks and validation history contain only chemical structure data (SMILES),
            scoped to an anonymous session identifier. This data is automatically purged after 30 days.
          </p>
          <p className="text-[var(--color-text-secondary)]">
            <strong>Your rights:</strong> You can view your stored data on the{' '}
            <a href="/bookmarks" className="text-[var(--color-primary)] hover:underline">Bookmarks</a> and{' '}
            <a href="/history" className="text-[var(--color-primary)] hover:underline">History</a> pages,
            or delete everything immediately using the button below.
          </p>
        </div>

        {/* Delete My Data */}
        <div className="card p-6 border border-red-500/20">
          <h2 className="text-lg font-semibold text-[var(--color-text-primary)] mb-3 flex items-center gap-2">
            <Trash2 className="w-5 h-5 text-red-500" />
            Delete All My Data
          </h2>
          <p className="text-[var(--color-text-secondary)] mb-4">
            This permanently deletes all bookmarks, validation history, and locally cached result snapshots
            associated with your current session. This action cannot be undone.
          </p>

          {!showConfirm && !purgeResult && (
            <button
              onClick={() => setShowConfirm(true)}
              className="px-4 py-2 rounded-lg bg-red-500/10 text-red-600 hover:bg-red-500/20 font-medium transition-colors"
            >
              Delete All My Data
            </button>
          )}

          {showConfirm && (
            <div className="flex items-center gap-3 p-4 rounded-lg bg-red-500/10">
              <AlertTriangle className="w-5 h-5 text-red-500 flex-shrink-0" />
              <span className="text-[var(--color-text-secondary)]">Are you sure? This cannot be undone.</span>
              <button
                onClick={handlePurge}
                disabled={purging}
                className="px-4 py-2 rounded-lg bg-red-600 text-white hover:bg-red-700 font-medium transition-colors disabled:opacity-50"
              >
                {purging ? 'Deleting...' : 'Yes, delete everything'}
              </button>
              <button
                onClick={() => setShowConfirm(false)}
                disabled={purging}
                className="px-4 py-2 rounded-lg bg-[var(--color-surface-sunken)] text-[var(--color-text-secondary)] hover:bg-[var(--color-surface-elevated)] font-medium transition-colors disabled:opacity-50"
              >
                Cancel
              </button>
            </div>
          )}

          {purgeResult && (
            <div className="p-4 rounded-lg bg-green-500/10 text-green-700">
              Data deleted successfully: {purgeResult.bookmarks} bookmark{purgeResult.bookmarks !== 1 ? 's' : ''} and {purgeResult.history} history entr{purgeResult.history !== 1 ? 'ies' : 'y'} removed.
              Local snapshots cleared.
            </div>
          )}

          {purgeError && (
            <div className="p-4 rounded-lg bg-red-500/10 text-red-600">
              {purgeError}
            </div>
          )}
        </div>

        {/* Contact */}
        <div className="card p-6 bg-gradient-to-br from-[var(--color-primary)]/5 to-[var(--color-accent)]/5">
          <h2 className="text-lg font-semibold text-[var(--color-text-primary)] mb-3">
            Questions?
          </h2>
          <p className="text-[var(--color-text-secondary)]">
            If you have any questions about this privacy policy, please open an issue on our{' '}
            <a
              href="https://github.com/Kohulan/ChemAudit"
              target="_blank"
              rel="noopener noreferrer"
              className="text-[var(--color-primary)] hover:underline font-medium"
            >
              GitHub repository
            </a>.
          </p>
        </div>

        {/* Last Updated */}
        <p className="text-center text-sm text-[var(--color-text-muted)]">
          Last updated: February 25, 2026
        </p>
      </motion.div>
    </div>
  );
}
