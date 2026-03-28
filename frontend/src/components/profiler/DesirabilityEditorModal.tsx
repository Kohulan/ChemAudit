import { useState, useEffect, useRef, useCallback, useId } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { X, Plus, Trash2 } from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';
import { useMPOProfiles, type MPOProfile } from '../../hooks/useMPOProfiles';
import type { MPOProperty } from '../../types/profiler';
import { cn } from '../../lib/utils';

// =============================================================================
// Available descriptors for MPO property dropdown (per D-15)
// =============================================================================
const AVAILABLE_DESCRIPTORS = [
  'MW', 'LogP', 'TPSA', 'HBD', 'HBA', 'RotBonds', 'Fsp3', 'QED', 'SA', 'PFI',
];

// =============================================================================
// Desirability function (from 07-RESEARCH.md)
// =============================================================================
function desirabilityValue(
  value: number,
  low: number,
  high: number,
  shape: 'sigmoid' | 'ramp' | 'step'
): number {
  if (high === low) return value >= high ? 1 : 0;
  const x = (value - low) / (high - low);
  const clamped = Math.max(0, Math.min(1, x));
  if (shape === 'sigmoid') return 1 / (1 + Math.exp(-12 * (clamped - 0.5)));
  if (shape === 'ramp') return clamped;
  return value >= high ? 1 : 0; // step
}

// =============================================================================
// SVG curve preview
// =============================================================================
interface CurvePreviewProps {
  low: number;
  high: number;
  shape: 'sigmoid' | 'ramp' | 'step';
}

const SVG_W = 200;
const SVG_H = 100;
const SAMPLES = 50;
const PAD = 10;

/**
 * Small SVG curve preview showing the desirability function shape.
 * Samples 50 points across the [low, high] range and renders as a polyline.
 * Updates on every prop change (no debounce — pure math per UI-SPEC).
 */
function CurvePreview({ low, high, shape }: CurvePreviewProps) {
  const rangeSpan = high - low || 1;
  const extendedLow = low - rangeSpan * 0.2;
  const extendedHigh = high + rangeSpan * 0.2;

  const points = Array.from({ length: SAMPLES }, (_, i) => {
    const t = i / (SAMPLES - 1);
    const xVal = extendedLow + t * (extendedHigh - extendedLow);
    const d = desirabilityValue(xVal, low, high, shape);
    const svgX = PAD + t * (SVG_W - 2 * PAD);
    const svgY = PAD + (1 - d) * (SVG_H - 2 * PAD);
    return `${svgX.toFixed(1)},${svgY.toFixed(1)}`;
  }).join(' ');

  return (
    <svg
      viewBox={`0 0 ${SVG_W} ${SVG_H}`}
      className="w-full h-[80px] border border-[var(--color-border)] rounded-lg bg-surface-sunken"
      aria-label={`Desirability curve preview: ${shape} shape from ${low} to ${high}`}
    >
      {/* Axes */}
      <line x1={PAD} y1={SVG_H - PAD} x2={SVG_W - PAD} y2={SVG_H - PAD} stroke="var(--color-border)" strokeWidth="1" />
      <line x1={PAD} y1={PAD} x2={PAD} y2={SVG_H - PAD} stroke="var(--color-border)" strokeWidth="1" />
      {/* Curve */}
      <polyline
        points={points}
        fill="none"
        stroke="var(--color-primary, #c41e3a)"
        strokeWidth="2"
        strokeLinecap="round"
        strokeLinejoin="round"
      />
      {/* Y=0 label */}
      <text x={PAD - 2} y={SVG_H - PAD} textAnchor="end" fontSize="8" fill="var(--color-text-muted)" dominantBaseline="middle">0</text>
      {/* Y=1 label */}
      <text x={PAD - 2} y={PAD} textAnchor="end" fontSize="8" fill="var(--color-text-muted)" dominantBaseline="middle">1</text>
    </svg>
  );
}

// =============================================================================
// Property row
// =============================================================================
interface PropertyRowProps {
  prop: MPOProperty;
  index: number;
  onChange: (index: number, updated: MPOProperty) => void;
  onRemove: (index: number) => void;
  activePropIndex: number;
  onActivate: (index: number) => void;
}

function PropertyRow({ prop, index, onChange, onRemove, activePropIndex, onActivate }: PropertyRowProps) {
  const isActive = activePropIndex === index;

  const update = (partial: Partial<MPOProperty>) =>
    onChange(index, { ...prop, ...partial });

  return (
    <div
      className={cn(
        'border rounded-xl p-3 space-y-2 cursor-pointer transition-colors',
        isActive
          ? 'border-[var(--color-primary)]/40 bg-surface-elevated'
          : 'border-[var(--color-border)] hover:border-[var(--color-border)]/80'
      )}
      onClick={() => onActivate(index)}
    >
      <div className="flex items-center gap-2 flex-wrap">
        {/* Property selector */}
        <select
          value={prop.property}
          onChange={(e) => update({ property: e.target.value })}
          onClick={(e) => e.stopPropagation()}
          className="text-sm bg-surface-elevated border border-[var(--color-border)] rounded-lg px-2 py-1 text-text-primary"
        >
          {AVAILABLE_DESCRIPTORS.map((d) => (
            <option key={d} value={d}>{d}</option>
          ))}
        </select>

        {/* Low */}
        <label className="flex items-center gap-1 text-xs text-text-muted">
          Low
          <input
            type="number"
            value={prop.low}
            onChange={(e) => update({ low: parseFloat(e.target.value) || 0 })}
            onClick={(e) => e.stopPropagation()}
            className="w-16 text-sm bg-surface-elevated border border-[var(--color-border)] rounded-lg px-2 py-1 text-text-primary"
          />
        </label>

        {/* High */}
        <label className="flex items-center gap-1 text-xs text-text-muted">
          High
          <input
            type="number"
            value={prop.high}
            onChange={(e) => update({ high: parseFloat(e.target.value) || 0 })}
            onClick={(e) => e.stopPropagation()}
            className="w-16 text-sm bg-surface-elevated border border-[var(--color-border)] rounded-lg px-2 py-1 text-text-primary"
          />
        </label>

        {/* Weight */}
        <label className="flex items-center gap-1 text-xs text-text-muted">
          Weight
          <input
            type="number"
            value={prop.weight}
            min={0}
            step={0.1}
            onChange={(e) => update({ weight: parseFloat(e.target.value) || 1 })}
            onClick={(e) => e.stopPropagation()}
            className="w-14 text-sm bg-surface-elevated border border-[var(--color-border)] rounded-lg px-2 py-1 text-text-primary"
          />
        </label>

        {/* Shape */}
        <select
          value={prop.shape}
          onChange={(e) => update({ shape: e.target.value as MPOProperty['shape'] })}
          onClick={(e) => e.stopPropagation()}
          className="text-sm bg-surface-elevated border border-[var(--color-border)] rounded-lg px-2 py-1 text-text-primary"
        >
          <option value="sigmoid">Sigmoid</option>
          <option value="ramp">Ramp</option>
          <option value="step">Step</option>
        </select>

        {/* Remove */}
        <button
          type="button"
          onClick={(e) => { e.stopPropagation(); onRemove(index); }}
          className="ml-auto text-text-muted hover:text-status-error transition-colors"
          aria-label={`Remove ${prop.property} property`}
        >
          <Trash2 className="w-4 h-4" />
        </button>
      </div>

      {/* Live curve preview — only shown for active (clicked) row */}
      {isActive && (
        <CurvePreview low={prop.low} high={prop.high} shape={prop.shape} />
      )}
    </div>
  );
}

// =============================================================================
// Delete confirmation modal
// =============================================================================
interface DeleteConfirmProps {
  profileName: string;
  onConfirm: () => void;
  onCancel: () => void;
}

function DeleteConfirmModal({ profileName, onConfirm, onCancel }: DeleteConfirmProps) {
  return (
    <div className="fixed inset-0 z-60 flex items-center justify-center">
      <div className="absolute inset-0 bg-surface-overlay backdrop-blur-sm" onClick={onCancel} />
      <motion.div
        className="relative bg-surface-elevated rounded-2xl p-6 max-w-sm w-full mx-4 shadow-2xl"
        initial={{ scale: 0.94, opacity: 0 }}
        animate={{ scale: 1, opacity: 1 }}
        transition={{ duration: 0.2, ease: 'easeOut' }}
        role="dialog"
        aria-modal="true"
      >
        <p className="text-sm text-text-primary mb-4">
          Delete &ldquo;{profileName}&rdquo;? This cannot be undone.
        </p>
        <div className="flex gap-2 justify-end">
          <ClayButton variant="ghost" size="sm" onClick={onCancel}>Cancel</ClayButton>
          <ClayButton variant="danger" size="sm" onClick={onConfirm}>Delete</ClayButton>
        </div>
      </motion.div>
    </div>
  );
}

// =============================================================================
// Main modal
// =============================================================================
interface DesirabilityEditorModalProps {
  isOpen: boolean;
  onClose: () => void;
  profile: MPOProfile;
  onSave: (profile: MPOProfile) => void;
}

/**
 * Full-screen modal for editing a custom MPO profile.
 *
 * Per D-14: live SVG curve preview (no debounce — pure math).
 * Accessibility: role="dialog", aria-modal, focus trap, aria-labelledby.
 * Saves to localStorage via useMPOProfiles().saveProfile().
 */
export function DesirabilityEditorModal({
  isOpen,
  onClose,
  profile,
  onSave,
}: DesirabilityEditorModalProps) {
  const titleId = useId();
  const { saveProfile, deleteProfile, isPreset } = useMPOProfiles();
  const [editedProfile, setEditedProfile] = useState<MPOProfile>(profile);
  const [activePropertyIndex, setActivePropertyIndex] = useState<number>(0);
  const [showDeleteConfirm, setShowDeleteConfirm] = useState(false);
  const modalRef = useRef<HTMLDivElement>(null);
  const closeBtnRef = useRef<HTMLButtonElement>(null);
  const previousFocusRef = useRef<HTMLElement | null>(null);

  // Sync when incoming profile changes
  useEffect(() => {
    setEditedProfile(profile);
    setActivePropertyIndex(0);
  }, [profile]);

  // Focus trap + return focus
  useEffect(() => {
    if (!isOpen) return;

    previousFocusRef.current = document.activeElement as HTMLElement;
    // Focus the close button when modal opens
    setTimeout(() => closeBtnRef.current?.focus(), 50);

    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === 'Escape') {
        onClose();
        return;
      }
      if (e.key !== 'Tab') return;

      const modal = modalRef.current;
      if (!modal) return;
      const focusable = modal.querySelectorAll<HTMLElement>(
        'button, [href], input, select, textarea, [tabindex]:not([tabindex="-1"])'
      );
      const first = focusable[0];
      const last = focusable[focusable.length - 1];
      if (e.shiftKey) {
        if (document.activeElement === first) {
          e.preventDefault();
          last?.focus();
        }
      } else {
        if (document.activeElement === last) {
          e.preventDefault();
          first?.focus();
        }
      }
    };

    document.addEventListener('keydown', handleKeyDown);
    return () => {
      document.removeEventListener('keydown', handleKeyDown);
      previousFocusRef.current?.focus();
    };
  }, [isOpen, onClose]);

  const handlePropertyChange = useCallback((index: number, updated: MPOProperty) => {
    setEditedProfile((prev) => {
      const props = [...prev.properties];
      props[index] = updated;
      return { ...prev, properties: props };
    });
  }, []);

  const handleRemoveProperty = useCallback((index: number) => {
    setEditedProfile((prev) => {
      const props = prev.properties.filter((_, i) => i !== index);
      return { ...prev, properties: props };
    });
    setActivePropertyIndex((prev) => Math.max(0, prev - 1));
  }, []);

  const handleAddProperty = useCallback(() => {
    const newProp: MPOProperty = {
      property: 'MW',
      low: 200,
      high: 500,
      weight: 1,
      shape: 'sigmoid',
    };
    setEditedProfile((prev) => ({
      ...prev,
      properties: [...prev.properties, newProp],
    }));
    setActivePropertyIndex(editedProfile.properties.length);
  }, [editedProfile.properties.length]);

  const handleSave = useCallback(() => {
    const toSave: MPOProfile = {
      ...editedProfile,
      id: editedProfile.id || String(Date.now()),
      createdAt: editedProfile.createdAt || Date.now(),
    };
    saveProfile(toSave);
    onSave(toSave);
    onClose();
  }, [editedProfile, saveProfile, onSave, onClose]);

  const handleDeleteConfirm = useCallback(() => {
    deleteProfile(editedProfile.id);
    setShowDeleteConfirm(false);
    onClose();
  }, [editedProfile.id, deleteProfile, onClose]);

  const profileIsPreset = isPreset(editedProfile.id);

  return (
    <AnimatePresence>
      {isOpen && (
        <div className="fixed inset-0 z-50 flex items-center justify-center">
          {/* Backdrop */}
          <motion.div
            className="absolute inset-0 bg-surface-overlay backdrop-blur-sm"
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
            onClick={onClose}
          />

          {/* Modal panel */}
          <motion.div
            ref={modalRef}
            role="dialog"
            aria-modal="true"
            aria-labelledby={titleId}
            className="relative max-w-2xl w-full mx-4 bg-surface-elevated rounded-2xl shadow-2xl overflow-hidden"
            initial={{ scale: 0.94, opacity: 0 }}
            animate={{ scale: 1, opacity: 1 }}
            exit={{ scale: 0.94, opacity: 0 }}
            transition={{ duration: 0.4, ease: 'easeOut' }}
          >
            {/* Header */}
            <div className="flex items-center justify-between px-6 py-4 border-b border-[var(--color-border)]">
              <h2 id={titleId} className="text-lg font-semibold text-text-primary font-display">
                {profileIsPreset ? `View Preset: ${editedProfile.name}` : 'Edit MPO Profile'}
              </h2>
              <button
                ref={closeBtnRef}
                type="button"
                onClick={onClose}
                aria-label="Close MPO editor"
                className="text-text-muted hover:text-text-primary transition-colors p-1 rounded-lg"
              >
                <X className="w-5 h-5" />
              </button>
            </div>

            {/* Content — scrollable */}
            <div className="px-6 py-4 max-h-[70vh] overflow-y-auto space-y-4">
              {/* Profile name */}
              {!profileIsPreset && (
                <div>
                  <label className="block text-sm font-medium text-text-secondary mb-1">
                    Profile name
                  </label>
                  <input
                    type="text"
                    value={editedProfile.name}
                    onChange={(e) => setEditedProfile((prev) => ({ ...prev, name: e.target.value }))}
                    className="w-full text-sm bg-surface-elevated border border-[var(--color-border)] rounded-xl px-3 py-2 text-text-primary"
                    placeholder="My custom MPO profile"
                  />
                </div>
              )}

              {/* Properties */}
              <div>
                <div className="flex items-center justify-between mb-2">
                  <p className="text-sm font-medium text-text-secondary">
                    Properties ({editedProfile.properties.length})
                  </p>
                  {!profileIsPreset && (
                    <ClayButton
                      variant="ghost"
                      size="sm"
                      onClick={handleAddProperty}
                      leftIcon={<Plus className="w-3.5 h-3.5" />}
                    >
                      Add Property
                    </ClayButton>
                  )}
                </div>

                {editedProfile.properties.length === 0 ? (
                  <p className="text-sm text-text-muted text-center py-4 border border-dashed border-[var(--color-border)] rounded-xl">
                    No properties configured. Add at least one property to compute the MPO score.
                  </p>
                ) : (
                  <div className="space-y-2">
                    {editedProfile.properties.map((prop, i) => (
                      profileIsPreset ? (
                        /* Read-only view for presets */
                        <div key={i} className="border border-[var(--color-border)] rounded-xl p-3">
                          <div className="flex items-center gap-2 flex-wrap text-sm text-text-secondary">
                            <span className="font-medium text-text-primary">{prop.property}</span>
                            <span>Low: {prop.low}</span>
                            <span>High: {prop.high}</span>
                            <span>Weight: {prop.weight}</span>
                            <span className="capitalize">{prop.shape}</span>
                          </div>
                          <CurvePreview low={prop.low} high={prop.high} shape={prop.shape} />
                        </div>
                      ) : (
                        <PropertyRow
                          key={i}
                          prop={prop}
                          index={i}
                          onChange={handlePropertyChange}
                          onRemove={handleRemoveProperty}
                          activePropIndex={activePropertyIndex}
                          onActivate={setActivePropertyIndex}
                        />
                      )
                    ))}
                  </div>
                )}
              </div>
            </div>

            {/* Footer */}
            <div className={cn(
              'flex items-center px-6 py-4 border-t border-[var(--color-border)]',
              profileIsPreset ? 'justify-end' : 'justify-between'
            )}>
              {!profileIsPreset && (
                <ClayButton
                  variant="danger"
                  size="sm"
                  onClick={() => setShowDeleteConfirm(true)}
                  leftIcon={<Trash2 className="w-3.5 h-3.5" />}
                >
                  Delete Profile
                </ClayButton>
              )}
              <div className="flex gap-2">
                <ClayButton variant="ghost" size="sm" onClick={onClose}>
                  Cancel
                </ClayButton>
                {!profileIsPreset && (
                  <ClayButton
                    variant="accent"
                    size="sm"
                    onClick={handleSave}
                    disabled={editedProfile.properties.length === 0 || !editedProfile.name.trim()}
                  >
                    Save Profile
                  </ClayButton>
                )}
              </div>
            </div>
          </motion.div>

          {/* Delete confirmation sub-modal */}
          {showDeleteConfirm && (
            <DeleteConfirmModal
              profileName={editedProfile.name}
              onConfirm={handleDeleteConfirm}
              onCancel={() => setShowDeleteConfirm(false)}
            />
          )}
        </div>
      )}
    </AnimatePresence>
  );
}
