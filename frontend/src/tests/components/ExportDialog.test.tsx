/**
 * ExportDialog Component Tests
 *
 * Tests enhanced export dialog with 9 formats and PDF section selection.
 */
import React, { forwardRef, createElement } from 'react';
import { describe, it, expect, vi } from 'vitest';
import { render, screen, fireEvent } from '../setup';
import { ExportDialog } from '../../components/batch/ExportDialog';

// Mock framer-motion: create forwarded-ref component factory for any HTML tag
function makeMockMotionComponent(tag: string) {
  return forwardRef(function MockMotion(props: Record<string, unknown>, ref: React.Ref<unknown>) {
    // Strip framer-motion props, pass only DOM-safe props
    const motionProps = [
      'initial', 'animate', 'exit', 'transition', 'whileHover', 'whileTap',
      'whileInView', 'layoutId', 'layout', 'onAnimationComplete', 'variants',
      'drag', 'dragConstraints',
    ];
    const domProps: Record<string, unknown> = {};
    for (const [key, val] of Object.entries(props)) {
      if (!motionProps.includes(key)) domProps[key] = val;
    }
    return createElement(tag, { ...domProps, ref });
  });
}

vi.mock('framer-motion', () => ({
  motion: new Proxy({}, {
    get(_target: unknown, tag: string) {
      return makeMockMotionComponent(tag);
    },
  }),
  AnimatePresence: ({ children }: { children: React.ReactNode }) => <>{children}</>,
}));

describe('ExportDialog', () => {
  const defaultProps = {
    jobId: 'test-job-123',
    isOpen: true,
    onClose: vi.fn(),
  };

  describe('Format Rendering', () => {
    it('renders all 9 export format options', () => {
      render(<ExportDialog {...defaultProps} />);

      expect(screen.getByText('CSV')).toBeInTheDocument();
      expect(screen.getByText('Excel')).toBeInTheDocument();
      expect(screen.getByText('SDF')).toBeInTheDocument();
      expect(screen.getByText('JSON')).toBeInTheDocument();
      expect(screen.getByText('PDF Report')).toBeInTheDocument();
      expect(screen.getByText('Fingerprints')).toBeInTheDocument();
      expect(screen.getByText('Deduplicated')).toBeInTheDocument();
      expect(screen.getByText('Scaffold-Grouped')).toBeInTheDocument();
      expect(screen.getByText('Property Matrix')).toBeInTheDocument();
    });

    it('shows "New" badge for advanced formats', () => {
      render(<ExportDialog {...defaultProps} />);

      const newBadges = screen.getAllByText('New');
      expect(newBadges.length).toBeGreaterThanOrEqual(4);
    });
  });

  describe('PDF Section Checkboxes', () => {
    it('shows section checkboxes when PDF format is selected', () => {
      render(<ExportDialog {...defaultProps} />);

      // Select PDF format
      const pdfLabel = screen.getByText('PDF Report');
      fireEvent.click(pdfLabel.closest('label')!);

      // Check for section checkboxes
      expect(screen.getByText('Validation Summary')).toBeInTheDocument();
      expect(screen.getByText('Score Distribution')).toBeInTheDocument();
      expect(screen.getByText('Alert Frequency')).toBeInTheDocument();
      expect(screen.getByText('Chemical Space')).toBeInTheDocument();
      expect(screen.getByText('Scaffold Treemap')).toBeInTheDocument();
      expect(screen.getByText('Statistics')).toBeInTheDocument();
      expect(screen.getByText('Correlation Matrix')).toBeInTheDocument();
      expect(screen.getByText('MMP Pairs')).toBeInTheDocument();
    });
  });

  describe('Fingerprint Format', () => {
    it('shows zip extension for formats that produce ZIP archives', () => {
      render(<ExportDialog {...defaultProps} />);

      // Multiple formats have .zip extension (fingerprint, dedup, property_matrix)
      const zipLabels = screen.getAllByText('.zip');
      expect(zipLabels.length).toBeGreaterThanOrEqual(2);
    });
  });

  describe('File Preview', () => {
    it('shows filename preview with job ID', () => {
      render(<ExportDialog {...defaultProps} />);

      // Default format is CSV
      const preview = screen.getByText(/chemaudit_csv_.*_test-job/);
      expect(preview).toBeInTheDocument();
    });
  });

  describe('Dialog Controls', () => {
    it('does not render when isOpen is false', () => {
      render(<ExportDialog {...defaultProps} isOpen={false} />);

      expect(screen.queryByText('Export Results')).not.toBeInTheDocument();
    });

    it('shows selected count when indices provided', () => {
      render(<ExportDialog {...defaultProps} selectedIndices={new Set([0, 1, 2])} />);

      expect(screen.getByText(/3 selected molecules/)).toBeInTheDocument();
    });

    it('calls onClose when Cancel is clicked', () => {
      const onClose = vi.fn();
      render(<ExportDialog {...defaultProps} onClose={onClose} />);

      fireEvent.click(screen.getByText('Cancel'));
      expect(onClose).toHaveBeenCalled();
    });
  });
});
