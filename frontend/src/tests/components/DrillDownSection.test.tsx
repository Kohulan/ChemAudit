/**
 * DrillDownSection Component Tests
 *
 * Tests the shared accordion component used by Profiler, Safety,
 * and Diagnostics sections in SingleValidation page.
 */
import { describe, it, expect, vi } from 'vitest';
import { render, screen, fireEvent } from '../setup';
import { DrillDownSection } from '../../components/ui/DrillDownSection';

// Mock framer-motion to render synchronously in tests
vi.mock('framer-motion', () => ({
  motion: {
    div: ({
      children,
      ...props
    }: React.HTMLAttributes<HTMLDivElement> & Record<string, unknown>) => {
      // Filter out framer-motion-specific props
      const {
        initial: _initial,
        animate: _animate,
        exit: _exit,
        transition: _transition,
        variants: _variants,
        ...htmlProps
      } = props;
      return <div {...htmlProps}>{children}</div>;
    },
  },
  AnimatePresence: ({ children }: { children: React.ReactNode }) => <>{children}</>,
}));

describe('DrillDownSection', () => {
  describe('Default Rendering', () => {
    it('renders collapsed by default with title visible', () => {
      render(
        <DrillDownSection title="Molecular Profile">
          <p>Inner content</p>
        </DrillDownSection>
      );

      expect(screen.getByText('Molecular Profile')).toBeInTheDocument();
    });

    it('shows summary badge when collapsed', () => {
      render(
        <DrillDownSection title="Safety" summary={<span>3 alerts</span>}>
          <p>Inner content</p>
        </DrillDownSection>
      );

      expect(screen.getByText('3 alerts')).toBeInTheDocument();
    });

    it('hides children when collapsed', () => {
      render(
        <DrillDownSection title="Safety">
          <p>Hidden content</p>
        </DrillDownSection>
      );

      expect(screen.queryByText('Hidden content')).not.toBeInTheDocument();
    });
  });

  describe('Expand/Collapse', () => {
    it('expands on click and shows children', () => {
      render(
        <DrillDownSection title="Diagnostics">
          <p>Expanded content</p>
        </DrillDownSection>
      );

      const button = screen.getByRole('button', { name: /diagnostics/i });
      fireEvent.click(button);

      expect(screen.getByText('Expanded content')).toBeInTheDocument();
    });

    it('hides summary badge when expanded', () => {
      render(
        <DrillDownSection title="Safety" summary={<span>3 alerts</span>}>
          <p>Content</p>
        </DrillDownSection>
      );

      // Summary visible when collapsed
      expect(screen.getByText('3 alerts')).toBeInTheDocument();

      // Click to expand
      const button = screen.getByRole('button', { name: /safety/i });
      fireEvent.click(button);

      // Summary should be hidden when expanded
      expect(screen.queryByText('3 alerts')).not.toBeInTheDocument();
    });

    it('collapses on second click', () => {
      render(
        <DrillDownSection title="Diagnostics">
          <p>Toggle content</p>
        </DrillDownSection>
      );

      const button = screen.getByRole('button', { name: /diagnostics/i });

      // Expand
      fireEvent.click(button);
      expect(screen.getByText('Toggle content')).toBeInTheDocument();

      // Collapse
      fireEvent.click(button);
      expect(screen.queryByText('Toggle content')).not.toBeInTheDocument();
    });
  });

  describe('defaultOpen Prop', () => {
    it('starts expanded when defaultOpen is true', () => {
      render(
        <DrillDownSection title="Profile" defaultOpen>
          <p>Already visible</p>
        </DrillDownSection>
      );

      expect(screen.getByText('Already visible')).toBeInTheDocument();
    });

    it('starts collapsed when defaultOpen is false', () => {
      render(
        <DrillDownSection title="Profile" defaultOpen={false}>
          <p>Not visible</p>
        </DrillDownSection>
      );

      expect(screen.queryByText('Not visible')).not.toBeInTheDocument();
    });
  });

  describe('onToggle Callback', () => {
    it('calls onToggle with true when expanding', () => {
      const onToggle = vi.fn();

      render(
        <DrillDownSection title="Safety" onToggle={onToggle}>
          <p>Content</p>
        </DrillDownSection>
      );

      const button = screen.getByRole('button', { name: /safety/i });
      fireEvent.click(button);

      expect(onToggle).toHaveBeenCalledTimes(1);
      expect(onToggle).toHaveBeenCalledWith(true);
    });

    it('calls onToggle with false when collapsing', () => {
      const onToggle = vi.fn();

      render(
        <DrillDownSection title="Safety" defaultOpen onToggle={onToggle}>
          <p>Content</p>
        </DrillDownSection>
      );

      const button = screen.getByRole('button', { name: /safety/i });
      fireEvent.click(button);

      expect(onToggle).toHaveBeenCalledTimes(1);
      expect(onToggle).toHaveBeenCalledWith(false);
    });
  });

  describe('Accessibility', () => {
    it('has correct aria-expanded=false when collapsed', () => {
      render(
        <DrillDownSection title="Profile">
          <p>Content</p>
        </DrillDownSection>
      );

      const button = screen.getByRole('button', { name: /profile/i });
      expect(button).toHaveAttribute('aria-expanded', 'false');
    });

    it('has correct aria-expanded=true when expanded', () => {
      render(
        <DrillDownSection title="Profile" defaultOpen>
          <p>Content</p>
        </DrillDownSection>
      );

      const button = screen.getByRole('button', { name: /profile/i });
      expect(button).toHaveAttribute('aria-expanded', 'true');
    });

    it('has aria-controls linking button to panel', () => {
      render(
        <DrillDownSection title="Profile" defaultOpen>
          <p>Content</p>
        </DrillDownSection>
      );

      const button = screen.getByRole('button', { name: /profile/i });
      const controlsId = button.getAttribute('aria-controls');
      expect(controlsId).toBeTruthy();

      // The panel element should have the matching id
      const panel = document.getElementById(controlsId!);
      expect(panel).toBeInTheDocument();
    });

    it('toggles on Enter key', () => {
      render(
        <DrillDownSection title="Profile">
          <p>Keyboard content</p>
        </DrillDownSection>
      );

      const button = screen.getByRole('button', { name: /profile/i });
      fireEvent.keyDown(button, { key: 'Enter' });

      expect(screen.getByText('Keyboard content')).toBeInTheDocument();
    });

    it('toggles on Space key', () => {
      render(
        <DrillDownSection title="Profile">
          <p>Keyboard content</p>
        </DrillDownSection>
      );

      const button = screen.getByRole('button', { name: /profile/i });
      fireEvent.keyDown(button, { key: ' ' });

      expect(screen.getByText('Keyboard content')).toBeInTheDocument();
    });
  });

  describe('Optional Props', () => {
    it('renders icon when provided', () => {
      render(
        <DrillDownSection
          title="Profile"
          icon={<span data-testid="custom-icon">*</span>}
        >
          <p>Content</p>
        </DrillDownSection>
      );

      expect(screen.getByTestId('custom-icon')).toBeInTheDocument();
    });

    it('applies custom className', () => {
      const { container } = render(
        <DrillDownSection title="Profile" className="custom-class">
          <p>Content</p>
        </DrillDownSection>
      );

      expect(container.firstChild).toHaveClass('custom-class');
    });
  });
});
