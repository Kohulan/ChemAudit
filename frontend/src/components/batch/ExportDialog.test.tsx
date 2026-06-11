import { render, screen } from '@testing-library/react';
import { describe, it, expect, vi } from 'vitest';
import { ExportDialog } from './ExportDialog';

vi.mock('../../services/api', () => ({
  api: { get: vi.fn(), post: vi.fn() },
}));

describe('ExportDialog accessibility', () => {
  it('exposes dialog role named by its heading', () => {
    render(<ExportDialog jobId="abcd1234" isOpen onClose={() => {}} />);
    expect(
      screen.getByRole('dialog', { name: /export results/i })
    ).toBeInTheDocument();
  });

  it('has an accessible name on the close button', () => {
    render(<ExportDialog jobId="abcd1234" isOpen onClose={() => {}} />);
    expect(
      screen.getByRole('button', { name: /close export dialog/i })
    ).toBeInTheDocument();
  });
});
