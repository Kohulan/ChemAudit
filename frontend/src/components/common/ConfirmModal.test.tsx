import { describe, it, expect, vi, beforeEach } from 'vitest';
import userEvent from '@testing-library/user-event';

import { render, screen } from '../../tests/setup';
import { ConfirmModal } from './ConfirmModal';

describe('ConfirmModal', () => {
  const onConfirm = vi.fn();
  const onCancel = vi.fn();

  beforeEach(() => vi.clearAllMocks());

  it('renders nothing when closed', () => {
    render(
      <ConfirmModal isOpen={false} message="Delete it?" onConfirm={onConfirm} onCancel={onCancel} />,
    );
    expect(screen.queryByRole('dialog')).not.toBeInTheDocument();
  });

  it('shows the message and is an accessible dialog when open', () => {
    render(
      <ConfirmModal isOpen message="Delete 'foo'?" onConfirm={onConfirm} onCancel={onCancel} />,
    );
    const dialog = screen.getByRole('dialog');
    expect(dialog).toHaveAttribute('aria-modal', 'true');
    expect(screen.getByText("Delete 'foo'?")).toBeInTheDocument();
  });

  it('calls onConfirm when the confirm button is clicked', async () => {
    const user = userEvent.setup();
    render(
      <ConfirmModal isOpen message="Delete?" confirmLabel="Delete" onConfirm={onConfirm} onCancel={onCancel} />,
    );
    await user.click(screen.getByRole('button', { name: 'Delete' }));
    expect(onConfirm).toHaveBeenCalledTimes(1);
    expect(onCancel).not.toHaveBeenCalled();
  });

  it('calls onCancel on cancel click and on Escape', async () => {
    const user = userEvent.setup();
    render(<ConfirmModal isOpen message="Delete?" onConfirm={onConfirm} onCancel={onCancel} />);
    await user.click(screen.getByRole('button', { name: 'Cancel' }));
    await user.keyboard('{Escape}');
    expect(onCancel).toHaveBeenCalledTimes(2);
    expect(onConfirm).not.toHaveBeenCalled();
  });
});
