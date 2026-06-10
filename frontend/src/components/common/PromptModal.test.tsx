import { describe, it, expect, vi, beforeEach } from 'vitest';
import userEvent from '@testing-library/user-event';

import { render, screen } from '../../tests/setup';
import { PromptModal } from './PromptModal';

describe('PromptModal', () => {
  const onConfirm = vi.fn();
  const onCancel = vi.fn();

  beforeEach(() => vi.clearAllMocks());

  it('renders nothing when closed', () => {
    render(
      <PromptModal isOpen={false} title="Name it" onConfirm={onConfirm} onCancel={onCancel} />,
    );
    expect(screen.queryByRole('dialog')).not.toBeInTheDocument();
  });

  it('disables confirm until a non-empty value is entered', async () => {
    const user = userEvent.setup();
    render(<PromptModal isOpen title="Name your profile" onConfirm={onConfirm} onCancel={onCancel} />);

    const confirm = screen.getByRole('button', { name: 'Save' });
    expect(confirm).toBeDisabled();

    await user.type(screen.getByRole('textbox'), '  my profile  ');
    expect(confirm).toBeEnabled();
  });

  it('confirms with the trimmed value', async () => {
    const user = userEvent.setup();
    render(<PromptModal isOpen title="Name it" onConfirm={onConfirm} onCancel={onCancel} />);
    await user.type(screen.getByRole('textbox'), '  Acme  ');
    await user.click(screen.getByRole('button', { name: 'Save' }));
    expect(onConfirm).toHaveBeenCalledWith('Acme');
  });

  it('submits on Enter and cancels on Escape', async () => {
    const user = userEvent.setup();
    render(<PromptModal isOpen title="Name it" onConfirm={onConfirm} onCancel={onCancel} />);
    const input = screen.getByRole('textbox');
    await user.type(input, 'foo{Enter}');
    expect(onConfirm).toHaveBeenCalledWith('foo');

    await user.keyboard('{Escape}');
    expect(onCancel).toHaveBeenCalledTimes(1);
  });
});
