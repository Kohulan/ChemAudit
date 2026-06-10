import { render } from '@testing-library/react';
import { describe, it, expect } from 'vitest';
import { Badge, CountBadge } from './Badge';

describe('Badge status differentiation', () => {
  it('adds a distinguishing icon to success badges', () => {
    const { container } = render(<Badge variant="success">Valid</Badge>);
    expect(container.querySelector('svg')).not.toBeNull();
  });

  it('adds a distinguishing icon to warning badges', () => {
    const { container } = render(<Badge variant="warning">Check</Badge>);
    expect(container.querySelector('svg')).not.toBeNull();
  });

  it('respects an explicit icon override', () => {
    const { container } = render(
      <Badge variant="success" icon={<span data-testid="custom" />}>Valid</Badge>
    );
    expect(container.querySelector('svg')).toBeNull();
  });

  it('keeps numeric CountBadge icon-free', () => {
    const { container } = render(<CountBadge count={5} variant="success" />);
    expect(container.querySelector('svg')).toBeNull();
  });
});
