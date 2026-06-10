import { describe, it, expect } from 'vitest';
import { MemoryRouter } from 'react-router-dom';

import { render, screen } from '../setup';
import { SingleValidationPage } from '../../pages/SingleValidation';
import { ConfigProvider } from '../../context/ConfigContext';
import { ValidationCacheProvider } from '../../contexts/ValidationCacheContext';

function renderPage(initialEntries = ['/']) {
  return render(
    <ConfigProvider>
      <ValidationCacheProvider>
        <MemoryRouter initialEntries={initialEntries}>
          <SingleValidationPage />
        </MemoryRouter>
      </ValidationCacheProvider>
    </ConfigProvider>,
  );
}

describe('SingleValidationPage', () => {
  it('renders the page heading and primary controls without crashing', () => {
    renderPage();
    expect(screen.getByRole('heading', { name: /molecule validation/i })).toBeInTheDocument();
    expect(screen.getByRole('button', { name: /^validate$/i })).toBeInTheDocument();
  });

  it('shows the primary validation tab labels', () => {
    renderPage();
    // Tab bar (TAB_ROW_1 / TAB_ROW_2) labels — present in the rendered tree.
    expect(screen.getAllByText('Safety').length).toBeGreaterThan(0);
    expect(screen.getAllByText('Standardize').length).toBeGreaterThan(0);
    expect(screen.getAllByText('Deep Validation').length).toBeGreaterThan(0);
  });

  it('exposes both primary actions (Validate and Score)', () => {
    renderPage();
    expect(screen.getByRole('button', { name: /^validate$/i })).toBeInTheDocument();
    expect(screen.getByRole('button', { name: /^score$/i })).toBeInTheDocument();
  });
});
