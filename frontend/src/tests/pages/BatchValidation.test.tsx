import { describe, it, expect } from 'vitest';
import { MemoryRouter } from 'react-router-dom';

import { render, screen } from '../setup';
import { BatchValidationPage } from '../../pages/BatchValidation';
import { ConfigProvider } from '../../context/ConfigContext';
import { BatchCacheProvider } from '../../contexts/BatchCacheContext';

function renderPage() {
  return render(
    <ConfigProvider>
      <BatchCacheProvider>
        <MemoryRouter>
          <BatchValidationPage />
        </MemoryRouter>
      </BatchCacheProvider>
    </ConfigProvider>,
  );
}

describe('BatchValidationPage', () => {
  it('renders the page heading without crashing', () => {
    renderPage();
    expect(screen.getByRole('heading', { name: /batch validation/i })).toBeInTheDocument();
  });

  it('starts in the upload state', () => {
    renderPage();
    // The upload state shows the feature cards; the results/processing UI is absent.
    expect(screen.getByText(/ML-Readiness Scoring/i)).toBeInTheDocument();
    expect(screen.queryByText(/Detailed Results/i)).not.toBeInTheDocument();
  });
});
