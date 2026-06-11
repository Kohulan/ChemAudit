import { useEffect } from 'react';
import { render, waitFor } from '@testing-library/react';
import { describe, it, expect, vi, beforeEach } from 'vitest';
import { MemoryRouter, useNavigate, useSearchParams } from 'react-router-dom';
import { ScrollToTop } from './ScrollToTop';

const scrollToMock = vi.fn();

beforeEach(() => {
  scrollToMock.mockClear();
  window.scrollTo = scrollToMock as unknown as typeof window.scrollTo;
});

function NavigateOnce({ to }: { to: string }) {
  const navigate = useNavigate();
  useEffect(() => {
    navigate(to);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);
  return null;
}

function ChangeSearchOnce() {
  const [, setSearchParams] = useSearchParams();
  useEffect(() => {
    setSearchParams({ smiles: 'CCO' });
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);
  return null;
}

describe('ScrollToTop', () => {
  it('scrolls to top when the pathname changes', async () => {
    render(
      <MemoryRouter initialEntries={['/about']}>
        <ScrollToTop />
        <NavigateOnce to="/" />
      </MemoryRouter>
    );
    await waitFor(() =>
      expect(scrollToMock).toHaveBeenCalledWith({ top: 0, left: 0, behavior: 'instant' })
    );
  });

  it('scrolls to top when the same path is pushed again (header re-click)', async () => {
    render(
      <MemoryRouter initialEntries={['/batch']}>
        <ScrollToTop />
        <NavigateOnce to="/batch" />
      </MemoryRouter>
    );
    await waitFor(() => expect(scrollToMock).toHaveBeenCalled());
  });

  it('does not scroll on an in-page query-param update', async () => {
    render(
      <MemoryRouter initialEntries={['/']}>
        <ScrollToTop />
        <ChangeSearchOnce />
      </MemoryRouter>
    );
    // allow effects to flush
    await new Promise((r) => setTimeout(r, 50));
    expect(scrollToMock).not.toHaveBeenCalled();
  });

  it('does not scroll on initial load (lets the browser restore)', async () => {
    render(
      <MemoryRouter initialEntries={['/about']}>
        <ScrollToTop />
      </MemoryRouter>
    );
    await new Promise((r) => setTimeout(r, 50));
    expect(scrollToMock).not.toHaveBeenCalled();
  });
});
