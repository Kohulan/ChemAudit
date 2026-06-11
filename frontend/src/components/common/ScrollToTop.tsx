import { useEffect, useRef } from 'react';
import { useLocation, useNavigationType } from 'react-router-dom';

/**
 * Resets window scroll on forward navigation (link clicks), app-wide.
 *
 * React Router does not touch scroll position, so navigating from a page
 * footer leaves the next page scrolled to the same offset. Rules:
 * - PUSH/REPLACE to a different pathname: scroll to top.
 * - PUSH to the same pathname with the same search (re-clicking a nav
 *   link for the page you are on): scroll to top.
 * - Same pathname with a changed search string (in-page param syncing,
 *   e.g. the validator writing ?smiles=...): leave scroll alone.
 * - POP (browser back/forward, initial load): leave scroll alone so the
 *   browser's native scroll restoration keeps working.
 *
 * Mounted once inside the Router in App.tsx; renders nothing.
 */
export function ScrollToTop() {
  const location = useLocation();
  const navigationType = useNavigationType();
  const prevLocation = useRef(location);

  useEffect(() => {
    const samePath = prevLocation.current.pathname === location.pathname;
    const sameSearch = prevLocation.current.search === location.search;
    prevLocation.current = location;

    if (navigationType === 'POP') return;
    if (samePath && !sameSearch) return;

    window.scrollTo({ top: 0, left: 0, behavior: 'instant' });
  }, [location, navigationType]);

  return null;
}
