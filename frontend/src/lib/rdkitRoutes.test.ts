import { describe, it, expect } from 'vitest';
import { routeNeedsRDKit } from './rdkitRoutes';

describe('routeNeedsRDKit', () => {
  it.each(['/about', '/privacy', '/diagnostics', '/about/'])(
    'returns false for RDKit-free route %s',
    (path) => {
      expect(routeNeedsRDKit(path)).toBe(false);
    }
  );

  // History and Bookmarks render molecule thumbnails, so they DO need RDKit.
  it.each(['/', '/batch', '/history', '/bookmarks', '/qsar-ready', '/structure-filter', '/dataset-audit', '/report/abc123'])(
    'returns true for structure route %s',
    (path) => {
      expect(routeNeedsRDKit(path)).toBe(true);
    }
  );
});
