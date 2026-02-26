# Bookmark Result Snapshots

**Date**: 2026-02-25
**Status**: Approved

## Problem

Bookmarks currently store only SMILES + metadata. Clicking a bookmark does nothing — no navigation, no result restore. Users expect to bookmark a complete validation result and return to it later.

## Decisions

- **Save scope**: All loaded tabs at bookmark time (validation, scoring, alerts, standardization, database lookups)
- **Storage**: Browser IndexedDB via `idb` package (per-device, no backend migration)
- **Batch restore**: Re-fetch by job_id (navigate to `/batch?job={id}`)
- **Restore UX**: Navigate to SingleValidation and instantly display saved results

## Architecture

### Storage Layer — `frontend/src/lib/bookmarkStore.ts`

IndexedDB database `chemaudit-bookmarks`, object store `snapshots`, keyed by bookmark `id` (number).

```ts
interface BookmarkSnapshot {
  bookmarkId: number;
  source: 'single_validation' | 'batch';
  molecule: string;
  savedAt: string;
  validationResult?: ValidationResponse;
  alertResult?: AlertScreenResponse;
  scoringResult?: ScoringResponse;
  standardizationResult?: StandardizeResponse;
  databaseResults?: { pubchem, chembl, coconut };
  jobId?: string;
}
```

Three exports: `saveSnapshot()`, `getSnapshot()`, `deleteSnapshot()`.

### Save Flow

1. `BookmarkButton` gets `onAfterBookmark?: (bookmarkId: number) => void`
2. After API creates bookmark, calls `onAfterBookmark(id)`
3. `SingleValidation` gathers loaded tab state → `saveSnapshot()`
4. Batch bookmarks save `{ source: 'batch', jobId }` only

### Restore Flow

1. `BookmarkList` gets `onOpenBookmark?: (bookmark: Bookmark) => void`
2. Rows get clickable Open button
3. `Bookmarks.tsx` handler:
   - Batch + job_id → navigate `/batch?job={job_id}`
   - Single → navigate `/` with `location.state = { bookmarkId, smiles }`
4. `SingleValidation` on mount checks `location.state.bookmarkId`:
   - Snapshot found → hydrate all tabs instantly
   - Not found → set SMILES, auto-trigger validation

### Visual Indicator

BookmarkList checks IndexedDB for which IDs have snapshots. Shows badge on rows with saved results.

### Cleanup

Delete snapshot from IndexedDB when bookmark is deleted (single or bulk).

## Files to Modify

- `frontend/src/lib/bookmarkStore.ts` (new)
- `frontend/src/components/bookmarks/BookmarkButton.tsx`
- `frontend/src/components/bookmarks/BookmarkList.tsx`
- `frontend/src/pages/Bookmarks.tsx`
- `frontend/src/pages/SingleValidation.tsx`
- `frontend/src/hooks/useValidation.ts` (restore/hydrate method)
- `package.json` (add `idb` dependency)
