/**
 * IndexedDB storage for bookmark result snapshots.
 *
 * Stores full validation tab results alongside bookmarks so users
 * can restore a complete validation view without re-running the pipeline.
 */

import { openDB, type IDBPDatabase } from 'idb';
import type { ValidationResponse } from '../types/validation';
import type { AlertScreenResponse } from '../types/alerts';
import type { ScoringResponse } from '../types/scoring';
import type { StandardizeResponse } from '../types/standardization';
import type { PubChemResult, ChEMBLResult, COCONUTResult } from '../types/integrations';

const DB_NAME = 'chemaudit-bookmarks';
const DB_VERSION = 1;
const STORE_NAME = 'snapshots';

export interface BookmarkSnapshot {
  bookmarkId: number;
  source: 'single_validation' | 'batch' | string;
  molecule: string;
  savedAt: string;
  // Single validation tab results (only what was loaded at bookmark time)
  validationResult?: ValidationResponse | null;
  alertResult?: AlertScreenResponse | null;
  scoringResult?: ScoringResponse | null;
  standardizationResult?: StandardizeResponse | null;
  databaseResults?: {
    pubchem: PubChemResult | null;
    chembl: ChEMBLResult | null;
    coconut: COCONUTResult | null;
  } | null;
  // Batch reference
  jobId?: string | null;
}

let dbPromise: Promise<IDBPDatabase> | null = null;

function getDB(): Promise<IDBPDatabase> {
  if (!dbPromise) {
    dbPromise = openDB(DB_NAME, DB_VERSION, {
      upgrade(db) {
        if (!db.objectStoreNames.contains(STORE_NAME)) {
          db.createObjectStore(STORE_NAME, { keyPath: 'bookmarkId' });
        }
      },
    });
  }
  return dbPromise;
}

export async function saveSnapshot(snapshot: BookmarkSnapshot): Promise<void> {
  const db = await getDB();
  await db.put(STORE_NAME, snapshot);
}

export async function getSnapshot(bookmarkId: number): Promise<BookmarkSnapshot | undefined> {
  const db = await getDB();
  return db.get(STORE_NAME, bookmarkId);
}

export async function deleteSnapshot(bookmarkId: number): Promise<void> {
  const db = await getDB();
  await db.delete(STORE_NAME, bookmarkId);
}

export async function getSnapshotIds(): Promise<Set<number>> {
  const db = await getDB();
  const keys = await db.getAllKeys(STORE_NAME);
  return new Set(keys as number[]);
}

export async function clearAllSnapshots(): Promise<void> {
  const db = await getDB();
  await db.clear(STORE_NAME);
}
