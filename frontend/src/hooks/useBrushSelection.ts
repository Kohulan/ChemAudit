/**
 * useBrushSelection Hook
 *
 * useReducer-based shared selection state for linked brushing across charts.
 * Selection is a Set<number> of molecule indices.
 */

import { useReducer, Dispatch } from 'react';

// ---------------------------------------------------------------------------
// Action types
// ---------------------------------------------------------------------------

export type SelectionAction =
  | { type: 'SET'; indices: Set<number> }
  | { type: 'TOGGLE'; index: number }
  | { type: 'ADD_RANGE'; indices: number[] }
  | { type: 'CLEAR' };

// ---------------------------------------------------------------------------
// Reducer
// ---------------------------------------------------------------------------

export function selectionReducer(
  state: Set<number>,
  action: SelectionAction
): Set<number> {
  switch (action.type) {
    case 'SET':
      return new Set(action.indices);
    case 'TOGGLE': {
      const next = new Set(state);
      if (next.has(action.index)) {
        next.delete(action.index);
      } else {
        next.add(action.index);
      }
      return next;
    }
    case 'ADD_RANGE': {
      const next = new Set(state);
      for (const idx of action.indices) {
        next.add(idx);
      }
      return next;
    }
    case 'CLEAR':
      return new Set();
    default:
      return state;
  }
}

// ---------------------------------------------------------------------------
// Action creators
// ---------------------------------------------------------------------------

export function setSelection(indices: Set<number>): SelectionAction {
  return { type: 'SET', indices };
}

export function toggleIndex(index: number): SelectionAction {
  return { type: 'TOGGLE', index };
}

export function addRange(indices: number[]): SelectionAction {
  return { type: 'ADD_RANGE', indices };
}

export function clearSelection(): SelectionAction {
  return { type: 'CLEAR' };
}

// ---------------------------------------------------------------------------
// Hook
// ---------------------------------------------------------------------------

export function useBrushSelection(): [Set<number>, Dispatch<SelectionAction>] {
  return useReducer(selectionReducer, new Set<number>());
}
