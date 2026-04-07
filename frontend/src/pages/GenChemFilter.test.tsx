/**
 * Vitest tests for the GenChem Filter page and types.
 *
 * Covers:
 * - PRESET_CONFIGS: exactly 4 presets matching PRESET_IDS (D-12)
 * - PRESET_CONFIGS: unique weight vectors per D-15
 * - PRESET_CONFIGS: weight vectors sum to 1.0
 * - PRESET_CONFIGS: per-preset property threshold assertions
 * - STAGE_COLORS: all 6 core stages + novelty stage defined (D-06)
 * - GenChemFilter component: renders page heading and empty state
 */
import { describe, it, expect } from 'vitest';
import { render, screen } from '@testing-library/react';
import { BrowserRouter } from 'react-router-dom';
import { PRESET_CONFIGS, PRESET_IDS, STAGE_COLORS } from '../types/genchem';

// =============================================================================
// PRESET_CONFIGS tests (D-12, D-15)
// =============================================================================

describe('PRESET_CONFIGS', () => {
  it('has exactly 4 presets matching PRESET_IDS', () => {
    const keys = Object.keys(PRESET_CONFIGS);
    expect(keys).toHaveLength(4);
    keys.forEach((k) => expect(PRESET_IDS.has(k)).toBe(true));
  });

  it('each preset has unique weight vector per D-15', () => {
    const weightVectors = Object.values(PRESET_CONFIGS).map(
      (c) => `${c.weight_validity}/${c.weight_qed}/${c.weight_alert_free}/${c.weight_sa}`,
    );
    const uniqueVectors = new Set(weightVectors);
    // All 4 presets must have different weight vectors
    expect(uniqueVectors.size).toBe(4);
  });

  it('all preset weight vectors sum to 1.0', () => {
    Object.entries(PRESET_CONFIGS).forEach(([_name, config]) => {
      const sum =
        config.weight_validity +
        config.weight_qed +
        config.weight_alert_free +
        config.weight_sa;
      expect(sum).toBeCloseTo(1.0, 5);
    });
  });

  it('drug_like preset has correct property thresholds', () => {
    const dl = PRESET_CONFIGS.drug_like;
    expect(dl.max_mw).toBe(500);
    expect(dl.max_rot_bonds).toBe(10);
    expect(dl.use_pains).toBe(true);
    expect(dl.use_nibr).toBe(false);
    // D-15: balanced weights
    expect(dl.weight_qed).toBeCloseTo(0.3, 5);
    expect(dl.weight_validity).toBeCloseTo(0.3, 5);
  });

  it('fragment_like preset has max_rings=3 and max_mw=300', () => {
    const fl = PRESET_CONFIGS.fragment_like;
    expect(fl.max_mw).toBe(300);
    expect(fl.max_rings).toBe(3);
    expect(fl.use_brenk).toBe(false);
    // D-15: emphasize alert-free and SA
    expect(fl.weight_alert_free).toBeCloseTo(0.3, 5);
    expect(fl.weight_sa).toBeCloseTo(0.3, 5);
  });

  it('lead_like preset emphasizes QED (D-15 weight_qed=0.4)', () => {
    const ll = PRESET_CONFIGS.lead_like;
    expect(ll.weight_qed).toBeCloseTo(0.4, 5);
    expect(ll.use_nibr).toBe(true);
    expect(ll.max_mw).toBe(350);
  });

  it('permissive preset disables all alert catalogs', () => {
    const pm = PRESET_CONFIGS.permissive;
    expect(pm.use_pains).toBe(false);
    expect(pm.use_brenk).toBe(false);
    expect(pm.use_kazius).toBe(false);
    expect(pm.use_nibr).toBe(false);
    // D-15: emphasize validity
    expect(pm.weight_validity).toBeCloseTo(0.4, 5);
    expect(pm.weight_alert_free).toBeCloseTo(0.1, 5);
  });
});

// =============================================================================
// STAGE_COLORS tests (D-06)
// =============================================================================

describe('STAGE_COLORS', () => {
  it('has entries for all 6 core stages plus novelty', () => {
    const expectedStages = [
      'parse',
      'valence',
      'alerts',
      'property',
      'sa_score',
      'dedup',
      'novelty',
    ];
    expectedStages.forEach((stage) => {
      expect(STAGE_COLORS[stage]).toBeDefined();
      expect(STAGE_COLORS[stage]).toMatch(/^#[0-9a-f]{6}$/i);
    });
  });
});

// =============================================================================
// GenChemFilter component render tests
// =============================================================================

describe('GenChemFilter page', () => {
  it('renders page heading', async () => {
    const GenChemFilter = (await import('./GenChemFilter')).default;
    render(
      <BrowserRouter>
        <GenChemFilter />
      </BrowserRouter>,
    );
    expect(screen.getByText('Generative Chemistry Filter')).toBeInTheDocument();
  });

  it('renders empty state when idle', async () => {
    const GenChemFilter = (await import('./GenChemFilter')).default;
    render(
      <BrowserRouter>
        <GenChemFilter />
      </BrowserRouter>,
    );
    expect(screen.getByText('No molecules filtered yet')).toBeInTheDocument();
  });

  it('renders Run Filter button', async () => {
    const GenChemFilter = (await import('./GenChemFilter')).default;
    render(
      <BrowserRouter>
        <GenChemFilter />
      </BrowserRouter>,
    );
    expect(screen.getByRole('button', { name: /run filter/i })).toBeInTheDocument();
  });
});
