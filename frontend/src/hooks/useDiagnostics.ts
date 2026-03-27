import { useState, useCallback } from 'react';
import { diagnosticsApi } from '../services/api';
import type {
  SMILESDiagnosticsResponse,
  InChIDiffResponse,
  RoundTripResponse,
  CrossPipelineResponse,
  FilePreValidationResponse,
} from '../types/diagnostics';

/**
 * Hook for Structure Quality Diagnostics API calls.
 *
 * Exposes independent per-tool state (result, loading, error) and action
 * callbacks. Each tool's state is isolated so a failure in one does not
 * affect the others (per UI-SPEC D-01 through D-05).
 *
 * Uses diagnosticsApi from services/api.ts (which uses the project's axios
 * instance with CSRF tokens, auth headers, and error interceptors).
 */
export function useDiagnostics() {
  // ── Result state — one per diagnostic tool ──
  const [smilesResult, setSmilesResult] = useState<SMILESDiagnosticsResponse | null>(null);
  const [inchiDiffResult, setInchiDiffResult] = useState<InChIDiffResponse | null>(null);
  const [roundtripResult, setRoundtripResult] = useState<RoundTripResponse | null>(null);
  const [crossPipelineResult, setCrossPipelineResult] = useState<CrossPipelineResponse | null>(null);
  const [fileResult, setFileResult] = useState<FilePreValidationResponse | null>(null);

  // ── Loading state — one per tool ──
  const [smilesLoading, setSmilesLoading] = useState(false);
  const [inchiDiffLoading, setInchiDiffLoading] = useState(false);
  const [roundtripLoading, setRoundtripLoading] = useState(false);
  const [crossPipelineLoading, setCrossPipelineLoading] = useState(false);
  const [fileLoading, setFileLoading] = useState(false);

  // ── Error state — one per tool ──
  const [smilesError, setSmilesError] = useState<string | null>(null);
  const [inchiDiffError, setInchiDiffError] = useState<string | null>(null);
  const [roundtripError, setRoundtripError] = useState<string | null>(null);
  const [crossPipelineError, setCrossPipelineError] = useState<string | null>(null);
  const [fileError, setFileError] = useState<string | null>(null);

  // Current SMILES for re-use across tools (e.g. round-trip, cross-pipeline)
  const [currentSmiles, setCurrentSmiles] = useState<string | null>(null);

  /**
   * Run SMILES diagnostics (DIAG-01): position-specific errors + fix suggestions.
   * Sets currentSmiles for downstream tools.
   */
  const analyzeMolecule = useCallback(async (smiles: string) => {
    setCurrentSmiles(smiles);
    setSmilesLoading(true);
    setSmilesError(null);
    try {
      const result = await diagnosticsApi.smiles(smiles);
      setSmilesResult(result);
    } catch (e: unknown) {
      const err = e as { error?: string; detail?: string };
      setSmilesError(err?.error ?? err?.detail ?? 'SMILES diagnostics failed');
    } finally {
      setSmilesLoading(false);
    }
  }, []);

  /**
   * Compare two InChI strings layer-by-layer (DIAG-02).
   */
  const compareInchi = useCallback(async (inchiA: string, inchiB: string) => {
    setInchiDiffLoading(true);
    setInchiDiffError(null);
    try {
      const result = await diagnosticsApi.inchiDiff(inchiA, inchiB);
      setInchiDiffResult(result);
    } catch (e: unknown) {
      const err = e as { error?: string; detail?: string };
      setInchiDiffError(err?.error ?? err?.detail ?? 'InChI comparison failed');
    } finally {
      setInchiDiffLoading(false);
    }
  }, []);

  /**
   * Check format round-trip lossiness (DIAG-03).
   * Defaults to smiles_inchi_smiles route.
   */
  const checkRoundtrip = useCallback(
    async (smiles: string, route: string = 'smiles_inchi_smiles') => {
      setRoundtripLoading(true);
      setRoundtripError(null);
      try {
        const result = await diagnosticsApi.roundtrip(smiles, route);
        setRoundtripResult(result);
      } catch (e: unknown) {
        const err = e as { error?: string; detail?: string };
        setRoundtripError(err?.error ?? err?.detail ?? 'Round-trip check failed');
      } finally {
        setRoundtripLoading(false);
      }
    },
    [],
  );

  /**
   * Compare standardization output across 3 pipelines (DIAG-04).
   */
  const comparePipelines = useCallback(async (molecule: string) => {
    setCrossPipelineLoading(true);
    setCrossPipelineError(null);
    try {
      const result = await diagnosticsApi.crossPipeline(molecule);
      setCrossPipelineResult(result);
    } catch (e: unknown) {
      const err = e as { error?: string; detail?: string };
      setCrossPipelineError(err?.error ?? err?.detail ?? 'Pipeline comparison failed');
    } finally {
      setCrossPipelineLoading(false);
    }
  }, []);

  /**
   * Pre-validate an SDF or CSV file (DIAG-05).
   */
  const prevalidateFile = useCallback(async (file: File) => {
    setFileLoading(true);
    setFileError(null);
    try {
      const result = await diagnosticsApi.filePrevalidate(file);
      setFileResult(result);
    } catch (e: unknown) {
      const err = e as { error?: string; detail?: string };
      setFileError(err?.error ?? err?.detail ?? 'File pre-validation failed');
    } finally {
      setFileLoading(false);
    }
  }, []);

  /**
   * Clear file result and error state (e.g. when user removes the file).
   */
  const clearFileResult = useCallback(() => {
    setFileResult(null);
    setFileError(null);
  }, []);

  return {
    // Results
    smilesResult,
    inchiDiffResult,
    roundtripResult,
    crossPipelineResult,
    fileResult,
    // Loading flags
    smilesLoading,
    inchiDiffLoading,
    roundtripLoading,
    crossPipelineLoading,
    fileLoading,
    // Errors
    smilesError,
    inchiDiffError,
    roundtripError,
    crossPipelineError,
    fileError,
    // Current molecule state
    currentSmiles,
    // Actions
    analyzeMolecule,
    compareInchi,
    checkRoundtrip,
    comparePipelines,
    prevalidateFile,
    clearFileResult,
  };
}
