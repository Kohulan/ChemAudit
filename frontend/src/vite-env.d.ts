/// <reference types="vite/client" />

interface ImportMetaEnv {
  /**
   * API Base URL
   * - Development: http://localhost:8000/api/v1
   * - Production: /api/v1 (relative, same-origin)
   */
  readonly VITE_API_URL: string;

  /**
   * API Documentation URL
   * - Development: http://localhost:8000/docs
   * - Production: /api/v1/docs
   */
  readonly VITE_API_DOCS_URL: string;

  /**
   * Debug mode flag
   * - Development: 'true'
   * - Production: 'false'
   */
  readonly VITE_DEBUG: string;

  /**
   * Vite mode (development | production)
   */
  readonly MODE: string;
}

interface ImportMeta {
  readonly env: ImportMetaEnv
}

interface Window {
  RDKit?: any
  initRDKitModule?: () => Promise<void>
}
