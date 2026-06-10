import { defineConfig, devices } from '@playwright/test';

/**
 * End-to-end tests for ChemAudit.
 *
 * Orchestrates the full stack: the FastAPI backend on :8001 and the Vite dev
 * server on :3002 (which proxies /api and /ws to the backend). Requires Redis
 * and Postgres reachable at their default URLs (see docker-compose).
 *
 * The backend launch command is overridable via E2E_BACKEND_CMD so CI / local
 * environments can point at the correct Python interpreter, e.g.:
 *   E2E_BACKEND_CMD="/path/to/venv/bin/python -m uvicorn app.main:app --port 8001"
 */
export default defineConfig({
  testDir: './e2e',
  timeout: 60_000,
  expect: { timeout: 15_000 },
  fullyParallel: false,
  retries: process.env.CI ? 1 : 0,
  reporter: 'list',
  use: {
    baseURL: 'http://localhost:3002',
    headless: true,
    trace: 'on-first-retry',
  },
  projects: [{ name: 'chromium', use: { ...devices['Desktop Chrome'] } }],
  webServer: [
    {
      command:
        process.env.E2E_BACKEND_CMD ||
        'python -m uvicorn app.main:app --port 8001',
      cwd: '../backend',
      port: 8001,
      timeout: 120_000,
      reuseExistingServer: !process.env.CI,
      env: { DEBUG: 'true' },
    },
    {
      command: 'npm run dev',
      port: 3002,
      timeout: 120_000,
      reuseExistingServer: !process.env.CI,
    },
  ],
});
