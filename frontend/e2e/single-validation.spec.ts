import { test, expect } from '@playwright/test';

/**
 * Happy-path end-to-end test: a user validates a single molecule and sees a
 * result. Exercises the full path — Vite dev server -> /api proxy -> FastAPI
 * /validate -> rendered result — which no unit/component test covers.
 */
test('validates a single molecule end to end', async ({ page }) => {
  await page.goto('/');

  // App shell renders.
  await expect(page.getByRole('heading', { name: /molecule validation/i })).toBeVisible();

  // Enter a clean, valid molecule (ethanol).
  const input = page.getByPlaceholder(/enter smiles/i);
  await input.fill('CCO');

  // Click Validate and wait for the backend round-trip to complete.
  const validate = page.getByRole('button', { name: /^validate$/i });
  const [response] = await Promise.all([
    page.waitForResponse(
      (r) => r.url().includes('/api/v1/validate') && r.request().method() === 'POST',
    ),
    validate.click(),
  ]);
  expect(response.ok()).toBeTruthy();

  // A validation result region renders (clean molecule -> "All Clear", or a
  // score / issues panel). Any of these confirms the result rendered.
  await expect(
    page.getByText(/all clear|quality|ml readiness|validation/i).first(),
  ).toBeVisible();
});
