import { expect, test } from "@playwright/test";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

test("opens a .dna file and shows its map", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".vt-doctab.active")).toContainText("pFA-KanMX4");
  await expect(page.locator(".veCircularView")).toBeVisible();
  await expect(page.getByText(/3941/).first()).toBeVisible();
  await page.screenshot({ path: "test-results/open-dna.png" });
});

test("API rejects requests without the session token", async ({ request }) => {
  const res = await request.get("/api/health");
  expect(res.status()).toBe(403);
  const ok = await request.get("/api/health", { headers: { "X-Vettore-Token": TOKEN } });
  expect(ok.status()).toBe(200);
});
