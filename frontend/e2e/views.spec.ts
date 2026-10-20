import { expect, test } from "@playwright/test";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

test.beforeEach(async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".vt-doctab.active")).toContainText("pFA-KanMX4");
});

test("view tabs switch the document view", async ({ page }) => {
  await expect(page.locator(".veCircularView")).toBeVisible();
  await page.screenshot({ path: "test-results/view-map.png" });

  await page.getByTestId("view-sequence").click();
  await expect(page.locator(".veRowView")).toBeVisible();
  await page.screenshot({ path: "test-results/view-sequence.png" });

  await page.getByTestId("view-features").click();
  await expect(page.locator(".ve-propertiesPanel")).toBeVisible();
  await expect(page.getByText("kanMX").first()).toBeVisible();
  await page.screenshot({ path: "test-results/view-features.png" });

  await page.getByTestId("view-enzymes").click();
  await expect(page.getByTestId("enzymes-view")).toBeVisible();
  // pFA-KanMX4: NdeI is a unique 6-cutter in the default set.
  await expect(page.getByTestId("enzymes-view").getByRole("cell", { name: "NdeI", exact: true })).toBeVisible();
  // Methylation is a column of its own; which enzymes it flags is covered by the unit tests.
  await expect(page.getByTestId("enzymes-view")).toContainText("Methylation");
  await page.screenshot({ path: "test-results/view-enzymes.png" });
  await expect(page.locator(".ve-menu-bar, .tg-menu-bar").first()).toBeVisible();

  await page.getByTestId("view-primers").click();
  await page.screenshot({ path: "test-results/view-primers.png" });

  await page.getByTestId("view-history").click();
  await expect(page.locator(".vt-history")).toContainText("Opened");

  await page.keyboard.press("Control+1");
  await expect(page.locator(".veCircularView")).toBeVisible();
});

test("status bar shows length, topology and enzyme set", async ({ page }) => {
  await expect(page.getByTestId("status-line")).toContainText("3941 bp · circular");
  await expect(page.locator(".vt-status-set")).toHaveText("Unique 6+ Cutters");
});
