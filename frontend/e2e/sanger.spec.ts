import { expect, test, type Page } from "@playwright/test";
import fs from "node:fs";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const DATA = path.resolve(import.meta.dirname, "../../tests/data");

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

/** Reference built from the real read by scripts/make_sanger_fixture.py. */
function referenceAroundTheRead(): string {
  return fs.readFileSync(path.resolve(import.meta.dirname, "fixtures/sanger-reference.txt"), "utf-8").trim();
}

test("align a Sanger read to the reference and show the chromatogram", async ({ page }) => {
  const reference = referenceAroundTheRead();
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("new").click();
  await page.getByLabel("Name").fill("pExpected");
  await page.getByTestId("new-seq-text").fill(reference);
  await page.getByTestId("new-seq-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText("pExpected");

  await menu(page, "Tools", "Align to Reference DNA Sequence…");
  await expect(page.getByTestId("sanger-dialog")).toBeVisible();
  await page.getByTestId("sanger-files").setInputFiles(path.join(DATA, "3730.ab1"));

  const summary = page.getByTestId("sanger-summary");
  await expect(summary).toBeVisible({ timeout: 30_000 });
  await expect(summary).toContainText("forward");
  await expect(summary).toContainText("100.00 %");
  await expect(summary).toContainText("none");

  // The read is stacked under the reference, with its chromatogram, in one view.
  await expect(page.getByTestId("sanger-alignment")).toBeVisible();
  await expect(page.locator(".scaleChromatogramButtonUp").first()).toBeVisible({ timeout: 20_000 });
  await expect(page.locator(".vt-alignment svg path").first()).toBeVisible({ timeout: 20_000 });
  await page.screenshot({ path: "test-results/sanger.png" });
});
