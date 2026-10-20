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

async function newSequence(page: Page, name: string, sequence: string) {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("new").click();
  await page.getByLabel("Name").fill(name);
  await page.getByTestId("new-seq-text").fill(sequence);
  await page.getByTestId("new-seq-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText(name);
}

test("align a Sanger read to the reference and show the chromatogram", async ({ page }) => {
  await newSequence(page, "pExpected", referenceAroundTheRead());

  await menu(page, "Tools", "Align to Reference DNA Sequence…");
  await expect(page.getByTestId("sanger-dialog")).toBeVisible();
  // Before any file is chosen the dialog says what to do.
  await expect(page.getByTestId("sanger-steps")).toContainText("Choose the .ab1 files");
  await page.getByTestId("sanger-files").setInputFiles(path.join(DATA, "3730.ab1"));

  const summary = page.getByTestId("sanger-summary");
  await expect(summary).toBeVisible({ timeout: 30_000 });
  await expect(page.getByTestId("sanger-verdict")).toContainText("No differences: the read matches pExpected");
  await expect(summary).toContainText("forward");
  await expect(summary).toContainText("100.00 %");
  await expect(summary).toContainText("none");
  await expect(summary).toContainText(/Q\d\d/);

  // The read is stacked under the reference, with its chromatogram, in one view.
  await expect(page.getByTestId("sanger-alignment")).toBeVisible();
  await expect(page.locator(".scaleChromatogramButtonUp").first()).toBeVisible({ timeout: 20_000 });
  await expect(page.locator(".vt-alignment svg path").first()).toBeVisible({ timeout: 20_000 });
  await page.screenshot({ path: "test-results/sanger.png" });
});

test("a planted mismatch is reported with the read quality behind it", async ({ page }) => {
  // The fixture is 4 × 74 bases of filler, then the read's core: base 497 lies inside the read.
  const reference = referenceAroundTheRead();
  const index = 4 * 74 + 200;
  const replacement = reference[index] === "A" ? "G" : "A";
  await newSequence(page, "pMutant", reference.slice(0, index) + replacement + reference.slice(index + 1));

  await menu(page, "Tools", "Align to Reference DNA Sequence…");
  await page.getByTestId("sanger-files").setInputFiles(path.join(DATA, "3730.ab1"));
  await expect(page.getByTestId("sanger-verdict")).toContainText("1 difference in 1 read", { timeout: 30_000 });
  const differences = page.getByTestId("sanger-differences");
  await expect(differences).toContainText(`${index + 1}: ${replacement} → ${reference[index]}`);
  await expect(differences).toContainText(/\(Q\d+/);
});
