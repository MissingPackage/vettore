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

test("a mutation seen by two reads is counted once", async ({ page }) => {
  // One base changed in the reference and the same trace loaded twice, as a forward and a
  // reverse read over the same insert would be: one mutated position, reported by both reads.
  const reference = referenceAroundTheRead();
  const index = 4 * 74 + 200;
  const replacement = reference[index] === "A" ? "G" : "A";
  await newSequence(page, "pMutant", reference.slice(0, index) + replacement + reference.slice(index + 1));

  await menu(page, "Tools", "Align to Reference DNA Sequence…");
  const trace = fs.readFileSync(path.join(DATA, "3730.ab1"));
  await page.getByTestId("sanger-files").setInputFiles([
    { name: "read_F.ab1", mimeType: "application/octet-stream", buffer: trace },
    { name: "read_R.ab1", mimeType: "application/octet-stream", buffer: trace },
  ]);
  // The verdict comes from the consensus of both reads: one difference, confirmed by two reads.
  const verdict = page.getByTestId("sanger-verdict");
  await expect(verdict).toContainText("1 difference between the consensus of 2 reads and pMutant", { timeout: 30_000 });
  await expect(verdict).toContainText("It is confirmed");
  // The per-read table still lists the difference under each read.
  const differences = page.getByTestId("sanger-differences");
  await expect(differences).toHaveCount(2);
  await expect(differences.nth(1)).toContainText(`${index + 1}: ${replacement} → ${reference[index]}`);
});

test("two reads give one consensus that opens as a new document", async ({ page }) => {
  const reference = referenceAroundTheRead();
  const index = 4 * 74 + 200;
  const replacement = reference[index] === "A" ? "G" : "A";
  await newSequence(page, "pMutant", reference.slice(0, index) + replacement + reference.slice(index + 1));

  await menu(page, "Tools", "Align to Reference DNA Sequence…");
  const trace = fs.readFileSync(path.join(DATA, "3730.ab1"));
  await page.getByTestId("sanger-files").setInputFiles([
    { name: "read_1.ab1", mimeType: "application/octet-stream", buffer: trace },
    { name: "read_2.ab1", mimeType: "application/octet-stream", buffer: trace },
  ]);
  const contig = page.getByTestId("sanger-contig");
  await expect(contig).toBeVisible({ timeout: 30_000 });
  // Both reads are the same forward trace: the covered stretch is read on one strand only.
  await expect(page.getByTestId("sanger-coverage")).toContainText("covered on one strand only");
  await expect(page.getByTestId("sanger-coverage")).toContainText("no read over 1–");
  await expect(page.getByTestId("sanger-coverage-strip")).toBeVisible();
  const differences = page.getByTestId("sanger-contig-differences");
  await expect(differences).toContainText(`${(index + 1).toLocaleString("en-US")}: ${replacement} → ${reference[index]}`);
  await expect(differences).toContainText("seen in 2 reads");
  await expect(differences).toContainText("confirmed");
  await page.screenshot({ path: "test-results/sanger-contig.png" });

  await page.getByTestId("sanger-save-consensus").click();
  await page.getByRole("button", { name: "Close" }).last().click();
  await expect(page.locator(".vt-doctab.active")).toContainText("pMutant consensus");
});
