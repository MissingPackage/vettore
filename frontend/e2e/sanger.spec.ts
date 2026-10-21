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

/** The two public contig reads agree on this reference (scripts/make_contig_fixture.py). */
function contigReference(): string {
  return fs.readFileSync(path.resolve(import.meta.dirname, "fixtures/contig-reference.txt"), "utf-8").trim();
}

const PAIR = ["contig-ACHLO006-09_1_F.ab1", "contig-ACHLO006-09_2_R.ab1"].map((f) => path.join(DATA, f));

test("the same trace chosen twice is counted once", async ({ page }) => {
  await newSequence(page, "pExpected", referenceAroundTheRead());
  await menu(page, "Tools", "Align to Reference DNA Sequence…");
  const trace = fs.readFileSync(path.join(DATA, "3730.ab1"));
  await page.getByTestId("sanger-files").setInputFiles([
    { name: "read.ab1", mimeType: "application/octet-stream", buffer: trace },
    { name: "read (1).ab1", mimeType: "application/octet-stream", buffer: trace },
  ]);
  await expect(page.getByTestId("sanger-warning")).toContainText("is the same trace as “read.ab1”: it is counted once",
    { timeout: 30_000 });
  await expect(page.getByTestId("sanger-verdict")).toContainText("No differences: the read matches pExpected");
  await expect(page.getByTestId("sanger-contig")).toHaveCount(0);
});

test("a forward and a reverse read give one consensus that opens as a new document", async ({ page }) => {
  // One base changed in the reference where both reads cover it: the clone differs there.
  const reference = contigReference();
  const index = 3 * 74 + 300;
  const replacement = reference[index] === "A" ? "G" : "A";
  await newSequence(page, "pCOI", reference.slice(0, index) + replacement + reference.slice(index + 1));

  await menu(page, "Tools", "Align to Reference DNA Sequence…");
  await page.getByTestId("sanger-files").setInputFiles(PAIR);
  const verdict = page.getByTestId("sanger-verdict");
  await expect(verdict).toContainText("1 difference between the consensus of 2 reads and pCOI", { timeout: 30_000 });
  await expect(verdict).toContainText("It is confirmed");
  // The forward read makes two N calls where it is alone: those bases are kept, unverified.
  await expect(verdict).toContainText("2 covered positions have no clear call in any read");

  await expect(page.getByTestId("sanger-contig")).toBeVisible();
  const coverage = page.getByTestId("sanger-coverage");
  await expect(coverage).toContainText("223–861 covered on both strands");
  await expect(coverage).toContainText("862–903 on one strand only");
  await expect(coverage).toContainText("no read over 1–222 and 904–1,125");
  await expect(page.getByTestId("sanger-coverage-strip")).toBeVisible();
  const differences = page.getByTestId("sanger-contig-differences");
  await expect(differences).toContainText(`${index + 1}: ${replacement} → ${reference[index]}`);
  await expect(differences).toContainText("seen in 2 reads");
  await expect(differences).toContainText("confirmed");
  await page.screenshot({ path: "test-results/sanger-contig.png" });

  await page.getByTestId("sanger-save-consensus").click();
  await page.getByRole("button", { name: "Close" }).last().click();
  await expect(page.locator(".vt-doctab.active")).toContainText("pCOI consensus");
});
