import { expect, test, type Page } from "@playwright/test";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

test("align a capsid gene of another serotype to the AAV2 genome", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  // NCBI records (see tests/test_alignment.py): AAV2 genome and AAV9 VP1.
  await page.getByTestId("file-input").setInputFiles([
    data("analysis-NC_001401.2.gb"),
    data("analysis-AY530579.1.gb"),
  ]);
  await expect(page.getByTestId("doctab")).toHaveCount(2);

  await menu(page, "Tools", "Align Multiple Sequences…");
  const dialog = page.getByTestId("msa-dialog");
  await expect(dialog).toBeVisible();
  // The active tab (the last file opened) is the default reference; switch to the genome.
  await dialog.getByTestId("msa-reference").selectOption({ index: 0 });
  await expect(dialog.getByTestId("msa-choice")).toHaveCount(1);
  await expect(dialog.getByTestId("msa-choice")).toBeChecked();

  await dialog.getByTestId("msa-align").click();
  const summary = dialog.getByTestId("msa-summary");
  await expect(summary).toBeVisible({ timeout: 30_000 });
  // VP1 of AAV2 is annotated at 2203..4410 in NC_001401.2; AAV9 VP1 lands there.
  await expect(summary).toContainText("2,203–4,410");
  await expect(summary).toContainText("forward");
  await expect(summary).toContainText("mismatches");
  await expect(dialog.getByTestId("msa-differences").locator("li").first()).toBeVisible();

  await expect(dialog.getByTestId("msa-alignment")).toBeVisible();
  await expect(page.locator(".vt-msa-view svg").first()).toBeVisible({ timeout: 20_000 });
  await page.screenshot({ path: "test-results/alignment.png" });
});

async function newSequence(page: Page, name: string, sequence: string) {
  const welcome = page.getByTestId("new");
  if (await welcome.isVisible().catch(() => false)) await welcome.click();
  else await menu(page, "File", "New DNA File…");
  await page.getByLabel("Name").fill(name);
  await page.getByTestId("new-seq-text").fill(sequence);
  await page.getByTestId("new-seq-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText(name);
}

test("with several files open nothing is ticked until the user chooses", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await newSequence(page, "pA", "ACGTTGCAAGGCTTACCGATGGCATTCAGGACTTGACCATG".repeat(5));
  await newSequence(page, "pB", "TTGACCGGATACGATCGGCTAAGCTTGCATGCCTGCAGGTC".repeat(5));
  await newSequence(page, "pC", "GGCCTTAAGGCTAGCATCGATTACGGATCCAAGCTTGGTAC".repeat(5));

  await menu(page, "Tools", "Align Multiple Sequences…");
  const dialog = page.getByTestId("msa-dialog");
  await expect(dialog.getByTestId("msa-choice")).toHaveCount(2);
  for (const box of await dialog.getByTestId("msa-choice").all()) await expect(box).not.toBeChecked();
  // A neutral prompt, not an error, says what to do.
  await expect(dialog.getByTestId("msa-choice-prompt")).toContainText("Tick the sequences to compare");
  await expect(dialog.getByTestId("msa-choice-problem")).toHaveCount(0);
  await expect(dialog.getByTestId("msa-align")).toBeDisabled();

  await dialog.getByTestId("msa-choice").first().check();
  await expect(dialog.getByTestId("msa-choice-prompt")).toHaveCount(0);
  await expect(dialog.getByTestId("msa-align")).toBeEnabled();
});

/** A sequence without repeats, so that alignments have one answer. */
function randomDna(length: number, seed: number): string {
  let state = seed;
  let out = "";
  for (let i = 0; i < length; i++) {
    state = (state * 1103515245 + 12345) % 2147483648;
    out += "ACGT"[Math.floor(state / 65536) % 4];
  }
  return out;
}

test("files with the same name are told apart by their tab number", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  const reference = randomDna(300, 11);
  const swap = (i: number) => reference.slice(0, i) + (reference[i] === "A" ? "C" : "A") + reference.slice(i + 1);
  await newSequence(page, "pRef", reference);
  await newSequence(page, "Pasted sequence", swap(50));
  await newSequence(page, "Pasted sequence", swap(200));
  await page.getByTestId("doctab").first().click();

  await menu(page, "Tools", "Align Multiple Sequences…");
  const dialog = page.getByTestId("msa-dialog");
  const choices = dialog.getByTestId("msa-choices");
  await expect(choices).toContainText("Pasted sequence (tab 2)");
  await expect(choices).toContainText("Pasted sequence (tab 3)");
  await expect(dialog.getByTestId("msa-reference")).toContainText("pRef (300 bp)");
  for (const box of await dialog.getByTestId("msa-choice").all()) await box.check();
  await dialog.getByTestId("msa-align").click();

  const rows = dialog.getByTestId("msa-summary").locator("tbody tr");
  await expect(rows).toHaveCount(2, { timeout: 30_000 });
  // Each variant carries its own single mismatch, worked out from the swaps above.
  await expect(rows.filter({ hasText: "Pasted sequence (tab 2)" })).toContainText("1 mismatch");
  await rows.filter({ hasText: "Pasted sequence (tab 3)" }).click();
  await expect(dialog.getByTestId("msa-details")).toContainText("Pasted sequence (tab 3)");
  await expect(dialog.getByTestId("msa-differences")).toContainText(`201: ${reference[200]} → `);
});

test("explains what to do when only one sequence is open", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("analysis-AY530579.1.gb"));
  await expect(page.getByTestId("doctab")).toHaveCount(1);
  await menu(page, "Tools", "Align Multiple Sequences…");
  await expect(page.getByTestId("msa-empty")).toContainText("Open at least two sequences");
});
