// Copy buttons in dialogs and views: they confirm only what really reached the clipboard, and
// when the system refuses it they show the text, selected, to copy by hand.
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

async function newSequence(page: Page, name: string, sequence: string, circular: boolean) {
  const welcome = page.getByTestId("new");
  if (await welcome.isVisible().catch(() => false)) {
    await welcome.click();
  } else {
    await menu(page, "File", "New DNA File…");
  }
  await page.getByLabel("Name").fill(name);
  await page.getByTestId("new-seq-text").fill(sequence);
  if (circular) {
    await page.locator("label.bp3-radio", { hasText: "Circular" }).click();
  }
  await page.getByTestId("new-seq-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText(name);
}

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

const BACKBONE = "ATTGCCTACGGTACATGCAACGTTGACCATGGCATTACGGACCATTGACCAGTTACGGATCAGCTTAGGCATAC".repeat(14);
const VECTOR = BACKBONE.slice(0, 500) + "GAATTC" + BACKBONE.slice(500, 900) + "GGATCC" + BACKBONE.slice(900, 1400);
const GFP = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCC";

async function gibson(page: Page) {
  await page.goto(`/?token=${TOKEN}`);
  await newSequence(page, "pVector", VECTOR, true);
  await newSequence(page, "GFP", GFP, false);
  await menu(page, "Actions", "Gibson Assembly…");
  const dialog = page.getByTestId("cloning-dialog");
  await dialog.getByTestId("vector-doc").selectOption({ index: 1 });
  await dialog.getByTestId("vector-enzymes").fill("EcoRI");
  await dialog.getByTestId("insert-doc-0").selectOption({ index: 2 });
  await expect(page.getByTestId("cloning-oligos").locator("tbody tr")).toHaveCount(2);
  return dialog;
}

/** Reference and a copy with base 101 changed from its own letter: one mismatch to list. */
async function alignTwo(page: Page) {
  await page.goto(`/?token=${TOKEN}`);
  const reference = randomDna(400, 7);
  const changed = reference[100] === "A" ? "C" : "A";
  await newSequence(page, "pRef", reference, false);
  await newSequence(page, "pVariant", reference.slice(0, 100) + changed + reference.slice(101), false);
  await menu(page, "Tools", "Align Multiple Sequences…");
  const dialog = page.getByTestId("msa-dialog");
  await dialog.getByTestId("msa-reference").selectOption({ index: 0 });
  await dialog.getByTestId("msa-choice").check();
  await dialog.getByTestId("msa-align").click();
  await expect(dialog.getByTestId("msa-differences")).toContainText(`101: ${reference[100]} → ${changed}`, { timeout: 30_000 });
  return dialog;
}

const copyBox = (page: Page) => page.getByTestId("copy-text-dialog");

test.describe("when the system refuses the clipboard", () => {
  test.beforeEach(async ({ page }) => {
    await page.addInitScript(() => {
      Object.defineProperty(navigator, "clipboard", { value: undefined });
      document.execCommand = () => false;
    });
  });

  test("Copy all in a cloning dialog shows the primers to copy by hand", async ({ page }) => {
    const dialog = await gibson(page);
    await dialog.getByTestId("cloning-copy").click();
    await expect(copyBox(page)).toContainText("did not let Vettore put this on the clipboard");
    await expect(page.getByTestId("copy-text")).toHaveValue(/^\S+\t[ACGT]+\n\S+\t[ACGT]+$/);
    await expect(page.getByTestId("cloning-copied")).toHaveCount(0);
    // Escape closes the box only: the cloning dialog and its inputs stay.
    await page.keyboard.press("Escape");
    await expect(copyBox(page)).toBeHidden();
    await expect(dialog).toBeVisible();
  });

  test("Copy for Ordering and the History copy show the text instead of claiming a copy", async ({ page }) => {
    await page.goto(`/?token=${TOKEN}`);
    await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
    await expect(page.locator(".veCircularView")).toBeVisible();
    await menu(page, "Primers", "Design PCR Primers…");
    const dialog = page.getByTestId("primer-design-dialog");
    await dialog.getByTestId("pd-feature").selectOption({ label: "KanR (459–1,268)" });
    await expect(dialog.getByTestId("pd-pair")).toHaveCount(5, { timeout: 30_000 });

    await page.getByTestId("pd-copy").click();
    await expect(page.getByTestId("copy-text")).toHaveValue(/^KanR_F\t[ACGT]{18,24}\nKanR_R\t[ACGT]{18,24}$/);
    await expect(page.getByTestId("pd-notice")).toHaveCount(0);
    await page.keyboard.press("Escape");
    await expect(copyBox(page)).toBeHidden();

    await page.getByTestId("pd-add").click();
    await page.getByTestId("pd-close").click();
    await page.getByTestId("view-history").click();
    await page.getByTestId("history-copy").first().click();
    await expect(page.getByTestId("copy-text")).toHaveValue(/^Added primers \(primer3\): KanR_F [ACGT]+, KanR_R [ACGT]+; product/);
    await expect(page.getByTestId("history-copy").first()).toHaveText("Copy");
  });

  test("Copy list in Align Multiple Sequences shows the list to copy by hand", async ({ page }) => {
    const dialog = await alignTwo(page);
    await dialog.getByTestId("msa-copy").click();
    await expect(page.getByTestId("copy-text")).toHaveValue(/^sequence\tkind\tstart\tend.*\npVariant\tmismatch\t101\t101\t/s);
    await expect(dialog.getByTestId("msa-copy")).toHaveText("Copy list (for a spreadsheet)");
  });
});

test.describe("with clipboard access", () => {
  test.use({ permissions: ["clipboard-read", "clipboard-write"] });

  test("Copy all confirms what it put on the clipboard", async ({ page }) => {
    const dialog = await gibson(page);
    await dialog.getByTestId("cloning-copy").click();
    await expect(page.getByTestId("cloning-copied")).toHaveText("Copied 2 primers: paste them into the order form.");
    const names = await dialog.getByTestId("cloning-oligos").locator("tbody tr td:first-child").allTextContents();
    const clipboard = await page.evaluate(() => navigator.clipboard.readText());
    expect(clipboard.split("\n").map((line) => line.split("\t")[0])).toEqual(names);
    await expect(copyBox(page)).toHaveCount(0);
  });

  test("Copy list confirms the copy on the button", async ({ page }) => {
    const dialog = await alignTwo(page);
    await dialog.getByTestId("msa-copy").click();
    await expect(dialog.getByTestId("msa-copy")).toHaveText("Copied: paste it into a spreadsheet");
    const clipboard = await page.evaluate(() => navigator.clipboard.readText());
    expect(clipboard).toMatch(/\npVariant\tmismatch\t101\t101\t/);
  });
});
