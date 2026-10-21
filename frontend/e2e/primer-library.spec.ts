// Primers → Primer Library (S5.3): add a primer to the shared library, import pasted rows,
// and see a library primer drawn on the open sequence (Show Library Primers).
import { expect, test, type Page } from "@playwright/test";
import { readFileSync } from "node:fs";
import path from "node:path";
import { TOKEN } from "../playwright.config";

// Both tests write the same shared library file: run them one after the other, or the second
// save rightly reports that the file changed on disk.
test.describe.configure({ mode: "serial" });

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

// A primer whose binding site is spelled out below, unique so the test does not depend on what
// else is in the shared library file.
const PRIMER = "GCATTCGAATGACCTGATTAACGT";
const SEQUENCE = "AAAACCCCGGGGTTTTACAC" + PRIMER + "TTGGCCAATTGGCCAATTGG";
// A fresh oligo each run, so an earlier run's import does not count as a duplicate.
const IMPORTED = Array.from({ length: 22 }, () => "ACGT"[Math.floor(Math.random() * 4)]).join("");

async function openLibrary(page: Page) {
  await menu(page, "Primers", "Primer Library…");
  const dialog = page.getByTestId("primer-library-dialog");
  await expect(dialog).toBeVisible();
  return dialog;
}

test("add a primer, import pasted rows, and keep it in the shared file", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  // A document so the dialog's "on this sequence" panel has something to work on. The menu bar
  // appears only once an editor is open, so the first document comes from the welcome screen.
  await page.getByTestId("new").click();
  await page.getByLabel("Name").fill("library test");
  await page.getByTestId("new-seq-text").fill(SEQUENCE);
  await page.getByTestId("new-seq-create").click();

  const dialog = await openLibrary(page);

  // Clean any leftover of this primer from an earlier run, so the add below is deterministic.
  const existing = dialog.getByTestId("pl-row").filter({ hasText: PRIMER });
  if (await existing.count()) {
    await existing.first().getByTestId("pl-remove").click();
    await page.getByTestId("pl-save").click();
    await expect(dialog.getByTestId("pl-notice")).toContainText("Saved");
  }

  await dialog.getByTestId("pl-add-name").fill("TestSeqP");
  await dialog.getByTestId("pl-add-seq").fill(PRIMER);
  await dialog.getByTestId("pl-add-location").fill("box 9 Z1");
  await dialog.getByTestId("pl-add-button").click();
  await expect(dialog.getByTestId("pl-table")).toContainText("TestSeqP");
  // Not saved yet, and already found on the open document.
  await expect(dialog.getByTestId("pl-bound")).toContainText("TestSeqP");
  await expect(dialog.getByTestId("pl-bound")).toContainText("forward");

  // A pasted spreadsheet block: one good row and one that is not DNA (left out). The import
  // panel renders alongside the main dialog body, so it is queried from the page.
  await dialog.getByTestId("pl-import").click();
  await page.getByTestId("pl-import-text").fill(`Imported1\t${IMPORTED}\nBadRow\tACGXYZ`);
  await page.getByTestId("pl-import-run").click();
  await expect(dialog.getByTestId("pl-notice")).toContainText("Imported 1");
  await expect(dialog.getByTestId("pl-notice")).toContainText("skipped");

  // Save writes both to the shared file; a reload reads them back.
  await page.getByTestId("pl-save").click();
  await expect(dialog.getByTestId("pl-notice")).toContainText("Saved");
  await dialog.getByTestId("pl-reload").click();
  await expect(dialog.getByTestId("pl-table")).toContainText("TestSeqP");
  await expect(dialog.getByTestId("pl-table")).toContainText("Imported1");

  // The search box filters to the primer that binds the open sequence.
  await dialog.getByTestId("pl-search").fill("TestSeqP");
  await expect(dialog.getByTestId("pl-count")).toContainText("1 of");

  await page.getByTestId("pl-close").click();

  // Show Library Primers draws it on the map as a primer.
  await menu(page, "Primers", "Show Library Primers");
  await menu(page, "Primers", "Primers View");
  await expect(page.locator(".veEditor").getByText("TestSeqP", { exact: false }).first()).toBeVisible();
});

test("showing library primers on an opened file leaves it unchanged and out of Undo", async ({ page }) => {
  // A 22-nt stretch of pBR322 itself, so the primer binds the file as opened.
  const embl = readFileSync(data("exchange-pBR322.embl"), "utf8");
  const sq = embl.slice(embl.indexOf("\nSQ")).split("\n").slice(1).join("").replace(/[^acgt]/gi, "").toUpperCase();
  const primer = sq.slice(1200, 1222);

  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("exchange-pBR322.embl"));
  await expect(page.locator(".veCircularView")).toBeVisible();
  const tab = page.locator(".vt-doctab.active");
  await expect(tab).not.toHaveAttribute("data-unsaved", "true");

  const dialog = await openLibrary(page);
  if (!(await dialog.getByTestId("pl-row").filter({ hasText: primer }).count())) {
    await dialog.getByTestId("pl-add-name").fill("pBR_1201F");
    await dialog.getByTestId("pl-add-seq").fill(primer);
    await dialog.getByTestId("pl-add-button").click();
    await page.getByTestId("pl-save").click();
    await expect(dialog.getByTestId("pl-notice")).toContainText("Saved");
  }
  await page.getByTestId("pl-close").click();

  await menu(page, "Primers", "Show Library Primers");
  await menu(page, "Primers", "Primers View");
  await expect(page.locator(".veEditor").getByText("pBR_1201F", { exact: false }).first()).toBeVisible();
  // A view, not an edit: the file still has no unsaved changes, and Undo has nothing to undo.
  await expect(tab).not.toHaveAttribute("data-unsaved", "true");
  await page.keyboard.press("Control+z");
  await expect(page.locator(".veEditor").getByText("pBR_1201F", { exact: false }).first()).toBeVisible();

  // Turning it off takes the primers away again, and the file is still unchanged.
  await menu(page, "Primers", "Show Library Primers");
  await expect(page.locator(".veEditor").getByText("pBR_1201F", { exact: false })).toHaveCount(0);
  await expect(tab).not.toHaveAttribute("data-unsaved", "true");
});
