// Saving in the desktop window, with the native dialogs stubbed: the engine writes real files.
import { expect, test, type Page } from "@playwright/test";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);
const DNA_MAGIC = Buffer.from([0x09, 0x00, 0x00, 0x00, 0x0e]);

let dir: string;

/** A copy of a fixture in a fresh folder, so saving never touches tests/data. */
function copyOf(fixture: string, name = fixture): string {
  const target = path.join(dir, name);
  fs.copyFileSync(data(fixture), target);
  return target;
}

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

/** Open a file from the welcome screen the way the desktop Open dialog does: by its path. */
async function openByPath(page: Page, file: string) {
  await page.evaluate((p) => { window.__openPaths = [p]; }, file);
  await page.getByTestId("open").click();
}

/** The save dialog will answer `file` (or Cancel when null). */
async function nextSaveAnswer(page: Page, file: string | null) {
  await page.evaluate((p) => { window.__nextSavePath = p; }, file);
}

const saveCalls = (page: Page) => page.evaluate(() => window.__saveCalls);

declare global {
  interface Window {
    __openPaths: string[];
    __nextSavePath: string | null;
    __saveCalls: [string, string, string][];
    __documentDrops: number;
  }
}

test.beforeEach(async ({ page }) => {
  dir = fs.mkdtempSync(path.join(os.tmpdir(), "vettore-e2e-files-"));
  // The desktop bridge as pywebview injects it, with the native dialogs answered by the test.
  await page.addInitScript(() => {
    window.__openPaths = [];
    window.__nextSavePath = null;
    window.__saveCalls = [];
    window.__documentDrops = 0;
    (window as any).pywebview = {
      api: {
        open_dialog: async () => window.__openPaths,
        save_dialog: async (name: string, directory: string, kind: string) => {
          window.__saveCalls.push([name, directory, kind]);
          return window.__nextSavePath;
        },
      },
    };
    window.__vettoreNativeDrop = true;
    // pywebview listens for drops on the document: the page must let them through.
    document.addEventListener("drop", () => { window.__documentDrops += 1; });
  });
  await page.goto(`/?token=${TOKEN}`);
});

test.afterEach(() => fs.rmSync(dir, { recursive: true, force: true }));

test("a dropped file opens from its path, and Save writes it back in place", async ({ page }) => {
  const file = copyOf("MT612431.1.gb", "colleague.gb");
  const before = fs.readFileSync(file, "utf-8");
  const transfer = await page.evaluateHandle(() => {
    const dt = new DataTransfer();
    dt.items.add(new File(["placeholder"], "colleague.gb"));
    return dt;
  });
  await page.dispatchEvent(".vt-shell", "drop", { dataTransfer: transfer });
  // What src/vettore/app.py sends once pywebview has resolved the full path.
  await page.evaluate((p) => window.dispatchEvent(new CustomEvent("vettore:dropped-paths", {
    detail: { paths: [p], unresolved: [] },
  })), file);
  const active = page.locator(".vt-doctab.active");
  await expect(active).toContainText("colleague");
  await expect(active).toHaveAttribute("title", file);
  expect(await page.evaluate(() => window.__documentDrops)).toBe(1);

  await menu(page, "Actions", "Reverse Complement Entire Sequence");
  await expect(active).toContainText("•");
  await menu(page, "File", "Save");
  await expect(active).not.toContainText("•");
  expect(await saveCalls(page)).toEqual([]); // no dialog: written in place
  const after = fs.readFileSync(file, "utf-8");
  expect(after).toMatch(/^LOCUS/);
  expect(after).not.toEqual(before);
  await expect(page.getByTestId("notice").filter({ hasText: "Saved “colleague.gb”" })).toBeVisible();
});

test("a dropped file whose path is unknown is opened by its content", async ({ page }) => {
  const transfer = await page.evaluateHandle((b64) => {
    const dt = new DataTransfer();
    dt.items.add(new File([Uint8Array.from(atob(b64), (c) => c.charCodeAt(0))], "pBR322.embl"));
    return dt;
  }, fs.readFileSync(data("exchange-pBR322.embl")).toString("base64"));
  await page.dispatchEvent(".vt-shell", "drop", { dataTransfer: transfer });
  await page.evaluate(() => window.dispatchEvent(new CustomEvent("vettore:dropped-paths", {
    detail: { paths: [], unresolved: ["pBR322.embl"] },
  })));
  await expect(page.locator(".vt-doctab.active")).toContainText("pBR322");
  await expect(page.getByTestId("status-line")).toContainText("4361 bp · circular");
});

test("Save on a .dna file explains the GenBank copy once and leaves the original untouched", async ({ page }) => {
  const original = copyOf("pFA-KanMX4.dna");
  const originalBytes = fs.readFileSync(original);
  await openByPath(page, original);
  const active = page.locator(".vt-doctab.active");
  await expect(active).toContainText("pFA-KanMX4");
  await menu(page, "Actions", "Reverse Complement Entire Sequence");

  // Cancelled the first time: nothing is written.
  await menu(page, "File", "Save");
  const dialog = page.getByTestId("save-copy-dialog");
  await expect(dialog).toContainText("“pFA-KanMX4.dna” is a .dna file");
  await expect(dialog).toContainText("leaves the original file untouched");
  await page.screenshot({ path: "test-results/files-save-copy.png" });
  await page.getByRole("button", { name: "Cancel" }).click();
  await expect(dialog).toBeHidden();
  expect(await saveCalls(page)).toEqual([]);

  // The explanation is not repeated; the save dialog suggests a .gb next to the original.
  await nextSaveAnswer(page, path.join(dir, "edited")); // no extension typed
  await menu(page, "File", "Save");
  await expect(dialog).toBeHidden();
  await expect(active).toContainText("edited");
  expect(await saveCalls(page)).toEqual([["pFA-KanMX4.gb", dir, "genbank"]]);
  const copy = path.join(dir, "edited.gb");
  expect(fs.readFileSync(copy, "utf-8")).toMatch(/^LOCUS/);
  expect(fs.readFileSync(original).equals(originalBytes)).toBe(true);

  // From now on Save writes the copy in place.
  await menu(page, "Actions", "Reverse Complement Entire Sequence");
  await menu(page, "File", "Save");
  await expect(active).not.toContainText("•");
  expect(await saveCalls(page)).toHaveLength(1);
});

test("one sequence of a multi-record file is saved as a separate file", async ({ page }) => {
  await openByPath(page, copyOf("exchange-plasmids.fasta"));
  await expect(page.getByTestId("doctab")).toHaveCount(3);
  await menu(page, "File", "Save");
  await expect(page.getByTestId("save-copy-dialog"))
    .toContainText("is one of 3 sequences in “exchange-plasmids.fasta”; saving over that file would replace the other 2.");
  await page.getByTestId("save-copy-confirm").click();
  await expect.poll(() => saveCalls(page)).toEqual([["X06403.1.gb", dir, "genbank"]]);
});

test("Save As writes .dna or FASTA by the extension, and FASTA warns what it leaves out", async ({ page }) => {
  await openByPath(page, copyOf("pFA-KanMX4.dna"));
  await expect(page.locator(".vt-doctab.active")).toContainText("pFA-KanMX4");

  await nextSaveAnswer(page, path.join(dir, "for-a-colleague.dna"));
  await menu(page, "File", "Save As…");
  await expect(page.locator(".vt-doctab.active")).toContainText("for-a-colleague");
  expect(fs.readFileSync(path.join(dir, "for-a-colleague.dna")).subarray(0, 5).equals(DNA_MAGIC)).toBe(true);

  await nextSaveAnswer(page, path.join(dir, "sequence-only.fasta"));
  await menu(page, "File", "Save As…");
  // The record is named after the sequence (the name inside the file), not the file.
  expect(fs.readFileSync(path.join(dir, "sequence-only.fasta"), "utf-8")).toMatch(/^>pFA-KanMX4 /);
  await expect(page.locator('[data-testid="notice"][data-kind="warning"]'))
    .toContainText("was saved as FASTA, which keeps only the sequence: its 9 features are not in that file");
});

test("Export writes a .dna copy without changing the document's file", async ({ page }) => {
  const file = copyOf("MT612431.1.gb");
  await openByPath(page, file);
  await nextSaveAnswer(page, path.join(dir, "copy.dna"));
  await menu(page, "File", "Export", "DNA File (.dna)…");
  await expect(page.getByTestId("notice").filter({ hasText: "Exported “copy.dna” (.dna)" })).toBeVisible();
  expect(await saveCalls(page)).toEqual([["MT612431.1.dna", dir, "dna"]]);
  expect(fs.readFileSync(path.join(dir, "copy.dna")).subarray(0, 5).equals(DNA_MAGIC)).toBe(true);
  await expect(page.locator(".vt-doctab.active")).toHaveAttribute("title", file);
});

test("a save the system refuses is explained", async ({ page }) => {
  await openByPath(page, copyOf("MT612431.1.gb"));
  await nextSaveAnswer(page, path.join(dir, "no-such-folder", "x.gb"));
  await menu(page, "File", "Save As…");
  await expect(page.locator('[data-testid="notice"][data-kind="error"]'))
    .toContainText("“x.gb” could not be saved in");
  await expect(page.locator('[data-testid="notice"][data-kind="error"]')).toContainText("the folder does not exist any more.");
});
