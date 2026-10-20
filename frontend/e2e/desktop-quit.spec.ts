// Closing the desktop window with unsaved work asks in Vettore's own dialog. The window's side
// (src/vettore/desktop/close_guard.py) is tested in tests/test_desktop_close.py; here the page
// gets the same request the window sends and answers through the desktop bridge, stubbed with
// the native dialogs answered by the test. The engine writes real files.
import { expect, test, type Page } from "@playwright/test";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

let dir: string;

declare global {
  interface Window {
    __openPaths: string[];
    __nextSavePath: string | null;
    __saveCalls: [string, string, string][];
    __quitCalls: number;
  }
}

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

/** Fails the test if the browser's own confirm/alert boxes appear. */
function forbidNativeDialogs(page: Page) {
  const seen: string[] = [];
  page.on("dialog", async (dialog) => {
    seen.push(dialog.message());
    await dialog.dismiss();
  });
  return seen;
}

/** What the desktop window runs when it is closed (close_guard.quit_script): true = handled. */
const closeWindow = (page: Page) =>
  page.evaluate(() => !window.dispatchEvent(new CustomEvent("vettore:quit-requested", { cancelable: true })));

const quitCalls = (page: Page) => page.evaluate(() => window.__quitCalls);
const saveCalls = (page: Page) => page.evaluate(() => window.__saveCalls);

async function nextSaveAnswer(page: Page, file: string | null) {
  await page.evaluate((p) => { window.__nextSavePath = p; }, file);
}

/** Open a file the way the desktop Open dialog does: by its path. */
async function openByPath(page: Page, file: string) {
  await page.evaluate((p) => { window.__openPaths = [p]; }, file);
  const welcome = page.getByTestId("open");
  if (await welcome.isVisible().catch(() => false)) await welcome.click();
  else await menu(page, "File", "Open…");
  await expect(page.locator(".vt-doctab.active")).toContainText(path.parse(file).name);
}

async function newSequence(page: Page, name: string, sequence: string) {
  const welcome = page.getByTestId("new");
  if (await welcome.isVisible().catch(() => false)) await welcome.click();
  else await menu(page, "File", "New DNA File…");
  await page.getByLabel("Name").fill(name);
  await page.getByTestId("new-seq-text").fill(sequence);
  await page.getByTestId("new-seq-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText(name);
}

/** An edited GenBank file, a new sequence never saved, and an untouched sample. */
async function unsavedWork(page: Page): Promise<string> {
  const file = path.join(dir, "colleague.gb");
  fs.copyFileSync(data("MT612431.1.gb"), file);
  await openByPath(page, file);
  await menu(page, "Actions", "Reverse Complement Entire Sequence");
  await expect(page.locator(".vt-doctab.active")).toContainText("• colleague");
  // pBR322 (J01749.1) bases 1-60.
  await newSequence(page, "pTest", "TTCTCATGTTTGACAGCTTATCATCGATAAGCTTTAATGCGGTAGTTTATCACAGTTAAA");
  await menu(page, "File", "Open Sample Construct", "pUC19");
  await expect(page.getByTestId("doctab")).toHaveCount(3);
  return file;
}

test.beforeEach(async ({ page }) => {
  dir = fs.mkdtempSync(path.join(os.tmpdir(), "vettore-e2e-quit-"));
  await page.addInitScript(() => {
    window.__openPaths = [];
    window.__nextSavePath = null;
    window.__saveCalls = [];
    window.__quitCalls = 0;
    (window as any).pywebview = {
      api: {
        open_dialog: async () => window.__openPaths,
        save_dialog: async (name: string, directory: string, kind: string) => {
          window.__saveCalls.push([name, directory, kind]);
          return window.__nextSavePath;
        },
        confirm_quit: async () => { window.__quitCalls += 1; },
      },
    };
  });
  await page.goto(`/?token=${TOKEN}`);
});

test.afterEach(() => fs.rmSync(dir, { recursive: true, force: true }));

test("closing the window lists the unsaved documents in Vettore's own dialog; Cancel keeps them",
  async ({ page }) => {
    const native = forbidNativeDialogs(page);
    await unsavedWork(page);

    expect(await closeWindow(page)).toBe(true);    // the window waits for Vettore's answer
    const dialog = page.getByTestId("quit-dialog");
    await expect(dialog).toContainText("2 documents have changes that are not saved:");
    const listed = dialog.getByTestId("quit-docs").locator("li");
    await expect(listed).toHaveText([
      "colleague — saved to colleague.gb",
      "pTest — never saved to a file: you choose where",
    ]);
    await expect(page.getByTestId("quit-save-all")).toHaveText("Save All");
    await expect(page.getByTestId("quit-save-all")).toBeFocused();
    await page.screenshot({ path: "test-results/desktop-quit-dialog.png", animations: "disabled" });

    await page.keyboard.press("Escape");
    await expect(dialog).toBeHidden();
    expect(await closeWindow(page)).toBe(true);
    await page.getByTestId("quit-cancel").click();
    await expect(dialog).toBeHidden();
    await expect(page.getByTestId("doctab")).toHaveCount(3);
    expect(await quitCalls(page)).toBe(0);
    expect(native).toEqual([]);
  });

test("Save All saves each document the way Save does, then closes the window", async ({ page }) => {
  const file = await unsavedWork(page);
  const before = fs.readFileSync(file, "utf-8");
  await nextSaveAnswer(page, path.join(dir, "pTest-final"));   // asked for the new sequence only

  await closeWindow(page);
  await page.getByTestId("quit-save-all").click();
  await expect.poll(() => quitCalls(page)).toBe(1);
  await expect(page.getByTestId("quit-dialog")).toBeHidden();

  const after = fs.readFileSync(file, "utf-8");
  expect(after).toMatch(/^LOCUS/);
  expect(after).not.toEqual(before);                            // written back in place
  expect(fs.readFileSync(path.join(dir, "pTest-final.gb"), "utf-8")).toMatch(/^LOCUS +pTest/);
  expect(await saveCalls(page)).toEqual([["pTest.gb", "", "genbank"]]);
  await expect(page.locator('[data-testid="doctab"][data-unsaved="true"]')).toHaveCount(0);
});

test("Don't Save closes the window and leaves the files as they were", async ({ page }) => {
  const file = await unsavedWork(page);
  const before = fs.readFileSync(file, "utf-8");
  await closeWindow(page);
  await page.getByTestId("quit-discard").click();
  await expect.poll(() => quitCalls(page)).toBe(1);
  expect(fs.readFileSync(file, "utf-8")).toEqual(before);
  expect(await saveCalls(page)).toEqual([]);
});

test("a save cancelled during Save All keeps the window open and says why", async ({ page }) => {
  await newSequence(page, "pTest", "ATGGTGAGCAAGGGCGAGGAGCTGTTCACC");
  await nextSaveAnswer(page, null);                             // Cancel in the save dialog
  await closeWindow(page);
  const dialog = page.getByTestId("quit-dialog");
  await expect(dialog).toContainText("“pTest” has changes that are not saved.");
  await expect(dialog).toContainText("never saved to a file: you choose where");
  await page.getByTestId("quit-save-all").click();               // labelled Save for one document
  await expect(dialog).toBeHidden();
  await expect(page.locator('[data-testid="notice"][data-kind="warning"]'))
    .toContainText("Vettore is still open because “pTest” was not saved.");
  expect(await quitCalls(page)).toBe(0);
  await expect(page.locator(".vt-doctab.active")).toContainText("• pTest");
});

test("with nothing unsaved the window closes at once", async ({ page }) => {
  await page.getByTestId("sample-puc19").click();
  await expect(page.locator(".vt-doctab.active")).toContainText("pUC19");
  expect(await closeWindow(page)).toBe(true);
  await expect.poll(() => quitCalls(page)).toBe(1);
  await expect(page.getByTestId("quit-dialog")).toHaveCount(0);
});

test("a tab's own question gives way to the window's", async ({ page }) => {
  await newSequence(page, "pTest", "ATGGTGAGCAAGGGCGAGGAGCTGTTCACC");
  await page.keyboard.press("Control+w");
  await expect(page.getByTestId("unsaved-dialog")).toBeVisible();
  await closeWindow(page);
  await expect(page.getByTestId("unsaved-dialog")).toBeHidden();
  await expect(page.getByTestId("quit-dialog")).toBeVisible();
  await page.getByTestId("quit-cancel").click();
  await expect(page.getByTestId("doctab")).toHaveCount(1);
});
