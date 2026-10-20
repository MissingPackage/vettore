// Exchanging sequences by copy and paste, as colleagues do in emails and chats.
import { expect, test, type Page } from "@playwright/test";
import fs from "node:fs";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

const readClipboard = (page: Page) => page.evaluate(() => navigator.clipboard.readText());

// Two records from NCBI's FASTA of pBR322 (J01749.1), bases 1-60 and 61-120.
const TWO_RECORDS = [
  ">part1 pBR322 bases 1-60",
  "TTCTCATGTTTGACAGCTTATCATCGATAAGCTTTAATGCGGTAGTTTATCACAGTTAAA",
  ">part2 pBR322 bases 61-120",
  "TTGCTAACGCAGTCAGGCACCGTGTATGAAATCTAACAATGCGCTCATCGTCATCCTCGG",
  "",
].join("\n");

test.describe("with clipboard access", () => {
  test.use({ permissions: ["clipboard-read", "clipboard-write"] });

  test.beforeEach(async ({ page }) => {
    await page.goto(`/?token=${TOKEN}`);
    await page.getByTestId("file-input").setInputFiles(data("sample-hybridization-params.dna"));
    await expect(page.locator(".vt-doctab.active")).toContainText("sample-hybridization-params");
  });

  test("New File from Clipboard opens what was copied", async ({ page }) => {
    await page.evaluate((text) => navigator.clipboard.writeText(text), TWO_RECORDS);
    await menu(page, "File", "New File from Clipboard");
    await expect(page.getByTestId("doctab").filter({ hasText: "part1" })).toBeVisible();
    await expect(page.locator(".vt-doctab.active")).toContainText("part2");
    await expect(page.getByTestId("notice").filter({ hasText: "2 sequences opened from the clipboard" })).toBeVisible();
    await expect(page.getByTestId("status-line")).toContainText("60 bp · linear");
  });

  test("Copy as FASTA copies the whole sequence, or the selection with its positions", async ({ page }) => {
    await menu(page, "Edit", "Copy as FASTA");
    await expect(page.getByTestId("notice").filter({ hasText: "Copied the whole sequence" })).toBeVisible();
    const whole = await readClipboard(page);
    expect(whole.startsWith(">sample-hybridization-params\n")).toBe(true);
    expect(whole.split("\n").slice(1).join("")).toHaveLength(2414);

    // The forward primer's binding site: bases 498-516 of the template.
    await page.getByTestId("view-primers").click();
    await page.locator(".ve-propertiesPanel").getByText("XhoI-hht2(US)-Fwd").click();
    await expect(page.getByTestId("status-line")).toContainText("19 bp");
    await menu(page, "Edit", "Copy as FASTA");
    expect((await readClipboard(page)).toUpperCase()).toBe(">SAMPLE-HYBRIDIZATION-PARAMS:498-516\nCTCGAGGAAAAGCTTCAAC\n");
  });

  test("the right-click Copy menu offers FASTA instead of the editor's own GenBank copy", async ({ page }) => {
    await page.getByTestId("view-sequence").click();
    const rows = page.locator(".veRowView").first();
    await rows.click({ position: { x: 200, y: 40 } });
    await page.keyboard.press("Control+a");
    await rows.click({ button: "right", position: { x: 300, y: 40 } });
    await page.locator(".bp3-menu-item").filter({ hasText: /^Copy/ }).first().hover();
    const submenu = page.locator(".bp3-menu-item:visible");
    await expect(submenu.filter({ hasText: "Copy as FASTA" })).toBeVisible();
    await expect(submenu.filter({ hasText: /Genbank/i })).toHaveCount(0);
    await submenu.filter({ hasText: "Copy as FASTA" }).click();
    expect(await readClipboard(page)).toMatch(/^>sample-hybridization-params:1-2414\n/);
  });
});

test("pasting on the start screen opens what was copied", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await expect(page.getByTestId("open")).toBeVisible();
  await page.evaluate((text) => {
    const transfer = new DataTransfer();
    transfer.setData("text/plain", text);
    document.body.dispatchEvent(new ClipboardEvent("paste", { clipboardData: transfer, bubbles: true }));
  }, TWO_RECORDS);
  await expect(page.getByTestId("doctab")).toHaveCount(2);
  await expect(page.locator(".vt-doctab.active")).toContainText("part2");
  await page.getByTestId("view-history").click();
  await expect(page.locator(".vt-history")).toContainText("Created from clipboard");
});

test.describe("when the system refuses the clipboard", () => {
  test.beforeEach(async ({ page }) => {
    await page.addInitScript(() => {
      Object.defineProperty(navigator, "clipboard", {
        value: {
          readText: () => Promise.reject(new DOMException("Read permission denied.", "NotAllowedError")),
          writeText: () => Promise.reject(new DOMException("Write permission denied.", "NotAllowedError")),
        },
      });
      document.execCommand = () => false;
    });
    await page.goto(`/?token=${TOKEN}`);
    await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
    await expect(page.locator(".vt-doctab.active")).toContainText("pFA-KanMX4");
  });

  test("New File from Clipboard opens a paste box that checks the text as you type", async ({ page }) => {
    await menu(page, "File", "New File from Clipboard");
    const dialog = page.getByTestId("paste-dialog");
    await expect(dialog).toContainText("could not read the clipboard directly");
    const open = page.getByTestId("paste-open");
    await expect(open).toBeDisabled();

    await page.getByTestId("paste-text").fill("TTCTCATGTT TGACAGCTTA\nTCATCGXTAA");
    await expect(page.getByTestId("paste-error")).toContainText("line 2, column 7 has “X”, which is not a DNA letter");
    await expect(open).toBeDisabled();
    await page.screenshot({ path: "test-results/files-paste-error.png" });

    await page.getByTestId("paste-text").fill(fs.readFileSync(data("MT612431.1.gb"), "utf-8"));
    await expect(page.getByTestId("paste-status")).toContainText("Recognised as GenBank: 7,204 bp, circular");
    await open.click();
    await expect(dialog).toBeHidden();
    await expect(page.locator(".vt-doctab.active")).toContainText("MT612431");
    await expect(page.getByTestId("status-line")).toContainText("7,204 bp · circular");
  });

  test("a plain sequence pasted in the box takes the name and topology chosen", async ({ page }) => {
    await menu(page, "File", "New File from Clipboard");
    await page.getByTestId("paste-text").fill("1 ttctcatgtt tgacagctta tcatcgataa\n31 gctttaatgc");
    await page.getByTestId("paste-name").fill("pBR322 start");
    await page.getByTestId("paste-circular").check({ force: true });
    await expect(page.getByTestId("paste-status")).toContainText("Recognised as plain sequence text: 40 bp");
    await page.getByTestId("paste-open").click();
    await expect(page.locator(".vt-doctab.active")).toContainText("pBR322 start");
    await expect(page.getByTestId("status-line")).toContainText("40 bp · circular");
  });

  test("Ctrl+Shift+V opens the paste box without touching the open sequence", async ({ page }) => {
    await page.locator(".veCircularView").click({ position: { x: 20, y: 20 } });
    await page.keyboard.press("Control+Shift+V");
    await expect(page.getByTestId("paste-dialog")).toBeVisible();
    await expect(page.locator(".vt-doctab.active")).not.toContainText("•");
  });

  test("Copy as FASTA shows the text to copy by hand", async ({ page }) => {
    await menu(page, "Edit", "Copy as FASTA");
    const dialog = page.getByTestId("copy-text-dialog");
    await expect(dialog).toContainText("did not let Vettore put this on the clipboard");
    await expect(page.getByTestId("copy-text")).toHaveValue(/^>pFA-KanMX4\nGAACGCGGCCGCCAGCTGAAGCTTCGTACGCTGCAGGTCGACGGATCCCCGGGTTAATTA\n/);
  });
});
