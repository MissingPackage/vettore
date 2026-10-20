// Closing a tab never loses work silently: edited files, and documents that exist only in
// Vettore (products, pastes, new files), are marked unsaved and closing them asks first.
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

/** Fails the test if the browser's own confirm/alert boxes appear. */
function forbidNativeDialogs(page: Page) {
  const seen: string[] = [];
  page.on("dialog", async (dialog) => {
    seen.push(dialog.message());
    await dialog.dismiss();
  });
  return seen;
}

/** The count of unsaved documents last reported to the engine (it guards quitting with it). */
function unsavedReports(page: Page) {
  const counts: number[] = [];
  page.on("request", (request) => {
    if (request.url().endsWith("/api/unsaved")) counts.push(request.postDataJSON().count);
  });
  return counts;
}

async function newSequence(page: Page, name: string, sequence: string, circular = false) {
  const welcome = page.getByTestId("new");
  if (await welcome.isVisible().catch(() => false)) await welcome.click();
  else await menu(page, "File", "New DNA File…");
  await page.getByLabel("Name").fill(name);
  await page.getByTestId("new-seq-text").fill(sequence);
  if (circular) await page.locator("label.bp3-radio", { hasText: "Circular" }).click();
  await page.getByTestId("new-seq-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText(name);
}

const tab = (page: Page, name: string) => page.getByTestId("doctab").filter({ hasText: name });

test("closing an edited file asks in Vettore's own dialog: Escape keeps it, Don't save closes it",
  async ({ page }) => {
    const native = forbidNativeDialogs(page);
    await page.goto(`/?token=${TOKEN}`);
    await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
    await expect(tab(page, "pFA-KanMX4")).not.toHaveAttribute("data-unsaved", "true");

    await menu(page, "Actions", "Reverse Complement Entire Sequence");
    await expect(page.locator(".vt-doctab.active")).toContainText("• pFA-KanMX4");

    await page.keyboard.press("Control+w");
    const dialog = page.getByTestId("unsaved-dialog");
    await expect(dialog).toContainText("“pFA-KanMX4” has changes that are not saved.");
    await expect(page.getByTestId("unsaved-save")).toBeFocused();
    await page.screenshot({ path: "test-results/shell-unsaved-dialog.png", animations: "disabled" });
    await page.keyboard.press("Escape");
    await expect(dialog).toBeHidden();
    await expect(tab(page, "pFA-KanMX4")).toBeVisible();

    await page.getByTestId("doctab").locator(".vt-doctab-close").click();
    await page.getByTestId("unsaved-cancel").click();
    await expect(tab(page, "pFA-KanMX4")).toBeVisible();

    await page.keyboard.press("Control+w");
    await page.getByTestId("unsaved-discard").click();
    await expect(page.getByTestId("welcome")).toBeVisible();
    expect(native).toEqual([]);
  });

test("a new file counts as unsaved until saved, and Enter saves it before the tab closes",
  async ({ page }) => {
    const native = forbidNativeDialogs(page);
    const reported = unsavedReports(page);
    await page.goto(`/?token=${TOKEN}`);
    await newSequence(page, "pTest", "ATGGTGAGCAAGGGCGAGGAGCTGTTCACC");
    await expect(page.locator(".vt-doctab.active")).toContainText("• pTest");
    await expect.poll(() => reported.at(-1)).toBe(1);

    await page.keyboard.press("Control+w");
    await expect(page.getByTestId("unsaved-dialog")).toContainText("never been saved to a file");
    // In the browser, saving downloads a GenBank file.
    const download = page.waitForEvent("download");
    await page.keyboard.press("Enter");
    expect((await download).suggestedFilename()).toBe("pTest.gb");
    await expect(page.getByTestId("welcome")).toBeVisible();
    await expect.poll(() => reported.at(-1)).toBe(0);
    expect(native).toEqual([]);
  });

test("an empty new file and an untouched sample close without asking", async ({ page }) => {
  forbidNativeDialogs(page);
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("sample-puc19").click();
  await expect(page.locator(".vt-doctab.active")).toContainText("pUC19");
  await expect(page.locator(".vt-doctab.active")).not.toHaveAttribute("data-unsaved", "true");
  await page.keyboard.press("Control+w");
  await expect(page.getByTestId("welcome")).toBeVisible();

  await page.getByTestId("new").click();
  await page.getByTestId("new-seq-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText("New DNA");
  await page.keyboard.press("Control+w");
  await expect(page.getByTestId("unsaved-dialog")).toHaveCount(0);
  await expect(page.getByTestId("welcome")).toBeVisible();
});

test.describe("with clipboard access", () => {
  test.use({ permissions: ["clipboard-read", "clipboard-write"] });

  test("pasted sequences are unsaved and get distinct names; the quit guard counts them",
    async ({ page }) => {
      const reported = unsavedReports(page);
      await page.goto(`/?token=${TOKEN}`);
      await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
      await expect(page.locator(".vt-doctab.active")).toContainText("pFA-KanMX4");

      // pBR322 (J01749.1) bases 1-60, then 61-120: plain sequence text, no name of its own.
      const pastes = [
        "TTCTCATGTTTGACAGCTTATCATCGATAAGCTTTAATGCGGTAGTTTATCACAGTTAAA",
        "TTGCTAACGCAGTCAGGCACCGTGTATGAAATCTAACAATGCGCTCATCGTCATCCTCGG",
      ];
      for (const text of pastes) {
        await page.evaluate((t) => navigator.clipboard.writeText(t), text);
        await menu(page, "File", "New File from Clipboard");
      }
      await expect(page.getByTestId("doctab")).toHaveCount(3);
      const names = await page.getByTestId("doctab").locator("span").allTextContents();
      expect(names).toEqual(["pFA-KanMX4", "• Pasted sequence", "• Pasted sequence 2"]);
      await expect.poll(() => reported.at(-1)).toBe(2);
    });
});

test("a cloning product is unsaved: closing it asks instead of dropping it", async ({ page }) => {
  const native = forbidNativeDialogs(page);
  await page.goto(`/?token=${TOKEN}`);
  const backbone = "ATTGCCTACGGTACATGCAACGTTGACCATGGCATTACGGACCATTGACCAGTTACGGATCAGCTTAGGCATAC".repeat(14);
  const vector = backbone.slice(0, 500) + "GAATTC" + backbone.slice(500, 900) + "GGATCC" + backbone.slice(900, 1400);
  const insert = "CCGAATTCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCC"
    + "ACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCGGATCCCC";
  await newSequence(page, "pVector", vector, true);
  await newSequence(page, "Insert", insert);

  await menu(page, "Actions", "Restriction and Insertion Cloning", "Insert Fragment…");
  const dialog = page.getByTestId("cloning-dialog");
  await dialog.getByTestId("vector-doc").selectOption({ index: 1 });
  await dialog.getByTestId("vector-enzymes").fill("EcoRI");
  await dialog.getByTestId("vector-enzymes").fill("BamHI");
  await dialog.getByTestId("insert-doc").selectOption({ index: 2 });
  await dialog.getByTestId("insert-enzymes").fill("EcoRI");
  await dialog.getByTestId("insert-enzymes").fill("BamHI");
  await expect(page.getByTestId("cloning-product")).toContainText("pVector + Insert");
  await page.getByTestId("cloning-create").click();

  await expect(page.locator(".vt-doctab.active")).toContainText("• pVector + Insert");
  await expect(page.locator(".vt-doctab.active")).toHaveAttribute("title", /not saved to a file yet/);
  await page.keyboard.press("Control+w");
  await expect(page.getByTestId("unsaved-dialog")).toContainText("“pVector + Insert (restriction cloning)” has changes that are not saved.");
  await page.getByTestId("unsaved-cancel").click();
  await expect(page.locator(".vt-doctab.active")).toContainText("pVector + Insert");
  expect(native).toEqual([]);
});
