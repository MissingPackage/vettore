// S6.1: the Lineage tab, and the History travelling inside the saved file.
import { expect, test, type Page } from "@playwright/test";
import fs from "node:fs";
import path from "node:path";
import { TOKEN } from "../playwright.config";

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

async function newSequence(page: Page, name: string, sequence: string, circular: boolean) {
  const welcome = page.getByTestId("new");
  if (await welcome.isVisible().catch(() => false)) await welcome.click();
  else await menu(page, "File", "New DNA File…");
  await page.getByLabel("Name").fill(name);
  await page.getByTestId("new-seq-text").fill(sequence);
  if (circular) {
    await page.locator("label.bp3-radio", { hasText: "Circular" }).click();
    await expect(page.getByTestId("topology-circular")).toBeChecked();
  }
  await page.getByTestId("new-seq-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText(name);
}

const BACKBONE = "ATTGCCTACGGTACATGCAACGTTGACCATGGCATTACGGACCATTGACCAGTTACGGATCAGCTTAGGCATAC".repeat(14);
const VECTOR = BACKBONE.slice(0, 500) + "GAATTC" + BACKBONE.slice(500, 900) + "GGATCC"
  + BACKBONE.slice(900, 1400);
const INSERT = "CC" + "GAATTC"
  + "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCC"
  + "GGATCC" + "CC";

/** Clone an insert into a vector, both made here, and land on the product's tab. */
async function cloneIntoVector(page: Page) {
  await page.goto(`/?token=${TOKEN}`);
  await newSequence(page, "pVector", VECTOR, true);
  await newSequence(page, "Insert", INSERT, false);
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
  await expect(page.locator(".vt-doctab.active")).toContainText("pVector + Insert");
}

test("the Lineage tab shows what the construct was made from", async ({ page }) => {
  await cloneIntoVector(page);

  await page.getByTestId("view-lineage").click();
  const view = page.getByTestId("lineage-view");
  await expect(view).toContainText("Lineage of “pVector + Insert");
  // Each input says which step was run on it, whatever it came out of itself.
  await expect(view).toContainText("used by Restriction cloning");
  // The document, then the two inputs the step ran on.
  await expect(view.getByTestId("lineage-node")).toHaveCount(3);
  await expect(view).toContainText("pVector");
  await expect(view).toContainText("Insert");
  await expect(view.getByTestId("lineage-open-snapshot")).toHaveCount(2);
  await page.screenshot({ path: "test-results/view-lineage.png" });

  // A snapshot opens as its own document, exactly as the input was when the step ran.
  await view.getByTestId("lineage-open-snapshot").first().click();
  await expect(page.locator(".vt-doctab.active")).toContainText("pVector (from the Lineage)");
});

test("a document with no inputs says so instead of showing a tree", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input")
    .setInputFiles(path.resolve(import.meta.dirname, "../../tests/data", "pFA-KanMX4.dna"));
  await expect(page.locator(".vt-doctab.active")).toContainText("pFA-KanMX4");
  await page.getByTestId("view-lineage").click();
  await expect(page.getByTestId("lineage-view"))
    .toContainText("was not made from other documents");
  await expect(page.getByTestId("lineage-node")).toHaveCount(1);
});

test("the saved file carries the History, and the preference leaves it out", async ({ page }) => {
  await cloneIntoVector(page);

  // Saved in the browser as a download: the file is what another computer would receive.
  const download = page.waitForEvent("download");
  await menu(page, "File", "Save");
  const file = await download;
  const saved = path.join(path.dirname(await file.path()), "pVector-Insert.gb");
  fs.copyFileSync(await file.path(), saved);
  expect(fs.readFileSync(saved, "utf-8")).toContain("##vettore-lineage-v1");

  // Opened again: the step and the snapshots of its inputs are back.
  await page.getByTestId("file-input").setInputFiles(saved);
  await expect(page.locator(".vt-doctab.active")).toContainText("pVector + Insert");
  await page.getByTestId("view-lineage").click();
  await expect(page.getByTestId("lineage-view")).toContainText("used by Restriction cloning");
  await expect(page.getByTestId("lineage-node")).toHaveCount(3);
  await page.getByTestId("view-history").click();
  await expect(page.locator(".vt-history")).toContainText("Restriction cloning");
  await expect(page.getByTestId("history-reopen")).toHaveCount(2);

  // Turned off, a saved file carries the sequence and its features only.
  await menu(page, "File", "Preferences…");
  const dialog = page.getByTestId("preferences-dialog");
  // Blueprint covers the checkbox with its own indicator: the label is what is clicked.
  await page.getByText("Keep the History and the Lineage inside the files you save").click();
  await expect(dialog.getByTestId("pref-save-history")).not.toBeChecked();
  await expect(dialog.getByTestId("pref-saved")).toContainText("leave the History out");
  await page.getByTestId("pref-close").click();   // the footer sits outside the body
  const plain = page.waitForEvent("download");
  await menu(page, "File", "Save");
  expect(fs.readFileSync(await (await plain).path(), "utf-8")).not.toContain("vettore-lineage");
});
