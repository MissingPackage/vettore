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

test("new DNA file, edit, unsaved marker, save as GenBank download", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("new").click();
  await page.getByLabel("Name").fill("test-insert");
  await page.getByTestId("new-seq-text").fill(">x\nATGGCTAGCAAGGAGGAGCTGTTCACCGGG\n123 TAA");
  await page.getByTestId("new-seq-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText("test-insert");
  await expect(page.getByTestId("status-line")).toContainText("33 bp · linear");
  // A new file exists only in Vettore until it is saved: unsaved from the start.
  await expect(page.locator(".vt-doctab.active")).toContainText("• test-insert");

  await menu(page, "Actions", "Reverse Complement Entire Sequence");
  await expect(page.locator(".vt-doctab.active")).toContainText("• test-insert");

  const download = page.waitForEvent("download");
  await menu(page, "File", "Save");
  const file = await download;
  expect(file.suggestedFilename()).toBe("test-insert.gb");
  const text = fs.readFileSync(await file.path(), "utf-8");
  expect(text).toMatch(/^LOCUS/);
  // Reverse complement of ATGGCTAGCAAGGAGGAGCTGTTCACCGGGTAA, in GenBank 10-base blocks.
  expect(text).toContain("1 ttacccggtg aacagctcct ccttgctagc cat");
  await expect(page.locator(".vt-doctab.active")).not.toContainText("•");
});

test("export the map as SVG", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".veCircularView")).toBeVisible();
  const download = page.waitForEvent("download");
  await menu(page, "File", "Export", "Map as SVG…");
  const file = await download;
  expect(file.suggestedFilename()).toBe("pFA-KanMX4-map.svg");
  const svg = fs.readFileSync(await file.path(), "utf-8");
  expect(svg).toContain("<svg");
  expect(svg).toContain("kanMX");
});

test("selecting a primer shows its Tm with the method", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("sample-hybridization-params.dna"));
  await page.getByTestId("view-primers").click();
  await page.locator(".ve-propertiesPanel").getByText("XhoI-hht2(US)-Fwd").click();
  await expect(page.getByTestId("status-line")).toContainText("19 bp");
  const tm = page.getByTestId("selection-tm");
  await expect(tm).toContainText(/Tm \d+\.\d °C/);
  await expect(tm).toHaveAttribute("title", /NEB Q5/);
});
