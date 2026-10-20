import { expect, test, type Page } from "@playwright/test";
import { TOKEN } from "../playwright.config";

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

// Same construction as tests/test_cloning_oligos.py: two BbsI sites facing a dropout, which
// leaves the backbone with the overhangs the gRNA protocol expects.
const FILLER = "ATTGCCTACGGTACATGCAACGTTGACCATAGCATTACGGACCATTGACCAGTTACGGATCAGCTTAGGCATAC";
const VECTOR = FILLER.repeat(8) + "CACC" + "TT" + "GTCTTC" + FILLER.slice(0, 60) + "GAAGAC" + "TT" + "GTTT";
const PROTOSPACER = "GGTCGCCCTCGAACTTCACCT";

test("design the two oligos that clone a gRNA into a cut vector", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("new").click();
  await page.getByLabel("Name").fill("pGuide");
  await page.getByTestId("new-seq-text").fill(VECTOR);
  await page.locator("label.bp3-radio", { hasText: "Circular" }).click();
  await page.getByTestId("new-seq-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText("pGuide");

  await menu(page, "Actions", "Anneal Oligos and Insert…");
  const dialog = page.getByTestId("cloning-dialog");
  await dialog.getByTestId("vector-doc").selectOption({ index: 1 });
  await dialog.getByTestId("vector-enzymes").fill("BbsI");
  await dialog.getByTestId("anneal-insert").fill(PROTOSPACER);

  const oligos = page.getByTestId("cloning-oligos");
  await expect(oligos).toBeVisible({ timeout: 20_000 });
  await expect(oligos).toContainText(`CACC${PROTOSPACER}`);
  await expect(page.getByTestId("cloning-warnings")).toContainText("T4 PNK");
  await page.screenshot({ path: "test-results/anneal-oligos.png" });

  await page.getByTestId("cloning-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText("pGuide + gRNA1_top");
});
