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

test("explains what to do when only one sequence is open", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("analysis-AY530579.1.gb"));
  await expect(page.getByTestId("doctab")).toHaveCount(1);
  await menu(page, "Tools", "Align Multiple Sequences…");
  await expect(page.getByTestId("msa-empty")).toContainText("Open at least two sequences");
});
