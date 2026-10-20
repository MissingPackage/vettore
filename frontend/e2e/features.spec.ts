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

test("detect common features on a real AAV plasmid and add them to the map", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("PZ629150.1.gb"));
  await expect(page.locator(".veCircularView")).toBeVisible();

  await menu(page, "Features", "Detect Common Features…");
  const table = page.getByTestId("detect-table");
  await expect(table).toBeVisible({ timeout: 30_000 });
  for (const name of ["SaCas9", "U6 promoter", "AmpR", "AAV2 ITR (truncated)", "ori"]) {
    await expect(table).toContainText(name);
  }
  await expect(table).toContainText("PZ629150.1:581..3925");   // provenance is shown
  await page.screenshot({ path: "test-results/detect-features.png" });

  await page.getByTestId("detect-add").click();
  await expect(page.getByTestId("detect-dialog")).toBeHidden();

  // The features are now on the map and in the Features view.
  await page.getByTestId("view-features").click();
  const features = page.locator(".ve-propertiesPanel");
  await expect(features.getByText("SaCas9", { exact: true }).first()).toBeVisible();
  await expect(features.getByText("U6 promoter", { exact: true }).first()).toBeVisible();

  await page.getByTestId("view-map").click();
  await expect(page.locator(".veCircularView")).toContainText("SaCas9");
  await page.screenshot({ path: "test-results/detect-map.png" });
});
