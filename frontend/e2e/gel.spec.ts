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

test("simulate an agarose gel of a digest and of the uncut plasmid", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".veCircularView")).toBeVisible();

  await menu(page, "Tools", "Simulate Agarose Gel…");
  const dialog = page.getByTestId("gel-dialog");
  await expect(dialog).toBeVisible();
  await dialog.getByTestId("gel-doc-0").selectOption({ index: 1 });
  await expect(page.getByTestId("gel-image")).toBeVisible();
  // Uncut circular plasmid: one supercoiled band.
  await expect(dialog).toContainText("pFA-KanMX4 / uncut");

  await dialog.getByTestId("gel-add-lane").click();
  await dialog.getByTestId("gel-doc-1").selectOption({ index: 1 });
  await dialog.getByTestId("gel-enzymes-1").fill("NdeI");
  await dialog.getByTestId("gel-enzymes-1").fill("NcoI");
  await expect(dialog).toContainText("pFA-KanMX4 / NdeI + NcoI");
  await page.screenshot({ path: "test-results/gel.png" });

  const bands = page.getByTestId("gel-image").locator("rect");
  expect(await bands.count()).toBeGreaterThan(20);
});
