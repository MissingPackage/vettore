// One way in, one way out: the editor's own import/export paths are not offered.
import { expect, test } from "@playwright/test";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

test("the toolbar has Save but none of the editor's own import, export or upload buttons", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".veCircularView")).toBeVisible();
  const toolbar = page.locator(".veToolbar");
  await expect(toolbar.locator('[data-test="saveTool"]')).toBeVisible();
  for (const tool of ["veDownloadTool", "veImportTool", "alignmentTool"]) {
    await expect(toolbar.locator(`[data-test="${tool}"]`)).toHaveCount(0);
  }
  await page.locator(".tg-menu-bar").getByText("File", { exact: true }).click();
  const fileMenu = page.locator(".bp3-popover .bp3-menu").first();
  await expect(fileMenu).toContainText("New File from Clipboard");
  await expect(fileMenu).not.toContainText("Import");
});
