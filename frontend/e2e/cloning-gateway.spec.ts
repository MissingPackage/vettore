import { expect, test, type Page } from "@playwright/test";
import path from "node:path";
import { TOKEN } from "../playwright.config";

// The three public records of tests/test_cloning_gateway.py (provenance in that file).
const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

test("LR reaction between an entry clone and a destination vector", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles([
    data("cloning-gateway-entry.gb"),
    data("cloning-gateway-destination.gb"),
  ]);
  await expect(page.locator(".vt-doctab")).toHaveCount(2);

  await menu(page, "Actions", "Gateway LR Reaction…");
  const dialog = page.getByTestId("cloning-dialog");
  await dialog.getByTestId("insert-doc").selectOption({ index: 1 });
  await dialog.getByTestId("vector-doc").selectOption({ index: 2 });

  const product = page.getByTestId("cloning-product");
  await expect(product).toContainText("expression clone", { timeout: 30_000 });
  await expect(page.getByTestId("cloning-products")).toContainText("Expression clone");
  await expect(page.getByTestId("cloning-products")).toContainText("attP by-product");
  // pAGRIKOLA has two Gateway cassettes, so one round leaves the other one with its ccdB.
  await expect(page.getByTestId("cloning-warnings")).toContainText("ccdB");
  await page.screenshot({ path: "test-results/cloning-gateway-lr.png" });

  await page.getByTestId("cloning-create").click();
  await page.getByTestId("view-history").click();
  await expect(page.locator(".vt-history")).toContainText("Gateway LR reaction");
  await expect(page.getByTestId("history-reopen")).toHaveCount(2);
});

test("a molecule without att sites is refused with a message that names it", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles([
    data("cloning-gateway-entry.gb"),
    data("pFA-KanMX4.dna"),
  ]);
  await expect(page.locator(".vt-doctab")).toHaveCount(2);

  await menu(page, "Actions", "Gateway LR Reaction…");
  const dialog = page.getByTestId("cloning-dialog");
  await dialog.getByTestId("insert-doc").selectOption({ index: 1 });
  await dialog.getByTestId("vector-doc").selectOption({ index: 2 });
  await expect(dialog.getByTestId("cloning-error")).toContainText("pFA-KanMX4", { timeout: 30_000 });
  await expect(dialog.getByTestId("cloning-error")).toContainText("att site");
  await expect(page.getByTestId("cloning-create")).toBeDisabled();
});
