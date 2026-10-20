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

test("LR picks the entry clone and the destination vector by their att sites", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles([
    data("pFA-KanMX4.dna"),
    data("cloning-gateway-destination.gb"),
    data("cloning-gateway-entry.gb"),
  ]);
  await expect(page.locator(".vt-doctab")).toHaveCount(3);
  // The destination tab is active: it still goes to the destination role, not the entry one.
  await page.locator(".vt-doctab", { hasText: "destination" }).click();

  await menu(page, "Actions", "Gateway LR Reaction…");
  const dialog = page.getByTestId("cloning-dialog");
  await expect(dialog.getByTestId("insert-doc").locator("option:checked")).toContainText("entry");
  await expect(dialog.getByTestId("insert-doc").locator("option:checked")).toContainText("attL1, attL2");
  await expect(dialog.getByTestId("vector-doc").locator("option:checked")).toContainText("destination");
  await expect(page.getByTestId("cloning-product")).toContainText("expression clone", { timeout: 30_000 });
  await expect(dialog.getByTestId("cloning-error")).toHaveCount(0);
});

test("when no open file has the att sites, the dialog asks for them instead of guessing", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles([data("cloning-gateway-entry.gb"), data("pFA-KanMX4.dna")]);
  await expect(page.locator(".vt-doctab")).toHaveCount(2);
  await menu(page, "Actions", "Gateway LR Reaction…");
  const dialog = page.getByTestId("cloning-dialog");
  await expect(dialog.getByTestId("insert-doc").locator("option:checked")).toContainText("entry");
  await expect(dialog.getByTestId("vector-doc")).toHaveValue("");
  await expect(dialog.getByTestId("cloning-missing"))
    .toContainText("Choose the destination vector (attR sites). None of the open files has attR sites");
  await expect(dialog.getByTestId("cloning-error")).toHaveCount(0);
});

test("if the att sites cannot be read, both inputs are left for the user to choose", async ({ page }) => {
  await page.route("**/api/cloning/gateway-sites*", (route) => route.abort());
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles([data("cloning-gateway-entry.gb"), data("cloning-gateway-destination.gb")]);
  await expect(page.locator(".vt-doctab")).toHaveCount(2);
  await menu(page, "Actions", "Gateway LR Reaction…");
  const dialog = page.getByTestId("cloning-dialog");
  await expect(dialog.getByTestId("cloning-missing")).toHaveText("Choose the entry clone (attL sites)");
  await expect(dialog.getByTestId("insert-doc")).toHaveValue("");
  await expect(dialog.getByTestId("vector-doc")).toHaveValue("");
});

test("Golden Gate opened on unrelated files says what to choose, without an error", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles([data("MT612431.1.gb"), data("pFA-KanMX4.dna")]);
  await expect(page.locator(".vt-doctab")).toHaveCount(2);
  await menu(page, "Actions", "Golden Gate Assembly…");
  const dialog = page.getByTestId("cloning-dialog");
  await expect(dialog.getByTestId("cloning-missing")).toContainText("Choose the inputs on the left", { timeout: 30_000 });
  await expect(dialog.getByTestId("cloning-error")).toHaveCount(0);
  // Once the user changes an input, a failure is an error again.
  await dialog.locator("select").first().selectOption("BsmBI");
  await expect(dialog.getByTestId("cloning-error")).toBeVisible({ timeout: 30_000 });
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
