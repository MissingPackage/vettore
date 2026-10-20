import { expect, test, type Page } from "@playwright/test";
import path from "node:path";
import { TOKEN } from "../playwright.config";

// Public records; provenance in tests/test_cloning_multisite.py and tests/test_cloning_gateway.py.
const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

test("a three-fragment LR reaction is set up from the open files and creates both products", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles([
    data("cloning-multisite-destination.gb"),
    data("cloning-multisite-3prime.gb"),
    data("cloning-gateway-entry.gb"),
    data("cloning-multisite-5prime.gb"),
  ]);
  await expect(page.locator(".vt-doctab")).toHaveCount(4);

  await menu(page, "Actions", "Gateway LR Reaction…");
  const dialog = page.getByTestId("cloning-dialog");
  // The entry clones are chained through the att sites, whatever the order of the tabs.
  await expect(dialog.getByTestId("insert-doc").locator("option:checked")).toContainText("5prime");
  await expect(dialog.getByTestId("insert-doc").locator("option:checked")).toContainText("attL4, attR1");
  await expect(dialog.getByTestId("entry-doc-2").locator("option:checked")).toContainText("gateway-entry");
  await expect(dialog.getByTestId("entry-doc-3").locator("option:checked")).toContainText("3prime");
  await expect(dialog.getByTestId("vector-doc").locator("option:checked")).toContainText("destination");

  await expect(page.getByTestId("cloning-product")).toContainText("11974 bp", { timeout: 30_000 });
  await expect(page.getByTestId("cloning-product")).toContainText("3 entry clones");
  await expect(page.getByTestId("cloning-products")).toContainText("Expression clone");
  await expect(page.getByTestId("cloning-products")).toContainText("attP by-product");
  await expect(dialog.getByTestId("cloning-error")).toHaveCount(0);
  await page.screenshot({ path: "test-results/cloning-multisite-lr.png" });

  // Create All opens the expression clone and the by-product, each in its own tab.
  await page.getByTestId("cloning-create-all").click();
  await expect(page.locator(".vt-doctab")).toHaveCount(6);
  await page.getByTestId("view-history").click();
  await expect(page.locator(".vt-history")).toContainText("Gateway LR reaction");
  await expect(page.getByTestId("history-reopen")).toHaveCount(4);   // the four parents
});

test("an entry clone can be added and removed, and the dialog says what is still missing", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles([
    data("cloning-gateway-entry.gb"),
    data("cloning-gateway-destination.gb"),
  ]);
  await expect(page.locator(".vt-doctab")).toHaveCount(2);
  await menu(page, "Actions", "Gateway LR Reaction…");
  const dialog = page.getByTestId("cloning-dialog");
  await expect(page.getByTestId("cloning-product")).toContainText("expression clone", { timeout: 30_000 });

  await dialog.getByTestId("add-entry").click();
  await expect(dialog.getByTestId("entry-doc-2")).toHaveValue("");
  await expect(dialog.getByTestId("cloning-missing"))
    .toContainText("Choose entry clone 2, or remove it from the reaction");
  await dialog.getByTestId("remove-entry-2").click();
  await expect(dialog.getByTestId("entry-doc-2")).toHaveCount(0);
  await expect(page.getByTestId("cloning-product")).toContainText("expression clone", { timeout: 30_000 });
});

test("a BP reaction shows the att pair of the donor and why mismatched sites cannot react", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles([
    data("cloning-gateway-expression.gb"),
    data("cloning-multisite-donor.gb"),
  ]);
  await expect(page.locator(".vt-doctab")).toHaveCount(2);
  await menu(page, "Actions", "Gateway BP Reaction…");
  const dialog = page.getByTestId("cloning-dialog");
  await expect(dialog.getByTestId("donor-pair")).toHaveText("attP1 and attP4");
  // attB1 pairs with attP1, but nothing pairs with the attB2 of the expression clone.
  await expect(dialog.getByTestId("cloning-missing")).toContainText("do not pair", { timeout: 30_000 });
  await expect(dialog.getByTestId("cloning-missing")).toContainText("attB1, attB2");
  await expect(dialog.getByTestId("cloning-missing")).toContainText("attP1, attP4");
  await expect(dialog.getByTestId("cloning-missing")).toContainText("has no attP2");
  await expect(page.getByTestId("cloning-create")).toBeDisabled();
});

test("a phiC31 attB vector is not offered as a Gateway input and says why", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles([
    data("cloning-phic31-attb.gb"),
    data("cloning-multisite-donor.gb"),
  ]);
  await expect(page.locator(".vt-doctab")).toHaveCount(2);
  await menu(page, "Actions", "Gateway BP Reaction…");
  const dialog = page.getByTestId("cloning-dialog");
  // Nothing is chosen for the attB side: the phiC31 vector does not qualify.
  await expect(dialog.getByTestId("insert-doc")).toHaveValue("", { timeout: 30_000 });
  await expect(dialog.getByTestId("insert-doc"))
    .toContainText("phiC31 attB — not a Gateway site");

  await dialog.getByTestId("insert-doc").selectOption({ index: 1 });
  await expect(dialog.getByTestId("cloning-missing"))
    .toContainText("carries a phiC31 attB site, which Gateway clonase does not recombine");
  await expect(page.getByTestId("cloning-create")).toBeDisabled();
});
