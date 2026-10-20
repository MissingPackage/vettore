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

test("preferences change the defaults, apply at once and survive a restart", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".veCircularView")).toBeVisible();
  await expect(page.locator(".vt-status-set")).toHaveText("Unique 6+ Cutters");

  await menu(page, "File", "Preferences…");
  const dialog = page.getByTestId("preferences-dialog");
  await expect(dialog).toBeVisible();
  await dialog.getByTestId("pref-enzyme-set").selectOption("uniqueDual");
  await expect(dialog.getByTestId("pref-saved")).toBeVisible();
  await dialog.getByTestId("pref-gel-ladder").selectOption("generuler-1kb");
  await dialog.getByTestId("pref-identity").selectOption("90");
  await page.screenshot({ path: "test-results/preferences.png" });
  await page.getByTestId("pref-close").click();   // the footer sits outside the body

  // A document opened now starts from the new enzyme set.
  await page.getByTestId("file-input").setInputFiles(data("sample-d.dna"));
  await expect(page.locator(".vt-status-set")).toHaveText("Unique & Dual Cutters");

  // The gel dialog opens on the ladder chosen in preferences.
  await menu(page, "Tools", "Simulate Agarose Gel…");
  await expect(page.getByTestId("gel-ladder")).toHaveValue("generuler-1kb");
  await page.keyboard.press("Escape");

  // And the choices are still there after the window is reopened.
  await page.reload();
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".vt-status-set")).toHaveText("Unique & Dual Cutters");
  await menu(page, "File", "Preferences…");
  await expect(page.getByTestId("pref-identity")).toHaveValue("90");

  // Restore defaults so the next test starts from a known state.
  await page.getByTestId("pref-reset").click();
  await expect(page.getByTestId("pref-enzyme-set")).toHaveValue("unique6");
});

test("when the ladder list cannot be read, Preferences says so under the field", async ({ page }) => {
  await page.route("**/api/gel/ladders", (route) =>
    route.fulfill({ status: 500, contentType: "text/plain", body: "Internal Server Error" }));
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".veCircularView")).toBeVisible();
  await menu(page, "File", "Preferences…");
  const error = page.getByTestId("pref-ladder-error");
  await expect(error).toContainText("The list of ladders could not be read");
  await expect(error).toContainText("Something went wrong inside Vettore (error 500)");
  await expect(error).toContainText("reopen Preferences");
});
