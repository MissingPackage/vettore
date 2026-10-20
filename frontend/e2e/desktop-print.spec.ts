// File → Print… opens the editor's print dialog; the browser (and the Windows window) print
// through the system's print preview. The Linux window drops print requests, and says so
// (src/vettore/desktop/printing.py; the real window is checked in tests/test_desktop_qt.py).
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

test.beforeEach(async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".vt-doctab.active")).toContainText("pFA-KanMX4");
});

test("File → Print… opens the print dialog with the map", async ({ page }) => {
  await menu(page, "File", "Print…");
  const dialog = page.locator(".bp3-dialog").filter({ hasText: "Print" });
  await expect(dialog.getByRole("button", { name: "Print", exact: true })).toBeVisible();
  await expect(dialog.locator(".veCircularView, .veLinearView").first()).toBeVisible();
});

test("a print request the Linux window drops is explained, with what to do instead", async ({ page }) => {
  // What the Linux window sends when Qt WebEngine reports a print request it will not serve.
  await page.evaluate(() => window.dispatchEvent(new CustomEvent("vettore:print-unavailable")));
  const notice = page.locator('[data-testid="notice"][data-kind="warning"]');
  await expect(notice).toContainText("Vettore's window cannot print on Linux yet, so nothing was sent to the printer.");
  await expect(notice).toContainText("File → Export → Map as PNG…");
});
