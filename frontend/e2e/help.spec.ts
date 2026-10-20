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
  await expect(page.locator(".veCircularView")).toBeVisible();
});

test("the shortcut list is in the Help menu, and the shortcuts it lists work",
  async ({ page }) => {
    await menu(page, "Help", "Keyboard Shortcuts…");
    const dialog = page.getByTestId("shortcuts-dialog");
    await expect(dialog).toBeVisible();
    for (const keys of ["Ctrl+O", "Ctrl+2", "Ctrl+F", "Ctrl+Z", "Ctrl+Shift+D"]) {
      await expect(dialog.getByText(keys, { exact: true })).toBeVisible();
    }
    await expect(dialog.getByText("Ctrl+Shift+G", { exact: true })).toBeVisible();
    await page.screenshot({ path: "test-results/shortcuts.png" });
    await page.keyboard.press("Escape");

    await page.keyboard.press("Control+Shift+G");
    await expect(page.getByTestId("gel-dialog")).toBeVisible();
    await page.keyboard.press("Escape");
    await expect(page.getByTestId("gel-dialog")).toBeHidden();

    // Spot-check three of the listed shortcuts in the running app.
    await page.keyboard.press("Control+2");
    await expect(page.locator(".veRowView")).toBeVisible();
    await page.keyboard.press("Control+f");
    await expect(page.locator(".tg-find-tool-input").first()).toBeVisible();
    await page.keyboard.press("Escape");
    await page.keyboard.press("Control+1");
    await expect(page.locator(".veCircularView")).toBeVisible();
  });

test("Getting Started walks through the day's work", async ({ page }) => {
  await menu(page, "Help", "Getting Started…");
  const dialog = page.getByTestId("getting-started-dialog");
  await expect(dialog).toContainText("Open a construct");
  await expect(dialog).toContainText("Read the map");
  await expect(dialog).toContainText("Simulate the cloning");
  await expect(dialog).toContainText("Check the clone");
  await expect(dialog).toContainText("sends nothing anywhere");
  // The tools of this release are all named.
  for (const tool of ["Insert Multiple Fragments", "Gateway", "TOPO", "Align Multiple Sequences",
    "AAV Packaging Check", "New File from Clipboard", "Open Sample Construct", "Design PCR Primers"]) {
    await expect(dialog).toContainText(tool);
  }
  await page.screenshot({ path: "test-results/getting-started.png" });

  // Every "Menu → Item…" the guide names is in that menu, under that name.
  const text = await dialog.innerText();
  const paths = [...text.matchAll(/\b(File|Edit|View|Enzymes|Features|Primers|Actions|Tools|Help) → ([^.,;:()→]+?…)/g)]
    .map((m) => [m[1], m[2].replace(/\s+/g, " ").trim()]);
  expect(paths.length).toBeGreaterThanOrEqual(7);
  await page.keyboard.press("Escape");
  for (const [menuName, item] of paths) {
    await page.locator(".tg-menu-bar").getByText(menuName, { exact: true }).click();
    await expect(page.locator(".bp3-popover .bp3-menu").first().getByText(item, { exact: true }),
      `${menuName} → ${item}`).toBeVisible();
    await page.locator(".tg-menu-bar").getByText(menuName, { exact: true }).click();
  }
});

test("About names the version, the licence and the log file", async ({ page }) => {
  await menu(page, "Help", "About Vettore…");
  const dialog = page.getByTestId("about-dialog");
  await expect(dialog.getByTestId("about-version")).toHaveText(/\d+\.\d+/);
  await expect(dialog).toContainText("GPL-3.0-or-later");
  await expect(dialog).toContainText(".vettore");
  await expect(dialog).toContainText("vettore.log");
  await expect(dialog).toContainText("primer3-py");
  await page.screenshot({ path: "test-results/about.png" });
});
