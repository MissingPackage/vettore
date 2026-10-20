// The menu bar: each command once, under one name, with one ellipsis character and shortcuts
// written as Windows and Linux users read them.
import { expect, test, type Page } from "@playwright/test";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

const MENUS = ["File", "Edit", "View", "Enzymes", "Features", "Primers", "Actions", "Tools", "Help"];

async function openMenu(page: Page, name: string) {
  await page.locator(".tg-menu-bar").getByText(name, { exact: true }).click();
  const menu = page.locator(".bp3-popover .bp3-menu").first();
  await expect(menu).toBeVisible();
  return menu;
}

async function closeMenu(page: Page, name: string) {
  await page.locator(".tg-menu-bar").getByText(name, { exact: true }).click();
  await expect(page.locator(".bp3-popover .bp3-menu")).toHaveCount(0);
}

/** The item labels of an open menu (top level only). */
const itemTexts = (menu: ReturnType<Page["locator"]>) =>
  menu.locator(":scope > li .bp3-menu-item > .bp3-text-overflow-ellipsis, :scope > .bp3-popover-wrapper .bp3-menu-item > .bp3-text-overflow-ellipsis")
    .allInnerTexts()
    .then((texts) => texts.map((t) => t.trim()).filter(Boolean));

test.beforeEach(async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".veCircularView")).toBeVisible();
});

test("Help holds the menu search and Vettore's own help, not the embedded editor's", async ({ page }) => {
  const help = await openMenu(page, "Help");
  await expect(help.locator("input")).toBeVisible();   // the menu search
  expect(await itemTexts(help)).toEqual(["Keyboard Shortcuts…", "Getting Started…", "About Vettore…"]);
  await expect(help).not.toContainText("OVE Version");
  await expect(help).not.toContainText("View Editor Hotkeys");
  await page.screenshot({ path: "test-results/shell-menu-help.png" });
});

test("every menu uses one ellipsis character and lists each command once", async ({ page }) => {
  const all: Record<string, string[]> = {};
  for (const name of MENUS) {
    const menu = await openMenu(page, name);
    all[name] = await itemTexts(menu);
    for (const text of all[name]) expect(text, `${name} → ${text}`).not.toContain("...");
    await closeMenu(page, name);
  }
  expect(all.Edit).toEqual(expect.arrayContaining(["Find…", "Go To…", "Select…"]));
  // Each of these lives in Actions only, under one name.
  for (const moved of ["Change Circular/Linear", "Reverse Complement Entire Sequence",
    "Reverse Complement Selection", "Rotate To Caret Position"]) {
    expect(all.Edit, `Edit → ${moved}`).not.toContain(moved);
  }
  expect(all.Actions).toEqual(expect.arrayContaining(["Reverse Complement Entire Sequence",
    "Set Origin at Cursor", "Change Topology"]));
  expect(all.Tools).not.toContain("Find…");
  expect(all.Tools).not.toContain("Go To…");
});

test("menu shortcuts read Ctrl+Shift+S, not ^ and ⇧ symbols", async ({ page }) => {
  const file = await openMenu(page, "File");
  const saveAs = file.locator(".bp3-menu-item", { hasText: "Save As…" });
  const keys = await saveAs.locator(".bp3-key-combo > .bp3-icon").evaluateAll((icons) => icons.map((icon) => ({
    written: getComputedStyle(icon, "::after").content,
    symbol: getComputedStyle(icon.querySelector("svg")!).display,
  })));
  expect(keys).toEqual([{ written: '"Ctrl+"', symbol: "none" }, { written: '"Shift+"', symbol: "none" }]);
  await expect(saveAs.locator(".bp3-key-combo")).toContainText("S");
  await page.screenshot({ path: "test-results/shell-menu-file.png" });
});

test("a sample construct opens from the File menu while another file is open", async ({ page }) => {
  await openMenu(page, "File");
  await page.locator(".bp3-popover .bp3-menu").getByText("Open Sample Construct", { exact: true }).click();
  await page.locator(".bp3-popover .bp3-menu").getByText("pUC19", { exact: true }).click();
  await expect(page.locator(".vt-doctab.active")).toContainText("pUC19");
  await expect(page.getByTestId("doctab")).toHaveCount(2);
  await expect(page.getByTestId("status-line")).toContainText("2,686 bp · circular");
});

test("the Primers view's Simulate PCR button opens Vettore's PCR dialog", async ({ page }) => {
  await page.getByTestId("file-input").setInputFiles(data("sample-hybridization-params.dna"));
  await expect(page.locator(".vt-doctab.active")).toContainText("sample-hybridization-params");
  await page.getByTestId("view-primers").click();
  await page.locator(".ve-propertiesPanel").getByRole("button", { name: "Simulate PCR" }).click();
  await expect(page.getByTestId("cloning-dialog")).toBeVisible();
  await expect(page.locator(".bp3-dialog-header")).toContainText("PCR");
  // The primer list is still what the Primers view shows.
  await page.keyboard.press("Escape");
  await expect(page.locator(".ve-propertiesPanel")).toContainText("XhoI-hht2(US)-Fwd");
});
