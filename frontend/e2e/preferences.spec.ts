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

// ---- Open sequence files with Vettore (packaged application only) ------------------------

const EXTENSIONS = [".dna", ".gb", ".gbk", ".gbff", ".ape", ".embl", ".fa", ".fasta", ".xdna", ".gck", ".ab1"];
const packaged = (patch: object = {}) => ({
  available: true, enabled: false, current: false, registered_executable: null,
  extensions: EXTENSIONS, overridden: [], platform: "win32", ...patch,
});

/** The menu bar comes with the first open document. */
async function openPreferences(page: Page) {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".vt-doctab.active")).toContainText("pFA-KanMX4");
  await menu(page, "File", "Preferences…");
}

test("a copy run from source code does not offer to open sequence files", async ({ page }) => {
  const answered = page.waitForResponse("**/api/desktop/file-associations");
  await openPreferences(page);
  expect((await (await answered).json()).available).toBe(false);   // the real engine says so
  await expect(page.getByTestId("preferences-dialog")).toBeVisible();
  await expect(page.getByTestId("pref-open-with-field")).toHaveCount(0);
});

test("the packaged app offers to open sequence files on double-click, and to undo it", async ({ page }) => {
  const posted: boolean[] = [];
  let current = packaged();
  let answer = () => {};
  await page.route("**/api/desktop/file-associations", async (route) => {
    if (route.request().method() === "POST") {
      const { enabled } = route.request().postDataJSON();
      posted.push(enabled);
      await new Promise<void>((resolve) => { answer = resolve; });   // the test decides when
      current = packaged(enabled
        ? { enabled: true, current: true, registered_executable: "C:\\Vettore\\Vettore.exe", overridden: [".dna"] }
        : {});
    }
    await route.fulfill({ json: current });
  });
  await openPreferences(page);
  const box = page.getByTestId("pref-open-with");
  await expect(page.getByTestId("pref-open-with-field"))
    .toContainText("Open sequence files with Vettore when you double-click them");
  await expect(page.getByTestId("pref-open-with-field")).toContainText(
    "Applies to .dna, .gb, .gbk, .gbff, .ape, .embl, .fa, .fasta, .xdna, .gck and .ab1 files, "
    + "for your Windows account only; untick it to undo, and they open as they did before.");
  await expect(box).not.toBeChecked();

  await page.getByText("Open sequence files with Vettore when you double-click them").click();
  await expect(page.getByTestId("pref-open-with-busy")).toHaveText("Applying…");
  await expect(box).toBeDisabled();
  answer();
  await expect(page.getByTestId("pref-open-with-done"))
    .toHaveText("Sequence files now open with Vettore when you double-click them.");
  await expect(box).toBeChecked();
  // Windows keeps .dna for the program the user chose with "Always use this app": say how to change it.
  await expect(page.getByTestId("pref-open-with-overridden")).toContainText(
    "Windows was told to always open .dna files with another program, so it still opens there.");
  await page.screenshot({ path: "test-results/preferences-open-with.png" });

  // The box follows what the system did, so it changes only once the answer is in.
  await page.getByText("Open sequence files with Vettore when you double-click them").click();
  await expect(page.getByTestId("pref-open-with-busy")).toBeVisible();
  answer();
  await expect(page.getByTestId("pref-open-with-done")).toHaveText("Sequence files open as they did before.");
  await expect(box).not.toBeChecked();
  await expect(page.getByTestId("pref-open-with-overridden")).toHaveCount(0);
  expect(posted).toEqual([true, false]);
});

test("a refusal to change the file associations is explained under the option", async ({ page }) => {
  await page.route("**/api/desktop/file-associations", (route) => route.request().method() === "POST"
    ? route.fulfill({ status: 422, json: { detail: "Vettore could not change which program opens "
      + "sequence files: Access is denied. Nothing was changed." } })
    : route.fulfill({ json: packaged({ platform: "linux", registered_executable: "/opt/old/Vettore" }) }));
  await openPreferences(page);
  await expect(page.getByTestId("pref-open-with-field")).toContainText("for your user account only");
  await expect(page.getByTestId("pref-open-with-moved")).toHaveText(
    "Sequence files open with the Vettore at /opt/old/Vettore, not with this copy. Tick the box to use this copy instead.");
  await page.getByText("Open sequence files with Vettore when you double-click them").click();
  await expect(page.getByTestId("pref-open-with-error"))
    .toHaveText("Vettore could not change which program opens sequence files: Access is denied. Nothing was changed.");
  await expect(page.getByTestId("pref-open-with")).not.toBeChecked();
  await expect(page.getByTestId("pref-open-with")).toBeEnabled();
});
