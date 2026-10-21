// File → Open from NCBI… (S6.3): a GenBank record by its accession, the plain message when the
// computer is offline, and the preference that hides the menu item. NCBI is replaced by the
// engine reading the same record from the test data, so the suite runs offline.
import { expect, test, type Page, type Route } from "@playwright/test";
import { readFileSync } from "node:fs";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

test.use({ viewport: { width: 1366, height: 768 } });

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

/** A document open, so the editor's menu bar is on screen. */
async function openFile(page: Page) {
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".veCircularView")).toBeVisible();
}

/** What the engine answers after downloading MT612431.1: the same record, opened from disk. */
async function answerWithRecord(page: Page, route: Route) {
  const opened = await page.request.post("/api/open-content", {
    headers: { "X-Vettore-Token": TOKEN },
    data: { filename: "MT612431.1.gb", content_b64: readFileSync(data("MT612431.1.gb")).toString("base64") },
  });
  await route.fulfill({ json: { path: null, accession: "MT612431.1", ...(await opened.json()) } });
}

test("an accession opens as an unsaved document, and only the accession is sent", async ({ page }) => {
  const sent: unknown[] = [];
  await page.route("**/api/ncbi/fetch", async (route) => {
    sent.push(route.request().postDataJSON());
    await answerWithRecord(page, route);
  });
  await page.goto(`/?token=${TOKEN}`);

  // The welcome screen offers it too: the menu bar comes with the first document.
  await page.getByTestId("welcome-ncbi").click();
  const dialog = page.getByTestId("ncbi-dialog");
  await expect(dialog).toBeVisible();
  await expect(dialog.getByTestId("ncbi-privacy")).toContainText("Only the accession number is sent");
  await expect(dialog.getByTestId("ncbi-privacy")).toContainText("No account is needed");
  await expect(page.getByTestId("ncbi-open")).toBeDisabled();
  await page.screenshot({ path: "test-results/ncbi-dialog-1366.png" });

  await page.getByTestId("ncbi-accession").fill("MT612431.1");
  await page.keyboard.press("Enter");
  await expect(dialog).toBeHidden();
  expect(sent).toEqual([{ accession: "MT612431.1" }]);

  await expect(page.locator(".vt-doctab.active")).toContainText("MT612431.1");
  await expect(page.getByTestId("status-line")).toContainText("7,204 bp · circular");
  // Downloaded, not on disk: closing it would lose it, so it counts as unsaved.
  await expect(page.locator(".vt-doctab.active")).toContainText("•");
  await menu(page, "View", "History");
  await expect(page.getByTestId("history-view")).toContainText("Opened from NCBI");
  await expect(page.getByTestId("history-view")).toContainText("MT612431.1, 7204 bp");
});

test("offline, the dialog says so and stays open to try again", async ({ page }) => {
  await page.route("**/api/ncbi/fetch", (route) => route.fulfill({
    status: 503,
    json: { detail: "Vettore could not reach NCBI: this computer seems to be offline, or the network "
      + "does not let Vettore through. Check the internet connection and try again; everything else "
      + "in Vettore works without it." },
  }));
  await page.goto(`/?token=${TOKEN}`);
  await openFile(page);
  await menu(page, "File", "Open from NCBI…");
  await page.getByTestId("ncbi-accession").fill("MT612431.1");
  await page.getByTestId("ncbi-open").click();
  await expect(page.getByTestId("ncbi-error")).toContainText("could not reach NCBI");
  await expect(page.getByTestId("ncbi-dialog")).toBeVisible();
  await expect(page.getByTestId("ncbi-accession")).toHaveValue("MT612431.1");
  await expect(page.locator(".vt-doctab")).toHaveCount(1);
});

test("a name that is not an accession is explained by the engine, nothing is downloaded", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("welcome-ncbi").click();
  await page.getByTestId("ncbi-accession").fill("pUC19");
  await page.getByTestId("ncbi-open").click();
  await expect(page.getByTestId("ncbi-error")).toContainText("is not a GenBank accession number");
});

test("turned off in Preferences, the menu item is hidden", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await openFile(page);
  await menu(page, "File", "Preferences…");
  const prefs = page.getByTestId("preferences-dialog");
  await expect(prefs.getByTestId("pref-ncbi")).toBeChecked();
  await prefs.getByText("Allow Open from NCBI").click();
  await expect(prefs.getByTestId("pref-ncbi")).not.toBeChecked();
  await page.getByTestId("pref-close").click();

  await page.locator(".tg-menu-bar").getByText("File", { exact: true }).click();
  await expect(page.locator(".bp3-popover .bp3-menu").getByText("Open…", { exact: true })).toBeVisible();
  await expect(page.locator(".bp3-popover .bp3-menu").getByText("Open from NCBI…")).toHaveCount(0);
  await page.keyboard.press("Escape");

  // Kept after a restart, until turned on again; the welcome screen leaves it out too.
  await page.reload();
  await expect(page.getByTestId("welcome")).toBeVisible();
  await expect(page.getByTestId("welcome-ncbi")).toHaveCount(0);
});
