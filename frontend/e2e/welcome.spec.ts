import { expect, test } from "@playwright/test";
import { TOKEN } from "../playwright.config";

test("the welcome screen says what Vettore does and how to get in", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  const welcome = page.getByTestId("welcome");
  await expect(welcome).toBeVisible();
  await expect(welcome).toContainText("plasmid maps");
  await expect(welcome.getByTestId("open")).toContainText("Ctrl+O");
  await expect(welcome.getByTestId("new")).toContainText("Ctrl+N");
  await expect(welcome.getByTestId("drop-hint")).toContainText(".ab1");
  await page.screenshot({ path: "test-results/welcome.png" });
});

test("a sample construct opens as an unsaved document named after the construct",
  async ({ page }) => {
    await page.goto(`/?token=${TOKEN}`);
    const sample = page.getByTestId("sample-puc19");
    await expect(sample).toBeVisible();
    await sample.click();

    await expect(page.locator(".vt-doctab.active")).toContainText("pUC19");
    await expect(page.locator(".veCircularView")).toBeVisible();
    await expect(page.getByTestId("status-line")).toContainText("2,686 bp · circular");

    // Unsaved: saving asks where to put it instead of writing inside the application.
    const download = page.waitForEvent("download");
    await page.locator(".tg-menu-bar").getByText("File", { exact: true }).click();
    await page.locator(".bp3-popover .bp3-menu").getByText("Save", { exact: true }).last().click();
    expect((await download).suggestedFilename()).toBe("pUC19.gb");
  });

test("the AAV sample opens with its ITRs on the map", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("sample-paav-ctnt-sacas9").click();
  await expect(page.locator(".vt-doctab.active")).toContainText("pAAV-cTNT-SaCas9-U6");
  await page.getByTestId("view-features").click();
  await expect(page.locator(".ve-propertiesPanel").getByText(/inverted terminal repeat/).first())
    .toBeVisible();
});

test("help is reachable before any document is open", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("welcome-getting-started").click();
  await expect(page.getByTestId("getting-started-dialog")).toContainText("Open a construct");
  // A dialog still opening can miss the key: press again until it has closed.
  await expect(async () => {
    await page.keyboard.press("Escape");
    await expect(page.getByTestId("getting-started-dialog")).toBeHidden({ timeout: 1000 });
  }).toPass({ timeout: 10_000 });
  await page.getByTestId("welcome-shortcuts").click();
  await expect(page.getByTestId("shortcuts-dialog")).toContainText("Ctrl+O");
});
