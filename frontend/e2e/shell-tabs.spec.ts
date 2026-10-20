// The document tab strip keeps the open document's tab in sight, however many files are open.
import { expect, test, type Page } from "@playwright/test";
import { TOKEN } from "../playwright.config";

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

/** The active tab lies wholly inside the visible part of the strip. */
async function expectActiveTabInSight(page: Page) {
  const { tab, strip } = await page.locator(".vt-doctab.active").evaluate((el) => {
    const t = el.getBoundingClientRect();
    const s = el.parentElement!.getBoundingClientRect();
    return { tab: { left: t.left, right: t.right }, strip: { left: s.left, right: s.right } };
  });
  expect(tab.left, "active tab starts inside the strip").toBeGreaterThanOrEqual(strip.left - 1);
  expect(tab.right, "active tab ends inside the strip").toBeLessThanOrEqual(strip.right + 1);
}

test.use({ viewport: { width: 1366, height: 768 } });

test("a new document's tab is scrolled into sight, and so is a tab chosen at the start", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  // Long names, as cloning products get: eight of them overflow a 1366-pixel window.
  for (let i = 1; i <= 8; i++) {
    const name = `cloning-gateway-destination PCR product ${i}`;
    if (i === 1) await page.getByTestId("new").click();
    else await menu(page, "File", "New DNA File…");
    await page.getByLabel("Name").fill(name);
    await page.getByTestId("new-seq-create").click();
    await expect(page.locator(".vt-doctab.active")).toContainText(name);
    await expectActiveTabInSight(page);
  }
  const overflow = await page.locator(".vt-doctabs").evaluate((el) => el.scrollWidth > el.clientWidth);
  expect(overflow, "the eight tabs are wider than the window").toBe(true);
  await page.screenshot({ path: "test-results/shell-tabs-many.png" });

  // Switching documents draws the other document's strip from its start: the tab chosen at
  // the far end must still be in sight once its document is shown.
  await page.getByTestId("doctab").first().click();
  await expect(page.locator(".vt-doctab.active")).toContainText("product 1");
  await expectActiveTabInSight(page);
  await page.getByTestId("doctab").last().click();
  await expect(page.locator(".vt-doctab.active")).toContainText("product 8");
  await expectActiveTabInSight(page);
});
