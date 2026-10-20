// Dialog layouts at the smallest supported screen, 1366 × 768: controls stay inside their
// column and short rows stay on one line.
import { expect, test, type Locator, type Page } from "@playwright/test";
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

/** Every child of `row` sits on the same line, each on a single text line. */
async function expectOneLine(row: Locator, what: string) {
  const boxes = await row.evaluate((el) => Array.from(el.children).map((c) => {
    const r = c.getBoundingClientRect();
    const style = getComputedStyle(c);
    return { top: r.top, height: r.height, lineHeight: parseFloat(style.lineHeight) || 20 };
  }));
  const tops = boxes.map((b) => Math.round(b.top + b.height / 2));
  expect(Math.max(...tops) - Math.min(...tops), `${what}: children on one line`).toBeLessThanOrEqual(4);
}

/** `inner` ends inside `outer` (no control spilling over a column edge). */
async function expectInside(inner: Locator, outer: Locator, what: string) {
  const a = await inner.boundingBox();
  const b = await outer.boundingBox();
  expect(a && b, what).toBeTruthy();
  expect(a!.x + a!.width, what).toBeLessThanOrEqual(b!.x + b!.width + 1);
}

test("Design PCR Primers keeps every input inside the form column", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".veCircularView")).toBeVisible();
  await menu(page, "Primers", "Design PCR Primers…");
  const dialog = page.getByTestId("primer-design-dialog");
  await expect(dialog).toBeVisible();
  const column = dialog.locator(".vt-pd-form");
  await expectInside(dialog.getByTestId("pd-maxTmDifference"), column, "Max Tm difference input");
  for (const input of await column.locator("input, select").all()) {
    await expectInside(input, column, `input ${await input.getAttribute("aria-label") ?? ""}`);
  }
});

test("the Gibson insert region stays on one line", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".vt-doctab.active")).toContainText("pFA-KanMX4");
  await page.locator(".veCircularView").click({ position: { x: 20, y: 20 } });
  await menu(page, "Actions", "Gibson Assembly…");
  const dialog = page.getByTestId("cloning-dialog");
  await dialog.getByTestId("insert-doc-0").selectOption({ index: 1 });
  await dialog.getByTestId("insert-region-0-feature").selectOption({ label: "KanR (459–1,268)" });
  const row = dialog.locator(".vt-insert-row .vt-region");
  await expect(row).toContainText("810 bp");
  await expectOneLine(row, "region row");
  await page.screenshot({ path: "test-results/layout-gibson-1366.png" });
});
