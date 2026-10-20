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

test("simulate an agarose gel of a digest and of the uncut plasmid", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".veCircularView")).toBeVisible();

  await menu(page, "Tools", "Simulate Agarose Gel…");
  const dialog = page.getByTestId("gel-dialog");
  await expect(dialog).toBeVisible();
  await expect(dialog).toContainText("Leave a lane without enzymes");
  await dialog.getByTestId("gel-doc-0").selectOption({ index: 1 });
  await expect(page.getByTestId("gel-image")).toBeVisible();
  // Uncut circular plasmid: one supercoiled band.
  await expect(dialog).toContainText("pFA-KanMX4 / uncut");

  await dialog.getByTestId("gel-add-lane").click();
  await dialog.getByTestId("gel-doc-1").selectOption({ index: 1 });
  await dialog.getByTestId("gel-enzymes-1").fill("NdeI");
  await dialog.getByTestId("gel-enzymes-1").fill("NcoI");
  await expect(dialog).toContainText("pFA-KanMX4 / NdeI + NcoI");
  await page.screenshot({ path: "test-results/gel.png" });

  const bands = page.getByTestId("gel-image").locator("rect");
  expect(await bands.count()).toBeGreaterThan(20);
});

async function openGel(page: Page) {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".veCircularView")).toBeVisible();
  await menu(page, "Tools", "Simulate Agarose Gel…");
  const dialog = page.getByTestId("gel-dialog");
  await dialog.getByTestId("gel-doc-0").selectOption({ index: 1 });
  await expect(page.getByTestId("gel-image")).toBeVisible();
  return dialog;
}

test("ladder sizes never overprint and the ladder's name is readable", async ({ page }) => {
  const dialog = await openGel(page);
  await dialog.getByTestId("gel-ladder").selectOption({ label: "NEB 1 kb Plus (N3200)" });
  // 517 and 500 bp run a few pixels apart: one of them is written, the 500 bp reference band.
  await expect(page.getByTestId("gel-ladder-label").filter({ hasText: /^500$/ })).toHaveCount(1);
  await expect(page.getByTestId("gel-ladder-label").filter({ hasText: /^517$/ })).toHaveCount(0);
  const ys = (await page.getByTestId("gel-ladder-label").evaluateAll(
    (els) => els.map((e) => Number(e.getAttribute("y"))))).sort((a, b) => a - b);
  for (let i = 1; i < ys.length; i++) expect(ys[i] - ys[i - 1]).toBeGreaterThanOrEqual(10);

  const header = page.getByTestId("gel-lane-header").first();
  await expect(header).toContainText("NEB 1 kb Plus");
  await expect(header).toContainText("(N3200)");
  await expect(header).not.toContainText("…");
});

test("a lane whose enzyme does not cut is drawn uncut, with a note under the gel", async ({ page }) => {
  // The engine draws such a lane uncut and explains it in lanes[i].note: modelled here on the
  // engine's uncut lane, so that the dialog's side is tested on its own.
  await page.route("**/api/gel", async (route) => {
    const body = route.request().postDataJSON();
    const idle = body.lanes.map((l: { enzymes: string[] }) => l.enzymes.includes("SwaI"));
    const response = await route.fetch({ postData: JSON.stringify({
      ...body,
      lanes: body.lanes.map((l: { enzymes: string[] }, i: number) => (idle[i] ? { ...l, enzymes: [] } : l)),
    }) });
    const json = await response.json();
    json.lanes.forEach((lane: { name: string; note?: string }, i: number) => {
      if (!idle[i]) return;
      lane.name = lane.name.replace("uncut", "SwaI");
      lane.note = "SwaI does not cut “pFA-KanMX4”, so the lane shows the uncut plasmid.";
    });
    await route.fulfill({ response, json });
  });
  const dialog = await openGel(page);
  await dialog.getByTestId("gel-add-lane").click();
  await dialog.getByTestId("gel-doc-1").selectOption({ index: 1 });
  // pFA-KanMX4 has no SwaI site (ATTTAAAT): the plasmid stays whole, 3,941 bp.
  await dialog.getByTestId("gel-enzymes-1").fill("SwaI");
  await expect(dialog.getByTestId("gel-note")).toHaveText(
    "Lane 2: SwaI does not cut “pFA-KanMX4”, so the lane shows the uncut plasmid.");
  await expect(page.getByTestId("gel-image")).toBeVisible();
  await expect(dialog.getByTestId("gel-error")).toHaveCount(0);
  await expect(dialog).toContainText("pFA-KanMX4 / SwaI: 3,941 bp");
});
