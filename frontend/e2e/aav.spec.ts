import { expect, test, type Page } from "@playwright/test";
import fs from "node:fs";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

/** The wild-type AAV2 ITR as the engine's catalogue holds it (NC_001401.2:1..145). */
function catalogueItr(): string {
  const file = path.resolve(import.meta.dirname,
                            "../../src/vettore/data/common_features.json");
  const catalogue = JSON.parse(fs.readFileSync(file, "utf-8")) as {
    features: { name: string; sequence: string }[];
  };
  return catalogue.features.find((f) => f.name === "AAV2 ITR")!.sequence;
}

/** Filler DNA with no catalogue element in it, built the same way on every run. */
function filler(length: number): string {
  let seed = 7;
  let out = "";
  while (out.length < length) {
    seed = (seed * 1103515245 + 12345) % 2147483648;
    out += "ACGT"[seed % 4];
  }
  return out;
}

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

test("the cargo of an AAV transfer plasmid is measured ITR to ITR", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("sample-paav-ctnt-sacas9").click();
  await expect(page.locator(".veCircularView")).toBeVisible();

  await menu(page, "Tools", "AAV Packaging Check…");
  const dialog = page.getByTestId("aav-dialog");
  await expect(dialog).toBeVisible();

  // PZ629150.1 annotates its ITRs at 1..131 and 4343..4483: 4483 bp ITR to ITR, which fits.
  await expect(dialog.getByTestId("aav-total")).toHaveText("4,483 bp");
  await expect(dialog.getByTestId("aav-insert")).toHaveText("4,211 bp");
  await expect(dialog.getByTestId("aav-verdict")).toContainText("Fits a single-stranded AAV");
  await expect(dialog.getByTestId("aav-bar")).toBeVisible();
  await page.waitForTimeout(400);                       // let the dialog finish opening
  await page.screenshot({ path: "test-results/aav-packaging.png" });

  // The same cassette is far too long for a self-complementary vector.
  await dialog.getByTestId("aav-kind").selectOption("sc");
  await expect(dialog.getByTestId("aav-verdict")).toContainText("Too long");
  await expect(dialog.getByTestId("aav-verdict")).toContainText("self-complementary");

  // The other way round the circle is the backbone, and it says so.
  await dialog.getByTestId("aav-side").selectOption("1");
  await expect(dialog.getByTestId("aav-total")).toHaveText("2,868 bp");
  await expect(dialog.getByTestId("aav-sizes")).toContainText("across the origin");
  await expect(dialog.getByTestId("aav-backbone-warning")).toContainText("backbone");
});

test("the limit caption sits under the limit marker", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("sample-paav-ctnt-sacas9").click();
  await expect(page.locator(".veCircularView")).toBeVisible();
  await menu(page, "Tools", "AAV Packaging Check…");
  const dialog = page.getByTestId("aav-dialog");
  for (const kind of ["ss", "sc"]) {
    await dialog.getByTestId("aav-kind").selectOption(kind);
    await expect(dialog.getByTestId("aav-bar")).toBeVisible();
    // The dialog opens with a zoom: measure once it is over.
    await expect(page.locator(".bp3-overlay-appear-active, .bp3-overlay-enter-active")).toHaveCount(0);
    const bar = await dialog.getByTestId("aav-bar").boundingBox();
    const marker = await dialog.locator(".vt-aav-limit").boundingBox();
    const caption = await dialog.getByTestId("aav-limit-caption").boundingBox();
    expect(bar && marker && caption).toBeTruthy();
    // The marker is above the caption, at the same share of the caption as of the bar.
    const at = (marker!.x + marker!.width / 2 - bar!.x) / bar!.width;
    expect(Math.abs(caption!.x + at * caption!.width - (marker!.x + marker!.width / 2)),
           `${kind}: caption under the marker`).toBeLessThanOrEqual(3);
    expect(caption!.x, `${kind}: caption inside the bar`).toBeGreaterThanOrEqual(bar!.x - 1);
    expect(caption!.x + caption!.width, `${kind}: caption inside the bar`).toBeLessThanOrEqual(bar!.x + bar!.width + 1);
    // On its own line: it does not run into the scale's end label.
    const end = await dialog.getByTestId("aav-scale").locator("span").nth(1).boundingBox();
    expect(caption!.y, `${kind}: caption below the end label`).toBeGreaterThanOrEqual(end!.y + end!.height - 1);
  }
});

test("ITRs annotated only by a note are found too", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("sample-paav-ef1a-dio").click();
  await expect(page.locator(".veCircularView")).toBeVisible();
  await menu(page, "Tools", "AAV Packaging Check…");
  // MT612431.1 carries "5' ITR" and "3' ITR" notes at 20..131 and 4477..4588.
  await expect(page.getByTestId("aav-total")).toHaveText("4,569 bp");
  await expect(page.getByTestId("aav-verdict")).toContainText("Fits");
});

test("ITRs that are there but not annotated are found with the catalogue", async ({ page }) => {
  const itr = catalogueItr();
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("new").click();
  await page.getByLabel("Name").fill("unannotated-aav");
  await page.getByTestId("new-seq-text").fill(itr + filler(2000) + itr);
  await page.getByTestId("new-seq-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText("unannotated-aav");

  await menu(page, "Tools", "AAV Packaging Check…");
  const dialog = page.getByTestId("aav-dialog");
  await expect(dialog.getByTestId("aav-no-itrs")).toBeVisible();
  await dialog.getByTestId("aav-detect").click();
  // 145 bp ITR + 2000 bp of filler + 145 bp ITR, ITRs included.
  await expect(dialog.getByTestId("aav-total")).toHaveText("2,290 bp", { timeout: 30_000 });
  await expect(dialog.getByTestId("aav-insert")).toHaveText("2,000 bp");
  await expect(dialog.getByTestId("aav-sizes")).toContainText("found by the catalogue");
  await expect(dialog).toContainText("not annotated in the document yet");
});

test("a plasmid without ITRs says so and offers to search for them", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".veCircularView")).toBeVisible();

  await menu(page, "Tools", "AAV Packaging Check…");
  const dialog = page.getByTestId("aav-dialog");
  await expect(dialog.getByTestId("aav-no-itrs")).toContainText("No ITR is annotated");
  await dialog.getByTestId("aav-detect").click();
  await expect(dialog.getByTestId("aav-error")).toContainText("annotate them by hand",
                                                              { timeout: 30_000 });
  await page.screenshot({ path: "test-results/aav-no-itr.png" });
});
