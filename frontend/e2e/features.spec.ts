import { expect, test, type Page } from "@playwright/test";
import fs from "node:fs";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

/** The bases of a GenBank record, from its ORIGIN section. */
function genbankSequence(name: string): string {
  const text = fs.readFileSync(data(name), "utf-8");
  return text.split(/\nORIGIN[^\n]*\n/)[1].split("\n//")[0].replace(/[^a-z]/gi, "").toUpperCase();
}

/** A new circular document with these bases and no features. */
async function newPlasmid(page: Page, name: string, sequence: string) {
  await page.getByTestId("new").click();
  await page.getByLabel("Name").fill(name);
  await page.getByTestId("new-seq-text").fill(sequence);
  await page.locator("label.bp3-radio", { hasText: "Circular" }).click();
  await page.getByTestId("new-seq-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText(name);
}

test("detect common features on a real AAV plasmid and add them to the map", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  // The published record already annotates what the catalogue finds; its bases alone do not.
  await newPlasmid(page, "pAAV-cTNT-SaCas9-U6", genbankSequence("PZ629150.1.gb"));
  await expect(page.locator(".veCircularView")).toBeVisible();

  await menu(page, "Features", "Detect Common Features…");
  await expect(page.getByTestId("detect-busy")).toBeVisible();
  const table = page.getByTestId("detect-table");
  await expect(table).toBeVisible({ timeout: 30_000 });
  await expect(page.getByTestId("detect-identity")).toHaveValue("95");
  for (const name of ["SaCas9", "U6 promoter (human)", "AmpR", "AAV2 ITR (truncated)", "ori"]) {
    await expect(table).toContainText(name);
  }
  await expect(table).toContainText("PZ629150.1:581..3925");   // provenance is shown
  await page.screenshot({ path: "test-results/detect-features.png" });

  await page.getByTestId("detect-add").click();
  await expect(page.getByTestId("detect-dialog")).toBeHidden();

  // The features are now on the map and in the Features view.
  await page.getByTestId("view-features").click();
  const features = page.locator(".ve-propertiesPanel");
  await expect(features.getByText("SaCas9", { exact: true }).first()).toBeVisible();
  await expect(features.getByText("U6 promoter (human)", { exact: true }).first()).toBeVisible();

  await page.getByTestId("view-map").click();
  await expect(page.locator(".veCircularView")).toContainText("SaCas9");
  await page.screenshot({ path: "test-results/detect-map.png" });
});

type Row = { name: string; start: number; end: number; isNew: boolean; annotatedAs?: string | null };
type Existing = { name: string; start: number; end: number };

test("an element already annotated under another name is shown as such and not ticked", async ({ page }) => {
  // The engine names the existing feature that covers a match (annotatedAs). An engine that
  // does not yet is given the same answer here, from the features the dialog sent.
  await page.route("**/api/features/detect", async (route) => {
    const sent = route.request().postDataJSON().sequenceData.features as Existing[] | Record<string, Existing>;
    const existing = Array.isArray(sent) ? sent : Object.values(sent);
    const response = await route.fetch();
    const json = await response.json() as { features: Row[]; new_count: number };
    if (json.features.every((f) => !("annotatedAs" in f))) {
      for (const f of json.features) {
        const cover = existing.find((e) => e.start <= f.start && f.end <= e.end && f.end - f.start >= (e.end - e.start) / 2);
        f.annotatedAs = cover?.name ?? null;
        f.isNew = !cover;
      }
      json.new_count = json.features.filter((f) => f.isNew).length;
    }
    await route.fulfill({ response, json });
  });
  await page.goto(`/?token=${TOKEN}`);
  // nptI (kanamycin resistance) is annotated at complement(820..1635), which the catalogue calls KanR.
  await page.getByTestId("file-input").setInputFiles(data("cloning-gateway-destination.gb"));
  await expect(page.locator(".veCircularView")).toBeVisible();
  await menu(page, "Features", "Detect Common Features…");
  const row = page.getByTestId("detect-row").filter({ hasText: "KanR" });
  await expect(row).toBeVisible({ timeout: 30_000 });
  await expect(row.getByTestId("detect-annotated")).toHaveText("· annotated as “nptI”");
  await expect(row.locator("input[type=checkbox]")).not.toBeChecked();
});

// tests/data/fidelity-*-cmv.gb: a CMV-EGFP plasmid built from catalogue sequences, with one
// vendor-style block annotation over the CMV enhancer and promoter (see tests/data/README.md).

test("elements inside a block named after them are already annotated", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("fidelity-block-cmv.gb"));
  await expect(page.locator(".veCircularView")).toBeVisible();
  await menu(page, "Features", "Detect Common Features…");
  await expect(page.getByTestId("detect-table")).toBeVisible({ timeout: 30_000 });

  // The 204 bp CMV promoter sits in the 584 bp "CMV promoter" block: not offered again.
  await expect(page.getByTestId("detect-row")).toHaveCount(7);
  const cmv = page.getByTestId("detect-row").filter({ hasText: /^CMV promoter/ });
  await expect(cmv.getByTestId("detect-annotated")).toHaveText("· already annotated");
  await expect(page.getByTestId("detect-row").filter({ hasText: "CMV enhancer" })
    .getByTestId("detect-annotated")).toHaveText("· annotated as “CMV promoter”");
  await expect(page.locator("[data-testid=detect-row] input[type=checkbox]:checked")).toHaveCount(0);
  await expect(page.getByTestId("detect-label-summary")).toHaveCount(0);
  await expect(page.getByTestId("detect-add")).toBeDisabled();
});

test("a label that names another element is flagged and not ticked", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("fidelity-mislabelled-cmv.gb"));
  await expect(page.locator(".veCircularView")).toBeVisible();
  await menu(page, "Features", "Detect Common Features…");
  await expect(page.getByTestId("detect-table")).toBeVisible({ timeout: 30_000 });

  const summary = page.getByTestId("detect-label-summary");
  await expect(summary).toContainText("1 label may be wrong");
  const enhancer = page.getByTestId("detect-row").filter({ hasText: "CMV enhancer" });
  await expect(enhancer.getByTestId("detect-label-warning")).toHaveText(
    "Annotated as “PGK promoter”, but the sequence is the CMV enhancer (100 % identity). Check the label.");
  const promoter = page.getByTestId("detect-row").filter({ hasText: /^CMV promoter/ });
  await expect(promoter.getByTestId("detect-label-warning")).toContainText("the CMV promoter");
  await expect(page.getByTestId("detect-label-warning")).toHaveCount(2);
  await expect(enhancer.locator("input[type=checkbox]")).not.toBeChecked();
  await expect(promoter.locator("input[type=checkbox]")).not.toBeChecked();
  // Amber, like every "check this" message, and not only by colour: the text says it.
  const background = await enhancer.locator("td").first()
    .evaluate((td) => getComputedStyle(td).backgroundColor);
  expect(background).toBe("rgb(255, 250, 235)");      // --vt-warning-bg
  await page.waitForTimeout(400);                     // let the dialog finish opening
  await page.screenshot({ path: "test-results/detect-label-warning.png" });

  // Ticking one adds the element under its own name.
  await promoter.locator("input[type=checkbox]").check();
  await page.getByTestId("detect-add").click();
  await expect(page.getByTestId("detect-dialog")).toBeHidden();
  await page.getByTestId("view-features").click();
  await expect(page.locator(".ve-propertiesPanel").getByText("CMV promoter", { exact: true })).toBeVisible();
});

test.describe("on a 1366 × 768 screen", () => {
  test.use({ viewport: { width: 1366, height: 768 } });

  test("the minimum identity control stays inside the dialog", async ({ page }) => {
    await page.goto(`/?token=${TOKEN}`);
    await page.getByTestId("file-input").setInputFiles(data("cloning-gateway-destination.gb"));
    await expect(page.locator(".veCircularView")).toBeVisible();
    await menu(page, "Features", "Detect Common Features…");
    const dialog = page.getByTestId("detect-dialog");
    await expect(page.getByTestId("detect-table")).toBeVisible({ timeout: 30_000 });
    // The dialog opens with a zoom: measure once it is over.
    await expect(page.locator(".bp3-overlay-appear-active, .bp3-overlay-enter-active")).toHaveCount(0);
    await page.screenshot({ path: "test-results/detect-features-1366.png" });
    const body = await dialog.boundingBox();
    const select = await page.getByTestId("detect-identity").boundingBox();
    expect(select!.x + select!.width).toBeLessThanOrEqual(body!.x + body!.width + 1);
  });
});

test("a match that covers part of its catalogue entry says which bases", async ({ page }) => {
  // mCherry in pT-PGK-mCherry has GTC for the catalogue's ATG and stops with TAG, not TAA:
  // the row shows the 708 bp from the second codon to the stop, and that they are bases
  // 4–711 of the entry.
  const fasta = fs.readFileSync(data("features-MT891341.1.fasta"), "utf-8");
  const bases = fasta.split("\n").filter((line) => !line.startsWith(">")).join("").toUpperCase();
  await page.goto(`/?token=${TOKEN}`);
  await newPlasmid(page, "pT-PGK-mCherry", bases);
  await menu(page, "Features", "Detect Common Features…");
  const row = page.getByTestId("detect-row").filter({ hasText: "mCherry" });
  await expect(row).toBeVisible({ timeout: 30_000 });
  await expect(row).toContainText("708 bp");
  await expect(row.getByTestId("detect-part")).toHaveText("· bases 4–711 of 711");
  await expect(row).toContainText("99.8 %");
});

test("a coding sequence with a broken reading frame is flagged and not ticked", async ({ page }) => {
  const catalogue = JSON.parse(fs.readFileSync(
    path.resolve(import.meta.dirname, "../../src/vettore/data/common_features.json"), "utf-8"));
  const egfp: string = catalogue.features.find((f: { name: string }) => f.name === "EGFP").sequence;
  // Filler that matches nothing in the catalogue, the same on every run.
  let state = 7;
  const filler = Array.from({ length: 400 }, () => {
    state = (state * 1103515245 + 12345) % 2 ** 31;
    return "ACGT"[state % 4];
  }).join("");
  // One base missing in the middle of EGFP: still a 99.9 % match, but not EGFP after it.
  await page.goto(`/?token=${TOKEN}`);
  await newPlasmid(page, "EGFP-frameshift", filler + egfp.slice(0, 360) + egfp.slice(361) + filler);
  await menu(page, "Features", "Detect Common Features…");
  const row = page.getByTestId("detect-row").filter({ hasText: "EGFP" });
  await expect(row).toBeVisible({ timeout: 30_000 });
  await expect(row.getByTestId("detect-frameshift")).toContainText("breaks the reading frame");
  await expect(row.getByRole("checkbox")).not.toBeChecked();
});
