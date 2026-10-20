// Files written in vector vendors' notation open as their authors meant them: the colour kept
// in a note (/note="color: #63b349; direction: RIGHT") and non-standard feature keys (ITR,
// Promoter, Miscellaneous, Transcription_t). tests/data/fidelity-notes.gb is NCBI PZ629150.1
// annotated that way at the published record's own coordinates (see tests/data/README.md).
import { expect, test, type Page } from "@playwright/test";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

/** The fill colours of the features drawn on the plasmid map. */
function mapColours(page: Page): Promise<string[]> {
  return page.locator(".veCircularView .veCircularViewFeature path[fill]")
    .evaluateAll((paths) => paths.map((p) => (p.getAttribute("fill") ?? "").toLowerCase()));
}

test("a vendor's file shows the author's colours and the standard feature types", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("fidelity-notes.gb"));
  await expect(page.locator(".veCircularView")).toBeVisible();
  await expect(page.locator(".veCircularView .veCircularViewFeature path[fill]")).toHaveCount(12);

  const colours = await mapColours(page);
  // From the notes: "color: #63b349" (5' ITR), "color: Blue" (cTNT promoter), "#FF378A"
  // (the terminator, written in capitals), "#ce0084" (SaCas9) and the others.
  for (const colour of ["#63b349", "#0000ff", "#ff378a", "#ce0084", "#5739ec", "#1f36a9",
                        "#129405", "#bde475", "#f60a2d", "#946a00"]) {
    expect(colours, colour).toContain(colour);
  }
  // "color: None" on the 3' ITR, and no colour note on the AmpR promoter: the type's colour.
  expect(colours).toContain("#966363");                // repeat_region
  expect(colours).toContain("#31b440");                // promoter
  await page.screenshot({ path: "test-results/fidelity-vendor-colours.png" });

  // The Features view groups them under the standard types.
  await page.getByTestId("view-features").click();
  const panel = page.locator(".ve-propertiesPanel");
  const row = (name: string) => panel.locator(".rt-tr-group").filter({
    has: page.getByText(name, { exact: true }),
  });
  await expect(row("5' ITR")).toContainText("repeat_region");            // "ITR"
  await expect(row("Kozak")).toContainText("misc_feature");              // "Miscellaneous"
  await expect(row("terminator")).not.toContainText("misc_feature");     // "Transcription_t"
  await expect(row("f1 ori")).toContainText("rep_origin");               // "Rep_origin"

  // ITR features keep their names, so the AAV packaging check finds them.
  await page.getByTestId("view-map").click();
  await menu(page, "Tools", "AAV Packaging Check…");
  const dialog = page.getByTestId("aav-dialog");
  await expect(dialog.getByTestId("aav-total")).toHaveText("4,483 bp");
});

test.describe("saving", () => {
  let dir: string;

  test.beforeEach(async ({ page }) => {
    dir = fs.mkdtempSync(path.join(os.tmpdir(), "vettore-e2e-fidelity-"));
    // The desktop bridge as pywebview injects it, with the native dialogs answered here.
    await page.addInitScript(() => {
      (window as any).__openPaths = [];
      (window as any).__nextSavePath = null;
      (window as any).pywebview = {
        api: {
          open_dialog: async () => (window as any).__openPaths,
          save_dialog: async () => (window as any).__nextSavePath,
        },
      };
    });
    await page.goto(`/?token=${TOKEN}`);
  });

  test.afterEach(() => fs.rmSync(dir, { recursive: true, force: true }));

  test("a vendor's file is saved with each colour once and the original keys kept", async ({ page }) => {
    const original = path.join(dir, "vendor.gb");
    fs.copyFileSync(data("fidelity-notes.gb"), original);
    await page.evaluate((p) => { (window as any).__openPaths = [p]; }, original);
    await page.getByTestId("open").click();
    await expect(page.locator(".vt-doctab.active")).toContainText("vendor");

    const copy = path.join(dir, "saved.gb");
    await page.evaluate((p) => { (window as any).__nextSavePath = p; }, copy);
    await menu(page, "File", "Save As…");
    await expect(page.locator(".vt-doctab.active")).toContainText("saved");
    const text = fs.readFileSync(copy, "utf-8");

    // One feature block of the saved GenBank file, by its label.
    const block = (label: string) => {
      const blocks = text.split(/\n(?= {5}\S)/);
      return blocks.find((b) => b.includes(`/label="${label}"`)) ?? "";
    };
    const itr = block("5' ITR");
    expect(itr).toMatch(/^ {5}repeat_region {3}1\.\.131/);
    expect(itr).toContain('/ApEinfo_fwdcolor="#63b349"');
    expect(itr).toContain('/note="direction: RIGHT"');
    expect(itr).toContain('/original_type="ITR"');
    expect(block("terminator")).toMatch(/^ {5}terminator {6}4324\.\.4328/);
    expect(block("cTNT promoter")).toContain('/ApEinfo_fwdcolor="#0000ff"');
    // The colour is written once, as a colour qualifier, never again in a note.
    expect(text).not.toMatch(/\/note="colou?r:/i);
    // A key the map shows as misc_feature keeps the file's own key in a note.
    expect(block("Cd36 sgRNA")).toMatch(/^ {5}misc_feature/);
    expect(block("Cd36 sgRNA")).toContain('/original_type="sgRNA"');
  });
});
