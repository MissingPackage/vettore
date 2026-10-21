// Tools → Find Diagnostic Digest (S5.1): from a cloning product to the digest that tells it apart
// from its wrong clones, its gel, the protocol text, and the hand-over to the gel simulation.
import { expect, test, type Page } from "@playwright/test";
import { readFileSync } from "node:fs";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

/** 1,000 bp of pACYC184 (no BamHI site) between two BamHI sites, as a PCR product. */
function insertSequence(): string {
  const fasta = readFileSync(data("exchange-plasmids.fasta"), "utf8");
  const record = fasta.split(">").find((r) => r.startsWith("X06403"))!;
  const sequence = record.split("\n").slice(1).join("").trim().toUpperCase();
  return `CGCGGATCC${sequence.slice(2000, 3000)}GGATCCCGC`;
}

/** pBR322 and a BamHI insert cloned into its BamHI site: the product is the active tab. */
async function cloneInsert(page: Page) {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("exchange-pBR322.embl"));
  await expect(page.locator(".veCircularView")).toBeVisible();
  await menu(page, "File", "New DNA File…");
  await page.getByLabel("Name").fill("cat PCR");
  await page.getByTestId("new-seq-text").fill(insertSequence());
  await page.getByTestId("new-seq-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText("cat PCR");

  await menu(page, "Actions", "Restriction and Insertion Cloning", "Insert Fragment…");
  const cloning = page.getByTestId("cloning-dialog");
  await cloning.getByTestId("vector-doc").selectOption({ index: 1 });
  await cloning.getByTestId("vector-enzymes").fill("BamHI");
  await cloning.getByTestId("insert-doc").selectOption({ index: 2 });
  await cloning.getByTestId("insert-enzymes").fill("BamHI");
  await expect(page.getByTestId("cloning-product")).toContainText("5367 bp");
  await page.getByTestId("cloning-create").click();
  await expect(page.getByTestId("cloning-dialog")).toHaveCount(0);
}

test("a cloning product brings its wrong clones and gets a ranked digest with its gel", async ({ page, context }) => {
  await context.grantPermissions(["clipboard-read", "clipboard-write"]);
  await cloneInsert(page);

  await menu(page, "Tools", "Find Diagnostic Digest…");
  const dialog = page.getByTestId("digest-dialog");
  await expect(dialog).toBeVisible();
  await expect(dialog.getByTestId("digest-wrong")).toContainText("insert in the other orientation");
  await expect(dialog.getByTestId("digest-wrong")).toContainText("parent vector");
  await expect(dialog.getByTestId("digest-wrong")).toContainText("from the cloning simulation");
  // "My enzymes" is empty by default: every commercial enzyme is tried.
  await expect(dialog.getByTestId("digest-set-commercial")).toBeChecked();
  await expect(dialog.getByTestId("digest-set-mine")).toBeDisabled();
  await expect(dialog.getByTestId("digest-mine-hint")).toContainText("File → Preferences");

  const results = dialog.getByTestId("digest-results");
  await expect(results).toBeVisible();
  expect(await results.locator("li").count()).toBeLessThanOrEqual(10);
  await expect(results.locator("li").first()).toContainText("Expected clone:");
  await expect(results.locator("li").first()).toContainText("insert in the other orientation:");
  await expect(dialog.getByTestId("digest-gel")).toBeVisible();
  // Ladder plus the expected clone and its two wrong clones.
  await expect(dialog.getByTestId("digest-gel").getByTestId("gel-lane-header")).toHaveCount(4);
  await page.screenshot({ path: "test-results/digest-planner.png" });

  await dialog.getByTestId("digest-result").nth(1).check();
  await page.getByTestId("digest-copy").click();
  await expect(dialog.getByTestId("digest-copied")).toBeVisible();
  const text = await page.evaluate(() => navigator.clipboard.readText());
  expect(text).toMatch(/Diagnostic digest for “.*pBR322 \+ cat PCR.*”/);
  expect(text).toContain("parent vector: ");

  await page.getByTestId("digest-send-to-gel").click();
  const gel = page.getByTestId("gel-dialog");
  await expect(gel).toBeVisible();
  await expect(gel.getByTestId("gel-image").getByTestId("gel-lane-header")).toHaveCount(4);
  await expect(gel.getByTestId("gel-doc-2")).toHaveValue(/diagnostic-/);
});

test("my enzymes are tried first; when they cannot tell the clones apart the dialog says why", async ({ page }) => {
  await cloneInsert(page);
  await menu(page, "File", "Preferences…");
  const prefs = page.getByTestId("preferences-dialog");
  await prefs.getByTestId("pref-my-enzymes-list").fill("EcoRI, PstI; Hind3");
  await prefs.getByTestId("pref-my-enzymes-add").click();
  await expect(prefs.getByTestId("pref-my-enzymes-unknown")).toContainText("“Hind3”");
  await expect(prefs.getByTestId("pref-my-enzymes-field")).toContainText("EcoRI");
  await expect(prefs.getByTestId("pref-my-enzymes-field")).toContainText("PstI");
  await page.getByTestId("pref-close").click();

  await menu(page, "Tools", "Find Diagnostic Digest…");
  const dialog = page.getByTestId("digest-dialog");
  await expect(dialog.getByTestId("digest-set-mine")).toBeChecked();
  await expect(dialog.getByTestId("digest-none")).toContainText("same length as the expected clone");
  await dialog.getByTestId("digest-try-commercial").click();
  await expect(dialog.getByTestId("digest-set-commercial")).toBeChecked();
  await expect(dialog.getByTestId("digest-results")).toBeVisible();
});

test("an open file can be added and removed as a wrong clone", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("exchange-plasmids.fasta"));
  await expect(page.locator(".vt-doctab")).toHaveCount(3);
  await menu(page, "Tools", "Find Diagnostic Digest…");
  const dialog = page.getByTestId("digest-dialog");
  await expect(dialog.getByTestId("digest-nothing-to-compare")).toBeVisible();
  await dialog.getByTestId("digest-add").selectOption({ index: 1 });
  await expect(dialog.getByTestId("digest-wrong")).toContainText("open file");
  await expect(dialog.getByTestId("digest-results")).toBeVisible();
  await dialog.getByTestId("digest-wrong-remove").click();
  await expect(dialog.getByTestId("digest-nothing-to-compare")).toBeVisible();
});
