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

async function openDesigner(page: Page) {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".veCircularView")).toBeVisible();
  await menu(page, "Primers", "Design PCR Primers…");
  const dialog = page.getByTestId("primer-design-dialog");
  await expect(dialog).toBeVisible();
  return dialog;
}

test("design primers around a feature, copy them and add them to the plasmid", async ({ page, context }) => {
  await context.grantPermissions(["clipboard-read", "clipboard-write"]);
  const dialog = await openDesigner(page);
  // Nothing selected yet: the dialog says how to choose the region.
  await expect(dialog.getByTestId("pd-region-error")).toContainText("Pick a feature");

  await dialog.getByTestId("pd-feature").selectOption({ label: "KanR (459–1,268)" });
  await expect(dialog.getByTestId("pd-start")).toHaveValue("459");
  await expect(dialog.getByTestId("pd-end")).toHaveValue("1268");
  const pairs = dialog.getByTestId("pd-pair");
  await expect(pairs).toHaveCount(5, { timeout: 30_000 });
  await expect(dialog.getByTestId("pd-pairs")).toContainText("°C");
  await expect(dialog).toContainText("SantaLucia & Hicks 2004");
  await page.screenshot({ path: "test-results/primer-design.png" });

  await page.getByTestId("pd-copy").click();
  await expect(page.getByTestId("pd-notice")).toContainText("Copied KanR_F and KanR_R");
  const clipboard = await page.evaluate(() => navigator.clipboard.readText());
  expect(clipboard).toMatch(/^KanR_F\t[ACGT]{18,24}\nKanR_R\t[ACGT]{18,24}$/);

  await page.getByTestId("pd-add").click();
  await expect(page.getByTestId("pd-notice")).toContainText("Added KanR_F and KanR_R to pFA-KanMX4");
  await page.getByTestId("pd-close").click();

  await page.getByTestId("view-primers").click();
  await expect(page.getByText("KanR_F", { exact: true }).first()).toBeVisible();
  await expect(page.getByText("KanR_R", { exact: true }).first()).toBeVisible();
  await page.getByTestId("view-history").click();
  await expect(page.locator(".vt-history")).toContainText("Added primers (primer3)");
});

test("a pair is added once, and the chosen feature stays shown", async ({ page }) => {
  const dialog = await openDesigner(page);
  await dialog.getByTestId("pd-feature").selectOption({ label: "KanR (459–1,268)" });
  await expect(dialog.getByTestId("pd-pair")).toHaveCount(5, { timeout: 30_000 });
  await expect(dialog.getByTestId("pd-feature").locator("option:checked")).toHaveText("KanR (459–1,268)");

  const add = page.getByTestId("pd-add");
  await add.click();
  await expect(page.getByTestId("pd-notice")).toContainText("Added KanR_F and KanR_R to pFA-KanMX4");
  await expect(add).toHaveText("Added");
  await expect(add).toBeDisabled();
  // Another pair can still be added; going back to the first says it is already there.
  await dialog.getByTestId("pd-pair").nth(1).click();
  await expect(add).toHaveText("Add to Sequence");
  await page.getByTestId("pd-close").click();

  // Opened again later: the pair is recognised as already in the file.
  await menu(page, "Primers", "Design PCR Primers…");
  await dialog.getByTestId("pd-feature").selectOption({ label: "KanR (459–1,268)" });
  await expect(dialog.getByTestId("pd-pair")).toHaveCount(5, { timeout: 30_000 });
  await expect(dialog.getByTestId("pd-already")).toContainText("Pair 1 is already in pFA-KanMX4 as KanR_F and KanR_R");
  await expect(add).toBeDisabled();
});

test("a region with no room for primers says what to choose instead, in the result pane", async ({ page }) => {
  const dialog = await openDesigner(page);
  // The whole 3,941 bp plasmid: an 18 nt primer on each side needs 36 bp outside the region.
  await dialog.getByTestId("pd-start").fill("1");
  await dialog.getByTestId("pd-end").fill("3941");
  const message = dialog.getByTestId("pd-no-room");
  await expect(message).toContainText("The region covers 3,941 of the 3,941 bp of “pFA-KanMX4”");
  await expect(message).toContainText("Choose a region at least 36 bp shorter than the plasmid.");
  await expect(message).toContainText("Actions → Design Primers to Amplify a Region");
  // Nothing is reported under Product size, a field the user never typed in.
  await expect(dialog.getByTestId("pd-product-error")).toHaveCount(0);

  await dialog.getByTestId("pd-end").fill("3905");
  await expect(message).toHaveCount(0);
});

test("impossible settings explain what to change", async ({ page }) => {
  const dialog = await openDesigner(page);
  await dialog.getByTestId("pd-start").fill("459");
  await dialog.getByTestId("pd-end").fill("1268");
  await expect(dialog.getByTestId("pd-pair")).toHaveCount(5, { timeout: 30_000 });

  // Inline check next to the field, before anything is sent.
  await dialog.getByTestId("pd-product-max").fill("500");
  await expect(dialog.getByTestId("pd-product-error")).toContainText("At least 846 bp");
  await dialog.getByTestId("pd-product-max").fill("");

  // A Tm window no primer can reach: primer3's reasons, in plain words.
  await dialog.getByTestId("pd-tmMax").fill("85");
  await dialog.getByTestId("pd-tmOpt").fill("82");
  await dialog.getByTestId("pd-tmMin").fill("80");
  await expect(dialog.getByTestId("pd-error")).toContainText("Tm too low", { timeout: 30_000 });
  await expect(dialog.getByTestId("pd-error")).toContainText("Widen the Tm");
});
