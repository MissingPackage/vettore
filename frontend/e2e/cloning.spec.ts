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

async function newSequence(page: Page, name: string, sequence: string, circular: boolean) {
  const welcome = page.getByTestId("new");
  if (await welcome.isVisible().catch(() => false)) {
    await welcome.click();
  } else {
    await menu(page, "File", "New DNA File…");
  }
  await page.getByLabel("Name").fill(name);
  await page.getByTestId("new-seq-text").fill(sequence);
  if (circular) {
    // Click the label: Blueprint covers the radio input with its own indicator.
    await page.locator("label.bp3-radio", { hasText: "Circular" }).click();
    await expect(page.getByTestId("topology-circular")).toBeChecked();
  }
  await page.getByTestId("new-seq-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText(name);
}

const BACKBONE = "ATTGCCTACGGTACATGCAACGTTGACCATGGCATTACGGACCATTGACCAGTTACGGATCAGCTTAGGCATAC".repeat(14);
const VECTOR = BACKBONE.slice(0, 500) + "GAATTC" + BACKBONE.slice(500, 900) + "GGATCC" + BACKBONE.slice(900, 1400);
const INSERT = "CC" + "GAATTC" + "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCC" + "GGATCC" + "CC";

test("restriction cloning creates a product tab with a history entry", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await newSequence(page, "pVector", VECTOR, true);
  await newSequence(page, "Insert", INSERT, false);

  await menu(page, "Actions", "Restriction and Insertion Cloning", "Insert Fragment…");
  const dialog = page.getByTestId("cloning-dialog");
  await dialog.getByTestId("vector-doc").selectOption({ index: 1 });
  await dialog.getByTestId("vector-enzymes").fill("EcoRI");
  await dialog.getByTestId("vector-enzymes").fill("BamHI");
  await dialog.getByTestId("insert-doc").selectOption({ index: 2 });
  await dialog.getByTestId("insert-enzymes").fill("EcoRI");
  await dialog.getByTestId("insert-enzymes").fill("BamHI");

  await expect(page.getByTestId("cloning-product")).toContainText("pVector + Insert");
  await expect(page.getByTestId("cloning-preview")).toBeVisible();
  await page.screenshot({ path: "test-results/cloning-restriction.png" });

  await page.getByTestId("cloning-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText("pVector + Insert");
  await page.getByTestId("view-history").click();
  await expect(page.locator(".vt-history")).toContainText("Restriction cloning");
  await expect(page.locator(".vt-history")).toContainText("EcoRI–BamHI fragment of Insert");
});

test("Gibson assembly designs primers with homology arms", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await newSequence(page, "pVector", VECTOR, true);
  await newSequence(page, "GFP", INSERT.slice(8, -8), false);

  await menu(page, "Actions", "Gibson Assembly…");
  const dialog = page.getByTestId("cloning-dialog");
  await dialog.getByTestId("vector-doc").selectOption({ index: 1 });
  await dialog.getByTestId("vector-enzymes").fill("EcoRI");
  await dialog.getByTestId("insert-doc-0").selectOption({ index: 2 });

  await expect(page.getByTestId("cloning-oligos")).toBeVisible();
  const rows = page.getByTestId("cloning-oligos").locator("tbody tr");
  await expect(rows).toHaveCount(2);
  await expect(rows.first()).toContainText("°C");
  // The suggested name says how the product was made, so it is not confused with another
  // product of the same two files.
  await expect(dialog.getByTestId("cloning-name")).toHaveValue("pVector + GFP (Gibson)");
  await page.screenshot({ path: "test-results/cloning-gibson.png" });
});

test("a Gibson insert from a plasmid is picked by its feature name", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await newSequence(page, "pVector", VECTOR, true);
  await page.getByTestId("file-input").setInputFiles(path.resolve(import.meta.dirname, "../../tests/data/pFA-KanMX4.dna"));
  await expect(page.locator(".vt-doctab.active")).toContainText("pFA-KanMX4");

  await menu(page, "Actions", "Gibson Assembly…");
  const dialog = page.getByTestId("cloning-dialog");
  await dialog.getByTestId("vector-doc").selectOption({ index: 1 });
  await dialog.getByTestId("vector-enzymes").fill("EcoRI");
  await dialog.getByTestId("insert-doc-0").selectOption({ index: 2 });
  // KanR is annotated at 459..1268 in pFA-KanMX4: 1268 − 459 + 1 = 810 bp.
  await dialog.getByTestId("insert-region-0-feature").selectOption({ label: "KanR (459–1,268)" });
  await expect(dialog.getByTestId("insert-region-0-start")).toHaveValue("459");
  await expect(dialog.getByTestId("insert-region-0-end")).toHaveValue("1268");
  await expect(dialog.locator(".vt-insert-row")).toContainText("810 bp");
  // The chosen feature stays shown while the positions are its own.
  await expect(dialog.getByTestId("insert-region-0-feature").locator("option:checked")).toHaveText("KanR (459–1,268)");
  await expect(page.getByTestId("cloning-oligos").locator("tbody tr")).toHaveCount(2, { timeout: 30_000 });
  // Typing another position leaves the feature: the picker goes back to its prompt.
  await dialog.getByTestId("insert-region-0-end").fill("1200");
  await expect(dialog.getByTestId("insert-region-0-feature")).toHaveValue("");
});

test("a second product of the same files gets its own name", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await newSequence(page, "pVector", VECTOR, true);
  await newSequence(page, "Insert", INSERT, false);
  for (const expected of ["pVector + Insert (restriction cloning)", "pVector + Insert (restriction cloning) 2"]) {
    await page.locator(".vt-doctab", { hasText: "pVector" }).first().click();
    await menu(page, "Actions", "Restriction and Insertion Cloning", "Insert Fragment…");
    const dialog = page.getByTestId("cloning-dialog");
    await dialog.getByTestId("vector-doc").selectOption({ index: 1 });
    await dialog.getByTestId("vector-enzymes").fill("EcoRI");
    await dialog.getByTestId("vector-enzymes").fill("BamHI");
    await dialog.getByTestId("insert-doc").selectOption({ index: 2 });
    await dialog.getByTestId("insert-enzymes").fill("EcoRI");
    await dialog.getByTestId("insert-enzymes").fill("BamHI");
    await expect(dialog.getByTestId("cloning-name")).toHaveValue(expected);
    await page.getByTestId("cloning-create").click();
    await expect(page.locator(".vt-doctab.active")).toContainText(expected);
  }
});

test("Mutagenesis and In-Fusion say in plain words what they do", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await newSequence(page, "pVector", VECTOR, true);
  await menu(page, "Actions", "Mutagenesis…");
  await expect(page.getByTestId("cloning-help")).toContainText("KLD mix (kinase, ligase and DpnI");
  await page.getByRole("button", { name: "Cancel" }).click();
  await menu(page, "Actions", "In-Fusion Cloning…");
  await expect(page.getByTestId("cloning-help")).toContainText("repeat the last 15 bp of each vector end");
});
