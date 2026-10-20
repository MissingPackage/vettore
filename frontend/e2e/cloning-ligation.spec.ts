import { expect, test, type Page } from "@playwright/test";
import { TOKEN } from "../playwright.config";

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
    await page.locator("label.bp3-radio", { hasText: "Circular" }).click();
    await expect(page.getByTestId("topology-circular")).toBeChecked();
  }
  await page.getByTestId("new-seq-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText(name);
}

// No EcoRI, BamHI or XhoI site in either filler (same backbone as cloning.spec.ts).
const BACKBONE = "ATTGCCTACGGTACATGCAACGTTGACCATGGCATTACGGACCATTGACCAGTTACGGATCAGCTTAGGCATAC".repeat(14);
const GFP = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCC";
const VECTOR = BACKBONE.slice(0, 500) + "GAATTC" + "AAAACCCC" + "CTCGAG" + BACKBONE.slice(500, 1000);
const PROMOTER = "CC" + "GAATTC" + GFP.slice(0, 60) + "GGATCC" + "CC";
const CDS = "CC" + "GGATCC" + GFP.slice(60) + "CTCGAG" + "CC";
// G^AATTC, G^GATCC, C^TCGAG: the stuffer between EcoRI and XhoI is replaced by both inserts.
const PRODUCT_LENGTH = 500 + 6 + 60 + 6 + (GFP.length - 60) + 6 + 500;

test("insert two fragments into a vector in one ligation", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await newSequence(page, "pVector", VECTOR, true);
  await newSequence(page, "Promoter", PROMOTER, false);
  await newSequence(page, "CDS", CDS, false);

  await menu(page, "Actions", "Restriction and Insertion Cloning", "Insert Multiple Fragments…");
  const dialog = page.getByTestId("cloning-dialog");
  // The open plasmid is the vector by default; the dialog says what is still missing.
  await expect(dialog.getByTestId("vector-doc")).toHaveValue(/doc-/);
  await expect(dialog.getByTestId("cloning-missing")).toContainText("Choose the enzymes that cut “pVector”");

  await dialog.getByTestId("vector-enzymes").fill("EcoRI");
  await dialog.getByTestId("vector-enzymes").fill("XhoI");
  await expect(dialog.getByTestId("vector-fragment")).toBeVisible();
  await dialog.getByTestId("insert-0-doc").selectOption({ index: 2 });
  await dialog.getByTestId("insert-0-enzymes").fill("EcoRI");
  await dialog.getByTestId("insert-0-enzymes").fill("BamHI");
  await dialog.getByTestId("insert-1-doc").selectOption({ index: 3 });
  await dialog.getByTestId("insert-1-enzymes").fill("BamHI");
  await dialog.getByTestId("insert-1-enzymes").fill("XhoI");

  const product = page.getByTestId("cloning-product");
  await expect(product).toContainText("pVector + Promoter + CDS", { timeout: 20_000 });
  await expect(product).toContainText(`${PRODUCT_LENGTH} bp, circular`);
  await expect(product).toContainText("Order: pVector → Promoter → CDS");
  await expect(page.getByTestId("cloning-preview")).toBeVisible();
  await page.screenshot({ path: "test-results/cloning-insert-multiple.png" });

  await page.getByTestId("cloning-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText("pVector + Promoter + CDS");
  await page.getByTestId("view-history").click();
  await expect(page.locator(".vt-history")).toContainText("Insert multiple fragments");
  // Every input is listed as a parent and can be reopened as it was.
  await expect(page.getByTestId("history-reopen")).toHaveCount(3);
});

test("ligate two blunt PCR products and choose among the products", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await newSequence(page, "Left PCR", GFP.slice(0, 60), false);
  await newSequence(page, "Right PCR", GFP.slice(60), false);

  await menu(page, "Actions", "Ligate Fragments…");
  const dialog = page.getByTestId("cloning-dialog");
  // Both open files are used whole by default: a first result with no input at all.
  const choices = dialog.getByTestId("cloning-products");
  await expect(choices).toContainText("4 possible products", { timeout: 20_000 });
  await expect(dialog.getByTestId("cloning-warnings")).toContainText("T4 polynucleotide kinase");
  await expect(dialog.getByTestId("cloning-product")).toContainText(`${GFP.length} bp, linear`);

  await choices.locator("label").nth(1).click();
  await dialog.getByTestId("cloning-name").fill("Joined PCR");
  await page.screenshot({ path: "test-results/cloning-ligate.png" });
  await page.getByTestId("cloning-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText("Joined PCR");
  await page.getByTestId("view-history").click();
  await expect(page.locator(".vt-history")).toContainText("Ligate fragments");
});

test("incompatible ends are explained before anything is created", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await newSequence(page, "Promoter", PROMOTER, false);
  await newSequence(page, "CDS", CDS, false);

  await menu(page, "Actions", "Ligate Fragments…");
  const dialog = page.getByTestId("cloning-dialog");
  // A sticky piece cut from the promoter file cannot join the blunt CDS PCR product.
  await dialog.getByTestId("fragment-0-doc").selectOption({ index: 1 });
  await dialog.getByTestId("fragment-0-mode").selectOption("digest");
  await dialog.getByTestId("fragment-0-enzymes").fill("EcoRI");
  await dialog.getByTestId("fragment-0-enzymes").fill("BamHI");
  await dialog.getByTestId("fragment-1-doc").selectOption({ index: 2 });
  await expect(dialog.getByTestId("cloning-error")).toContainText("cannot ligate", { timeout: 20_000 });
  await expect(page.getByTestId("cloning-create")).toBeDisabled();
});
