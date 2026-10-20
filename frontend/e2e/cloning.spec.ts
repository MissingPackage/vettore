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
    const radio = page.getByTestId("topology-circular");
    await radio.waitFor();
    await radio.check({ force: true });
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
  await page.screenshot({ path: "test-results/cloning-gibson.png" });
});
