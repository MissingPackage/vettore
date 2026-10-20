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

const BACKBONE = "ATTGCCTACGGTACATGCAACGTTGACCATGGCATTACGGACCATTGACCAGTTACGGATCAGCTTAGGCATAC".repeat(14);
const VECTOR = BACKBONE.slice(0, 600) + "GATATC" + BACKBONE.slice(600, 1000);
const ORF = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCC";
const rc = (s: string) =>
  s.split("").reverse().map((c) => ({ A: "T", C: "G", G: "C", T: "A" })[c as "A" | "C" | "G" | "T"]).join("");

test("TA cloning shows both orientations of a PCR product in a T-vector", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await newSequence(page, "pT", VECTOR, true);
  await newSequence(page, "ORF PCR", ORF, false);

  await menu(page, "Actions", "TA or GC Cloning…");
  const dialog = page.getByTestId("cloning-dialog");
  await expect(dialog.getByTestId("cloning-help")).toContainText("T-vector");
  // Vector and PCR product are picked from the open files; only the cloning site is missing.
  await expect(dialog.getByTestId("cloning-missing")).toContainText("“pT” is circular: enter where it is opened");
  await dialog.getByTestId("vector-position").fill("603");   // between GAT and ATC

  const choices = dialog.getByTestId("cloning-products");
  await expect(choices).toContainText("2 possible products", { timeout: 20_000 });
  await expect(choices).toContainText("Insert forward");
  await expect(choices).toContainText("Insert reversed");
  // GAT + T + insert + A + ATC: one T from the vector, one A from the PCR product.
  await expect(dialog.getByTestId("cloning-product")).toContainText(`${VECTOR.length + ORF.length + 2} bp, circular`);
  await expect(dialog.getByTestId("cloning-warnings")).toContainText("either orientation");
  await page.screenshot({ path: "test-results/cloning-ta.png" });

  await page.getByTestId("cloning-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText("pT + ORF PCR");
  await page.getByTestId("view-history").click();
  await expect(page.locator(".vt-history")).toContainText("TA cloning");
});

test("directional TOPO asks for CACC on the forward primer and adds it on request", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await newSequence(page, "pENTR-like", VECTOR, true);
  await newSequence(page, "cDNA", BACKBONE.slice(0, 80) + ORF + BACKBONE.slice(80, 160), false);

  await menu(page, "Actions", "TOPO Cloning…");
  const dialog = page.getByTestId("cloning-dialog");
  await dialog.getByTestId("ta-kind").selectOption("directional-TOPO");
  await dialog.getByTestId("vector-position").fill("603");
  await dialog.getByTestId("insert-source").selectOption("pcr");
  await expect(dialog.getByTestId("template-doc")).toHaveValue(/doc-/);
  await dialog.getByTestId("forward-primer").fill(ORF.slice(0, 22));
  await dialog.getByTestId("reverse-primer").fill(rc(ORF.slice(-22)));

  await expect(dialog.getByTestId("cacc-missing")).toBeVisible();
  await expect(dialog.getByTestId("cloning-error")).toContainText("does not start with CACC", { timeout: 20_000 });
  await dialog.getByTestId("add-cacc").click();
  await expect(dialog.getByTestId("forward-primer")).toHaveValue(`CACC${ORF.slice(0, 22)}`);

  const product = dialog.getByTestId("cloning-product");
  await expect(product).toContainText(`${VECTOR.length + 4 + ORF.length} bp, circular`, { timeout: 20_000 });
  await expect(product).toContainText("insert forward");
  await expect(dialog.getByTestId("cloning-products")).toHaveCount(0);
  await expect(dialog.getByTestId("cloning-oligos")).toContainText(`CACC${ORF.slice(0, 22)}`);
  await page.screenshot({ path: "test-results/cloning-topo.png" });

  await page.getByTestId("cloning-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText("pENTR-like + cDNA PCR product");
});
