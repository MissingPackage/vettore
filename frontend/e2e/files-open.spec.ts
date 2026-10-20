// Opening the files colleagues send: other formats, wrong extensions, several records, quirks.
import { expect, test, type Page } from "@playwright/test";
import fs from "node:fs";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

const tab = (page: Page, name: string) => page.getByTestId("doctab").filter({ hasText: name });

async function dropFile(page: Page, name: string, content: Buffer) {
  const transfer = await page.evaluateHandle(({ name, b64 }) => {
    const dt = new DataTransfer();
    dt.items.add(new File([Uint8Array.from(atob(b64), (c) => c.charCodeAt(0))], name));
    return dt;
  }, { name, b64: content.toString("base64") });
  await page.dispatchEvent(".vt-shell", "drop", { dataTransfer: transfer });
}

test.beforeEach(async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
});

test("opens EMBL, .xdna, .gck and plain-sequence files", async ({ page }) => {
  await page.getByTestId("file-input").setInputFiles([
    data("exchange-pBR322.embl"),
    data("exchange-sample-b.xdna"),
    data("exchange-artificial.gck"),
    data("exchange-pUC19-numbered.seq"),
  ]);
  for (const name of ["exchange-pBR322", "exchange-sample-b", "exchange-artificial", "exchange-pUC19-numbered"]) {
    await expect(tab(page, name)).toBeVisible();
  }
  await tab(page, "exchange-pBR322").click();
  await expect(page.getByTestId("status-line")).toContainText("4,361 bp · circular");
  await tab(page, "exchange-pUC19-numbered").click();
  await expect(page.getByTestId("status-line")).toContainText("2,686 bp · linear");
});

test("a file with several sequences opens one tab each and says so", async ({ page }) => {
  await page.getByTestId("file-input").setInputFiles(data("exchange-plasmids.fasta"));
  for (const name of ["L09137.2", "J01749.1", "X06403.1"]) await expect(tab(page, name)).toBeVisible();
  const notice = page.getByTestId("notice").filter({ hasText: "3 sequences opened" });
  await expect(notice).toContainText("3 sequences opened from “exchange-plasmids.fasta”, one tab each");
  await expect(notice).toHaveAttribute("data-kind", "info");
});

test("a GenBank file saved as .txt opens by its content", async ({ page }) => {
  await page.getByTestId("file-input").setInputFiles({
    name: "from-a-colleague.txt",
    mimeType: "text/plain",
    buffer: fs.readFileSync(data("MT612431.1.gb")),
  });
  await expect(page.locator(".vt-doctab.active")).toContainText("from-a-colleague");
  await expect(page.getByTestId("status-line")).toContainText("7,204 bp · circular");
  await page.getByTestId("view-history").click();
  await expect(page.locator(".vt-history")).toContainText("(GenBank)");
});

test("errors show their whole message, stack up and can be dismissed", async ({ page }) => {
  await page.getByTestId("file-input").setInputFiles([
    { name: "sequence.txt", mimeType: "text/plain", buffer: Buffer.from("ACGTACGT\nACGTXACGT\n") },
    { name: "report.docx", mimeType: "application/octet-stream", buffer: Buffer.from([0x50, 0x4b, 3, 4, 0, 0, 0, 1]) },
  ]);
  const errors = page.locator('[data-testid="notice"][data-kind="error"]');
  await expect(errors).toHaveCount(2);
  await expect(errors.nth(0)).toContainText("“sequence.txt” could not be read as a DNA sequence: line 2, column 5 has “X”");
  // The end of a long message is visible too: nothing is cut off.
  await expect(errors.nth(1)).toContainText("and plain sequence text (.txt, .seq).");
  for (const i of [0, 1]) {
    const text = errors.nth(i).locator(".vt-toast-text");
    const { scrollWidth, clientWidth, scrollHeight, clientHeight } = await text.evaluate((el) => ({
      scrollWidth: el.scrollWidth, clientWidth: el.clientWidth, scrollHeight: el.scrollHeight, clientHeight: el.clientHeight,
    }));
    expect(scrollWidth).toBeLessThanOrEqual(clientWidth);
    expect(scrollHeight).toBeLessThanOrEqual(clientHeight);
  }
  await page.screenshot({ path: "test-results/files-errors.png" });
  await errors.nth(0).getByTestId("notice-dismiss").click();
  await expect(errors).toHaveCount(1);
  await expect(errors.first()).toContainText("report.docx");
  // The same error again is counted, not stacked.
  await page.getByTestId("file-input").setInputFiles({
    name: "report.docx", mimeType: "application/octet-stream", buffer: Buffer.from([0x50, 0x4b, 3, 4, 0, 0, 0, 1]),
  });
  await expect(errors).toHaveCount(1);
  await expect(errors.first()).toContainText("(×2)");
});

test("features that cannot be placed are reported by name", async ({ page }) => {
  await page.getByTestId("file-input").setInputFiles(data("exchange-quirks-features.gb"));
  await expect(page.locator(".vt-doctab.active")).toContainText("exchange-quirks-features");
  const warnings = page.locator('[data-testid="notice"][data-kind="warning"]');
  await expect(warnings.first()).toContainText("“test site between bases” (between two bases)");
  await expect(warnings.nth(1)).toContainText("“test region past the end”");
  await page.getByTestId("notice-dismiss-all").click();
  await expect(page.getByTestId("notice")).toHaveCount(0);
});

test("dropping files on the window opens them", async ({ page }) => {
  await dropFile(page, "exchange-pBR322.embl", fs.readFileSync(data("exchange-pBR322.embl")));
  await expect(page.locator(".vt-doctab.active")).toContainText("exchange-pBR322");
  await expect(page.getByTestId("status-line")).toContainText("4,361 bp · circular");
});

test("the file chooser offers every format the engine reads", async ({ page }) => {
  const accept = (await page.getByTestId("file-input").getAttribute("accept")) ?? "";
  for (const ext of [".dna", ".gb", ".gbk", ".ape", ".embl", ".fasta", ".xdna", ".gck", ".ab1", ".txt", ".seq"]) {
    expect(accept.split(",")).toContain(ext);
  }
});
