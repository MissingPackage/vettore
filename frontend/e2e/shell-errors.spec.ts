// Engine errors reach the user as sentences, whatever the engine sent back.
import { expect, test } from "@playwright/test";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

test("a crash while listing the samples is reported, not swallowed", async ({ page }) => {
  // A server crash answers with the plain text "Internal Server Error", not JSON.
  await page.route("**/api/app/samples", (route) =>
    route.fulfill({ status: 500, contentType: "text/plain", body: "Internal Server Error" }));
  await page.goto(`/?token=${TOKEN}`);
  const notice = page.getByTestId("notice").filter({ hasText: "Something went wrong inside Vettore (error 500)" });
  await expect(notice).toBeVisible();
  await expect(notice).toContainText("Help → About Vettore");
  await expect(page.getByTestId("notices")).not.toContainText("[object");
});

test("a request the engine rejects as malformed reads as a sentence", async ({ page }) => {
  await page.route("**/api/open-content", (route) => route.fulfill({
    status: 422,
    contentType: "application/json",
    body: JSON.stringify({
      detail: [{ type: "missing", loc: ["body", "content_b64"], msg: "Field required", input: {} }],
    }),
  }));
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  const notice = page.getByTestId("notice").filter({ hasText: "could not use the request" });
  await expect(notice).toBeVisible();
  await expect(notice).toContainText("content_b64: field required");
  await expect(page.getByTestId("notices")).not.toContainText("[object");
});

test("the status bar says when the engine cannot give a Tm, and why", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(page.locator(".vt-doctab.active")).toContainText("pFA-KanMX4");
  await page.route("**/api/tm", (route) => route.fulfill({
    status: 422,
    contentType: "application/json",
    body: JSON.stringify({ detail: "Tm needs an unambiguous DNA sequence (A, C, G, T)" }),
  }));
  // PacI's 8-base site is long enough for a Tm.
  await page.getByTestId("view-enzymes").click();
  await page.getByTestId("enzymes-view").getByRole("cell", { name: "PacI", exact: true }).click();
  const tm = page.getByTestId("selection-tm-error");
  await expect(tm).toHaveText("Tm: not available");
  await expect(tm).toHaveAttribute("title", "Tm needs an unambiguous DNA sequence (A, C, G, T)");
});
