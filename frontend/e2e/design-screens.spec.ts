// Design harness: screenshots of every screen and dialog at the two window sizes the app is
// used at, plus layout checks that catch overlaps, truncation and missing focus rings.
// Screenshots go to test-results/design/<width>x<height>/ (git-ignored) for visual review.
import { expect, test, type Locator, type Page } from "@playwright/test";
import fs from "node:fs";
import path from "node:path";
import { TOKEN } from "../playwright.config";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);

const SIZES = [
  { width: 1366, height: 768 },
  { width: 1920, height: 1080 },
];

// Same construction as cloning.spec.ts: a circular vector with EcoRI and BamHI sites.
const BACKBONE = "ATTGCCTACGGTACATGCAACGTTGACCATGGCATTACGGACCATTGACCAGTTACGGATCAGCTTAGGCATAC".repeat(14);
const VECTOR = BACKBONE.slice(0, 500) + "GAATTC" + BACKBONE.slice(500, 900) + "GGATCC" + BACKBONE.slice(900, 1400);
const INSERT = "CC" + "GAATTC" + "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCC" + "GGATCC" + "CC";

async function menu(page: Page, ...items: string[]) {
  await page.locator(".tg-menu-bar").getByText(items[0], { exact: true }).click();
  for (const item of items.slice(1)) {
    await page.locator(".bp3-popover .bp3-menu").getByText(item, { exact: true }).last().click();
  }
}

/** Open a fixture; `name` is the document name (GenBank files take it from the LOCUS line). */
async function open(page: Page, file: string, name = path.parse(file).name) {
  await page.getByTestId("file-input").setInputFiles(data(file));
  await expect(page.locator(".vt-doctab.active")).toContainText(name);
}

async function newSequence(page: Page, name: string, sequence: string, circular: boolean) {
  const welcome = page.getByTestId("new");
  if (await welcome.isVisible().catch(() => false)) await welcome.click();
  else await menu(page, "File", "New DNA File…");
  await page.getByLabel("Name").fill(name);
  await page.getByTestId("new-seq-text").fill(sequence);
  if (circular) await page.locator("label.bp3-radio", { hasText: "Circular" }).click();
  await page.getByTestId("new-seq-create").click();
  await expect(page.locator(".vt-doctab.active")).toContainText(name);
}

/** The element fits horizontally: nothing is clipped or pushed into a scrollbar. */
async function expectNoHorizontalOverflow(locator: Locator) {
  const { scroll, client } = await locator.evaluate((el) => ({ scroll: el.scrollWidth, client: el.clientWidth }));
  expect(scroll, "content wider than its container").toBeLessThanOrEqual(client + 1);
}

/** The open dialog lies inside the window, so its title and footer buttons are reachable. */
async function expectDialogInViewport(page: Page) {
  const dialog = page.locator(".bp3-dialog").last();
  await expect(dialog).toBeVisible();
  const box = (await dialog.boundingBox())!;
  const viewport = page.viewportSize()!;
  expect(box.x).toBeGreaterThanOrEqual(0);
  expect(box.x + box.width).toBeLessThanOrEqual(viewport.width + 1);
  // Tall dialogs scroll inside the Blueprint container; the footer must still be on screen.
  const footer = dialog.locator(".bp3-dialog-footer");
  if (await footer.count()) await expect(footer.first()).toBeInViewport();
}

async function closeDialog(page: Page) {
  await page.locator(".bp3-dialog-close-button").last().click();
  await expect(page.locator(".bp3-dialog")).toHaveCount(0);
}

/**
 * WCAG contrast ratio between an element's text and the background behind it, blending every
 * ancestor's background colour from the page down (enough for flat surfaces like these).
 */
async function contrastOf(locator: Locator): Promise<number> {
  return locator.evaluate((el) => {
    type Rgba = { r: number; g: number; b: number; a: number };
    const parse = (c: string): Rgba => {
      const [r, g, b, a = 1] = (c.match(/[\d.]+/g) ?? ["0", "0", "0", "0"]).map(Number);
      return { r, g, b, a };
    };
    const over = (top: Rgba, under: Rgba): Rgba => ({
      r: top.r * top.a + under.r * (1 - top.a),
      g: top.g * top.a + under.g * (1 - top.a),
      b: top.b * top.a + under.b * (1 - top.a),
      a: 1,
    });
    const chain: Element[] = [];
    for (let n: Element | null = el; n; n = n.parentElement) chain.unshift(n);
    let background: Rgba = { r: 255, g: 255, b: 255, a: 1 };
    for (const n of chain) {
      const c = parse(getComputedStyle(n).backgroundColor);
      if (c.a > 0) background = over(c, background);
    }
    const text = over(parse(getComputedStyle(el).color), background);
    const luminance = ({ r, g, b }: Rgba) => {
      const f = (v: number) => (v /= 255) <= 0.03928 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4;
      return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b);
    };
    const [hi, lo] = [luminance(text), luminance(background)].sort((x, y) => y - x);
    return (hi + 0.05) / (lo + 0.05);
  });
}

/** Text meant to be read passes WCAG AA (4.5:1). */
async function expectReadable(locator: Locator, what: string) {
  expect(await contrastOf(locator.first()), `${what} contrast`).toBeGreaterThanOrEqual(4.5);
}

for (const size of SIZES) {
  const dir = `test-results/design/${size.width}x${size.height}`;
  // Animations are finished first, so dialogs are captured fully open rather than fading in.
  const shot = (page: Page, name: string) =>
    page.screenshot({ path: `${dir}/${name}.png`, animations: "disabled" });

  test.describe(`design at ${size.width}x${size.height}`, () => {
    test.use({ viewport: size });

    test("welcome, document views and shell chrome", async ({ page }) => {
      await page.goto(`/?token=${TOKEN}`);
      await expect(page.locator(".vt-welcome")).toBeVisible();
      // Text never falls back to Blueprint's icon font (wide spaces, heavy "+"), and form
      // controls use the page font rather than the browser's own.
      const fonts = await page.evaluate(() => ({
        body: getComputedStyle(document.body).fontFamily,
        button: getComputedStyle(document.querySelector(".vt-welcome button")!).fontFamily,
      }));
      expect(fonts.body).not.toContain("Icons16");
      expect(fonts.button).toBe(fonts.body);
      if (await page.locator(".vt-recent-dir").count()) await expectReadable(page.locator(".vt-recent-dir"), "recent file folder");
      await shot(page, "01-welcome");

      // Keyboard users see where focus is on Vettore's own controls.
      await page.getByTestId("open").focus();
      await page.keyboard.press("Tab");
      const outline = await page.getByTestId("new").evaluate((el) => getComputedStyle(el).outlineStyle);
      expect(outline).not.toBe("none");
      await shot(page, "02-welcome-focus");

      await open(page, "pFA-KanMX4.dna");
      await expect(page.locator(".veCircularView")).toBeVisible();
      await expectNoHorizontalOverflow(page.locator(".vt-viewbar"));
      await expectNoHorizontalOverflow(page.locator(".vt-status-bar"));
      await expectReadable(page.locator(".vt-viewtab:not(.active)"), "inactive view tab");
      await expectReadable(page.locator(".vt-enzyme-set"), "enzyme set label");
      await expectReadable(page.locator(".vt-status-set"), "status bar enzyme set");
      // The editor's toolbar buttons stay transparent (only filled primary buttons are recoloured).
      const toolbarBackgrounds = await page.locator(".veToolbar .bp3-button.bp3-minimal:not(.bp3-active)")
        .evaluateAll((els) => els.map((e) => getComputedStyle(e).backgroundColor));
      expect(toolbarBackgrounds.length).toBeGreaterThan(0);
      for (const bg of toolbarBackgrounds) expect(bg).toBe("rgba(0, 0, 0, 0)");
      await shot(page, "03-map");
      // Keyboard focus on the view tabs is drawn with the focus colour.
      await page.getByTestId("view-map").focus();
      await page.keyboard.press("Tab");
      await expect(page.getByTestId("view-sequence")).toBeFocused();
      expect(await page.getByTestId("view-sequence").evaluate((el) => getComputedStyle(el).outlineStyle)).toBe("solid");
      await shot(page, "03-viewtab-focus");

      await page.getByTestId("view-sequence").click();
      await expect(page.locator(".veRowView")).toBeVisible();
      await shot(page, "04-sequence");

      await page.getByTestId("view-enzymes").click();
      await expect(page.getByTestId("enzymes-view")).toBeVisible();
      await expectReadable(page.getByTestId("enzymes-view").locator(".vt-table-toolbar .vt-muted"), "enzyme count");
      await shot(page, "05-enzymes");

      await page.getByTestId("view-features").click();
      await expect(page.locator(".ve-propertiesPanel")).toBeVisible();
      await shot(page, "06-features");

      await open(page, "sample-hybridization-params.dna");
      await page.getByTestId("view-primers").click();
      await expect(page.locator(".ve-propertiesPanel")).toBeVisible();
      await shot(page, "07-primers");

      await page.getByTestId("view-history").click();
      await expect(page.locator(".vt-history")).toContainText("Opened");
      await shot(page, "08-history");

      // Several documents: tabs stay readable and the view bar keeps its controls.
      await open(page, "MT612431.1.gb", "MT612431");
      await open(page, "PZ629150.1.gb", "PZ629150");
      await expectNoHorizontalOverflow(page.locator(".vt-viewbar"));
      await shot(page, "09-many-tabs");

      await expectReadable(page.locator(".vt-doctab:not(.active) span"), "inactive document tab");

      await page.getByTestId("view-map").click();
      await page.locator(".tg-menu-bar").getByText("Actions", { exact: true }).click();
      const actionsMenu = page.locator(".bp3-popover .bp3-menu").first();
      await expect(actionsMenu).toBeVisible();
      // One item per row, also for disabled items with a shortcut (the editor wraps them in an
      // inline tooltip target, which used to put two on one line).
      const tops = await actionsMenu.locator("a.bp3-menu-item")
        .evaluateAll((els) => els.map((e) => Math.round(e.getBoundingClientRect().top)));
      expect(new Set(tops).size).toBe(tops.length);
      await shot(page, "10-actions-menu");
      await page.locator(".tg-menu-bar").getByText("Actions", { exact: true }).click();
      await expect(actionsMenu).toBeHidden();

      // A file that cannot be read: the error is shown in the window, not lost.
      await page.getByTestId("file-input").setInputFiles({
        name: "broken.gb", mimeType: "text/plain", buffer: Buffer.from("not a sequence file"),
      });
      await expect(page.locator(".vt-error")).toBeVisible();
      await shot(page, "11-open-error");
    });

    test("cloning dialogs, empty and with a result", async ({ page }) => {
      await page.goto(`/?token=${TOKEN}`);
      await newSequence(page, "pVector", VECTOR, true);
      await newSequence(page, "Insert", INSERT, false);

      const dialogs: [string, string[]][] = [
        ["20-restriction", ["Actions", "Restriction and Insertion Cloning", "Insert Fragment…"]],
        ["21-delete-fragment", ["Actions", "Restriction and Insertion Cloning", "Delete Restriction Fragment…"]],
        ["22-pcr", ["Actions", "PCR…"]],
        ["23-design-primers", ["Actions", "Design Primers to Amplify a Region…"]],
        ["24-mutagenesis", ["Actions", "Mutagenesis…"]],
        ["25-anneal-oligos", ["Actions", "Anneal Oligos and Insert…"]],
        ["26-gibson", ["Actions", "Gibson Assembly…"]],
        ["27-hifi", ["Actions", "NEBuilder HiFi DNA Assembly…"]],
        ["28-infusion", ["Actions", "In-Fusion Cloning…"]],
        ["29-golden-gate", ["Actions", "Golden Gate Assembly…"]],
      ];
      const dialog = page.getByTestId("cloning-dialog");
      for (const [name, path] of dialogs) {
        await menu(page, ...path);
        await expect(dialog).toBeVisible();
        await expectDialogInViewport(page);
        await expectReadable(dialog.locator(".vt-cloning-form > p"), `${name} explanation`);
        if (await dialog.locator(".vt-field-hint").count()) await expectReadable(dialog.locator(".vt-field-hint"), `${name} hint`);
        // One control height across the form: no 24px inputs next to 30px selects.
        const heights = await dialog.locator(".vt-cloning-form .bp3-input, .vt-cloning-form select")
          // offsetHeight ignores the dialog's opening zoom, which scales bounding boxes.
          .evaluateAll((els) => els.map((e) => (e as HTMLElement).offsetHeight));
        expect(new Set(heights), `${name} control heights`).toEqual(new Set([30]));
        await shot(page, name);
        await closeDialog(page);
      }

      // Restriction cloning with a product: preview, product line and buttons.
      await menu(page, "Actions", "Restriction and Insertion Cloning", "Insert Fragment…");
      await dialog.getByTestId("vector-doc").selectOption({ index: 1 });
      await dialog.getByTestId("vector-enzymes").fill("EcoRI");
      await dialog.getByTestId("vector-enzymes").fill("BamHI");
      await dialog.getByTestId("insert-doc").selectOption({ index: 2 });
      await dialog.getByTestId("insert-enzymes").fill("EcoRI");
      await dialog.getByTestId("insert-enzymes").fill("BamHI");
      await expect(page.getByTestId("cloning-product")).toContainText("pVector + Insert");
      await expectDialogInViewport(page);
      await expectReadable(dialog.locator(".vt-enzyme-picker .bp3-tag"), "enzyme chip");
      await shot(page, "30-restriction-result");
      await closeDialog(page);

      // Gibson with designed primers, then a failing input to show the error callout.
      await menu(page, "Actions", "Gibson Assembly…");
      await dialog.getByTestId("vector-doc").selectOption({ index: 1 });
      await dialog.getByTestId("vector-enzymes").fill("EcoRI");
      await dialog.getByTestId("insert-doc-0").selectOption({ index: 2 });
      await expect(page.getByTestId("cloning-oligos")).toBeVisible();
      await expectDialogInViewport(page);
      await shot(page, "31-gibson-result");
      await closeDialog(page);

      await menu(page, "Actions", "PCR…");
      await dialog.getByTestId("template-doc").selectOption({ index: 1 });
      const inputs = dialog.locator(".vt-primer-input input.vt-mono");
      await inputs.nth(0).fill("GGGGGGGGGGGGGGGGGGGG");
      await inputs.nth(1).fill("CCCCCCCCCCCCCCCCCCCC");
      await expect(page.getByTestId("cloning-error")).toBeVisible();
      await expectReadable(page.getByTestId("cloning-error"), "error message");
      await shot(page, "32-pcr-error");
      await closeDialog(page);
    });

    test("tool dialogs", async ({ page }) => {
      await page.goto(`/?token=${TOKEN}`);
      await page.getByTestId("new").click();
      await expect(page.locator(".bp3-dialog")).toBeVisible();
      await expectDialogInViewport(page);
      await shot(page, "40-new-dna-file");
      await closeDialog(page);

      await open(page, "pFA-KanMX4.dna");
      await menu(page, "Tools", "Simulate Agarose Gel…");
      const gel = page.getByTestId("gel-dialog");
      await expect(gel).toBeVisible();
      await expectDialogInViewport(page);
      await shot(page, "41-gel-empty");
      await expectReadable(gel.locator(".vt-field-hint"), "gel lane hint");
      await gel.getByTestId("gel-doc-0").selectOption({ index: 1 });
      await gel.getByTestId("gel-add-lane").click();
      await gel.getByTestId("gel-doc-1").selectOption({ index: 1 });
      await gel.getByTestId("gel-enzymes-1").fill("NdeI");
      await gel.getByTestId("gel-enzymes-1").fill("NcoI");
      await expect(page.getByTestId("gel-image")).toBeVisible();
      await expectDialogInViewport(page);
      await shot(page, "42-gel");
      await closeDialog(page);

      // Sanger reads against the reference built around the real read (see sanger.spec.ts).
      const reference = fs.readFileSync(path.resolve(import.meta.dirname, "fixtures/sanger-reference.txt"), "utf-8").trim();
      await newSequence(page, "pExpected", reference, false);
      await menu(page, "Tools", "Align to Reference DNA Sequence…");
      await expect(page.getByTestId("sanger-dialog")).toBeVisible();
      await expectDialogInViewport(page);
      await shot(page, "43-sanger-empty");
      await page.getByTestId("sanger-files").setInputFiles(data("3730.ab1"));
      await expect(page.getByTestId("sanger-summary")).toBeVisible({ timeout: 30_000 });
      await expect(page.locator(".vt-alignment svg path").first()).toBeVisible({ timeout: 20_000 });
      await expectDialogInViewport(page);
      await shot(page, "44-sanger");
      await closeDialog(page);

      await open(page, "PZ629150.1.gb", "PZ629150");
      await menu(page, "Features", "Detect Common Features…");
      await expect(page.getByTestId("detect-table")).toBeVisible({ timeout: 30_000 });
      await expectDialogInViewport(page);
      await shot(page, "45-detect-features");
      await closeDialog(page);

      await menu(page, "Features", "Add Feature…");
      await expect(page.locator(".bp3-dialog")).toBeVisible();
      await shot(page, "46-add-feature");
      await closeDialog(page);

      await menu(page, "Actions", "Simulate Digestion…");
      await expect(page.getByText("Virtual Digest").first()).toBeVisible();
      await shot(page, "47-simulate-digestion");
    });
  });
}

// The shared classes from docs/DESIGN.md, rendered with the app's real stylesheet. New dialogs
// adopt them; this keeps their contrast and look under test before they do.
const SPECIMEN = `
  <div class="bp3-dialog-body" style="margin:0">
    <p class="vt-dialog-intro">Cut the vector and the insert with the same enzymes and ligate them.</p>
    <p class="vt-section-title">Inputs</p>
    <label class="vt-field"><span class="vt-field-label">Vector</span>
      <input class="bp3-input" value="pUC19"><span class="vt-field-hint">One enzyme, or two for directional cloning</span></label>
    <label class="vt-field vt-field-invalid"><span class="vt-field-label">Overlap length (bp)</span>
      <input class="bp3-input" value="4"><span class="vt-field-error">Use 15 to 40 bp.</span></label>
    <div class="vt-notice vt-notice-danger"><span class="vt-notice-title">No product</span>The primer does not bind.</div>
    <div class="vt-notice vt-notice-warning"><span class="vt-notice-title">Check before ordering</span>BsaI also cuts the insert.</div>
    <div class="vt-notice vt-notice-success"><span class="vt-notice-title">Matches</span>All reads agree with the reference.</div>
    <div class="vt-notice vt-notice-info"><span class="vt-notice-title">Note</span>Positions are approximate.</div>
    <p class="vt-error">The ladder file could not be read.</p>
    <div class="vt-empty-state"><p class="vt-empty-state-title">No product yet</p><p>Choose a vector on the left.</p></div>
    <p>Press <span class="vt-kbd">Ctrl</span>+<span class="vt-kbd">O</span> · <span class="vt-busy">Simulating…</span>
      · <span class="vt-num">59.7 °C</span> · <span class="vt-ok">100.00 %</span> · <span class="vt-warn">98.1 %</span>
      · <span class="vt-bad">failed</span> · <span class="vt-muted">89 of 1165 nt</span></p>
  </div>`;

test("shared classes render readable with the app stylesheet", async ({ page }) => {
  await page.goto(`/?token=${TOKEN}`);
  await expect(page.locator(".vt-welcome")).toBeVisible();
  await page.evaluate((html) => {
    const panel = document.createElement("div");
    panel.id = "vt-specimen";
    panel.className = "bp3-dialog";
    panel.style.cssText = "position:fixed;top:20px;left:20px;width:560px;z-index:100;padding:16px";
    panel.innerHTML = html;
    document.body.appendChild(panel);
  }, SPECIMEN);
  const specimen = page.locator("#vt-specimen");
  for (const selector of [
    ".vt-dialog-intro", ".vt-section-title", ".vt-field-hint", ".vt-field-error", ".vt-notice-title",
    ".vt-notice", ".vt-error", ".vt-empty-state p", ".vt-empty-state-title", ".vt-kbd", ".vt-busy",
    ".vt-ok", ".vt-warn", ".vt-bad", ".vt-muted",
  ]) {
    const all = specimen.locator(selector);
    for (let i = 0; i < await all.count(); i++) await expectReadable(all.nth(i), selector);
  }
  // The invalid field is outlined in the danger colour.
  const shadow = await specimen.locator(".vt-field-invalid .bp3-input").evaluate((el) => getComputedStyle(el).boxShadow);
  expect(shadow).toContain("rgb(180, 35, 24)");
  await expectNoHorizontalOverflow(specimen);
  await specimen.screenshot({ path: "test-results/design/specimen.png", animations: "disabled" });
});

test("the app icon is served next to the page", async ({ page, request }) => {
  await page.goto(`/?token=${TOKEN}`);
  const href = await page.locator('link[rel="icon"][type="image/svg+xml"]').getAttribute("href");
  expect(href).toBe("/favicon.svg");
  const svg = await request.get(href!);
  expect(svg.status()).toBe(200);
  expect(await svg.text()).toContain("<svg");
  const ico = await request.get("/vettore.ico");
  expect(ico.status()).toBe(200);
  // ICO header: reserved 0, type 1 (icon).
  expect([...(await ico.body()).subarray(0, 4)]).toEqual([0, 0, 1, 0]);
});
