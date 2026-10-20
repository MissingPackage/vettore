// Double-clicking a sequence file while Vettore runs opens it in the window already open.
// A double-click starts `Vettore <file>`; that process finds the running copy and hands the
// file over (src/vettore/desktop/instance.py). These tests run their own copy of Vettore, with
// its own settings folder, so no other test's page can take the files, and start the second
// process for real.
import { expect, test } from "@playwright/test";
import { spawn, type ChildProcess } from "node:child_process";
import fs from "node:fs";
import net from "node:net";
import os from "node:os";
import path from "node:path";

const data = (name: string) => path.resolve(import.meta.dirname, "../../tests/data", name);
const PYTHON = path.resolve(import.meta.dirname, "..", "..", ".venv",
  process.platform === "win32" ? "Scripts/python.exe" : "bin/python");
const TOKEN = "handover-token";

let home: string;
let dir: string;
let port: number;
let running: ChildProcess;

function freePort(): Promise<number> {
  return new Promise((resolve, reject) => {
    const server = net.createServer();
    server.listen(0, "127.0.0.1", () => {
      const { port: p } = server.address() as net.AddressInfo;
      server.close(() => resolve(p));
    });
    server.on("error", reject);
  });
}

/** `vettore` as a double-click starts it, without a window (the test plays the page). */
function vettore(...args: string[]): ChildProcess {
  return spawn(PYTHON, ["-m", "vettore.app", "--no-window", "--single-instance", ...args], {
    env: { ...process.env, VETTORE_HOME: home },
    stdio: "ignore",
  });
}

function exitCode(proc: ChildProcess): Promise<number | null> {
  return new Promise((resolve) => {
    if (proc.exitCode !== null) resolve(proc.exitCode);
    else proc.on("exit", (code) => resolve(code));
  });
}

function copyOf(fixture: string, name = fixture): string {
  const target = path.join(dir, name);
  fs.copyFileSync(data(fixture), target);
  return target;
}

test.beforeAll(async () => {
  home = fs.mkdtempSync(path.join(os.tmpdir(), "vettore-e2e-home-"));
  dir = fs.mkdtempSync(path.join(os.tmpdir(), "vettore-e2e-handover-"));
  port = await freePort();
  running = vettore("--port", String(port), "--token", TOKEN);
  const deadline = Date.now() + 60_000;
  for (;;) {
    try {
      const res = await fetch(`http://127.0.0.1:${port}/api/health`, { headers: { "X-Vettore-Token": TOKEN } });
      if (res.ok && fs.existsSync(path.join(home, "instance.json"))) break;
    } catch {
      // not up yet
    }
    if (Date.now() > deadline) throw new Error("the running copy of Vettore did not start");
    await new Promise((r) => setTimeout(r, 200));
  }
});

test.afterAll(() => {
  running?.kill();
  fs.rmSync(dir, { recursive: true, force: true });
  fs.rmSync(home, { recursive: true, force: true });
});

test("a file double-clicked while Vettore runs opens as a tab in the open window", async ({ page }) => {
  await page.goto(`http://127.0.0.1:${port}/?token=${TOKEN}`);
  await expect(page.getByTestId("welcome")).toBeVisible();

  const file = copyOf("MT612431.1.gb", "from-explorer.gb");
  const second = vettore(file);
  // The second start exits once the page has taken the file: no second copy keeps running.
  expect(await exitCode(second)).toBe(0);
  const active = page.locator(".vt-doctab.active");
  await expect(active).toContainText("from-explorer");
  await expect(active).toHaveAttribute("title", file);
  await expect(page.getByTestId("status-line")).toContainText("7,204 bp");

  // The same file again: its tab comes forward, it is not opened twice.
  await page.getByTestId("file-input").setInputFiles(data("pFA-KanMX4.dna"));
  await expect(active).toContainText("pFA-KanMX4");
  expect(await exitCode(vettore(file))).toBe(0);
  await expect(active).toContainText("from-explorer");
  await expect(page.getByTestId("doctab")).toHaveCount(2);

  // Several files at once (a multiple selection opened together) all arrive.
  const fasta = copyOf("exchange-plasmids.fasta");
  const dna = copyOf("pFA-KanMX4.dna", "yeast.dna");
  const both = [vettore(fasta), vettore(dna)];
  expect(await Promise.all(both.map(exitCode))).toEqual([0, 0]);
  await expect(page.getByTestId("doctab")).toHaveCount(6);   // the FASTA file holds 3 records
  await expect(page.getByTestId("doctab").filter({ hasText: "yeast" })).toHaveAttribute("title", dna);
});

test("a file that cannot be opened is explained in the open window", async ({ page }) => {
  await page.goto(`http://127.0.0.1:${port}/?token=${TOKEN}`);
  await expect(page.getByTestId("welcome")).toBeVisible();
  const missing = path.join(dir, "moved-away.gb");
  expect(await exitCode(vettore(missing))).toBe(0);
  await expect(page.locator('[data-testid="notice"][data-kind="error"]'))
    .toContainText("“moved-away.gb” was not found");
});
