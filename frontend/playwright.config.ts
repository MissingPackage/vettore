import { defineConfig } from "@playwright/test";

const PORT = 8765;
export const TOKEN = "e2e-token";

export default defineConfig({
  testDir: "e2e",
  timeout: 60_000,
  use: {
    baseURL: `http://127.0.0.1:${PORT}`,
    viewport: { width: 1400, height: 900 },
    // Use a pre-installed Chromium when the pinned browser build is not downloaded.
    launchOptions: process.env.PW_CHROMIUM ? { executablePath: process.env.PW_CHROMIUM } : {},
  },
  webServer: {
    command: `${process.platform === "win32" ? "..\\.venv\\Scripts\\python" : "../.venv/bin/python"} -m vettore.app --no-window --port ${PORT} --token ${TOKEN}`,
    url: `http://127.0.0.1:${PORT}/index.html`,
    reuseExistingServer: false,
    timeout: 60_000,
  },
});
