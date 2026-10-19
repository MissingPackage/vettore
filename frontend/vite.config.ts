import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  define: { global: "globalThis" },
  // Blueprint 3 CSS (used by OVE) is rejected by the lightningcss minifier.
  build: { outDir: "dist", chunkSizeWarningLimit: 8000, cssMinify: false },
  server: { proxy: { "/api": "http://127.0.0.1:8765" } },
});
