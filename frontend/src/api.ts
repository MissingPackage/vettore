// Bridge to the local Python engine (HTTP on 127.0.0.1, guarded by a session token).

export interface OveAnnotation {
  id?: string;
  name: string;
  type?: string;
  start: number;
  end: number;
  strand?: number;
  forward?: boolean;
  color?: string;
  locations?: { start: number; end: number }[];
  notes?: Record<string, string[]>;
  /** Full primer sequence 5'->3' (primers only), including 5' tails. */
  bases?: string;
}

export interface SequenceData {
  name: string;
  description?: string;
  sequence: string;
  circular: boolean;
  features: OveAnnotation[] | Record<string, OveAnnotation>;
  primers: OveAnnotation[] | Record<string, OveAnnotation>;
}

export interface OpenResult {
  path: string | null;
  documents: SequenceData[];
}

interface DesktopApi {
  open_dialog(): Promise<string[]>;
  save_dialog(suggestedName: string): Promise<string | null>;
}

declare global {
  interface Window {
    pywebview?: { api: DesktopApi };
  }
}

const params = new URLSearchParams(window.location.search);
const token = params.get("token") ?? sessionStorage.getItem("vettore-token") ?? "";
if (token) sessionStorage.setItem("vettore-token", token);

async function call<T>(path: string, body?: unknown): Promise<T> {
  const res = await fetch(`/api/${path}`, {
    method: body === undefined ? "GET" : "POST",
    headers: { "Content-Type": "application/json", "X-Vettore-Token": token },
    body: body === undefined ? undefined : JSON.stringify(body),
  });
  const payload = await res.json().catch(() => ({}));
  if (!res.ok) throw new Error(payload.detail ?? `${res.status} ${res.statusText}`);
  return payload as T;
}

export function desktop(): DesktopApi | null {
  return window.pywebview?.api ?? null;
}

/** Wait briefly for pywebview to inject its bridge (it arrives after page load). */
export function whenDesktopReady(timeoutMs = 1500): Promise<DesktopApi | null> {
  if (desktop()) return Promise.resolve(desktop());
  return new Promise((resolve) => {
    const done = () => resolve(desktop());
    window.addEventListener("pywebviewready", done, { once: true });
    setTimeout(done, timeoutMs);
  });
}

async function toBase64(blob: Blob): Promise<string> {
  const bytes = new Uint8Array(await blob.arrayBuffer());
  let binary = "";
  for (let i = 0; i < bytes.length; i += 0x8000) {
    binary += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
  }
  return btoa(binary);
}

export const api = {
  cloning: <T>(endpoint: string, body: object) => call<T>(`cloning/${endpoint}`, body),
  gel: <T>(body: object) => call<T>("gel", body),
  sangerAlign: <T>(body: object) => call<T>("sanger/align", body),
  gelLadders: () => call<{ ladders: { id: string; label: string }[]; percentages: number[];
    default_ladder: string; default_percent: number }>("gel/ladders"),
  recent: () => call<{ paths: string[] }>("recent"),
  unsaved: (count: number) => call<{ count: number }>("unsaved", { count }),
  tm: (sequence: string) =>
    call<{ tm: number; method: string; presets: { key: string; label: string; description: string; tm: number }[] }>(
      "tm",
      { sequence },
    ),
  saveBinary: async (path: string, blob: Blob) =>
    call<{ path: string }>("save-binary", { path, content_b64: await toBase64(blob) }),
  startupFiles: () => call<{ results: OpenResult[]; errors: string[] }>("startup-files"),
  openPath: (path: string) => call<OpenResult>("open-path", { path }),
  openContent: async (file: File) =>
    call<OpenResult>("open-content", { filename: file.name, content_b64: await toBase64(file) }),
  savePath: (path: string, sequenceData: SequenceData) =>
    call<{ path: string }>("save-path", { path, sequenceData }),
  exportText: (sequenceData: SequenceData, format: "genbank" | "fasta") =>
    call<{ content: string }>("export", { sequenceData, format }),
};
