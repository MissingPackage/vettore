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
  /** Engine format id ("genbank", "embl", "fasta", "raw", ...) and its label for people. */
  format?: string;
  format_label?: string;
  /** Things the user should know, e.g. features that could not be placed. */
  warnings?: string[];
  /** Save may overwrite the file: it holds exactly this one document. */
  writable?: boolean;
}

export interface SaveResult {
  /** Where the file was written (an extension is appended when none was typed). */
  path: string;
  format: string;
  format_label: string;
}

/** A sample construct shipped with Vettore (offered on the welcome screen). */
export interface Sample {
  id: string;
  label: string;
  description: string;
}

interface DesktopApi {
  open_dialog(): Promise<string[]>;
  /** `kind` (genbank, fasta, dna, svg, png) orders the file types; `directory` is where it starts. */
  save_dialog(suggestedName: string, directory?: string, kind?: string): Promise<string | null>;
}

declare global {
  interface Window {
    pywebview?: { api: DesktopApi };
  }
}

const params = new URLSearchParams(window.location.search);
const token = params.get("token") ?? sessionStorage.getItem("vettore-token") ?? "";
if (token) sessionStorage.setItem("vettore-token", token);

/** Call an engine endpoint under /api/; throws with the engine's message on error. */
export async function call<T>(path: string, body?: unknown): Promise<T> {
  let res: Response;
  try {
    res = await fetch(`/api/${path}`, {
      method: body === undefined ? "GET" : "POST",
      headers: { "Content-Type": "application/json", "X-Vettore-Token": token },
      body: body === undefined ? undefined : JSON.stringify(body),
    });
  } catch {
    throw new Error("Vettore's engine is not responding. Close Vettore and start it again; "
      + "your files on disk are not affected.");
  }
  const payload = await res.json().catch(() => ({}));
  if (!res.ok) {
    // The engine's messages are written for people; anything else gets a plain fallback.
    if (typeof payload.detail === "string") throw new Error(payload.detail);
    throw new Error(`Something went wrong inside Vettore (error ${res.status}). The details were `
      + "written to Vettore's log file; please send it if this happens again.");
  }
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

export async function toBase64(blob: Blob): Promise<string> {
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
  detectFeatures: <T>(body: object) => call<T>("features/detect", body),
  gelLadders: () => call<{ ladders: { id: string; label: string }[]; percentages: number[];
    default_ladder: string; default_percent: number }>("gel/ladders"),
  recent: () => call<{ paths: string[] }>("recent"),
  samples: () => call<{ samples: Sample[] }>("app/samples"),
  openSample: (id: string) => call<OpenResult>(`app/samples/${encodeURIComponent(id)}`),
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
  /** Pasted text (sequence, FASTA, GenBank or EMBL record), recognised like file content. */
  openText: (text: string, name = "") => call<OpenResult>("open-text", { text, name }),
  /** The engine writes the format of the path's extension; `format` applies when there is none. */
  savePath: (path: string, sequenceData: SequenceData, format?: "genbank" | "fasta" | "dna") =>
    call<SaveResult>("save-path", { path, sequenceData, format }),
  exportText: (sequenceData: SequenceData, format: "genbank" | "fasta") =>
    call<{ content: string }>("export", { sequenceData, format }),
  /** File content for a browser download: text for GenBank/FASTA, base64 for .dna. */
  exportFile: (sequenceData: SequenceData, format: "genbank" | "fasta" | "dna") =>
    call<{ content?: string; content_b64?: string }>("export", { sequenceData, format }),
};
