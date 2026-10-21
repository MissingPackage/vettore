// Pure helpers for the lab primer library (Primers → Primer Library…). The engine
// (src/vettore/primer_library.py) owns the binding rule, the shared .csv and the products;
// this module holds the view-side types, text and the small client-side checks the dialog needs
// so it can show what an import will keep before anything is sent.

export interface LibraryPrimer {
  name: string;
  sequence: string;
  notes: string;
  location: string;
  added: string;
  length: number;
  /** Tm (°C) of a plain-DNA primer; null for one with ambiguity codes. */
  tm: number | null;
  /** The sequence has IUPAC codes (a degenerate primer). */
  degenerate?: boolean;
  /** Long enough to be looked for on sequences (at least 15 nt). */
  searched?: boolean;
}

export interface RejectedRow {
  name: string;
  sequence: string;
  reason: string;
}

export interface Binding {
  start: number;
  end: number;
  strand: number;
  footprint: number;
  tail: number;
  mismatches: number;
}

/** A library primer together with where it binds the open document. */
export interface BoundPrimer extends LibraryPrimer {
  /** Where it binds; empty when it binds too many places to draw (see `sites`). */
  bindings: Binding[];
  /** How many places it binds. */
  sites: number;
}

/** Today as YYYY-MM-DD in the user's own time zone, for "date added". */
export function today(now = new Date()): string {
  const pad = (n: number) => String(n).padStart(2, "0");
  return `${now.getFullYear()}-${pad(now.getMonth() + 1)}-${pad(now.getDate())}`;
}

/** One sentence for the rows of a file or a paste that were left out, with why. */
export function rejectedText(rejected: RejectedRow[], where: "file" | "paste"): string {
  const shown = rejected.slice(0, 3).map((r) => `“${r.name || r.sequence || "a row"}”: ${r.reason}`);
  const more = rejected.length > 3 ? ` (and ${rejected.length - 3} more)` : "";
  const n = rejected.length;
  const rows = `${n} row${n === 1 ? "" : "s"}`;
  return where === "file"
    ? `${rows} of the library file could not be read and stay in the file as they are, so a colleague can correct them. ${shown.join(" ")}${more}`
    : `${rows} left out: ${shown.join(" ")}${more}`;
}

export interface LibraryFile {
  path: string;
  exists: boolean;
  mtime: number | null;
  /** The file's separator (";" when Excel is set to Italian), kept for Export. */
  delimiter: string;
  primers: LibraryPrimer[];
  rejected: RejectedRow[];
}

/** Note key the engine marks map-only overlay primers with (src/vettore/fileio.py). */
export const OVERLAY_NOTE = "vettore_overlay";

const DNA = /^[ACGTRYSWKMBDHVN]+$/;

/** A primer sequence cleaned for comparison, as the engine does (primer_library.clean_sequence):
 * upper case, U→T, no spaces or dashes, a 5′-…-3′ frame dropped. "" when not DNA. */
export function cleanSequence(raw: string): string {
  let seq = (raw || "").toUpperCase().replace(/U/g, "T").replace(/\s+/g, "")
    .replace(/[\u2019\u2032]/g, "'").replace(/-/g, "");
  if (seq.startsWith("5'") && seq.endsWith("3'")) seq = seq.slice(2, -2);
  return DNA.test(seq) ? seq : "";
}

/** The file's base name, for showing which file the library lives in. */
export function baseName(path: string): string {
  return path.split(/[\\/]/).pop() || path;
}

/** Library primers whose name, sequence, notes or location contain the query (case-insensitive).
 * An empty query keeps every primer. */
export function filterPrimers<T extends LibraryPrimer>(primers: T[], query: string): T[] {
  const q = query.trim().toLowerCase();
  if (!q) return primers;
  const seq = cleanSequence(q);
  return primers.filter((p) =>
    p.name.toLowerCase().includes(q) ||
    p.notes.toLowerCase().includes(q) ||
    p.location.toLowerCase().includes(q) ||
    (!!seq && p.sequence.includes(seq)));
}

/** Whether a sequence is already in the library (same oligo, 5'→3'); returns its primer. */
export function sequenceAlreadyIn(primers: LibraryPrimer[], sequence: string): LibraryPrimer | null {
  const seq = cleanSequence(sequence);
  return seq ? primers.find((p) => p.sequence === seq) ?? null : null;
}

/** Range text, 1-based inclusive like the editor, with a wrap note across the origin. */
export function range(start: number, end: number): string {
  const a = start + 1;
  const b = end + 1;
  return a <= b ? `${a.toLocaleString("en-US")}–${b.toLocaleString("en-US")}`
                : `${a.toLocaleString("en-US")}→${b.toLocaleString("en-US")} (across the origin)`;
}

/** One line describing a binding for the dialog: strand, position and any mismatches/tail. */
export function bindingText(b: Binding): string {
  const dir = b.strand === 1 ? "→ forward" : "← reverse";
  const extras: string[] = [];
  if (b.mismatches) extras.push(`${b.mismatches} mismatch${b.mismatches > 1 ? "es" : ""}`);
  if (b.tail) extras.push(`${b.tail} nt 5′ tail`);
  const note = extras.length ? ` (${extras.join(", ")})` : " (exact)";
  return `${dir}, ${range(b.start, b.end)}${note}`;
}

/** The .csv rows copied from the Copy for Ordering format: name⇥sequence per line. */
export function orderingText(primers: { name: string; sequence: string }[]): string {
  return primers.map((p) => `${p.name}\t${p.sequence}`).join("\n");
}

const CSV_COLUMNS = ["name", "sequence", "notes", "location", "added"] as const;

function csvCell(value: string, delimiter: string): string {
  return value.includes(delimiter) || /["\r\n]/.test(value) ? `"${value.replace(/"/g, '""')}"` : value;
}

/** The library as .csv text (header + one row per primer), for Export as .csv, separated by
 * `delimiter` (the library file's own). */
export function toCsv(primers: LibraryPrimer[], delimiter = ","): string {
  const rows = primers.map((p) => CSV_COLUMNS.map((c) => csvCell(String(p[c] ?? ""), delimiter)).join(delimiter));
  return [CSV_COLUMNS.join(delimiter), ...rows].join("\r\n") + "\r\n";
}

/** Add imported primers to the library, skipping any whose sequence is already there. Returns
 * the merged list and how many were added versus skipped as duplicates. */
export function mergeImport(existing: LibraryPrimer[], incoming: LibraryPrimer[]): {
  merged: LibraryPrimer[];
  added: number;
  duplicates: number;
} {
  const have = new Set(existing.map((p) => p.sequence));
  const merged = [...existing];
  let added = 0;
  let duplicates = 0;
  for (const p of incoming) {
    if (have.has(p.sequence)) {
      duplicates += 1;
      continue;
    }
    have.add(p.sequence);
    merged.push(p);
    added += 1;
  }
  return { merged, added, duplicates };
}

/** A primer added to the library now, with today's date (from the engine) and a location. */
export function newEntry(name: string, sequence: string, today: string,
                         notes = "", location = ""): LibraryPrimer {
  const seq = cleanSequence(sequence);
  return { name: name.trim() || "(unnamed)", sequence: seq, notes, location, added: today,
           length: seq.length, tm: null };
}
