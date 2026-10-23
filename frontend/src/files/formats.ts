// File formats as the user meets them: what Open accepts and what Save writes.

/**
 * Extensions the engine reads, for the browser's file chooser. Keep in sync with FORMATS in
 * src/vettore/fileio.py (tests/test_fileio_formats.py checks it). The engine recognises files
 * by their content, so a file with another extension still opens if the user picks it.
 */
export const OPEN_EXTENSIONS = [
  ".dna", ".gb", ".gbk", ".gbff", ".genbank", ".ape", ".embl", ".emb",
  ".fa", ".fasta", ".fna", ".fas", ".fsa", ".xdna", ".gck", ".ab1", ".abi", ".txt", ".seq",
];

export const ACCEPT = OPEN_EXTENSIONS.join(",");

/** Formats Vettore writes; the engine picks the format from the extension of the path. */
export type SaveKind = "genbank" | "fasta" | "dna";

export const SAVE_EXTENSION: Record<SaveKind, string> = {
  genbank: ".gb",
  fasta: ".fasta",
  dna: ".dna",
};

/** A document name made safe as a file name on Windows and Linux. */
export function safeFileName(name: string): string {
  const cleaned = name.replace(/[\\/:*?"<>|\u0000-\u001f]+/g, "_").replace(/^[\s.]+|[\s.]+$/g, "");
  return cleaned || "Untitled";
}

export function dirName(path: string): string {
  const i = Math.max(path.lastIndexOf("/"), path.lastIndexOf("\\"));
  return i > 0 ? path.slice(0, i) : "";
}

/** The file name of a path, on either platform's separator. */
export function baseName(path: string): string {
  return path.split(/[\\/]/).pop() ?? path;
}

/** The file name without its extension: the name a document opened from a file takes. */
export function stem(path: string): string {
  return baseName(path).replace(/\.[^.]+$/, "");
}
