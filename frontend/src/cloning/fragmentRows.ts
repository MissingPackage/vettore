// State of one ligation fragment in a form, and the engine body it becomes (FragmentIn in
// cloning/ligation.py). Kept free of React and of the HTTP bridge so it can be unit-tested.
import { docData, selectionOf } from "./defaults";
import type { OpenDoc } from "./types";

export type FragmentMode = "digest" | "pcr";

export interface FragmentRow {
  editorName: string | null;
  mode: FragmentMode;
  enzymes: string[];
  /** Digest fragment chosen by the user (engine numbering), or null = automatic. */
  fragmentIndex: number | null;
  start: number | null;
  end: number | null;
}

const byName = (docs: OpenDoc[], editorName: string | null) =>
  docs.find((d) => d.editorName === editorName) ?? null;

/** A new row for a document: circular files are cut, linear ones used as the form prefers
 * (as they are for ligations, cut for restriction cloning). */
export function newRow(docs: OpenDoc[], editorName: string | null, prefer: FragmentMode): FragmentRow {
  const d = byName(docs, editorName);
  const circular = !!d?.data()?.circular;
  const mode: FragmentMode = circular ? "digest" : prefer;
  const region = mode === "pcr" ? selectionOf(d) : { start: null, end: null };
  return { editorName, mode, enzymes: [], fragmentIndex: null, ...region };
}

/** The engine body for a row, or a sentence saying what is missing. */
export function rowBody(docs: OpenDoc[], row: FragmentRow, what: string): object | string {
  const d = byName(docs, row.editorName);
  if (!d) return `Choose the file for the ${what}`;
  if (row.mode === "digest") {
    if (!row.enzymes.length) return `Choose the enzymes that cut “${d.name}” (${what})`;
    return { doc: docData(d), enzymes: row.enzymes, fragment_index: row.fragmentIndex };
  }
  if ((row.start === null) !== (row.end === null)) return `Enter both ends of the region of “${d.name}”`;
  if (row.start === null && d.data().circular) {
    return `“${d.name}” is circular: choose the region to use as a PCR product, or cut it with enzymes`;
  }
  return { doc: docData(d), enzymes: [], start: row.start, end: row.end };
}
