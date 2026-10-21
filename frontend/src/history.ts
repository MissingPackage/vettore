// The History of a document: what it is made of, step by step, with a copy of everything each
// step was run on. Kept apart from `documents.ts` (and free of the editor) because a saved
// file carries the History (S6.1), so reading one is reading outside content.
import type { SequenceData } from "./api";
import type { WrongClone } from "./cloning/types";
import type { Doc } from "./documents";

/** An input document of a step, kept as it was when the step ran. */
export interface HistoryParent {
  name: string;
  /** The input exactly as it was; missing when a saved file could not keep the snapshot. */
  data?: SequenceData;
  /** The file the input was open from, so the Lineage view can offer to open the file itself. */
  path?: string | null;
  /** The input's own History, which is what makes the lineage a tree rather than one step. */
  history?: HistoryEntry[];
  /** The saved file left this snapshot out because the History did not fit (engine's flag). */
  snapshots_dropped?: boolean;
}

export interface HistoryEntry {
  time: string;
  action: string;
  detail?: string;
  /** Snapshots of the input documents of a cloning step, so they can be reopened later. */
  parents?: HistoryParent[];
  /** The other outcomes of that cloning step (other products, reverse insert, empty vector),
   * offered by Find Diagnostic Digest as the wrong clones to tell apart. */
  wrongClones?: WrongClone[];
}

/** How deep the Lineage view follows a tree: the spec's "readable for 20 steps". */
export const LINEAGE_DEPTH = 20;

/**
 * The History of a document as it was read from a file: a file is outside content, so only
 * entries that look like History are kept, and the tree is followed to `LINEAGE_DEPTH`.
 */
export function cleanHistory(stored: unknown, depth = LINEAGE_DEPTH): HistoryEntry[] {
  if (!Array.isArray(stored) || depth <= 0) return [];
  const out: HistoryEntry[] = [];
  for (const raw of stored) {
    if (!raw || typeof raw !== "object") continue;
    const step = raw as Record<string, unknown>;
    if (typeof step.action !== "string" || !step.action) continue;
    const entry: HistoryEntry = {
      time: typeof step.time === "string" ? step.time : "",
      action: step.action,
    };
    if (typeof step.detail === "string") entry.detail = step.detail;
    if (Array.isArray(step.wrongClones)) entry.wrongClones = step.wrongClones as WrongClone[];
    if (Array.isArray(step.parents)) {
      const parents: HistoryParent[] = [];
      for (const rawParent of step.parents) {
        if (!rawParent || typeof rawParent !== "object") continue;
        const parent = rawParent as Record<string, unknown>;
        const data = parent.data as SequenceData | undefined;
        parents.push({
          name: typeof parent.name === "string" ? parent.name : "an input",
          ...(data && typeof data === "object" && typeof data.sequence === "string"
            ? { data } : {}),
          ...(typeof parent.path === "string" ? { path: parent.path } : {}),
          ...(parent.snapshots_dropped ? { snapshots_dropped: true } : {}),
          history: cleanHistory(parent.history, depth - 1),
        });
      }
      if (parents.length) entry.parents = parents;
    }
    out.push(entry);
  }
  return out;
}

/**
 * A document just read from a file, with the History the file carried taken off it: the
 * History belongs to the Doc, so the editor never holds a second, ageing copy of it.
 */
export function takeHistory(data: SequenceData): {
  data: SequenceData;
  history: HistoryEntry[];
} {
  if (!data.history) return { data, history: [] };
  const { history, ...rest } = data;
  return { data: rest as SequenceData, history: cleanHistory(history) };
}

/** What to write: the document, with its History when the preference keeps it in the file. */
export function dataToSave(data: SequenceData, doc: Doc, keepHistory: boolean): SequenceData {
  return keepHistory && doc.history.length ? { ...data, history: doc.history } : data;
}

/** A document as an input of a step: its content, where it came from, and its own History. */
export function snapshotOf(doc: Doc, data: SequenceData): HistoryParent {
  return { name: doc.name, data, path: doc.path, history: doc.history };
}
