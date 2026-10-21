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

/** Limits on a History read from a file. Real work never meets them: a cloning has a handful
 * of inputs, a construct a few dozen steps, and a step's description is a sentence. A file
 * that does meet them is read in part rather than shown as it asks to be. */
const MAX_STEPS = 500;
const MAX_PARENTS = 100;
const MAX_TEXT = 4000;

/** A string from a file, cut to a length a view can draw. */
function text(value: unknown, fallback = ""): string {
  return typeof value === "string" ? value.slice(0, MAX_TEXT) : fallback;
}

/**
 * A snapshot kept only when it is shaped like a document the editor can open: a sequence, and
 * features and primers in the shape the editor reads. A snapshot that is not goes to the
 * editor's store the moment someone opens it, so it is dropped here instead.
 */
function snapshot(value: unknown): SequenceData | undefined {
  if (!value || typeof value !== "object") return undefined;
  const data = value as Record<string, unknown>;
  if (typeof data.sequence !== "string") return undefined;
  const annotations = (kind: unknown) =>
    Array.isArray(kind) || (!!kind && typeof kind === "object");
  if (!annotations(data.features) || !annotations(data.primers)) return undefined;
  if (typeof data.name !== "string") return undefined;
  return data as unknown as SequenceData;
}

/** The wrong clones of a step, kept only where each one is shaped as the planner reads it. */
function wrongClones(value: unknown): WrongClone[] | undefined {
  if (!Array.isArray(value) || !value.length || value.length > MAX_PARENTS) return undefined;
  const kept = value.filter((clone) => {
    if (!clone || typeof clone !== "object") return false;
    const { name, data } = clone as Record<string, unknown>;
    return typeof name === "string" && !!snapshot(data);
  });
  return kept.length === value.length ? (kept as WrongClone[]) : undefined;
}

/**
 * The History of a document as it was read from a file: a file is outside content, so only
 * entries that look like History are kept, each within the limits above, and the tree is
 * followed to `LINEAGE_DEPTH`.
 */
export function cleanHistory(stored: unknown, depth = LINEAGE_DEPTH): HistoryEntry[] {
  if (!Array.isArray(stored) || depth <= 0) return [];
  const out: HistoryEntry[] = [];
  for (const raw of stored.slice(0, MAX_STEPS)) {
    if (!raw || typeof raw !== "object") continue;
    const step = raw as Record<string, unknown>;
    if (typeof step.action !== "string" || !step.action) continue;
    const entry: HistoryEntry = { time: text(step.time), action: text(step.action) };
    if (typeof step.detail === "string") entry.detail = text(step.detail);
    const clones = wrongClones(step.wrongClones);
    if (clones) entry.wrongClones = clones;
    if (Array.isArray(step.parents)) {
      const parents: HistoryParent[] = [];
      for (const rawParent of step.parents.slice(0, MAX_PARENTS)) {
        if (!rawParent || typeof rawParent !== "object") continue;
        const parent = rawParent as Record<string, unknown>;
        const data = snapshot(parent.data);
        parents.push({
          name: text(parent.name, "an input"),
          ...(data ? { data } : {}),
          ...(typeof parent.path === "string" ? { path: text(parent.path) } : {}),
          ...(parent.snapshots_dropped || !data ? { snapshots_dropped: true } : {}),
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
