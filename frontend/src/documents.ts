// Open documents: one OVE editor instance (editorName) per document.
import { actions, updateEditor } from "@teselagen/ove";
import type { SequenceData } from "./api";
import { store } from "./store";
import type { EnzymeSetId } from "./enzymeSets";
import { loadPreferences } from "./preferences";
import { panelsFor, type ViewId } from "./views";

export interface HistoryEntry {
  time: string;
  action: string;
  detail?: string;
  /** Snapshots of the input documents of a cloning step, so they can be reopened later. */
  parents?: { name: string; data: SequenceData }[];
}

/** The file a document was opened from, as the engine described it. */
export interface DocSource {
  /** Format label, e.g. "GenBank", ".dna", "EMBL", "plain sequence text". */
  label: string;
  /** How many sequences the file held (one tab each). */
  records: number;
  fileName: string;
  /** Some of the file's content (features) could not be shown. */
  incomplete: boolean;
}

export interface Doc {
  editorName: string;
  name: string;
  /** File the document was opened from or last saved to. */
  path: string | null;
  /**
   * Save may overwrite `path`: Vettore wrote it, or it is a one-record GenBank/FASTA file that
   * the document holds completely. Otherwise Save makes a GenBank copy.
   */
  saveInPlace: boolean;
  source: DocSource | null;
  /** The "Save a GenBank copy" explanation was already shown for this document. */
  copyNoticeShown: boolean;
  /**
   * The document exists only in Vettore and has never been saved: a cloning product, a pasted
   * sequence, a new file with a sequence, a tool's result. It counts as unsaved until it is.
   */
  neverSaved: boolean;
  view: ViewId;
  enzymeSet: EnzymeSetId;
  history: HistoryEntry[];
}

let counter = 0;

export function canSaveInPlace(doc: Doc): boolean {
  return !!doc.path && doc.saveInPlace;
}

export function now(): string {
  return new Date().toLocaleString();
}

export function createDocument(
  sequenceData: SequenceData,
  path: string | null,
  firstEntry: HistoryEntry,
  options?: { saveInPlace?: boolean; source?: DocSource; neverSaved?: boolean },
): Doc {
  const editorName = `doc-${++counter}`;
  updateEditor(store, editorName, {
    sequenceData,
    sequenceDataHistory: {},
    readOnly: false,
    panelsShown: panelsFor(sequenceData.circular),
  });
  return {
    editorName,
    name: sequenceData.name,
    path,
    saveInPlace: options?.saveInPlace ?? false,
    source: options?.source ?? null,
    copyNoticeShown: false,
    neverSaved: options?.neverSaved ?? false,
    view: "map",
    enzymeSet: loadPreferences().enzymeSet,
    history: [firstEntry],
  };
}

export function removeDocument(editorName: string) {
  // Reset the editor slice so a reopened document starts clean.
  updateEditor(store, editorName, { sequenceData: { sequence: "" }, sequenceDataHistory: {} });
}

export function editorState(editorName: string) {
  return (store.getState() as { VectorEditor: Record<string, any> }).VectorEditor[editorName];
}

export { isDirty, needsSave, uniqueName } from "./files/unsaved";

export function markSaved(editorName: string) {
  const id = editorState(editorName)?.sequenceData?.stateTrackingId;
  store.dispatch(actions.lastSavedIdUpdate(id ?? "initialLoadId", { editorName }));
}

export function baseName(path: string): string {
  return path.split(/[\\/]/).pop() ?? path;
}

export function stem(path: string): string {
  return baseName(path).replace(/\.[^.]+$/, "");
}
