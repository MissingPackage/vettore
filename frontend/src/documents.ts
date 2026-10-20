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

export interface Doc {
  editorName: string;
  name: string;
  /** File the document was opened from or last saved to. */
  path: string | null;
  view: ViewId;
  enzymeSet: EnzymeSetId;
  history: HistoryEntry[];
}

const WRITABLE = /\.(gb|gbk|gbff|genbank|ape|fa|fasta|fna|fas)$/i;
let counter = 0;

export function canSaveInPlace(doc: Doc): boolean {
  return !!doc.path && WRITABLE.test(doc.path);
}

export function now(): string {
  return new Date().toLocaleString();
}

export function createDocument(
  sequenceData: SequenceData,
  path: string | null,
  firstEntry: HistoryEntry,
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

export function isDirty(editorState: any): boolean {
  const id = editorState?.sequenceData?.stateTrackingId;
  return !!id && id !== "initialLoadId" && id !== editorState?.lastSavedId;
}

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
