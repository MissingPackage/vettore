// Which documents would lose work if closed, and names that tell new documents apart.
// Pure functions (no editor store), re-exported by documents.ts.

/** The editor's sequence changed since it was opened or last saved. */
export function isDirty(editorState: any): boolean {
  const id = editorState?.sequenceData?.stateTrackingId;
  return !!id && id !== "initialLoadId" && id !== editorState?.lastSavedId;
}

/**
 * Closing the document would lose work: it was edited since it was opened or saved, or it has
 * never been saved at all (a product, a paste, a new file). An untouched sample or opened file
 * does not count.
 */
export function needsSave(doc: { neverSaved: boolean }, editorState: any): boolean {
  return doc.neverSaved || isDirty(editorState);
}

/** `name`, or `name 2`, `name 3`... — the first that no name in `taken` already uses. */
export function uniqueName(name: string, taken: Iterable<string>): string {
  const used = new Set(taken);
  if (!used.has(name)) return name;
  let n = 2;
  while (used.has(`${name} ${n}`)) n++;
  return `${name} ${n}`;
}
