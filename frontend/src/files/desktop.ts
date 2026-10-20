// Messages between the page and the desktop window (src/vettore/desktop/): closing the window
// with unsaved work, files handed over by a second start of Vettore, print requests the
// window cannot serve.

/**
 * Dispatched on `window` by the desktop window when it is closed while documents have unsaved
 * changes (src/vettore/desktop/close_guard.py). It is cancelable: the page cancels it to say it
 * will ask the user; a page that does not, within a few seconds, gets the system's own question.
 */
export const QUIT_EVENT = "vettore:quit-requested";

/** Dispatched on `window` when a second start handed files over (src/vettore/desktop/instance.py). */
export const FILES_EVENT = "vettore:files-handed-over";

/** How often the page asks the engine for files handed over by a second start (ms). */
export const STARTUP_POLL_MS = 500;

/**
 * Dispatched on `window` when the desktop window dropped a print request: the Linux window
 * (Qt WebEngine) does not print (src/vettore/desktop/printing.py).
 */
export const PRINT_EVENT = "vettore:print-unavailable";

/** What to do instead of printing from the Linux window. */
export const PRINT_UNAVAILABLE = "Vettore's window cannot print on Linux yet, so nothing was sent "
  + "to the printer. To print the map, use File → Export → Map as PNG… and print the image from "
  + "your image viewer, or save it as a PDF from there.";

/** What the window's unsaved-changes dialog says about each document. */
export interface QuitDoc {
  editorName: string;
  name: string;
  path: string | null;
  /** Save writes back to `path` (see Doc.saveInPlace). */
  saveInPlace: boolean;
}

/** The first line of the dialog. */
export function quitHeadline(docs: QuitDoc[]): string {
  return docs.length === 1
    ? `“${docs[0].name}” has changes that are not saved.`
    : `${docs.length} documents have changes that are not saved:`;
}

/** Where Save puts each document, in a few words: the same thing Save would do. */
export function quitDetail(doc: QuitDoc): string {
  const fileName = doc.path?.split(/[\\/]/).pop();
  if (!fileName) return "never saved to a file: you choose where";
  return doc.saveInPlace ? `saved to ${fileName}` : `saved as a GenBank copy of ${fileName}`;
}

/** The message when Save All stopped at a document, so the window stayed open. */
export function notSavedMessage(name: string): string {
  return `Vettore is still open because “${name}” was not saved. Save it with File → Save, `
    + "or close the window again and choose Don't Save to discard the changes.";
}
