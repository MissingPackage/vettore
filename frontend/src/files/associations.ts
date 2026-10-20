// What File → Preferences says about opening sequence files with Vettore on double-click
// (src/vettore/desktop/associations.py). Pure functions, tested in associations.test.ts.
import type { FileAssociations } from "../api";

/** ".dna, .gb and .ab1": the extensions as a sentence. */
export function extensionList(extensions: string[]): string {
  if (extensions.length <= 1) return extensions.join("");
  return `${extensions.slice(0, -1).join(", ")} and ${extensions[extensions.length - 1]}`;
}

/** One sentence: what the option does and how to undo it. */
export function associationHint(status: FileAssociations): string {
  const account = status.platform === "win32" ? "your Windows account" : "your user account";
  return `Applies to ${extensionList(status.extensions)} files, for ${account} only; `
    + "untick it to undo, and they open as they did before.";
}

/** The files are registered to a Vettore in another folder (this copy was moved or replaced). */
export function movedNote(status: FileAssociations): string | null {
  if (!status.registered_executable || status.current) return null;
  return `Sequence files open with the Vettore at ${status.registered_executable}, not with this `
    + "copy. Tick the box to use this copy instead.";
}

/** Extensions Windows keeps for a program the user chose with "Always use this app". */
export function overriddenNote(status: FileAssociations): string | null {
  if (!status.enabled || !status.overridden.length) return null;
  const list = extensionList(status.overridden);
  const they = status.overridden.length === 1 ? "it still opens" : "they still open";
  return `Windows was told to always open ${list} files with another program, so ${they} there. `
    + "To change that, right-click such a file, choose Open with → Choose another app, pick "
    + "Vettore and tick “Always use this app”.";
}
