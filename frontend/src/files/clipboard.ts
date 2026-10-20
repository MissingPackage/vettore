// Clipboard access with graceful failure: the desktop window or the browser may refuse it,
// in which case the caller shows a box the user can paste into or copy from by hand.

/** The clipboard text, or null when reading is not permitted or not available. */
export async function readClipboard(): Promise<string | null> {
  try {
    if (!navigator.clipboard?.readText) return null;
    return await navigator.clipboard.readText();
  } catch {
    return null;
  }
}

/** Put text on the clipboard; false when neither the Clipboard API nor the copy command works. */
export async function writeClipboard(text: string): Promise<boolean> {
  try {
    if (navigator.clipboard?.writeText) {
      await navigator.clipboard.writeText(text);
      return true;
    }
  } catch {
    // fall through to the copy command
  }
  const area = document.createElement("textarea");
  area.value = text;
  area.setAttribute("readonly", "");
  area.style.position = "fixed";
  area.style.opacity = "0";
  document.body.appendChild(area);
  area.select();
  try {
    return document.execCommand("copy");
  } catch {
    return false;
  } finally {
    area.remove();
  }
}
