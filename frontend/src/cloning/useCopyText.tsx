// Copy buttons in dialogs and views (primers to order, difference lists, history steps) go
// through the clipboard helper. When the system refuses the clipboard, the text opens selected
// in a box to copy by hand, so a button never claims to have copied what it did not.
import { useCallback, useState } from "react";
import { writeClipboard } from "../files/clipboard";
import CopyTextDialog from "../files/CopyTextDialog";

export function useCopyText() {
  const [fallback, setFallback] = useState<{ title: string; text: string } | null>(null);

  /** True only when the text is on the clipboard: confirm the copy then, and only then. */
  const copy = useCallback(async (text: string, title: string): Promise<boolean> => {
    if (await writeClipboard(text)) return true;
    setFallback({ title, text });
    return false;
  }, []);

  /** Render this next to the button: the copy-by-hand box, while it is open. */
  const fallbackDialog = fallback ? (
    // React passes events from a portal up the component tree: without this, Escape in the box
    // would also close the dialog that holds the copy button.
    <span onKeyDown={(e) => e.stopPropagation()}>
      <CopyTextDialog title={fallback.title} text={fallback.text} onClose={() => setFallback(null)} />
    </span>
  ) : null;

  return { copy, fallbackDialog };
}
