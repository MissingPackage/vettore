import { useLayoutEffect, useState, type ReactNode } from "react";
import { createPortal } from "react-dom";

type Placement = "afterToolbar" | "overlay";

const SLOT_CLASS: Record<Placement, string> = {
  afterToolbar: "vt-editor-slot",
  overlay: "vt-editor-overlay",
};

function findOrCreateSlot(editorName: string, placement: Placement): HTMLElement | null {
  const editor = document.querySelector<HTMLElement>(`.veEditor.${editorName}`);
  if (!editor) return null;
  const anchor = editor.querySelector<HTMLElement>(
    placement === "afterToolbar" ? ":scope > .veToolbar-outer" : ":scope > .tg-editor-container",
  );
  if (!anchor) return null;
  const parent = placement === "afterToolbar" ? editor : anchor;
  let slot = parent.querySelector<HTMLElement>(`:scope > .${SLOT_CLASS[placement]}`);
  if (!slot) {
    slot = document.createElement("div");
    slot.className = SLOT_CLASS[placement];
    if (placement === "afterToolbar") anchor.after(slot);
    else anchor.appendChild(slot);
  }
  return slot;
}

/**
 * Renders children inside the OVE editor.
 * - "afterToolbar": between OVE's menu/toolbar and its panels, so the window reads like
 *   a desktop editor (menu, toolbar, document tabs, view tabs, content). OVE lays the editor out as a
 *   flex column, so the panels shrink to make room.
 * - "overlay": on top of the panel area, for Vettore's own views (Enzymes, History) while the
 *   OVE menu and toolbar stay available.
 */
export default function EditorSlot({
  editorName,
  placement = "afterToolbar",
  children,
}: {
  editorName: string;
  placement?: Placement;
  children: ReactNode;
}) {
  const [host, setHost] = useState<HTMLElement | null>(null);

  useLayoutEffect(() => {
    let slot: HTMLElement | null = null;
    const attach = () => {
      slot = findOrCreateSlot(editorName, placement);
      if (slot) setHost(slot);
      return !!slot;
    };
    if (attach()) return () => slot?.remove();
    const observer = new MutationObserver(() => attach() && observer.disconnect());
    observer.observe(document.body, { childList: true, subtree: true });
    return () => {
      observer.disconnect();
      slot?.remove();
    };
  }, [editorName, placement]);

  return host ? createPortal(children, host) : null;
}
