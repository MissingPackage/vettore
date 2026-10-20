// Plug-in contract for dialogs reached from the menu bar (analysis tools, help, preferences).
// A module under tools/ exports `tools: Tool[]`; tools/index.ts collects them, the menu bar
// lists them and App renders the open one. Nothing else needs to change to add a tool.
import type { ComponentType } from "react";
import type { SequenceData } from "../api";
import type { OpenDoc } from "../cloning/types";

export type MenuName =
  | "File" | "Edit" | "View" | "Enzymes" | "Features" | "Primers" | "Actions" | "Tools" | "Help";

export interface ToolProps {
  /** Open documents, in tab order. Read their current content with `doc.data()`. */
  docs: OpenDoc[];
  /** The document in the active tab, or null on the welcome screen. */
  active: OpenDoc | null;
  onClose: () => void;
  /** Open a new document tab; `action`/`detail` become its first History entry. */
  openDocument: (data: SequenceData, action: string, detail?: string) => void;
  /** Add a History entry to an open document (e.g. after changing it). */
  logHistory: (editorName: string, action: string, detail?: string) => void;
}

export interface Tool {
  id: string;
  menu: MenuName;
  /** Menu text without the trailing ellipsis (added automatically). */
  label: string;
  /** OVE hotkey syntax, e.g. "mod+shift+a". */
  hotkey?: string;
  /** Disabled in the menu while no document is open. */
  requiresDocument?: boolean;
  Component: ComponentType<ToolProps>;
}
