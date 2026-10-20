// The editor component ships its own ways to open and save files (toolbar import/export
// buttons, a GenBank copy in the right-click menus, an alignment upload that needs a server).
// They bypass the engine, so Vettore removes or rewires them: File → Open… and File → Save are
// the only ways in and out, and copying text goes through Copy as FASTA.
import "./editorFiles.css";

/** The editor's toolbar without its import, export and alignment-upload buttons. */
export const EDITOR_TOOLS = [
  "saveTool",
  "undoTool",
  "redoTool",
  "cutsiteTool",
  "featureTool",
  "partTool",
  "oligoTool",
  "orfTool",
  "editTool",
  "findTool",
  "visibilityTool",
];

/** Right-click menus that carry the editor's copy options (selection, features, primers, ...). */
const MENUS_WITH_COPY = [
  "selectionLayerRightClicked",
  "searchLayerRightClicked",
  "featureRightClicked",
  "partRightClicked",
  "primerRightClicked",
  "orfRightClicked",
  "translationRightClicked",
  "warningRightClicked",
];

type MenuItem = Record<string, unknown> | string;

/** Replace the editor's own "Copy Genbank For Selection" with Copy as FASTA. */
export function rewireCopyItems(items: MenuItem[], copyAsFasta: () => void): MenuItem[] {
  return items.map((item) => {
    if (typeof item !== "object" || item === null) return item;
    if (item.className === "openVeCopyGenbankForSelection") {
      return { text: "Copy as FASTA", className: "vt-copy-fasta", onClick: copyAsFasta };
    }
    return Array.isArray(item.submenu)
      ? { ...item, submenu: rewireCopyItems(item.submenu as MenuItem[], copyAsFasta) }
      : item;
  });
}

/** `rightClickOverrides` for the editor: every copy submenu offers Copy as FASTA instead. */
export function rightClickOverrides(copyAsFasta: () => void) {
  const rewire = (items: MenuItem[]) => (Array.isArray(items) ? rewireCopyItems(items, copyAsFasta) : items);
  return Object.fromEntries(MENUS_WITH_COPY.map((key) => [key, rewire]));
}
