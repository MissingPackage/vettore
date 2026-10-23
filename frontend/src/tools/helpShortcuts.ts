// The keyboard shortcuts Help → Keyboard Shortcuts lists.
//
// Only shortcuts that are really bound in the running app are listed here:
// - Ctrl+O / Ctrl+N / Ctrl+W, Ctrl+Shift+V, Ctrl+Shift+G and Ctrl+1 … Ctrl+6 are bound by
//   Vettore itself (App.tsx);
// - the rest are Open Vector Editor commands, read from the editor's own hotkey registry
//   (its command definitions in @teselagen/ove) and checked in the running app (e2e/help.spec.ts).
// A shortcut the menus advertise but nothing binds must not be added here.

export interface Shortcut {
  keys: string;
  what: string;
}

export interface ShortcutGroup {
  title: string;
  shortcuts: Shortcut[];
}

export const SHORTCUT_GROUPS: ShortcutGroup[] = [
  {
    title: "Files and tabs",
    shortcuts: [
      { keys: "Ctrl+O", what: "Open a sequence file" },
      { keys: "Ctrl+N", what: "New DNA file" },
      { keys: "Ctrl+Shift+V", what: "New file from a sequence you copied" },
      { keys: "Ctrl+S", what: "Save (GenBank)" },
      { keys: "Ctrl+Shift+S", what: "Save As…" },
      { keys: "Ctrl+P", what: "Print the current view" },
      { keys: "Ctrl+W", what: "Close the open document" },
      { keys: "Ctrl+Shift+G", what: "Simulate an agarose gel" },
    ],
  },
  {
    title: "Views",
    shortcuts: [
      { keys: "Ctrl+1", what: "Map" },
      { keys: "Ctrl+2", what: "Sequence" },
      { keys: "Ctrl+3", what: "Enzymes" },
      { keys: "Ctrl+4", what: "Features" },
      { keys: "Ctrl+5", what: "Primers" },
      { keys: "Ctrl+6", what: "History" },
      { keys: "Ctrl+7", what: "Lineage" },
    ],
  },
  {
    title: "Finding and selecting",
    shortcuts: [
      { keys: "Ctrl+F", what: "Find DNA or protein, with mismatches" },
      { keys: "Ctrl+G", what: "Go to a position" },
      { keys: "Ctrl+A", what: "Select the whole sequence" },
      { keys: "Ctrl+I", what: "Select everything except the current selection" },
      { keys: "Ctrl+/", what: "Search the menu commands by name" },
    ],
  },
  {
    title: "Editing",
    shortcuts: [
      { keys: "Ctrl+Z", what: "Undo" },
      { keys: "Ctrl+Shift+Z", what: "Redo" },
      { keys: "Ctrl+E", what: "Reverse complement the selection" },
      { keys: "Ctrl+B", what: "Set the origin at the cursor" },
      { keys: "Ctrl+K", what: "Add a feature over the selection" },
    ],
  },
  {
    title: "The editor's own panels",
    shortcuts: [
      {
        keys: "Ctrl+Shift+D",
        what: "Digest panel: pick enzymes and see the fragments. A view tab brings the map back",
      },
      {
        keys: "Ctrl+Shift+P",
        what: "PCR panel of the editor. Vettore's own simulation is Actions → PCR…",
      },
    ],
  },
];

export const ALL_SHORTCUTS: Shortcut[] = SHORTCUT_GROUPS.flatMap((g) => g.shortcuts);
