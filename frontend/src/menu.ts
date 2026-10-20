// Menu bar built on top of OVE's command system.
// Items with `cmd` reuse OVE commands (and their hotkeys); items with `onClick` call Vettore.
import type { MutableRefObject } from "react";
import { baseName } from "./documents";
import { READY_ACTIONS } from "./cloning/registry";
import { ENZYME_SETS, type EnzymeSetId } from "./enzymeSets";
import { TOOLS } from "./tools";
import type { MenuName } from "./tools/types";
import type { ViewId } from "./views";

export interface MenuHandlers {
  newFile: () => void;
  /** Open what the clipboard holds (sequence, FASTA or GenBank text) as a new document. */
  newFromClipboard: () => void;
  open: () => void;
  openRecent: (path: string) => void;
  close: () => void;
  save: () => void;
  saveAs: () => void;
  exportGenbank: () => void;
  exportFasta: () => void;
  exportDna: () => void;
  /** Copy the selection, or the whole sequence, as FASTA text. */
  copyAsFasta: () => void;
  exportSvg: () => void;
  exportPng: () => void;
  view: (view: ViewId) => void;
  cloning: (actionId: string) => void;
  gel: () => void;
  sanger: () => void;
  detectFeatures: () => void;
  useEnzymeSet: (set: EnzymeSetId) => void;
  /** Open a dialog from the tool registry (tools/index.ts). */
  tool: (id: string) => void;
  recent: string[];
  hasDocument: boolean;
}

type MenuItem = Record<string, unknown> | string;

const findMenu = (menu: MenuItem[], text: string) =>
  menu.find((m) => typeof m === "object" && m.text === text) as Record<string, unknown> | undefined;

export function buildMenu(ref: MutableRefObject<MenuHandlers>) {
  const h = () => ref.current;
  return (defaultMenu: MenuItem[]): MenuItem[] => {
    const edit = findMenu(defaultMenu, "Edit");
    const view = findMenu(defaultMenu, "View");
    const help = findMenu(defaultMenu, "Help");
    const recent = h().recent;
    const hasDocument = h().hasDocument;

    // Registry tools of one menu, after a separator.
    const toolItems = (menu: MenuName): MenuItem[] => {
      const items = TOOLS.filter((t) => t.menu === menu).map((t) => ({
        text: `${t.label}…`,
        hotkey: t.hotkey,
        disabled: !!t.requiresDocument && !hasDocument,
        "data-test": `tool-${t.id}`,
        onClick: () => h().tool(t.id),
      }));
      return items.length ? ["--", ...items] : [];
    };
    const withTools = (menu: Record<string, unknown> | undefined, name: MenuName) =>
      menu && { ...menu, submenu: [...((menu.submenu as MenuItem[]) ?? []), ...toolItems(name)] };

    // Copy as FASTA sits right after the editor's own Copy.
    const editItems = [...((edit?.submenu as MenuItem[]) ?? [])];
    const copyAt = editItems.findIndex((m) => m === "copy" || (typeof m === "object" && m.cmd === "copy"));
    editItems.splice(copyAt >= 0 ? copyAt + 1 : editItems.length, 0, {
      text: "Copy as FASTA",
      disabled: !hasDocument,
      "data-test": "copy-as-fasta",
      onClick: () => h().copyAsFasta(),
    });
    const editMenu = edit && { ...edit, submenu: editItems };

    const file = {
      text: "File",
      "data-test": "file",
      submenu: [
        { text: "New DNA File…", hotkey: "mod+n", onClick: () => h().newFile() },
        { text: "New File from Clipboard", hotkey: "mod+shift+v", "data-test": "new-from-clipboard", onClick: () => h().newFromClipboard() },
        { text: "Open…", hotkey: "mod+o", onClick: () => h().open() },
        {
          text: "Open Recent",
          disabled: recent.length === 0,
          submenu: recent.map((path) => ({ text: baseName(path), label: path, onClick: () => h().openRecent(path) })),
        },
        "--",
        { text: "Close", hotkey: "mod+w", onClick: () => h().close() },
        { cmd: "saveSequence", text: "Save" },
        { cmd: "saveSequenceAs", text: "Save As…" },
        "--",
        {
          text: "Export",
          submenu: [
            { text: "GenBank (.gb)…", onClick: () => h().exportGenbank() },
            { text: "FASTA (.fasta)…", onClick: () => h().exportFasta() },
            { text: "DNA File (.dna)…", onClick: () => h().exportDna() },
            "--",
            { text: "Map as SVG…", onClick: () => h().exportSvg() },
            { text: "Map as PNG…", onClick: () => h().exportPng() },
          ],
        },
        "--",
        { text: "Print…", cmd: "print" },
        ...toolItems("File"),
      ],
    };

    const cloningItem = (id: string, label?: string) => {
      const action = READY_ACTIONS.find((a) => a.id === id);
      return action ? [{ text: `${label ?? action.title}…`, onClick: () => h().cloning(id) }] : [];
    };

    const actions = {
      text: "Actions",
      submenu: [
        {
          text: "Restriction and Insertion Cloning",
          submenu: [
            ...cloningItem("restriction", "Insert Fragment"),
            ...cloningItem("insert-multiple", "Insert Multiple Fragments"),
            ...cloningItem("delete-fragment", "Delete Restriction Fragment"),
          ],
        },
        ...cloningItem("pcr", "PCR"),
        ...cloningItem("design-primers", "Design Primers to Amplify a Region"),
        ...cloningItem("mutagenesis", "Mutagenesis"),
        ...cloningItem("anneal-oligos", "Anneal Oligos and Insert"),
        ...cloningItem("ligate", "Ligate Fragments"),
        "--",
        ...cloningItem("gibson", "Gibson Assembly"),
        ...cloningItem("hifi", "NEBuilder HiFi DNA Assembly"),
        ...cloningItem("infusion", "In-Fusion Cloning"),
        ...cloningItem("golden-gate", "Golden Gate Assembly"),
        ...cloningItem("gateway-bp", "Gateway BP Reaction"),
        ...cloningItem("gateway-lr", "Gateway LR Reaction"),
        ...cloningItem("ta", "TA or GC Cloning"),
        ...cloningItem("topo", "TOPO Cloning"),
        "--",
        "reverseComplementEntireSequence",
        "reverseComplementSelection",
        { text: "Set Origin at Cursor", cmd: "rotateToCaretPosition" },
        { text: "Change Topology", cmd: "changeCircularityCmd", submenu: ["toggleCircular", "toggleLinear"] },
        { text: "New File from Selection", cmd: "createNewFromSubsequence" },
        "--",
        { text: "Simulate Digestion…", cmd: "simulateDigestion" },
        ...toolItems("Actions"),
      ],
    };

    const features = {
      text: "Features",
      submenu: [
        { text: "Add Feature…", cmd: "newFeature" },
        { text: "Detect Common Features…", onClick: () => h().detectFeatures() },
        "--",
        { text: "Show Features", cmd: "toggleFeatures" },
        { text: "Feature Types…", cmd: "featureTypesCmd" },
        "--",
        { text: "Features View", hotkey: "mod+4", onClick: () => h().view("features") },
        ...toolItems("Features"),
      ],
    };

    const primers = {
      text: "Primers",
      submenu: [
        { text: "Add Primer…", cmd: "newPrimer" },
        "--",
        { text: "Show Primers", cmd: "togglePrimers" },
        "--",
        { text: "Primers View", hotkey: "mod+5", onClick: () => h().view("primers") },
        ...toolItems("Primers"),
      ],
    };

    const enzymes = {
      text: "Enzymes",
      submenu: [
        { text: "Choose Enzymes…", cmd: "openFilterCutsites" },
        {
          text: "Use Enzyme Set",
          submenu: ENZYME_SETS.map((set) => ({ text: set.label, onClick: () => h().useEnzymeSet(set.id) })),
        },
        { text: "Enzyme Manager…", cmd: "restrictionEnzymesManager" },
        { text: "New Custom Enzyme…", cmd: "openCreateCustomEnzyme" },
        "--",
        { text: "Show Enzymes", cmd: "toggleCutsites" },
        "--",
        { text: "Enzymes View", hotkey: "mod+3", onClick: () => h().view("enzymes") },
        ...toolItems("Enzymes"),
      ],
    };

    const tools = {
      text: "Tools",
      submenu: [
        { text: "Align to Reference DNA Sequence…", onClick: () => h().sanger() },
        { text: "Simulate Agarose Gel…", hotkey: "mod+shift+g", onClick: () => h().gel() },
        "--",
        "find",
        "goTo",
        "--",
        { text: "Minimum ORF Size…", cmd: "minOrfSizeCmd" },
        "toggleOrfs",
        ...toolItems("Tools"),
      ],
    };

    const views = {
      text: "View",
      submenu: [
        { text: "Map", hotkey: "mod+1", onClick: () => h().view("map") },
        { text: "Sequence", hotkey: "mod+2", onClick: () => h().view("sequence") },
        { text: "Enzymes", hotkey: "mod+3", onClick: () => h().view("enzymes") },
        { text: "Features", hotkey: "mod+4", onClick: () => h().view("features") },
        { text: "Primers", hotkey: "mod+5", onClick: () => h().view("primers") },
        { text: "History", hotkey: "mod+6", onClick: () => h().view("history") },
        "--",
        ...((view?.submenu as MenuItem[]) ?? []),
        ...toolItems("View"),
      ],
    };

    // Menu order: File, Edit, View, Enzymes, Features, Primers, Actions, Tools, Help.
    const helpMenu = withTools(help, "Help") ?? (toolItems("Help").length
      ? { text: "Help", submenu: toolItems("Help").slice(1) } : undefined);
    return [file, withTools(editMenu, "Edit"), views, enzymes, features, primers, actions, tools, helpMenu]
      .filter(Boolean) as MenuItem[];
  };
}
