// Menu bar built on top of OVE's command system.
// Items with `cmd` reuse OVE commands (and their hotkeys); items with `onClick` call Vettore.
import type { MutableRefObject } from "react";
import { baseName } from "./documents";
import { READY_ACTIONS } from "./cloning/registry";
import { ENZYME_SETS, type EnzymeSetId } from "./enzymeSets";
import type { ViewId } from "./views";

export interface MenuHandlers {
  newFile: () => void;
  open: () => void;
  openRecent: (path: string) => void;
  close: () => void;
  save: () => void;
  saveAs: () => void;
  exportGenbank: () => void;
  exportFasta: () => void;
  exportSvg: () => void;
  exportPng: () => void;
  view: (view: ViewId) => void;
  cloning: (actionId: string) => void;
  gel: () => void;
  useEnzymeSet: (set: EnzymeSetId) => void;
  recent: string[];
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

    const file = {
      text: "File",
      "data-test": "file",
      submenu: [
        { text: "New DNA File…", hotkey: "mod+n", onClick: () => h().newFile() },
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
            { text: "GenBank…", onClick: () => h().exportGenbank() },
            { text: "FASTA…", onClick: () => h().exportFasta() },
            "--",
            { text: "Map as SVG…", onClick: () => h().exportSvg() },
            { text: "Map as PNG…", onClick: () => h().exportPng() },
          ],
        },
        "--",
        { text: "Print…", cmd: "print" },
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
            ...cloningItem("delete-fragment", "Delete Restriction Fragment"),
          ],
        },
        ...cloningItem("pcr", "PCR"),
        ...cloningItem("design-primers", "Design Primers to Amplify a Region"),
        ...cloningItem("mutagenesis", "Mutagenesis"),
        "--",
        ...cloningItem("gibson", "Gibson Assembly"),
        ...cloningItem("hifi", "NEBuilder HiFi DNA Assembly"),
        ...cloningItem("infusion", "In-Fusion Cloning"),
        ...cloningItem("golden-gate", "Golden Gate Assembly"),
        "--",
        "reverseComplementEntireSequence",
        "reverseComplementSelection",
        { text: "Set Origin at Cursor", cmd: "rotateToCaretPosition" },
        { text: "Change Topology", cmd: "changeCircularityCmd", submenu: ["toggleCircular", "toggleLinear"] },
        { text: "New File from Selection", cmd: "createNewFromSubsequence" },
        "--",
        { text: "Simulate Digestion…", cmd: "simulateDigestion" },
      ],
    };

    const features = {
      text: "Features",
      submenu: [
        { text: "Add Feature…", cmd: "newFeature" },
        "--",
        { text: "Show Features", cmd: "toggleFeatures" },
        { text: "Feature Types…", cmd: "featureTypesCmd" },
        "--",
        { text: "Features View", hotkey: "mod+4", onClick: () => h().view("features") },
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
      ],
    };

    const tools = {
      text: "Tools",
      submenu: [
        { text: "Simulate Agarose Gel…", hotkey: "mod+shift+g", onClick: () => h().gel() },
        "--",
        "find",
        "goTo",
        "--",
        { text: "Minimum ORF Size…", cmd: "minOrfSizeCmd" },
        "toggleOrfs",
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
      ],
    };

    // Menu order: File, Edit, View, Enzymes, Features, Primers, Actions, Tools, Help.
    return [file, edit, views, enzymes, features, primers, actions, tools, help].filter(Boolean) as MenuItem[];
  };
}
