// Document view tabs mapped onto Open Vector Editor panels.
import { actions } from "@teselagen/ove";
import { store } from "./store";

export type ViewId = "map" | "sequence" | "enzymes" | "features" | "primers" | "history";

export const VIEWS: { id: ViewId; label: string; hotkey: string }[] = [
  { id: "map", label: "Map", hotkey: "1" },
  { id: "sequence", label: "Sequence", hotkey: "2" },
  { id: "enzymes", label: "Enzymes", hotkey: "3" },
  { id: "features", label: "Features", hotkey: "4" },
  { id: "primers", label: "Primers", hotkey: "5" },
  { id: "history", label: "History", hotkey: "6" },
];

const PROPERTIES_TAB: Partial<Record<ViewId, string>> = {
  features: "features",
  primers: "primers",
};

/** Views rendered by Vettore itself on top of the editor panels. */
export const OVERLAY_VIEWS: ViewId[] = ["enzymes", "history"];

/** A single full-width panel group: one view at a time fills the document area. */
export function panelsFor(circular: boolean, view: ViewId = "map") {
  const panelId = panelIdFor(view, circular);
  return [
    [
      { id: circular ? "circular" : "rail", name: "Map", active: panelId !== "sequence" && panelId !== "properties" },
      { id: "sequence", name: "Sequence", active: panelId === "sequence" },
      { id: "properties", name: "Properties", active: panelId === "properties" },
    ],
  ];
}

function panelIdFor(view: ViewId, circular: boolean): string {
  if (view === "sequence") return "sequence";
  if (view in PROPERTIES_TAB) return "properties";
  return circular ? "circular" : "rail";
}

export function showView(editorName: string, view: ViewId, circular: boolean) {
  if (OVERLAY_VIEWS.includes(view)) return;
  const meta = { editorName };
  store.dispatch(actions.panelsShownUpdate(panelsFor(circular, view), meta));
  const tab = PROPERTIES_TAB[view];
  if (tab) store.dispatch(actions.propertiesViewTabUpdate(tab, meta));
}
