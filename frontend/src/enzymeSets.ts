// Applies enzyme sets to an OVE editor (pure logic lives in enzymeLogic.ts).
import { actions } from "@teselagen/ove";
import { enzymesInSet, type EnzymeSetId } from "./enzymeLogic";
import { store } from "./store";

export { DEFAULT_ENZYME_SET, ENZYME_SETS, enzymesInSet, type EnzymeSetId } from "./enzymeLogic";

const NO_MATCH = { value: "__vettore_no_enzymes__", label: "(none)" };

export function applyEnzymeSet(editorName: string, set: EnzymeSetId) {
  const ed = (store.getState() as any).VectorEditor[editorName];
  const sequence: string = ed?.sequenceData?.sequence ?? "";
  const circular = !!ed?.sequenceData?.circular;
  const names = enzymesInSet(set, sequence, circular);
  // An empty filter means "show everything" in OVE, so an empty set needs a placeholder.
  const filter = names.length ? names.map((n) => ({ value: n.toLowerCase(), label: n })) : [NO_MATCH];
  store.dispatch(actions.filteredRestrictionEnzymesUpdate(filter, { editorName }));
}
