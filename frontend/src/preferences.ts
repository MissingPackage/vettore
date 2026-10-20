// User preferences kept in the window's local storage (persisted by the desktop shell).
import { DEFAULT_ENZYME_SET, ENZYME_SETS, type EnzymeSetId } from "./enzymeLogic";

const KEY = "vettore-preferences";

export interface Preferences {
  /** Enzyme set shown on the map of newly opened documents. */
  enzymeSet: EnzymeSetId;
}

const DEFAULTS: Preferences = { enzymeSet: DEFAULT_ENZYME_SET };

export function loadPreferences(): Preferences {
  try {
    const stored = JSON.parse(localStorage.getItem(KEY) ?? "{}") as Partial<Preferences>;
    const prefs = { ...DEFAULTS, ...stored };
    if (!ENZYME_SETS.some((s) => s.id === prefs.enzymeSet)) prefs.enzymeSet = DEFAULTS.enzymeSet;
    return prefs;
  } catch {
    return { ...DEFAULTS };
  }
}

export function savePreferences(patch: Partial<Preferences>): Preferences {
  const prefs = { ...loadPreferences(), ...patch };
  localStorage.setItem(KEY, JSON.stringify(prefs));
  return prefs;
}
