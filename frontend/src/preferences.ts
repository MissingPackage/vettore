// User preferences kept in the window's local storage (persisted by the desktop shell).
// Changes apply immediately: components read them through `usePreferences`, which re-renders
// every subscriber when something is saved.
import { useEffect, useState } from "react";
import { DEFAULT_ENZYME_SET, ENZYMES, ENZYME_SETS, type EnzymeSetId } from "./enzymeLogic";

const KEY = "vettore-preferences";

export interface Preferences {
  /** Enzyme set shown on the map of newly opened documents. */
  enzymeSet: EnzymeSetId;
  /** Ladder pre-selected in Simulate Agarose Gel (engine id, e.g. "neb-1kb-plus"). */
  gelLadder: string;
  /** Agarose percentage pre-selected in Simulate Agarose Gel. */
  gelPercent: number;
  /** Minimum identity (%) for Detect Common Features. */
  detectMinIdentity: number;
  /** The frame around the work (menus, tabs, status bar): light, dark, or as the system is set. */
  chrome: ChromeTheme;
  /** The enzymes in the lab's freezer: Find Diagnostic Digest tries these first. Empty until
   * the user fills it in, and then the commercial set is used. */
  myEnzymes: string[];
}

export type ChromeTheme = "light" | "dark" | "system";
export const CHROME_THEMES: { id: ChromeTheme; label: string }[] = [
  { id: "light", label: "Light" },
  { id: "dark", label: "Dark menus and tabs" },
  { id: "system", label: "As the system is set" },
];

/** Identity thresholds offered for feature detection; the engine's own default is 95. */
export const IDENTITY_CHOICES = [90, 95, 98, 100];

const DEFAULTS: Preferences = {
  enzymeSet: DEFAULT_ENZYME_SET,
  gelLadder: "neb-1kb-plus",
  gelPercent: 1.0,
  detectMinIdentity: 95,
  chrome: "light",
  myEnzymes: [],
};

function clean(stored: Partial<Preferences>): Preferences {
  const prefs = { ...DEFAULTS, ...stored };
  if (!ENZYME_SETS.some((s) => s.id === prefs.enzymeSet)) prefs.enzymeSet = DEFAULTS.enzymeSet;
  if (typeof prefs.gelLadder !== "string" || !prefs.gelLadder) prefs.gelLadder = DEFAULTS.gelLadder;
  if (!(prefs.gelPercent > 0 && prefs.gelPercent <= 3)) prefs.gelPercent = DEFAULTS.gelPercent;
  if (!IDENTITY_CHOICES.includes(prefs.detectMinIdentity)) {
    prefs.detectMinIdentity = DEFAULTS.detectMinIdentity;
  }
  if (!CHROME_THEMES.some((t) => t.id === prefs.chrome)) prefs.chrome = DEFAULTS.chrome;
  prefs.myEnzymes = knownEnzymes(prefs.myEnzymes);
  return prefs;
}

export function defaultPreferences(): Preferences {
  return { ...DEFAULTS, myEnzymes: [] };
}

/** Enzyme names as the database writes them, each once; anything else is dropped. */
export function knownEnzymes(names: unknown): string[] {
  if (!Array.isArray(names)) return [];
  const out: string[] = [];
  for (const raw of names) {
    const name = typeof raw === "string" ? ENZYMES[raw.trim().toLowerCase()]?.name : undefined;
    if (name && !out.includes(name)) out.push(name);
  }
  return out;
}

/**
 * Enzyme names typed or pasted as a list (commas, spaces, semicolons or new lines between
 * them, as copied from a spreadsheet): the names Vettore knows, and the words it does not.
 */
export function parseEnzymeList(text: string): { known: string[]; unknown: string[] } {
  const words = text.split(/[\s,;]+/).map((w) => w.trim()).filter(Boolean);
  const known = knownEnzymes(words);
  const unknown = words.filter((w) => !ENZYMES[w.toLowerCase()]);
  return { known, unknown: [...new Set(unknown)] };
}

export function loadPreferences(): Preferences {
  try {
    return clean(JSON.parse(localStorage.getItem(KEY) ?? "{}") as Partial<Preferences>);
  } catch {
    return defaultPreferences();
  }
}

const listeners = new Set<(prefs: Preferences) => void>();

export function savePreferences(patch: Partial<Preferences>): Preferences {
  const prefs = clean({ ...loadPreferences(), ...patch });
  try {
    localStorage.setItem(KEY, JSON.stringify(prefs));
  } catch {
    // A full or disabled storage must not break the dialog: the change still applies to
    // this session through the listeners below.
  }
  listeners.forEach((listener) => listener(prefs));
  return prefs;
}

/** Current preferences, re-read whenever any of them is saved. */
export function usePreferences(): Preferences {
  const [prefs, setPrefs] = useState(loadPreferences);
  useEffect(() => {
    listeners.add(setPrefs);
    setPrefs(loadPreferences());
    return () => {
      listeners.delete(setPrefs);
    };
  }, []);
  return prefs;
}

/** Puts data-chrome="light" or "dark" on <html> (theme/dark-chrome.css); "system" follows the
 * operating system, also when it changes while Vettore is open. Returns the undo. */
export function applyChrome(theme: ChromeTheme): () => void {
  const root = document.documentElement;
  if (theme !== "system") {
    root.dataset.chrome = theme;
    return () => {};
  }
  const query = window.matchMedia?.("(prefers-color-scheme: dark)");
  const set = () => { root.dataset.chrome = query?.matches ? "dark" : "light"; };
  set();
  query?.addEventListener?.("change", set);
  return () => query?.removeEventListener?.("change", set);
}
