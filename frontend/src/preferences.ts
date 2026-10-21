// User preferences kept in the window's local storage (persisted by the desktop shell).
// Changes apply immediately: components read them through `usePreferences`, which re-renders
// every subscriber when something is saved.
import { useEffect, useState } from "react";
import { DEFAULT_ENZYME_SET, ENZYME_SETS, type EnzymeSetId } from "./enzymeLogic";

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
  return prefs;
}

export function defaultPreferences(): Preferences {
  return { ...DEFAULTS };
}

export function loadPreferences(): Preferences {
  try {
    return clean(JSON.parse(localStorage.getItem(KEY) ?? "{}") as Partial<Preferences>);
  } catch {
    return { ...DEFAULTS };
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
