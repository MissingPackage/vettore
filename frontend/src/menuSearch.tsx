// Ranking for the menu search (Help ▸ search box, Ctrl+/ or the "Search commands" button).
// The editor's own search matched letters in order anywhere in a name, so "gel" listed
// "Change Topology" before "Simulate Agarose Gel…", and it ranked only by name length. Here a
// name that starts with what was typed comes first, then a word that starts with it, then the
// text anywhere, then a match on the menu it sits in ("tools" lists the Tools menu).
import type { ReactNode } from "react";

export interface SearchableItem {
  text?: unknown;
  onClick?: unknown;
  hidden?: boolean;
  hideFromMenuSearch?: boolean;
  hiddenButSearchableText?: string;
  /** The menus above the item; the editor can put React elements here, which are skipped. */
  path?: unknown[];
  [key: string]: unknown;
}

const MAX_RESULTS = 12;

/** 0 best; null when the item does not match. */
export function matchScore(query: string, text: string, path: unknown[] = [], extra = ""): number | null {
  const q = query.trim().toLowerCase();
  if (!q) return null;
  const t = text.toLowerCase();
  if (t.startsWith(q)) return 0;
  if (t.split(/[\s(/–-]+/).some((word) => word.startsWith(q))) return 1;
  if (t.includes(q)) return 2;
  if (path.some((p) => typeof p === "string" && p.toLowerCase().includes(q))) return 3;
  if (typeof extra === "string" && extra.toLowerCase().includes(q)) return 3;
  return null;
}

/** The name with the typed text in bold, where it appears. */
function highlighted(query: string, text: string): ReactNode {
  const at = text.toLowerCase().indexOf(query.trim().toLowerCase());
  if (at < 0) return text;
  const end = at + query.trim().length;
  return <>{text.slice(0, at)}<strong>{text.slice(at, end)}</strong>{text.slice(end)}</>;
}

/** Drop-in for the editor's itemListPredicate: (query, every menu item) → items to list. */
export function rankMenuItems<T extends SearchableItem>(query: string, items: T[]): T[] {
  const scored: { item: T; text: string; score: number }[] = [];
  for (const item of items) {
    if (typeof item.text !== "string" || !item.onClick || item.hidden || item.hideFromMenuSearch) continue;
    const score = matchScore(query, item.text, item.path, item.hiddenButSearchableText);
    if (score !== null) scored.push({ item, text: item.text, score });
  }
  scored.sort((a, b) => a.score - b.score || a.text.length - b.text.length);
  return scored.slice(0, MAX_RESULTS).map(({ item, text }) => ({
    ...item,
    isSimpleText: true,
    justText: highlighted(query, text),
  }));
}
