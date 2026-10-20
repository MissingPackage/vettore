// Starting values for the cloning forms, so that the first click already does something useful:
// the vector is the open plasmid, the inserts are the other open files, regions come from the
// editor selection. Pure functions of the open documents (unit-tested).
import type { SequenceData } from "../api";
import type { OpenDoc } from "./types";

const isCircular = (d: OpenDoc) => !!d.data()?.circular;

/** The document's content, named as its tab: products list their parents by this name. */
export function docData(d: OpenDoc): SequenceData {
  return { ...d.data(), name: d.name };
}

/** The vector: the active document if it is circular, else the first circular one. */
export function defaultVector(docs: OpenDoc[], active: OpenDoc | null): string | null {
  if (active && isCircular(active)) return active.editorName;
  return docs.find(isCircular)?.editorName ?? active?.editorName ?? null;
}

/**
 * `count` documents other than `exclude`, in tab order with linear ones (PCR products,
 * fragments) first; null when there are not enough documents.
 */
export function defaultInserts(docs: OpenDoc[], exclude: string | null, count: number): (string | null)[] {
  const others = docs.filter((d) => d.editorName !== exclude);
  const ordered = [...others.filter((d) => !isCircular(d)), ...others.filter(isCircular)];
  return Array.from({ length: count }, (_, i) => ordered[i]?.editorName ?? null);
}

/** A linear document to use as a PCR product: the active one if linear, else the last linear. */
export function defaultLinear(docs: OpenDoc[], active: OpenDoc | null, exclude: string | null): string | null {
  if (active && !isCircular(active) && active.editorName !== exclude) return active.editorName;
  const linear = docs.filter((d) => !isCircular(d) && d.editorName !== exclude);
  return linear[linear.length - 1]?.editorName ?? null;
}

/** The att sites each Gateway input carries: [insert side, vector side]. */
export const GATEWAY_SITE_KINDS = { BP: ["attB", "attP"], LR: ["attL", "attR"] } as const;

/**
 * Gateway inputs chosen by their att sites (`sites[i]` lists the sites of `docs[i]`, as
 * /api/cloning/gateway-sites answers): for LR the entry clone (attL) and the destination
 * vector (attR), for BP the attB insert and the donor vector (attP). The active document is
 * preferred when it fits, then the one with more sites of the kind, then (for the attB insert)
 * a linear PCR product. A role that no open file fits stays null, for the user to choose.
 */
export function gatewayDefaults(
  reaction: "BP" | "LR",
  docs: OpenDoc[],
  sites: string[][],
  active: OpenDoc | null,
): [string | null, string | null] {
  const [insertKind, vectorKind] = GATEWAY_SITE_KINDS[reaction];
  const pick = (kind: string, exclude: string | null, preferLinear: boolean) => {
    const fits = docs
      .map((d, i) => ({ d, i, n: (sites[i] ?? []).filter((s) => s.startsWith(kind)).length }))
      .filter((x) => x.n > 0 && x.d.editorName !== exclude);
    fits.sort((x, y) =>
      Number(y.d.editorName === active?.editorName) - Number(x.d.editorName === active?.editorName)
      || y.n - x.n
      || (preferLinear ? Number(isCircular(x.d)) - Number(isCircular(y.d)) : 0)
      || x.i - y.i);
    return fits[0]?.d.editorName ?? null;
  };
  const insert = pick(insertKind, null, reaction === "BP");
  return [insert, pick(vectorKind, insert, false)];
}

interface AttSite {
  kind: string;
  n: number;
}

const attSites = (names: string[]): AttSite[] =>
  names.filter((s) => /^att[BLPR][1-5]$/.test(s)).map((s) => ({ kind: s.slice(0, 4), n: Number(s.slice(4)) }));

const partner = (kind: string) => (kind === "attL" ? "attR" : "attL");

/** The most entry clones one MultiSite LR reaction takes. */
export const MAX_ENTRY_CLONES = 4;

/**
 * Inputs for an LR reaction, chosen by their att sites: the destination vector and the entry
 * clones that chain through it, in reaction order. A MultiSite reaction lines the molecules up
 * through the numbered sites (each attL pairs with the attR of the same number), so a
 * destination with attR4 and attR3 takes, for example, an attL4/attR1 clone, an attL1/attL2
 * clone and an attR2/attL3 clone. The longest chain the open files allow is offered; when no
 * chain is found, the single best entry clone and destination vector are (as before).
 */
export function gatewayChain(
  docs: OpenDoc[],
  sites: string[][],
  active: OpenDoc | null,
): { entries: (string | null)[]; vector: string | null } {
  const found = docs.map((_, i) => attSites(sites[i] ?? []).filter((s) => s.kind === "attL" || s.kind === "attR"));
  const vectors = docs
    .map((_, i) => i)
    .filter((i) => found[i].length >= 2 && found[i].every((s) => s.kind === "attR")
      && new Set(found[i].map((s) => s.n)).size === 2)
    .sort((x, y) => Number(docs[y].editorName === active?.editorName) - Number(docs[x].editorName === active?.editorName));

  for (const v of vectors) {
    // A MultiSite cassette reads 5 - 4 - 1 - 2 - 3, so the entry clones are listed starting
    // from whichever of the destination's two sites comes first in that order.
    const [first, second] = [...new Set(found[v].map((s) => s.n))]
      .sort((x, y) => [5, 4, 1, 2, 3].indexOf(x) - [5, 4, 1, 2, 3].indexOf(y));
    let best: number[] | null = null;
    for (const [from, to] of [[first, second], [second, first]]) {
      const walk = (need: AttSite, used: number[]): number[][] => {
        const chains: number[][] = [];
        docs.forEach((_, i) => {
          if (i === v || used.includes(i) || found[i].length !== 2) return;
          const match = found[i].find((s) => s.kind === need.kind && s.n === need.n);
          if (!match) return;
          const open = found[i].find((s) => s !== match)!;
          const next = { kind: partner(open.kind), n: open.n };
          const path = [...used, i];
          if (next.kind === "attR" && next.n === to) chains.push(path);
          if (path.length < MAX_ENTRY_CLONES) chains.push(...walk(next, path));
        });
        return chains;
      };
      for (const chain of walk({ kind: "attL", n: from }, [])) {
        if (!best || chain.length > best.length) best = chain;
      }
    }
    if (best) return { entries: best.map((i) => docs[i].editorName), vector: docs[v].editorName };
  }
  const [entry, vector] = gatewayDefaults("LR", docs, sites, active);
  return { entries: [entry], vector };
}

/** The att pairs a donor vector offers a BP reaction: [[1, 2]] for a donor with attP1/attP2. */
export function donorPairs(sites: string[] | undefined): number[][] {
  const numbers = [...new Set(attSites(sites ?? []).filter((s) => s.kind === "attP").map((s) => s.n))].sort();
  const pairs: number[][] = [];
  numbers.forEach((a, i) => numbers.slice(i + 1).forEach((b) => pairs.push([a, b])));
  return pairs;
}

/** The editor selection of a document as a region, or nulls. */
export function selectionOf(d: OpenDoc | null): { start: number | null; end: number | null } {
  const s = d?.selection() ?? null;
  return s ? { start: s.start, end: s.end } : { start: null, end: null };
}

/** Where a selection opens a circular vector: right after its last base (0..length). */
export function openingFromSelection(d: OpenDoc | null): number | null {
  const s = d?.selection() ?? null;
  if (!s) return null;
  return s.end + 1;
}
