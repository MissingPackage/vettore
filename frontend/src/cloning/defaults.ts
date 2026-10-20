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
