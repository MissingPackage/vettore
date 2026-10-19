// Enzyme sets (unique cutters, 6+ cutters, ...), computed with OVE's own cut-site engine (@teselagen/sequence-utils)
// so that what the map shows and what the set contains always agree. Enzyme data is the
// REBASE-derived, nonredundant commercial database in data/enzymes.json (see vettore.enzymes),
// cross-checked against Biopython in __tests__/enzymeDatabase.test.ts.
import { getCutsitesFromSequence } from "@teselagen/sequence-utils";
import database from "./data/enzymes.json";

export interface Enzyme {
  name: string;
  site: string;
  forwardRegex: string;
  reverseRegex: string;
  topSnipOffset: number;
  bottomSnipOffset: number;
  isType2S: boolean;
  suppliers: string[];
  isoschizomers: string[];
  goldenGate: boolean;
}

/** Enzyme database keyed by lower-case name (the shape OVE expects for `additionalEnzymes`). */
export const ENZYMES = database.enzymes as Record<string, Enzyme>;
export const ENZYME_SOURCE = database.source;

export type EnzymeSetId =
  | "unique6"
  | "unique"
  | "uniqueDual"
  | "sixPlus"
  | "type2s"
  | "goldenGate"
  | "commercial";

export const ENZYME_SETS: { id: EnzymeSetId; label: string }[] = [
  { id: "unique6", label: "Unique 6+ Cutters" },
  { id: "unique", label: "Unique Cutters" },
  { id: "uniqueDual", label: "Unique & Dual Cutters" },
  { id: "sixPlus", label: "6+ Cutters" },
  { id: "type2s", label: "Type IIS Enzymes" },
  { id: "goldenGate", label: "Golden Gate Enzymes" },
  { id: "commercial", label: "Nonredundant Commercial" },
];

export const DEFAULT_ENZYME_SET: EnzymeSetId = "unique6";

interface RawCutsite {
  topSnipPosition: number;
  bottomSnipPosition: number;
}

/**
 * Cut sites per enzyme, without OVE's phantom cuts: on a linear molecule a cut at either end
 * (e.g. a Type IIS site whose cleavage would fall outside the sequence) produces no fragment.
 */
export function cutsites<T extends RawCutsite>(
  sequence: string,
  circular: boolean,
  enzymes: Enzyme[] = Object.values(ENZYMES),
): Record<string, T[]> {
  const found: Record<string, T[]> = getCutsitesFromSequence(sequence, circular, enzymes);
  if (circular) return found;
  const inside = (p: number) => p > 0 && p < sequence.length;
  const out: Record<string, T[]> = {};
  for (const [name, sites] of Object.entries(found)) {
    const kept = sites.filter((s) => inside(s.topSnipPosition) && inside(s.bottomSnipPosition));
    if (kept.length) out[name] = kept;
  }
  return out;
}

/** Number of specified (non-N) bases in the recognition site. */
export function specifiedBases(site: string): number {
  return site.replace(/n/gi, "").length;
}

export function enzymesInSet(set: EnzymeSetId, sequence: string, circular: boolean): string[] {
  const cutsitesByName = cutsites(sequence, circular);
  // The database is keyed by lower-case name; cut sites are keyed by display name.
  const byName = (n: string) => ENZYMES[n.toLowerCase()];
  const cutters = Object.keys(cutsitesByName);
  const count = (n: string) => cutsitesByName[n].length;
  const six = (n: string) => specifiedBases(byName(n)?.site ?? "") >= 6;
  switch (set) {
    case "unique6":
      return cutters.filter((n) => count(n) === 1 && six(n));
    case "unique":
      return cutters.filter((n) => count(n) === 1);
    case "uniqueDual":
      return cutters.filter((n) => count(n) <= 2);
    case "sixPlus":
      return cutters.filter(six);
    case "type2s":
      return cutters.filter((n) => byName(n)?.isType2S);
    case "goldenGate":
      return cutters.filter((n) => byName(n)?.goldenGate);
    case "commercial":
      return cutters;
  }
}


export interface CutInfo {
  /** Cut location as bench scientists quote it: the base 5′ of the top-strand cut, 1-based. */
  location: number;
  /** 0-based inclusive recognition site range, for selecting it in the editor. */
  start: number;
  end: number;
}

export interface EnzymeRow {
  name: string;
  site: string;
  ends: "blunt" | "5′ overhang" | "3′ overhang";
  overhang: string;
  cuts: CutInfo[];
}

interface Cutsite extends RawCutsite {
  start: number;
  end: number;
  overhangSize: number;
  overhangBps?: string;
  topSnipBeforeBottom: boolean;
  recognitionSiteRange?: { start: number; end: number };
}

/** Rows for the Enzymes view: every enzyme of the set with its cut locations. */
export function enzymeRows(set: EnzymeSetId, sequence: string, circular: boolean): EnzymeRow[] {
  const names = new Set(enzymesInSet(set, sequence, circular).map((n) => n.toLowerCase()));
  const db = ENZYMES;
  const enzymes = [...names].map((n) => db[n]).filter(Boolean);
  const byName = cutsites<Cutsite>(sequence, circular, enzymes);
  return Object.entries(byName)
    .map(([name, sites]) => {
      const first = sites[0];
      const ends: EnzymeRow["ends"] = !first.overhangSize
        ? "blunt"
        : first.topSnipBeforeBottom
          ? "5′ overhang"
          : "3′ overhang";
      return {
        name,
        site: (db[name.toLowerCase()]?.site ?? "").toUpperCase(),
        ends,
        overhang: (first.overhangBps ?? "").toUpperCase(),
        cuts: sites
          .map((c) => ({
            location: c.topSnipPosition === 0 && circular ? sequence.length : c.topSnipPosition,
            start: c.recognitionSiteRange?.start ?? c.start,
            end: c.recognitionSiteRange?.end ?? c.end,
          }))
          .sort((a, b) => a.location - b.location),
      };
    })
    .sort((a, b) => a.name.localeCompare(b.name));
}
