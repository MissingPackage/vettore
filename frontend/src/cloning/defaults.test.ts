import { describe, expect, it } from "vitest";
import {
  defaultInserts, defaultLinear, defaultVector, docData, donorPairs, gatewayChain, gatewayDefaults,
  openingFromSelection, selectionOf,
} from "./defaults";
import { newRow, rowBody } from "./fragmentRows";
import type { OpenDoc } from "./types";

function doc(editorName: string, name: string, circular: boolean,
             selection: { start: number; end: number } | null = null): OpenDoc {
  return {
    editorName,
    name,
    data: () => ({ name: `${name} (stored)`, sequence: "ACGT".repeat(50), circular, features: [], primers: [] }),
    selection: () => selection,
  };
}

const plasmid = doc("d1", "pUC", true);
const pcr1 = doc("d2", "GFP PCR", false, { start: 9, end: 99 });
const plasmid2 = doc("d3", "pDonor", true);
const pcr2 = doc("d4", "WPRE PCR", false);
const docs = [plasmid, pcr1, plasmid2, pcr2];

describe("defaults", () => {
  it("names the data as the tab, so products list their parents by tab name", () => {
    expect(docData(pcr1).name).toBe("GFP PCR");
  });

  it("takes the active plasmid as the vector, else the first circular file", () => {
    expect(defaultVector(docs, plasmid2)).toBe("d3");
    expect(defaultVector(docs, pcr1)).toBe("d1");
    expect(defaultVector([pcr1], pcr1)).toBe("d2");
    expect(defaultVector([], null)).toBeNull();
  });

  it("fills inserts with the other files in tab order, linear ones first", () => {
    expect(defaultInserts(docs, "d1", 2)).toEqual(["d2", "d4"]);
    expect(defaultInserts(docs, "d1", 3)).toEqual(["d2", "d4", "d3"]);
    expect(defaultInserts(docs, null, 4)).toEqual(["d2", "d4", "d1", "d3"]);
    expect(defaultInserts([plasmid], "d1", 2)).toEqual([null, null]);
  });

  it("picks a linear file as the PCR product", () => {
    expect(defaultLinear(docs, pcr1, "d1")).toBe("d2");
    expect(defaultLinear(docs, plasmid, "d1")).toBe("d4");
    expect(defaultLinear([plasmid], plasmid, null)).toBeNull();
  });

  it("picks Gateway inputs by their att sites, whichever tab is active", () => {
    const entry = doc("e", "pENTR", true);
    const dest = doc("x", "pDEST", true);
    const other = doc("o", "pUC19", true);
    const attB = doc("b", "attB PCR", false);
    const expression = doc("c", "pEXP", true);
    const donor = doc("p", "pDONR", true);
    const open = [other, dest, entry, attB, expression, donor];
    const sites = [[], ["attR1", "attR2"], ["attL1", "attL2"], ["attB1", "attB2"], ["attB1", "attB2"], ["attP1", "attP2"]];
    // The destination tab is active: it still goes to the destination role.
    expect(gatewayDefaults("LR", open, sites, dest)).toEqual(["e", "x"]);
    expect(gatewayDefaults("LR", open, sites, entry)).toEqual(["e", "x"]);
    expect(gatewayDefaults("LR", open, sites, other)).toEqual(["e", "x"]);
    // BP: the linear attB product before the circular clone with attB sites, unless that is active.
    expect(gatewayDefaults("BP", open, sites, null)).toEqual(["b", "p"]);
    expect(gatewayDefaults("BP", open, sites, expression)).toEqual(["c", "p"]);
  });

  it("leaves a Gateway input empty when no open file has the sites it needs", () => {
    const entry = doc("e", "pENTR", true);
    const other = doc("o", "pUC19", true);
    expect(gatewayDefaults("LR", [entry, other], [["attL1", "attL2"], []], other)).toEqual(["e", null]);
    expect(gatewayDefaults("LR", [entry, other], [], entry)).toEqual([null, null]);
    expect(gatewayDefaults("BP", [entry, other], [["attL1", "attL2"], []], null)).toEqual([null, null]);
  });

  it("chains MultiSite entry clones through the att sites of the destination vector", () => {
    // A three-fragment reaction: attR4/attR3 destination, 5' entry (attL4-attR1), middle
    // entry (attL1-attL2), 3' entry (attR2-attL3), in any tab order.
    const dest = doc("x", "pDEST R4-R3", true);
    const five = doc("p5", "p5E", true);
    const mid = doc("pm", "pME", true);
    const three = doc("p3", "p3E", true);
    const other = doc("o", "pUC19", true);
    const open = [other, three, dest, five, mid];
    const sites = [[], ["attL3", "attR2"], ["attR3", "attR4"], ["attL4", "attR1"], ["attL1", "attL2"]];
    expect(gatewayChain(open, sites, null)).toEqual({ entries: ["p5", "pm", "p3"], vector: "x" });
    // Two entry clones into a destination for two fragments (attR4/attR2).
    const dest2 = doc("x2", "pDEST R4-R2", true);
    const two = [five, mid, dest2];
    expect(gatewayChain(two, [["attL4", "attR1"], ["attL1", "attL2"], ["attR2", "attR4"]], null))
      .toEqual({ entries: ["p5", "pm"], vector: "x2" });
  });

  it("falls back to one entry clone when the open files do not chain", () => {
    const dest = doc("x", "pDEST", true);
    const entry = doc("e", "pENTR", true);
    const stray = doc("s", "p5E", true);
    const open = [entry, stray, dest];
    const sites = [["attL1", "attL2"], ["attL4", "attR1"], ["attR1", "attR2"]];
    expect(gatewayChain(open, sites, null)).toEqual({ entries: ["e"], vector: "x" });
    expect(gatewayChain([entry], [["attL1", "attL2"]], null)).toEqual({ entries: ["e"], vector: null });
  });

  it("offers the att pairs a donor vector carries", () => {
    expect(donorPairs(["attP1", "attP2"])).toEqual([[1, 2]]);
    expect(donorPairs(["attP1", "attP4"])).toEqual([[1, 4]]);
    expect(donorPairs(["attP1", "attP2", "attP4"])).toEqual([[1, 2], [1, 4], [2, 4]]);
    expect(donorPairs(["attL1", "attL2"])).toEqual([]);
    expect(donorPairs(undefined)).toEqual([]);
  });

  it("reads regions and opening points from the editor selection", () => {
    expect(selectionOf(pcr1)).toEqual({ start: 9, end: 99 });
    expect(selectionOf(plasmid)).toEqual({ start: null, end: null });
    expect(openingFromSelection(pcr1)).toBe(100);
    expect(openingFromSelection(plasmid)).toBeNull();
  });
});

describe("fragment rows", () => {
  it("cuts circular files and uses linear ones as PCR products when the form prefers it", () => {
    expect(newRow(docs, "d1", "pcr").mode).toBe("digest");
    const row = newRow(docs, "d2", "pcr");
    expect(row).toMatchObject({ mode: "pcr", start: 9, end: 99 });
    expect(newRow(docs, "d2", "digest").mode).toBe("digest");
  });

  it("says what is missing before it builds a request", () => {
    expect(rowBody(docs, newRow(docs, null, "digest"), "insert 1")).toBe("Choose the file for the insert 1");
    expect(rowBody(docs, newRow(docs, "d1", "digest"), "vector")).toBe("Choose the enzymes that cut “pUC” (vector)");
    expect(rowBody(docs, { ...newRow(docs, "d3", "pcr"), mode: "pcr" }, "insert 2")).toMatch(/is circular: choose the region/);
    expect(rowBody(docs, { ...newRow(docs, "d4", "pcr"), start: 3 }, "x")).toMatch(/Enter both ends/);
    expect(rowBody(docs, { ...newRow(docs, "d1", "digest"), enzymes: ["EcoRI"] }, "vector")).toEqual({
      doc: expect.objectContaining({ name: "pUC" }),
      enzymes: ["EcoRI"],
      fragment_index: null,
    });
    expect(rowBody(docs, newRow(docs, "d4", "pcr"), "x")).toEqual({
      doc: expect.objectContaining({ name: "WPRE PCR" }), enzymes: [], start: null, end: null,
    });
  });
});
