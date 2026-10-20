import { describe, expect, it } from "vitest";
import { defaultInserts, defaultLinear, defaultVector, docData, openingFromSelection, selectionOf } from "./defaults";
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
