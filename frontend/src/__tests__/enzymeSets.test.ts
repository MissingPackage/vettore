import { describe, expect, it } from "vitest";
import { enzymeRows, enzymesInSet, specifiedBases } from "../enzymeLogic";

// pUC19-like test fragment with known single sites: EcoRI GAATTC, BamHI GGATCC, HindIII AAGCTT.
const SEQ =
  "TTTTGAATTCAAAAAAAAAAGGATCCTTTTTTTTTTAAGCTTCCCCCCCCCCGGATCCAAAAA" +
  "GCGCAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTT";

describe("specifiedBases", () => {
  it("ignores N positions", () => {
    expect(specifiedBases("gccnnnnnggc")).toBe(6);
    expect(specifiedBases("gaattc")).toBe(6);
  });
});

describe("enzymesInSet", () => {
  it("unique 6+ cutters include single 6-bp sites and exclude double cutters", () => {
    const names = enzymesInSet("unique6", SEQ, false);
    expect(names).toContain("EcoRI");
    expect(names).toContain("HindIII");
    expect(names).not.toContain("BamHI"); // cuts twice
  });

  it("unique & dual cutters include the double cutter", () => {
    expect(enzymesInSet("uniqueDual", SEQ, false)).toContain("BamHI");
  });

  it("unique 6+ never includes 4-cutters", () => {
    const names = enzymesInSet("unique6", SEQ, false);
    expect(names).not.toContain("HhaI"); // GCGC, single site but 4 bp
    expect(enzymesInSet("unique", SEQ, false)).toContain("HhaI");
  });
});

describe("enzymeRows", () => {
  it("reports cut locations, ends and overhangs", () => {
    // EcoRI G^AATTC: site at 1-based 5..10, top-strand cut after base 5.
    const rows = enzymeRows("unique6", SEQ, false);
    const ecori = rows.find((r) => r.name === "EcoRI")!;
    expect(ecori.site).toBe("GAATTC");
    expect(ecori.ends).toBe("5′ overhang");
    expect(ecori.overhang).toBe("AATT");
    expect(ecori.cuts).toEqual([{ location: 5, start: 4, end: 9 }]);
  });

  it("lists every site of a dual cutter in order", () => {
    const bamhi = enzymeRows("uniqueDual", SEQ, false).find((r) => r.name === "BamHI")!;
    expect(bamhi.cuts.map((c) => c.location)).toEqual([21, 53]);
  });
});
