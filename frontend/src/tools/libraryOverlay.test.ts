import { describe, expect, it } from "vitest";
import { OVERLAY_NOTE, type BoundPrimer } from "./primerLibraryText";
import { isOverlay, overlayAnnotations, removeOverlay, withOverlay } from "./libraryOverlay";

const bound: BoundPrimer[] = [
  {
    name: "M13F", sequence: "GTAAAACGACGGCCAGT", notes: "", location: "box 1", added: "", length: 17, tm: 55,
    sites: 1,
    bindings: [{ start: 10, end: 26, strand: 1, footprint: 17, tail: 0, mismatches: 0 }],
  },
  {
    name: "Rep", sequence: "ACGTACGTACGTACGTAC", notes: "", location: "", added: "", length: 18, tm: null,
    sites: 2,
    bindings: [
      { start: 30, end: 47, strand: 1, footprint: 18, tail: 0, mismatches: 0 },
      { start: 73, end: 90, strand: -1, footprint: 18, tail: 0, mismatches: 1 },
    ],
  },
];

describe("overlayAnnotations", () => {
  it("makes one tagged primer per binding and numbers repeats", () => {
    const anns = overlayAnnotations(bound);
    expect(anns).toHaveLength(3);
    expect(anns.every(isOverlay)).toBe(true);
    expect(anns[0].name).toBe("M13F");                 // a single binding keeps the bare name
    expect(anns.map((a) => a.name)).toContain("Rep (1)");
    expect(anns.map((a) => a.name)).toContain("Rep (2)");
    const rev = anns.find((a) => a.strand === -1)!;
    expect(rev.forward).toBe(false);
    // Drawn where it anneals: the lower coordinate first, as every primer in the editor.
    expect([rev.start, rev.end]).toEqual([73, 90]);
    expect((rev.notes as Record<string, string[]>)[OVERLAY_NOTE]).toEqual(["library"]);
  });
  it("carries the primer sequence and location in the note", () => {
    const [m13] = overlayAnnotations(bound);
    expect(m13.bases).toBe("GTAAAACGACGGCCAGT");
    expect((m13.notes as Record<string, string[]>).note[0]).toContain("box 1");
  });
});

describe("a tailed primer", () => {
  it("draws only the bases that anneal and names the whole primer in the note", () => {
    const tailed: BoundPrimer = {
      name: "BamF", sequence: "CGCGGATCC" + "GTAAAACGACGGCCAGT", notes: "", location: "", added: "",
      length: 26, tm: null, sites: 1,
      bindings: [{ start: 100, end: 116, strand: 1, footprint: 17, tail: 9, mismatches: 0 }],
    };
    const [a] = overlayAnnotations([tailed]);
    expect(a.bases).toBe("GTAAAACGACGGCCAGT");
    expect(a.bases!.length).toBe(a.end - a.start + 1);
    expect((a.notes as Record<string, string[]>).note[0]).toContain("CGCGGATCCGTAAAACGACGGCCAGT");
  });
});

describe("removeOverlay / withOverlay", () => {
  const own: any = [
    { name: "real-fwd", type: "primer_bind", start: 0, end: 15, strand: 1 },
    { name: "stale-overlay", type: "primer_bind", start: 5, end: 20, strand: 1, notes: { [OVERLAY_NOTE]: ["library"] } },
  ];
  it("keeps only the document's own primers", () => {
    const kept = removeOverlay(own);
    expect(kept.map((p) => p.name)).toEqual(["real-fwd"]);
  });
  it("replaces any previous overlay with a fresh one", () => {
    const merged = withOverlay(own, bound);
    const ownNames = merged.filter((p) => !isOverlay(p)).map((p) => p.name);
    expect(ownNames).toEqual(["real-fwd"]);              // stale overlay dropped, real one kept
    expect(merged.filter(isOverlay)).toHaveLength(3);
  });
  it("accepts a record of annotations as well as a list", () => {
    const asRecord = { a: own[0], b: own[1] };
    expect(removeOverlay(asRecord).map((p) => p.name)).toEqual(["real-fwd"]);
  });
});
