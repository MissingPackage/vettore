import { describe, expect, it } from "vitest";
import { parseEnzymeList } from "../preferences";
import { blockedByMethylation, protocolText, type PlannerResult } from "./digestPlannerText";

describe("blockedByMethylation", () => {
  it("names an enzyme whose site overlaps a Dam site that blocks it", () => {
    // ClaI ATCGAT preceded by G: GATCGAT carries a Dam GATC over the ClaI site (blocked).
    const blocked = blockedByMethylation([{ sequence: `${"A".repeat(30)}GATCGAT${"C".repeat(30)}`, circular: false }]);
    expect(blocked).toContain("ClaI");
  });

  it("leaves out enzymes that methylation does not stop", () => {
    // BamHI GGATCC always holds a GATC and is cut anyway.
    const blocked = blockedByMethylation([{ sequence: `${"A".repeat(30)}GGATCC${"T".repeat(30)}`, circular: false }]);
    expect(blocked).not.toContain("BamHI");
  });
});

describe("protocolText", () => {
  const result: PlannerResult = {
    enzymes: ["EcoRI", "PstI"],
    buffers: [null, null],
    mixed_buffers: false,
    same_pattern: [["XmnI", "PstI"]],
    sentence: "",
    clear: true,
    min_separation: 0.2,
    lanes: [
      { name: "Expected clone", sizes: [3236, 1125, 40], uncut: false, text: "", bands: [] },
      { name: "parent vector", sizes: [4361], uncut: true, text: "", bands: [] },
    ],
  };

  it("lists the enzymes, the gel and the bands of every outcome", () => {
    const text = protocolText("pBR322 + cat", result, { percent: 1, ladderLabel: "NEB 1 kb Plus (N3200)" }, [true, true]);
    expect(text).toContain("Diagnostic digest for “pBR322 + cat”");
    expect(text).toContain("Digest: EcoRI + PstI (double digest)");
    expect(text).toContain("Gel: 1.0 % agarose, ladder NEB 1 kb Plus (N3200)");
    expect(text).toContain("Expected clone: 3,236 bp + 1,125 bp + 40 bp (too small to see)");
    expect(text).toContain("parent vector: 4,361 bp, not cut (a supercoiled plasmid runs ahead of its size)");
    expect(text).toContain("Same bands with: XmnI + PstI");
    expect(text).not.toContain("Note:");
  });

  it("names each enzyme's buffer and says what to do when they differ", () => {
    const mixed = { ...result, buffers: ["rCutSmart", "NEBuffer r3.1"], mixed_buffers: true };
    const text = protocolText("pX", mixed, { percent: 1, ladderLabel: "L" }, [true, true]);
    expect(text).toContain("Digest: EcoRI (rCutSmart) + PstI (NEBuffer r3.1) (double digest)");
    expect(text).toContain("Note: EcoRI and PstI go in different buffers");
  });
});

describe("parseEnzymeList", () => {
  it("reads two columns copied from a spreadsheet as enzyme and buffer", () => {
    const parsed = parseEnzymeList("EcoRI\trCutSmart\nsali\tNEBuffer r3.1\nPstI");
    expect(parsed.known).toEqual(["EcoRI", "SalI", "PstI"]);
    expect(parsed.buffers).toEqual({ EcoRI: "rCutSmart", SalI: "NEBuffer r3.1" });
    expect(parsed.unknown).toEqual([]);
  });

  it("still reads a plain list and points out unknown names", () => {
    const parsed = parseEnzymeList("EcoRI, PstI; Hind3");
    expect(parsed.known).toEqual(["EcoRI", "PstI"]);
    expect(parsed.unknown).toEqual(["Hind3"]);
    expect(parsed.buffers).toEqual({});
  });
});
