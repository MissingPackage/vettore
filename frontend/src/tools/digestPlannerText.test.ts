import { describe, expect, it } from "vitest";
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
  });
});
