import { describe, expect, it } from "vitest";
import { enzymeRows, methylationWarnings } from "../enzymeLogic";

// ClaI (ATCGAT) is blocked by Dam when followed by the right context: GATCGAT contains GATC
// whose methylated adenine sits inside the ClaI site.
const FILLER = "TTGCCTACGGTACATGCAACGTTGACCATAGCATTACGGACCATTGACCAGTTACGGCATAC";
const WITH_DAM = FILLER + "G" + "ATCGAT" + FILLER;          // GATCGAT
// A ClaI site followed by C would itself create GATC — the classic Dam surprise — so the
// control puts a T there instead.
const WITHOUT_DAM = FILLER + "C" + "ATCGAT" + "T" + FILLER;
const WITH_DCM = FILLER + "CC" + "AGG" + FILLER;            // CCAGG, a Dcm site

function warnings(sequence: string) {
  const rows = enzymeRows("commercial", sequence, false);
  return methylationWarnings(sequence, false, rows);
}

describe("methylation warnings", () => {
  it("flags a ClaI site whose adenine is Dam methylated", () => {
    expect(warnings(WITH_DAM).ClaI).toEqual([
      { kind: "Dam", status: "blocked", blocking: true },
    ]);
  });

  it("flags ClaI when the base after the site creates a GATC", () => {
    expect(warnings(FILLER + "C" + "ATCGAT" + "C" + FILLER).ClaI?.[0].status).toBe("blocked");
  });

  it("stays quiet about an enzyme REBASE records as cut despite the overlap", () => {
    // BamHI's GGATCC always contains a Dam site and is cut anyway.
    const found = warnings(FILLER + "GGATCC" + FILLER);
    expect(found.BamHI).toBeUndefined();
  });

  it("does not flag the same enzyme in a context without GATC", () => {
    expect(warnings(WITHOUT_DAM).ClaI).toBeUndefined();
  });

  it("flags enzymes affected by an overlapping Dcm site", () => {
    const found = warnings(WITH_DCM);
    const flagged = Object.values(found).flat().filter((n) => n.kind === "Dcm");
    expect(flagged.length).toBeGreaterThan(0);
  });

  it("says so when REBASE records no sensitivity for an overlapping site", () => {
    const notes = Object.values(warnings(WITH_DAM)).flat();
    const unknown = notes.filter((n) => n.status === "not recorded");
    expect(unknown.every((n) => n.blocking === false)).toBe(true);
  });

  it("reports nothing when the sequence has neither motif", () => {
    expect(warnings("ACGTACGTACGTACGTACGTACGTACGTACGT")).toEqual({});
  });
});
