import { describe, expect, it } from "vitest";
import { cleanSequence } from "../NewSequenceDialog";

describe("cleanSequence", () => {
  it("strips FASTA headers, digits and whitespace and upper-cases", () => {
    expect(cleanSequence(">seq1\n1 acgt acgt\n11 nnAC\n")).toEqual({ sequence: "ACGTACGTNNAC", rejected: 0 });
  });
  it("counts and drops non-IUPAC characters", () => {
    expect(cleanSequence("ACGT*XZ")).toEqual({ sequence: "ACGT", rejected: 3 });
  });
});
