import { describe, expect, it } from "vitest";
import { type Contig, type ContigDifference, contigVerdict, coverageSentence, evidence, listed } from "./contigText";

function contig(over: Partial<Contig> = {}): Contig {
  return {
    reads: 2,
    consensus_length: 1125,
    coverage: [
      { start: 1, end: 222, strands: 0 },
      { start: 223, end: 861, strands: 2 },
      { start: 862, end: 903, strands: 1 },
      { start: 904, end: 1125, strands: 0 },
    ],
    covered: [{ start: 223, end: 903 }],
    verified: [{ start: 223, end: 903, reads: 1 }],
    unclear: 0,
    doubtful: 0,
    differences: [],
    disagreements: [],
    ...over,
  };
}

const mismatch = (over: Partial<ContigDifference> = {}): ContigDifference => ({
  kind: "mismatch", position: 523, reference: "G", consensus: "A", reads: 2, quality: 50, confirmed: true, ...over,
});

describe("coverageSentence", () => {
  it("names the stretches on both strands, on one, and without a read, in reference numbering", () => {
    expect(coverageSentence(contig())).toBe(
      "223–861 covered on both strands; 862–903 on one strand only; no read over 1–222 and 904–1,125.");
  });

  it("says so when every base is read on both strands", () => {
    expect(coverageSentence(contig({ coverage: [{ start: 1, end: 900, strands: 2 }] })))
      .toBe("Every base is covered on both strands.");
  });

  it("keeps a long list short", () => {
    expect(listed(["a", "b", "c", "d", "e", "f"])).toBe("a, b, c and 3 more stretches");
  });
});

describe("contigVerdict", () => {
  it("passes a clone whose consensus matches", () => {
    const v = contigVerdict(contig(), "pCOI", 61.4);
    expect(v.ok).toBe(true);
    expect(v.text).toBe("No differences: the consensus of 2 reads matches pCOI wherever they cover it (61.4 % of pCOI).");
  });

  it("separates confirmed differences from those to check in the trace", () => {
    const v = contigVerdict(contig({
      differences: [mismatch(), mismatch({ position: 880, reads: 1, quality: 29, confirmed: false })],
    }), "pCOI", 61.4);
    expect(v.ok).toBe(false);
    expect(v.text).toContain("2 differences between the consensus of 2 reads and pCOI");
    expect(v.text).toContain("1 is confirmed, 1 needs a look at the trace");
    expect(v.text).toContain("one read at Q30 or better, and no read at Q20 or better shows something else");
  });

  it("mentions positions where the reads disagreed, not those without a call", () => {
    const call = { read: "F", strand: 1, call: "A", quality: 40 };
    const v = contigVerdict(contig({
      disagreements: [
        { kind: "base", position: 300, reference: "A", consensus: "A", calls: [call] },
        { kind: "no call", position: 895, reference: "N", consensus: "N", calls: [call] },
      ],
    }), "pCOI", 61.4);
    expect(v.text).toContain("The reads disagree at 1 place (listed under Contig)");
  });

  it("does not pass a clone where good reads conflict or only weak calls back a base", () => {
    const v = contigVerdict(contig({ doubtful: 3 }), "pCOI", 61.4);
    expect(v.ok).toBe(false);
    expect(v.text).toContain("3 covered positions are not verified: good reads conflict there, or only calls under Q20");
  });

  it("does not pass a clone where some covered bases have no clear call", () => {
    const v = contigVerdict(contig({ unclear: 2 }), "pCOI", 61.4);
    expect(v.ok).toBe(false);
    expect(v.text).toContain("2 covered positions have no clear call in any read");
  });
});

describe("evidence", () => {
  it("says why a difference is trusted or not", () => {
    expect(evidence(mismatch())).toBe("seen in 2 reads, best Q50");
    expect(evidence(mismatch({ reads: 1, quality: 38 }))).toBe("1 read at Q38");
    expect(evidence(mismatch({ reads: 1, quality: 22, confirmed: false }))).toBe("1 read at Q22, under Q30");
    expect(evidence(mismatch({ reads: 2, quality: 15, confirmed: false })))
      .toBe("seen in 2 reads, best Q15; another read shows something else, or the reads are under Q20");
  });
});
