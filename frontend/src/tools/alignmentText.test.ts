import { describe, expect, it } from "vitest";
import {
  choiceProblem,
  DEFAULT_LIMITS,
  defaultChoice,
  describeDifference,
  differencesAsTsv,
  ineligible,
  summarizeCounts,
} from "./alignmentText";

describe("describeDifference", () => {
  it("writes each kind of difference in reference numbering", () => {
    expect(describeDifference({ kind: "mismatch", start: 1204, end: 1204, reference: "A", sequence: "G" }))
      .toBe("1,204: A → G");
    expect(describeDifference({ kind: "deletion", start: 601, end: 604, reference: "ACGT", sequence: "" }))
      .toBe("601–604: 4 bp missing (ACGT)");
    expect(describeDifference({ kind: "deletion", start: 9, end: 9, reference: "T", sequence: "" }))
      .toBe("9: 1 bp missing (T)");
    expect(describeDifference({ kind: "insertion", start: 1001, end: 1002, reference: "", sequence: "GATTACA" }))
      .toBe("between 1,001 and 1,002: 7 bp extra (GATTACA)");
  });

  it("shortens long stretches", () => {
    const d = { kind: "insertion" as const, start: 5, end: 6, reference: "", sequence: "A".repeat(40) };
    expect(describeDifference(d)).toBe("between 5 and 6: 40 bp extra (AAAAAAAAAAAA…)");
  });
});

describe("summarizeCounts", () => {
  it("says none when the sequences agree", () => {
    expect(summarizeCounts({ mismatches: 0, insertions: 0, deletions: 0, inserted_bases: 0, deleted_bases: 0 }))
      .toBe("none");
  });
  it("counts each kind with singular and plural", () => {
    expect(summarizeCounts({ mismatches: 3, insertions: 1, deletions: 2, inserted_bases: 7, deleted_bases: 1500 }))
      .toBe("3 mismatches · 1 insertion (7 bp) · 2 deletions (1,500 bp)");
    expect(summarizeCounts({ mismatches: 1, insertions: 0, deletions: 0, inserted_bases: 0, deleted_bases: 0 }))
      .toBe("1 mismatch");
  });
});

describe("choosing sequences", () => {
  const docs = [
    { editorName: "a", length: 5000 },
    { editorName: "b", length: 4000 },
    { editorName: "c", length: 0 },
    { editorName: "d", length: 20000 },
    { editorName: "e", length: 3000 },
  ];

  it("ticks every other usable document by default", () => {
    expect(defaultChoice(docs, "a", DEFAULT_LIMITS)).toEqual(["b", "e"]);
  });

  it("never ticks more than the limit allows", () => {
    const many = Array.from({ length: 15 }, (_, i) => ({ editorName: `d${i}`, length: 1000 }));
    expect(defaultChoice(many, "d0", DEFAULT_LIMITS)).toHaveLength(DEFAULT_LIMITS.max_sequences - 1);
  });

  it("explains why a document cannot be used", () => {
    expect(ineligible(0, DEFAULT_LIMITS)).toBe("empty");
    expect(ineligible(10, DEFAULT_LIMITS)).toBe("too short (under 15 bp)");
    expect(ineligible(20000, DEFAULT_LIMITS)).toBe("too long (20,000 bp; the limit is 15,000 bp)");
    expect(ineligible(5000, DEFAULT_LIMITS)).toBeNull();
  });

  it("asks for at least one and at most the limit", () => {
    expect(choiceProblem(0, DEFAULT_LIMITS)).toMatch(/at least one/);
    expect(choiceProblem(3, DEFAULT_LIMITS)).toBeNull();
    expect(choiceProblem(11, DEFAULT_LIMITS)).toBe("Too many: at most 9 besides the reference. Untick 2.");
  });
});

it("copies differences as tab-separated lines", () => {
  const text = differencesAsTsv("pX", [
    { kind: "mismatch", start: 10, end: 10, reference: "A", sequence: "C" },
  ]);
  expect(text.split("\n")).toEqual([
    "sequence\tkind\tstart\tend\treference bases\tsequence bases",
    "pX\tmismatch\t10\t10\tA\tC",
  ]);
});
