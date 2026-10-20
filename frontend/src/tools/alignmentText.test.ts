import { describe, expect, it } from "vitest";
import {
  choicePrompt,
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

  it("ticks nothing by default when several files could be compared", () => {
    expect(defaultChoice(docs, "a", DEFAULT_LIMITS)).toEqual([]);
    const many = Array.from({ length: 15 }, (_, i) => ({ editorName: `d${i}`, length: 1000 }));
    expect(defaultChoice(many, "d0", DEFAULT_LIMITS)).toEqual([]);
  });

  it("ticks the other file when it is the only one that can be compared", () => {
    expect(defaultChoice(docs.slice(0, 2), "a", DEFAULT_LIMITS)).toEqual(["b"]);
    // "c" is empty and "d" too long: "b" is the only usable one.
    expect(defaultChoice(docs.slice(0, 4), "a", DEFAULT_LIMITS)).toEqual(["b"]);
    expect(defaultChoice(docs.slice(0, 1), "a", DEFAULT_LIMITS)).toEqual([]);
  });

  it("explains why a document cannot be used", () => {
    expect(ineligible(0, DEFAULT_LIMITS)).toBe("empty");
    expect(ineligible(10, DEFAULT_LIMITS)).toBe("too short (under 15 bp)");
    expect(ineligible(20000, DEFAULT_LIMITS)).toBe("too long (20,000 bp; the limit is 15,000 bp)");
    expect(ineligible(5000, DEFAULT_LIMITS)).toBeNull();
  });

  it("prompts, without an error, while nothing is ticked", () => {
    expect(choicePrompt(0)).toMatch(/^Tick the sequences to compare/);
    expect(choicePrompt(1)).toBeNull();
    expect(choiceProblem(0, DEFAULT_LIMITS)).toBeNull();
  });

  it("allows at most the limit", () => {
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
