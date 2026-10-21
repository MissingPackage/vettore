import { describe, expect, it } from "vitest";
import {
  conservationBins,
  conservationLevel,
  DEFAULT_MULTIPLE_LIMITS,
  defaultMultipleChoice,
  describeBin,
  ineligibleForMultiple,
  multipleChoiceProblem,
  multiplePrompt,
  orientationNote,
  progressText,
  variableRegions,
  type MultipleRow,
} from "./multipleText";

const limits = DEFAULT_MULTIPLE_LIMITS;

describe("choice", () => {
  it("refuses empty, short and long documents with the reason", () => {
    expect(ineligibleForMultiple(0, limits)).toBe("empty");
    expect(ineligibleForMultiple(10, limits)).toBe("too short (under 15 bp)");
    expect(ineligibleForMultiple(20000, limits)).toBe("too long (20,000 bp; the limit is 15,000 bp)");
    expect(ineligibleForMultiple(5000, limits)).toBeNull();
  });

  it("ticks both when only two can take part, else only the active one", () => {
    const two = [{ editorName: "a", length: 500 }, { editorName: "b", length: 600 }, { editorName: "c", length: 0 }];
    expect(defaultMultipleChoice(two, "a", limits)).toEqual(["a", "b"]);
    const three = [...two, { editorName: "d", length: 700 }];
    expect(defaultMultipleChoice(three, "d", limits)).toEqual(["d"]);
    expect(defaultMultipleChoice(three, "c", limits)).toEqual([]);
  });

  it("explains too many sequences or too many bases", () => {
    expect(multipleChoiceProblem(Array(52).fill(100), limits)).toBe("Too many: at most 50. Untick 2.");
    expect(multipleChoiceProblem(Array(11).fill(15000), limits))
      .toBe("Together these are 165,000 bp; at most 150,000 bp can be aligned at once. "
        + "Untick some, or align the regions you need.");
    expect(multipleChoiceProblem(Array(10).fill(15000), limits)).toBeNull();
  });

  it("prompts until two are ticked", () => {
    expect(multiplePrompt(1)).toContain("Tick at least two");
    expect(multiplePrompt(2)).toBeNull();
  });
});

describe("conservation", () => {
  it("folds columns into bars with their average and lowest value", () => {
    const bins = conservationBins([100, 100, 50, 100, 20, 80], 3);
    expect(bins).toEqual([
      { start: 1, end: 2, mean: 100, min: 100 },
      { start: 3, end: 4, mean: 75, min: 50 },
      { start: 5, end: 6, mean: 50, min: 20 },
    ]);
    expect(conservationBins([90, 80], 400)).toHaveLength(2);
  });

  it("describes a bar in words", () => {
    expect(describeBin({ start: 1201, end: 1225, mean: 84.2, min: 50 }))
      .toBe("Columns 1,201–1,225: 84 % identical on average (lowest 50 %)");
    expect(describeBin({ start: 7, end: 7, mean: 100, min: 100 })).toBe("Column 7: 100 % identical on average");
  });

  it("grades bars", () => {
    expect(conservationLevel(95)).toBe("high");
    expect(conservationLevel(70)).toBe("mid");
    expect(conservationLevel(30)).toBe("low");
  });

  it("finds the longest variable stretches", () => {
    const c = [...Array(5).fill(100), ...Array(12).fill(40), ...Array(5).fill(100), ...Array(20).fill(50), 100];
    expect(variableRegions(c)).toEqual([{ start: 23, end: 42 }, { start: 6, end: 17 }]);
    expect(variableRegions(Array(50).fill(100))).toEqual([]);
  });
});

describe("wording", () => {
  it("says how far the alignment has got", () => {
    expect(progressText("Aligning: step 3 of 9", 0.416)).toBe("Aligning: step 3 of 9… 42 %");
  });

  it("describes strand and origin", () => {
    const row = { strand: -1, origin_moved_to: 1237 } as MultipleRow;
    expect(orientationNote(row)).toBe("reverse complement, origin moved to its base 1,237");
    expect(orientationNote({ strand: 1, origin_moved_to: null } as MultipleRow)).toBe("forward");
  });
});
