import { describe, expect, it } from "vitest";
import type { OveAnnotation } from "../api";
import {
  cargoOptions,
  chooseCargo,
  findItrs,
  looksLikeItr,
  SC_LIMIT_BP,
  spanLength,
  SS_LIMIT_BP,
  verdictFor,
} from "./aavPackaging";

// Coordinates below are the ones the GenBank records carry, converted to OVE's 0-based
// inclusive numbering: PZ629150.1 (7079 bp) annotates its ITRs at 1..131 and 4343..4483,
// MT612431.1 (7204 bp) at 20..131 and 4477..4588.
const itr = (name: string, start: number, end: number, extra: object = {}): OveAnnotation => ({
  name, type: "repeat_region", start, end, ...extra,
});

const PZ_FEATURES: OveAnnotation[] = [
  itr("AAV2 inverted terminal repeat (ITR); AAV2 ITR (5')", 0, 130,
      { notes: { rpt_type: ["inverted"] } }),
  { name: "SaCas9 fusion protein", type: "CDS", start: 580, end: 3924,
    notes: { translation: ["MAPKKKRKVGIHGVPAAKRNYILGLDIGITSVGYGIIDYETRDVIDAGVRLFKEANVITRDEN"] } },
  itr("AAV2 inverted terminal repeat (ITR); AAV2 ITR (3')", 4342, 4482,
      { notes: { rpt_type: ["inverted"] } }),
  { name: "f1 ori", type: "rep_origin", start: 4557, end: 5012, notes: {} },
  { name: "AmpR", type: "CDS", start: 5200, end: 6060, notes: {} },
  { name: "ori", type: "rep_origin", start: 6430, end: 7018, notes: {} },
];
const PZ_LENGTH = 7079;

describe("recognising ITRs", () => {
  it("accepts the two annotation styles the sample constructs use", () => {
    expect(looksLikeItr({ name: "5' ITR", type: "misc_feature", start: 0, end: 10, notes: {} }))
      .toBe(true);
    expect(looksLikeItr({ name: "AAV2 ITR (5')", type: "repeat_region", start: 0, end: 10,
                          notes: { rpt_type: ["inverted"] } })).toBe(true);
    expect(looksLikeItr({ name: "repeat", type: "repeat_region", start: 0, end: 10,
                          notes: { note: ["inverted terminal repeat"] } })).toBe(true);
  });

  it("does not take the letters ITR inside a protein translation for an ITR", () => {
    const found = findItrs(PZ_FEATURES, PZ_LENGTH);
    expect(found.map((f) => f.start)).toEqual([0, 4342]);
  });

  it("ignores a feature whose name merely contains the letters", () => {
    expect(looksLikeItr({ name: "mitral marker", type: "CDS", start: 0, end: 9, notes: {} }))
      .toBe(false);
  });

  it("measures each ITR, including one that crosses the origin", () => {
    const [first] = findItrs([itr("ITR", 7070, 130, {})], PZ_LENGTH);
    expect(first.length).toBe(140);     // 7070..7078 plus 0..130
  });
});

describe("spans", () => {
  it("counts both ends of a plain span", () => {
    expect(spanLength(0, 130, PZ_LENGTH)).toBe(131);
  });

  it("wraps through the origin when the end comes first", () => {
    expect(spanLength(4342, 130, PZ_LENGTH)).toBe(7079 - 4342 + 131);
  });
});

describe("choosing the packaged side", () => {
  it("takes the side without the origin and the resistance marker", () => {
    const [a, b] = findItrs(PZ_FEATURES, PZ_LENGTH);
    const options = cargoOptions([a, b], { circular: true, features: PZ_FEATURES }, PZ_LENGTH);
    const { cargo, reason } = chooseCargo(options);
    // ITR 1..131 through to the end of ITR 4343..4483: 4483 bp, ITRs included.
    expect(cargo.totalBp).toBe(4483);
    expect(cargo.insertBp).toBe(4483 - 131 - 141);
    expect(cargo.crossesOrigin).toBe(false);
    expect(cargo.backbone).toEqual([]);
    expect(reason).toContain("backbone");
  });

  it("offers the span across the origin as the alternative on a circular plasmid", () => {
    const [a, b] = findItrs(PZ_FEATURES, PZ_LENGTH);
    const options = cargoOptions([a, b], { circular: true, features: PZ_FEATURES }, PZ_LENGTH);
    expect(options).toHaveLength(2);
    const other = options[1];
    expect(other.crossesOrigin).toBe(true);
    expect(other.totalBp).toBe(7079 - 4342 + 131);
    expect(other.backbone).toContain("ori");
  });

  it("keeps a cassette that really does cross the origin", () => {
    // Same plasmid rotated so that the cassette spans position 1: ITRs at 6500..6630 and 3700.
    const features: OveAnnotation[] = [
      itr("5' ITR", 6499, 6629),
      itr("3' ITR", 3699, 3839),
      { name: "ori", type: "rep_origin", start: 4500, end: 5088, notes: {} },
    ];
    const [first, second] = findItrs(features, PZ_LENGTH);
    const options = cargoOptions([first, second], { circular: true, features }, PZ_LENGTH);
    const { cargo } = chooseCargo(options);
    expect(cargo.crossesOrigin).toBe(true);
    expect(cargo.from.start).toBe(6499);
    expect(cargo.totalBp).toBe(spanLength(6499, 3839, PZ_LENGTH));
    expect(cargo.backbone).toEqual([]);
  });

  it("gives a linear sequence a single reading", () => {
    const features: OveAnnotation[] = [itr("5' ITR", 19, 130), itr("3' ITR", 4476, 4587)];
    const [a, b] = findItrs(features, 7204);
    const options = cargoOptions([a, b], { circular: false, features }, 7204);
    expect(options).toHaveLength(1);
    expect(options[0].totalBp).toBe(4569);
  });
});

describe("packaging verdict", () => {
  it("passes a cassette that fits a single-stranded vector", () => {
    const verdict = verdictFor(4483, "ss");
    expect(verdict.level).toBe("ok");
    expect(verdict.detail).toContain((SS_LIMIT_BP - 4483).toLocaleString());
  });

  it("warns between the limit and the size where genomes come out truncated", () => {
    expect(verdictFor(4900, "ss").level).toBe("tight");
    expect(verdictFor(5200, "ss").level).toBe("over");
  });

  it("halves the room for a self-complementary vector", () => {
    expect(SC_LIMIT_BP).toBe(2350);
    expect(verdictFor(2300, "sc").level).toBe("ok");
    expect(verdictFor(4483, "sc").level).toBe("over");
    expect(verdictFor(4483, "sc").headline).toContain("self-complementary");
  });

  it("reports how full the capsid is", () => {
    expect(verdictFor(SS_LIMIT_BP, "ss").load).toBeCloseTo(1);
    expect(verdictFor(SS_LIMIT_BP / 2, "ss").load).toBeCloseTo(0.5);
  });
});
