import { describe, expect, it } from "vitest";
import {
  bindingText,
  cleanSequence,
  filterPrimers,
  mergeImport,
  newEntry,
  orderingText,
  range,
  sequenceAlreadyIn,
  toCsv,
  type LibraryPrimer,
} from "./primerLibraryText";

const primer = (over: Partial<LibraryPrimer>): LibraryPrimer => ({
  name: "P", sequence: "ACGTACGT", notes: "", location: "", added: "", length: 8, tm: null, ...over,
});

describe("cleanSequence", () => {
  it("upper-cases, reads U as T and drops spaces", () => {
    expect(cleanSequence(" acg uac gt ")).toBe("ACGTACGT");
  });
  it("returns empty for anything that is not DNA", () => {
    expect(cleanSequence("ACGT123")).toBe("");
    expect(cleanSequence("hello")).toBe("");
  });
  it("keeps IUPAC ambiguity codes", () => {
    expect(cleanSequence("acgryn")).toBe("ACGRYN");
  });
});

describe("filterPrimers", () => {
  const primers = [
    primer({ name: "M13F", sequence: "GTAAAACGACGGCCAGT", location: "box 3 A4" }),
    primer({ name: "KanR_R", sequence: "TTAGAAAAACTCATCGAGCA", notes: "kanamycin" }),
  ];
  it("keeps everything for an empty query", () => {
    expect(filterPrimers(primers, "  ")).toHaveLength(2);
  });
  it("matches name, notes and location case-insensitively", () => {
    expect(filterPrimers(primers, "m13")).toHaveLength(1);
    expect(filterPrimers(primers, "kanamycin")).toHaveLength(1);
    expect(filterPrimers(primers, "box 3")).toHaveLength(1);
  });
  it("matches a sequence fragment regardless of case and spaces", () => {
    expect(filterPrimers(primers, "gtaaa")[0].name).toBe("M13F");
    expect(filterPrimers(primers, "acgt")).toHaveLength(0);
  });
});

describe("sequenceAlreadyIn", () => {
  const primers = [primer({ name: "A", sequence: "ACGTACGTAC" })];
  it("finds a duplicate regardless of case and spaces", () => {
    expect(sequenceAlreadyIn(primers, "acgt acgt ac")?.name).toBe("A");
  });
  it("returns null when the sequence is new or invalid", () => {
    expect(sequenceAlreadyIn(primers, "TTTTTTTT")).toBeNull();
    expect(sequenceAlreadyIn(primers, "nope")).toBeNull();
  });
});

describe("range", () => {
  it("is 1-based inclusive", () => {
    expect(range(0, 15)).toBe("1–16");
  });
  it("notes a wrap across the origin", () => {
    expect(range(50, 5)).toContain("across the origin");
  });
});

describe("bindingText", () => {
  it("describes an exact forward binding", () => {
    expect(bindingText({ start: 0, end: 15, strand: 1, footprint: 16, tail: 0, mismatches: 0 }))
      .toBe("→ forward, 1–16 (exact)");
  });
  it("names mismatches and a tail on a reverse binding", () => {
    const text = bindingText({ start: 20, end: 5, strand: -1, footprint: 16, tail: 4, mismatches: 2 });
    expect(text).toContain("← reverse");
    expect(text).toContain("2 mismatches");
    expect(text).toContain("4 nt 5′ tail");
  });
});

describe("toCsv", () => {
  it("writes a header and quotes cells with commas or quotes", () => {
    const csv = toCsv([primer({ name: "M13F", sequence: "ACGT", notes: 'reads, "universal"', location: "box 1" })]);
    const [header, row] = csv.trim().split("\r\n");
    expect(header).toBe("name,sequence,notes,location,added");
    expect(row).toContain('"reads, ""universal"""');
    expect(row.startsWith("M13F,ACGT,")).toBe(true);
  });
});

describe("mergeImport", () => {
  const existing = [primer({ name: "A", sequence: "ACGTACGTAC" })];
  it("adds new primers and skips duplicate sequences", () => {
    const incoming = [primer({ name: "A-copy", sequence: "ACGTACGTAC" }), primer({ name: "B", sequence: "TTTTGGGGCC" })];
    const { merged, added, duplicates } = mergeImport(existing, incoming);
    expect(added).toBe(1);
    expect(duplicates).toBe(1);
    expect(merged.map((p) => p.name)).toEqual(["A", "B"]);
  });
});

describe("orderingText and newEntry", () => {
  it("copies name⇥sequence lines", () => {
    expect(orderingText([{ name: "F", sequence: "ACGT" }, { name: "R", sequence: "TTTT" }]))
      .toBe("F\tACGT\nR\tTTTT");
  });
  it("stamps a new entry with today's date and cleans the sequence", () => {
    const e = newEntry("  KanR_F ", "acg t", "2026-10-21", "note", "box 1");
    expect(e).toMatchObject({ name: "KanR_F", sequence: "ACGT", added: "2026-10-21", location: "box 1", length: 4 });
  });
});
