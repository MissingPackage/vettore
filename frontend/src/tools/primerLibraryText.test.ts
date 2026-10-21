import { describe, expect, it } from "vitest";
import {
  bindingText,
  cleanSequence,
  filterPrimers,
  mergeImport,
  newEntry,
  orderingText,
  range,
  heldText,
  rejectedText,
  saveChanges,
  sequenceAlreadyIn,
  today,
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
  it("reads a 5′-…-3′ frame as the engine does", () => {
    expect(cleanSequence("5'-ACG TAC-3'")).toBe("ACGTAC");
    expect(cleanSequence("5′-ACGTAC-3′")).toBe("ACGTAC");
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
  it("uses the library's separator, quoting cells that contain it", () => {
    const csv = toCsv([primer({ name: "M13F", sequence: "ACGT", notes: "5 µl; 4 °C, box 1" })], ";");
    const [header, row] = csv.trim().split("\r\n");
    expect(header).toBe("name;sequence;notes;location;added");
    expect(row).toBe('M13F;ACGT;"5 µl; 4 °C, box 1";;');
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
    expect(e).toMatchObject({ name: "KanR_F", sequence: "ACGT", added: "2026-10-21", location: "box 1", length: 4, fresh: true });
  });
});

describe("saveChanges", () => {
  const a = primer({ name: "A", sequence: "AAAACCCCGG", notes: "n" });
  const b = primer({ name: "B", sequence: "TTTTGGGGCC" });
  it("sends only what changed: additions, changed fields, removals", () => {
    const c = primer({ name: "C", sequence: "GGGGAAAATT" });
    const changes = saveChanges([a, b], [{ ...a, notes: "new" }, c]);
    expect(changes.added.map((p) => p.name)).toEqual(["C"]);
    expect(changes.edited).toEqual([{ sequence: "AAAACCCCGG", was: { name: "A", notes: "n", location: "", added: "" }, twins: 1, notes: "new" }]);
    expect(changes.removed).toEqual([{ sequence: "TTTTGGGGCC", was: { name: "B", notes: "", location: "", added: "" }, twins: 1 }]);
  });
  it("marks a primer removed and added again, so it comes back", () => {
    const again = { ...b, name: "B2", fresh: true };
    const [e] = saveChanges([a, b], [a, again]).edited;
    expect(e).toMatchObject({ sequence: "TTTTGGGGCC", readd: true, name: "B2", was: { name: "B" } });
  });
  it("counts a second tube of the oligo in the file", () => {
    const twin = { name: "B copy", sequence: "ttttggggcc", reason: "it is the same oligo as “B”." };
    expect(saveChanges([a, b], [a], [twin]).removed[0].twins).toBe(2);
  });
  it("sends nothing when nothing changed", () => {
    expect(saveChanges([a, b], [a, b])).toEqual({ added: [], edited: [], removed: [] });
  });
});

describe("bindingText in a repeat", () => {
  it("says the primer primes at several positions", () => {
    expect(bindingText({ start: 0, end: 19, strand: 1, footprint: 20, tail: 0, mismatches: 0, positions: 6 }))
      .toContain("primes at 6 overlapping positions");
  });
});

describe("heldText", () => {
  it("names the rows and says nothing in them is lost", () => {
    const text = heldText([{ name: "M13F", sequence: "GTAAAACGACGGCCAGT", line: 3 }]);
    expect(text).toContain("“M13F” (line 3)");
    expect(text).toContain("nothing is lost");
  });
});

describe("today and rejectedText", () => {
  it("writes the local date as YYYY-MM-DD", () => {
    expect(today(new Date(2026, 0, 5, 23, 30))).toBe("2026-01-05");
  });
  it("says why rows were left out, and that file rows are kept", () => {
    const rows = [{ name: "Ino", sequence: "ACGI", reason: "its sequence has characters that are not DNA ('I')." }];
    expect(rejectedText(rows, "paste")).toContain("not DNA");
    expect(rejectedText(rows, "file")).toContain("stay in the file");
    const dup = [{ name: "M13F copy", sequence: "GTAAAACG", reason: "it is the same oligo as “M13F”." }];
    expect(rejectedText(dup, "paste")).not.toContain("not DNA");
  });
});
