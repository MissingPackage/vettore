import { describe, expect, it } from "vitest";
import type { SequenceData } from "./api";
import { cleanHistory, dataToSave, takeHistory, type HistoryEntry } from "./history";
import type { Doc } from "./documents";
import { buildLineage, lineageDepth, lineageText, sizeOf } from "./lineage";

function seq(name: string, length: number, circular = false): SequenceData {
  return { name, sequence: "A".repeat(length), circular, features: [], primers: [] };
}

function doc(name: string, history: HistoryEntry[], path: string | null = null): Doc {
  return {
    editorName: `doc-${name}`, name, path, saveInPlace: !!path, source: null,
    copyNoticeShown: false, neverSaved: !path, view: "lineage", enzymeSet: "unique6",
    history,
  };
}

/** A product of one cloning step on a vector that was itself made from two fragments. */
const GRANDPARENT: HistoryEntry[] = [
  { time: "19/10/2026", action: "Opened", detail: "fragments.gb" },
  {
    time: "19/10/2026", action: "Gibson assembly", detail: "2 fragments",
    parents: [
      { name: "fragment A", path: "/lab/a.gb", data: seq("fragment A", 800), history: [] },
      { name: "fragment B", path: null, data: seq("fragment B", 400), history: [] },
    ],
  },
];

const HISTORY: HistoryEntry[] = [
  { time: "20/10/2026", action: "Opened", detail: "backbone.gb" },
  {
    time: "21/10/2026", action: "Restriction cloning", detail: "EcoRI + BamHI",
    parents: [
      { name: "vector", path: "/lab/vector.gb", data: seq("vector", 4000, true),
        history: GRANDPARENT },
      { name: "insert", path: null, data: seq("insert", 1200), history: [] },
    ],
  },
];

describe("buildLineage", () => {
  it("puts the document at the top and its inputs under it", () => {
    const tree = buildLineage(doc("pProduct", HISTORY), seq("pProduct", 5200, true));
    expect(tree.name).toBe("pProduct");
    expect(tree.size).toBe("5,200 bp, circular");
    expect(tree.operation).toBe("Restriction cloning");
    expect(tree.children.map((c) => c.name)).toEqual(["vector", "insert"]);
  });

  it("follows each input back to what it was made from", () => {
    const tree = buildLineage(doc("pProduct", HISTORY), seq("pProduct", 5200, true));
    const [vector, insert] = tree.children;
    // The vector was itself assembled, so the tree carries on above it.
    expect(vector.operation).toBe("Gibson assembly");
    expect(vector.size).toBe("4,000 bp, circular");
    expect(vector.children.map((c) => c.name)).toEqual(["fragment A", "fragment B"]);
    expect(insert.children).toEqual([]);
    expect(lineageDepth(tree)).toBe(3);
  });

  it("names the step that was run on each input, and says when its own origin is unknown", () => {
    const [vector, insert] = buildLineage(doc("p", HISTORY), seq("p", 10)).children;
    expect(vector.operation).toBe("Gibson assembly");
    expect(vector.usedBy).toBe("Restriction cloning");
    // The insert was opened, pasted or typed somewhere else: the History does not say.
    expect(insert.operation).toBe("where it came from is not recorded");
    expect(insert.usedBy).toBe("Restriction cloning");
  });

  it("keeps the path of an input that came from a file, and none for one that did not", () => {
    const tree = buildLineage(doc("p", HISTORY), seq("p", 10));
    expect(tree.children[0].path).toBe("/lab/vector.gb");
    expect(tree.children[1].path).toBeNull();
  });

  it("marks an input whose copy the file could not keep", () => {
    const history: HistoryEntry[] = [{
      time: "t", action: "Restriction cloning",
      parents: [{ name: "vector", snapshots_dropped: true, history: [] }],
    }];
    const [vector] = buildLineage(doc("p", history), seq("p", 10)).children;
    expect(vector.dropped).toBe(true);
    expect(vector.size).toBe("size unknown");
  });

  it("stands alone for a document that was only opened", () => {
    const tree = buildLineage(doc("pBR322", [{ time: "t", action: "Opened" }]), seq("x", 4361, true));
    expect(tree.children).toEqual([]);
    expect(lineageDepth(tree)).toBe(1);
  });

  it("stops following the tree at 20 steps", () => {
    // A chain of 40 steps, each made from the one before it.
    let history: HistoryEntry[] = [{ time: "t", action: "Opened" }];
    for (let i = 0; i < 40; i++) {
      history = [{
        time: "t", action: `step ${i}`,
        parents: [{ name: `input ${i}`, data: seq(`input ${i}`, 100), history }],
      }];
    }
    expect(lineageDepth(buildLineage(doc("p", history), seq("p", 100)))).toBe(20);
  });

  it("writes the tree out with its shape as indentation", () => {
    const text = lineageText(buildLineage(doc("pProduct", HISTORY), seq("pProduct", 5200, true)));
    expect(text.split("\n")).toEqual([
      "pProduct — Restriction cloning (5,200 bp, circular)",
      "  vector — Gibson assembly (4,000 bp, circular), used by Restriction cloning [vector.gb]",
      "    fragment A — where it came from is not recorded (800 bp, linear), "
        + "used by Gibson assembly [a.gb]",
      "    fragment B — where it came from is not recorded (400 bp, linear), "
        + "used by Gibson assembly",
      "  insert — where it came from is not recorded (1,200 bp, linear), "
        + "used by Restriction cloning",
    ]);
  });
});

describe("sizeOf", () => {
  it("names the size and the topology, and says so when there is no copy", () => {
    expect(sizeOf(seq("x", 1234, true))).toBe("1,234 bp, circular");
    expect(sizeOf(undefined)).toBe("size unknown");
  });
});

describe("the History a file carries", () => {
  it("is taken off the data the editor keeps", () => {
    const { data, history } = takeHistory({ ...seq("p", 10), history: HISTORY });
    expect("history" in data).toBe(false);
    expect(history).toHaveLength(2);
    expect(history[1].parents?.[0].name).toBe("vector");
  });

  it("keeps only what looks like History, because a file is outside content", () => {
    expect(cleanHistory("nonsense")).toEqual([]);
    expect(cleanHistory([{ action: "" }, 7, null, { time: 1, action: "Opened" }])).toEqual([
      { time: "", action: "Opened" },
    ]);
    const [step] = cleanHistory([{
      time: "t", action: "Cloning",
      parents: [{ name: "v", data: { sequence: "ACGT" } }, { data: "not a document" }, 3],
    }]);
    expect(step.parents?.map((parent) => parent.name)).toEqual(["v", "an input"]);
    expect(step.parents?.[1].data).toBeUndefined();
  });

  it("goes back into the data on save only when the preference keeps it", () => {
    const d = doc("p", HISTORY);
    expect(dataToSave(seq("p", 10), d, true).history).toBe(HISTORY);
    expect(dataToSave(seq("p", 10), d, false).history).toBeUndefined();
    expect(dataToSave(seq("p", 10), doc("p", []), true).history).toBeUndefined();
  });
});
