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

describe("a History from a file is outside content", () => {
  it("does not throw on more inputs than a function call takes arguments", () => {
    // 150,000 parents in one step used to take the whole window down: Math.max(...children).
    const parents = Array.from({ length: 150_000 }, (_, i) => ({ name: `a${i}` }));
    const history = cleanHistory([{ time: "t", action: "Gibson assembly", parents }]);
    expect(history[0].parents?.length).toBe(100);
    const tree = buildLineage(doc("p", history), seq("p", 10));
    expect(() => lineageDepth(tree)).not.toThrow();
    expect(lineageDepth(tree)).toBe(2);
  });

  it("cuts the number of steps and the length of what a step says", () => {
    const many = Array.from({ length: 900 }, (_, i) => ({ time: "t", action: `step ${i}` }));
    expect(cleanHistory(many)).toHaveLength(500);
    const [step] = cleanHistory([{ time: "t", action: "A".repeat(5_000_000), detail: "B".repeat(5_000_000) }]);
    expect(step.action.length).toBe(4000);
    expect(step.detail?.length).toBe(4000);
  });

  it("drops a snapshot that is not shaped like a document the editor can open", () => {
    const parents = [
      { name: "good", data: seq("good", 10) },
      { name: "no features", data: { name: "x", sequence: "ACGT" } },
      { name: "features are not annotations", data: { name: "x", sequence: "ACGT", features: "nope", primers: [] } },
      { name: "no sequence", data: { name: "x", features: [], primers: [] } },
    ];
    const [step] = cleanHistory([{ time: "t", action: "Cloning", parents }]);
    expect(step.parents?.map((p) => !!p.data)).toEqual([true, false, false, false]);
    // A missing copy is marked, so the view says so instead of offering a dead button.
    expect(step.parents?.slice(1).every((p) => p.snapshots_dropped)).toBe(true);
  });

  it("keeps the wrong clones of a step only when every one of them is a document", () => {
    const real = [{ name: "empty vector", data: seq("empty vector", 4000, true) }];
    expect(cleanHistory([{ time: "t", action: "Cloning", wrongClones: real }])[0].wrongClones)
      .toHaveLength(1);
    for (const junk of [[{ evil: "yes" }], [{ name: "x" }], [{ name: "x", data: { sequence: 1 } }]]) {
      expect(cleanHistory([{ time: "t", action: "Cloning", wrongClones: junk }])[0].wrongClones)
        .toBeUndefined();
    }
  });
});

describe("the step a node came out of", () => {
  it("is the cloning, not the save that followed it", () => {
    // Saving, exporting and printing are history entries too; none of them made the document.
    const history: HistoryEntry[] = [
      { time: "t1", action: "Opened", detail: "backbone.gb" },
      { time: "t2", action: "Restriction cloning", detail: "EcoRI + BamHI",
        parents: [{ name: "vector", data: seq("vector", 4000, true), history: [] }] },
      { time: "t3", action: "Saved", detail: "/lab/pProduct.gb (GenBank)" },
      { time: "t4", action: "Exported", detail: "/lab/pProduct.fasta (FASTA)" },
    ];
    const tree = buildLineage(doc("pProduct", history), seq("pProduct", 5200, true));
    expect(tree.operation).toBe("Restriction cloning");
    expect(tree.time).toBe("t2");
  });

  it("is how it came into being when no step was run on inputs", () => {
    const history: HistoryEntry[] = [
      { time: "t1", action: "Opened", detail: "pBR322.gb" },
      { time: "t2", action: "Saved", detail: "/lab/pBR322.gb" },
    ];
    expect(buildLineage(doc("pBR322", history), seq("x", 4361, true)).operation).toBe("Opened");
  });

  it("is the cloning for an input that was saved after it was made", () => {
    const vectorHistory: HistoryEntry[] = [
      { time: "t0", action: "Gibson assembly", detail: "2 fragments",
        parents: [{ name: "fragment", data: seq("fragment", 500), history: [] }] },
      { time: "t1", action: "Saved", detail: "/lab/vector.gb (GenBank)" },
    ];
    const history: HistoryEntry[] = [{
      time: "t2", action: "Restriction cloning",
      parents: [{ name: "vector", path: "/lab/vector.gb", data: seq("vector", 4000, true),
                  history: vectorHistory }],
    }];
    const [vector] = buildLineage(doc("p", history), seq("p", 4500, true)).children;
    expect(vector.operation).toBe("Gibson assembly");
    expect(vector.usedBy).toBe("Restriction cloning");
  });
});
