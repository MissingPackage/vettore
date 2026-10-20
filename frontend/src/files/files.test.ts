import { describe, expect, it } from "vitest";
import { fastaRecord, rangeLength, selectionRange, subsequence } from "./fasta";
import { addNotice, dismissNotice, messageOf } from "./notices";
import { ACCEPT, dirName, safeFileName } from "./formats";
import { copyReason } from "./saveCopy";
import { EDITOR_TOOLS, rewireCopyItems, rightClickOverrides } from "./editorFiles";

// 70 bp: pBR322 (NCBI J01749.1) bases 1-70.
const PBR322 = "TTCTCATGTTTGACAGCTTATCATCGATAAGCTTTAATGCGGTAGTTTATCACAGTTAAATTGCTAACGC";

describe("Copy as FASTA", () => {
  it("copies the whole sequence in 60-base lines when nothing is selected", () => {
    expect(selectionRange(70, { start: -1, end: -1 })).toBeNull();
    expect(selectionRange(70, undefined)).toBeNull();
    expect(fastaRecord("pBR322 start", PBR322, null)).toBe(
      `>pBR322_start\n${PBR322.slice(0, 60)}\n${PBR322.slice(60)}\n`,
    );
  });

  it("copies a selection with its 1-based positions", () => {
    const range = selectionRange(70, { start: 9, end: 18 });
    expect(range).toEqual({ start: 9, end: 18 });
    expect(fastaRecord("pBR322", PBR322, range)).toBe(">pBR322:10-19\nTTGACAGCTT\n");
  });

  it("follows a selection across the origin of a circular sequence", () => {
    // Bases 66-70 (AACGC) then 1-5 (TTCTC).
    const range = { start: 65, end: 4 };
    expect(rangeLength(range, 70)).toBe(10);
    expect(subsequence(PBR322, range)).toBe("AACGCTTCTC");
    expect(fastaRecord("p", PBR322, range)).toBe(">p:66-5\nAACGCTTCTC\n");
  });

  it("ignores a selection outside the sequence", () => {
    expect(selectionRange(70, { start: 10, end: 90 })).toBeNull();
  });
});

describe("notices", () => {
  it("keeps several messages instead of overwriting the last one", () => {
    let list = addNotice([], "error", "A could not be opened.");
    list = addNotice(list, "error", "B could not be opened.");
    expect(list.map((n) => n.text)).toEqual(["A could not be opened.", "B could not be opened."]);
  });

  it("counts a repeated message instead of stacking it", () => {
    let list = addNotice([], "error", "Same.");
    list = addNotice(list, "info", "Other.");
    list = addNotice(list, "error", "Same.");
    expect(list.map((n) => [n.text, n.count])).toEqual([["Other.", 1], ["Same.", 2]]);
  });

  it("drops the oldest confirmations first when too many pile up", () => {
    let list = addNotice([], "info", "saved 0");
    for (let i = 1; i <= 6; i++) list = addNotice(list, "error", `error ${i}`);
    expect(list).toHaveLength(6);
    expect(list.some((n) => n.kind === "info")).toBe(false);
    expect(dismissNotice(list, list[0].id)).toHaveLength(5);
  });

  it("shows the message of an error, not its type", () => {
    expect(messageOf(new Error("“x.gb” is empty."))).toBe("“x.gb” is empty.");
    expect(messageOf("plain")).toBe("plain");
  });
});

describe("file names", () => {
  it("makes document names safe as file names", () => {
    expect(safeFileName('pAAV: CMV/GFP "v2"')).toBe("pAAV_ CMV_GFP _v2_");
    expect(safeFileName("  ..  ")).toBe("Untitled");
    expect(safeFileName("MT612431.1")).toBe("MT612431.1");
  });

  it("finds the folder of Windows and Linux paths", () => {
    expect(dirName("C:\\Users\\me\\plasmids\\p.dna")).toBe("C:\\Users\\me\\plasmids");
    expect(dirName("/home/me/p.gb")).toBe("/home/me");
    expect(dirName("p.gb")).toBe("");
  });

  it("offers every readable extension in the browser's file chooser", () => {
    for (const ext of [".dna", ".gb", ".embl", ".fasta", ".xdna", ".gck", ".ab1", ".txt", ".seq"]) {
      expect(ACCEPT.split(",")).toContain(ext);
    }
  });
});

describe("the editor's own file paths", () => {
  it("offers Copy as FASTA where the editor offered its own GenBank copy", () => {
    const copy = () => undefined;
    const items = [
      "newFeature",
      { text: "Copy", submenu: [{ text: "Copy" }, { text: "Copy Genbank For Selection", className: "openVeCopyGenbankForSelection" }] },
    ];
    const [, rewired] = rewireCopyItems(items, copy) as [string, { submenu: { text: string; onClick?: () => void }[] }];
    expect(rewired.submenu.map((i) => i.text)).toEqual(["Copy", "Copy as FASTA"]);
    expect(rewired.submenu[1].onClick).toBe(copy);
    expect(Object.keys(rightClickOverrides(copy))).toContain("featureRightClicked");
  });

  it("keeps the toolbar free of import, export and alignment-upload buttons", () => {
    for (const tool of ["downloadTool", "importTool", "alignmentTool"]) expect(EDITOR_TOOLS).not.toContain(tool);
    expect(EDITOR_TOOLS).toContain("saveTool");
  });
});

describe("why Save makes a GenBank copy", () => {
  const source = { label: "GenBank", records: 1, fileName: "p.gb", incomplete: false };
  it("explains each case in one sentence", () => {
    expect(copyReason("p", { ...source, label: ".dna", fileName: "p.dna" })).toContain("“p.dna” is a .dna file");
    expect(copyReason("J01749.1", { ...source, label: "FASTA", records: 3, fileName: "all.fasta" }))
      .toBe("“J01749.1” is one of 3 sequences in “all.fasta”; saving over that file would replace the other 2.");
    expect(copyReason("p", { ...source, incomplete: true })).toContain("could not be shown");
    expect(copyReason("p", { ...source, fileName: "p.txt" })).toBe("“p.txt” holds a GenBank record under another file extension.");
    expect(copyReason("p", { ...source, label: "EMBL", fileName: "p.embl" })).toBe("“p.embl” is an EMBL file, which Vettore reads but does not write.");
    expect(copyReason("p", { ...source, label: "Sanger trace (.ab1)", fileName: "r.ab1" })).toContain("is a Sanger trace (.ab1) file");
  });
});
