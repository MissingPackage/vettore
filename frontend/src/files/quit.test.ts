import { describe, expect, it } from "vitest";
import { notSavedMessage, quitDetail, quitHeadline } from "./quit";

const doc = (name: string, path: string | null, saveInPlace = false) =>
  ({ editorName: name, name, path, saveInPlace });

describe("closing the window with unsaved work", () => {
  it("names a single document, and counts several", () => {
    expect(quitHeadline([doc("pUC19", null)])).toBe("“pUC19” has changes that are not saved.");
    expect(quitHeadline([doc("a", null), doc("b", null)])).toBe("2 documents have changes that are not saved:");
  });

  it("says where Save puts each document, as File → Save would", () => {
    expect(quitDetail(doc("pTest", null))).toBe("never saved to a file: you choose where");
    expect(quitDetail(doc("colleague", "/data/colleague.gb", true))).toBe("saved to colleague.gb");
    expect(quitDetail(doc("pFA", "C:\\Users\\me\\pFA-KanMX4.dna"))).toBe("saved as a GenBank copy of pFA-KanMX4.dna");
  });

  it("explains why the window stayed open and what to do", () => {
    const text = notSavedMessage("pTest");
    expect(text).toContain("“pTest” was not saved");
    expect(text).toContain("File → Save");
    expect(text).toContain("Don't Save");
  });
});
