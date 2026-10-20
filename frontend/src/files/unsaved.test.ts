import { describe, expect, it } from "vitest";
import { isDirty, needsSave, uniqueName } from "./unsaved";

describe("unsaved documents", () => {
  const untouched = { sequenceData: { stateTrackingId: "initialLoadId" }, lastSavedId: undefined };
  const edited = { sequenceData: { stateTrackingId: "abc" }, lastSavedId: "initialLoadId" };
  const savedAfterEdit = { sequenceData: { stateTrackingId: "abc" }, lastSavedId: "abc" };

  it("follows the editor's edits since opening or saving", () => {
    expect(isDirty(untouched)).toBe(false);
    expect(isDirty(edited)).toBe(true);
    expect(isDirty(savedAfterEdit)).toBe(false);
    expect(isDirty(undefined)).toBe(false);
  });

  it("counts a document never saved (product, paste, new file) even when untouched", () => {
    expect(needsSave({ neverSaved: true }, untouched)).toBe(true);
    expect(needsSave({ neverSaved: true }, savedAfterEdit)).toBe(true);
  });

  it("does not count an untouched opened file or sample, but counts it once edited", () => {
    expect(needsSave({ neverSaved: false }, untouched)).toBe(false);
    expect(needsSave({ neverSaved: false }, edited)).toBe(true);
  });
});

describe("names of new documents", () => {
  it("keeps a name nobody uses", () => {
    expect(uniqueName("Pasted sequence", ["pUC19"])).toBe("Pasted sequence");
  });

  it("numbers a name already open, from 2, skipping numbers in use", () => {
    expect(uniqueName("Pasted sequence", ["Pasted sequence"])).toBe("Pasted sequence 2");
    expect(uniqueName("Pasted sequence", ["Pasted sequence", "Pasted sequence 2"])).toBe("Pasted sequence 3");
    expect(uniqueName("pUC19", ["pUC19", "pUC19 3"])).toBe("pUC19 2");
  });
});
