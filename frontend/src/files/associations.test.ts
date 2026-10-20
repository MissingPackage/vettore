import { describe, expect, it } from "vitest";
import type { FileAssociations } from "../api";
import { associationHint, extensionList, movedNote, overriddenNote } from "./associations";

const status = (patch: Partial<FileAssociations> = {}): FileAssociations => ({
  available: true,
  enabled: true,
  current: true,
  registered_executable: "C:\\Vettore\\Vettore.exe",
  extensions: [".dna", ".gb", ".ab1"],
  overridden: [],
  platform: "win32",
  ...patch,
});

describe("opening sequence files with Vettore", () => {
  it("lists the extensions as a sentence", () => {
    expect(extensionList([".dna"])).toBe(".dna");
    expect(extensionList([".dna", ".gb", ".ab1"])).toBe(".dna, .gb and .ab1");
  });

  it("says in one sentence what the option does and how to undo it", () => {
    expect(associationHint(status())).toBe(
      "Applies to .dna, .gb and .ab1 files, for your Windows account only; untick it to undo, "
      + "and they open as they did before.");
    expect(associationHint(status({ platform: "linux" }))).toContain("for your user account only");
  });

  it("explains when the files are registered to a Vettore in another folder", () => {
    expect(movedNote(status())).toBeNull();
    expect(movedNote(status({ registered_executable: null, current: false, enabled: false }))).toBeNull();
    expect(movedNote(status({ current: false, enabled: false, registered_executable: "D:\\old\\Vettore.exe" })))
      .toBe("Sequence files open with the Vettore at D:\\old\\Vettore.exe, not with this copy. "
        + "Tick the box to use this copy instead.");
  });

  it("tells how to change an extension Windows keeps for another program", () => {
    expect(overriddenNote(status())).toBeNull();
    expect(overriddenNote(status({ overridden: [".dna"] })))
      .toContain("always open .dna files with another program, so it still opens there.");
    expect(overriddenNote(status({ overridden: [".dna", ".gb"] })))
      .toContain(".dna and .gb files with another program, so they still open there.");
    expect(overriddenNote(status({ overridden: [".dna"] }))).toContain("Open with → Choose another app");
    expect(overriddenNote(status({ overridden: [".dna"], enabled: false }))).toBeNull();
  });
});
