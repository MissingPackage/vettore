import { describe, expect, it } from "vitest";
import { ALL_SHORTCUTS, SHORTCUT_GROUPS } from "./helpShortcuts";

describe("keyboard shortcut list", () => {
  it("binds each key combination to one action only", () => {
    const keys = ALL_SHORTCUTS.map((s) => s.keys);
    expect(new Set(keys).size).toBe(keys.length);
  });

  it("describes every shortcut in words", () => {
    for (const shortcut of ALL_SHORTCUTS) {
      expect(shortcut.what.length).toBeGreaterThanOrEqual(3);
      expect(shortcut.keys).toMatch(/^(Ctrl|Alt)\+/);
    }
  });

  it("lists one view shortcut per view tab", () => {
    const views = SHORTCUT_GROUPS.find((g) => g.title === "Views");
    expect(views?.shortcuts.map((s) => s.keys)).toEqual([
      "Ctrl+1", "Ctrl+2", "Ctrl+3", "Ctrl+4", "Ctrl+5", "Ctrl+6",
    ]);
  });

  it("lists the shortcuts Vettore binds itself", () => {
    const keys = ALL_SHORTCUTS.map((s) => s.keys);
    expect(keys).toEqual(expect.arrayContaining(["Ctrl+Shift+G", "Ctrl+Shift+V"]));
  });
});
