import { describe, expect, it } from "vitest";
import { matchScore, rankMenuItems } from "./menuSearch";

const item = (text: string, path: string[] = []) => ({ text, path, onClick: () => {} });

describe("menu search ranking", () => {
  it("puts a name containing the word before letters scattered in order", () => {
    const items = [item("Change Topology", ["Actions"]), item("Simulate Agarose Gel…", ["Tools"])];
    const names = rankMenuItems("gel", items).map((i) => i.text);
    expect(names).toEqual(["Simulate Agarose Gel…"]);
  });

  it("ranks a name start, then a word start, then anywhere, then the menu path", () => {
    expect(matchScore("gib", "Gibson Assembly…")).toBe(0);
    expect(matchScore("ass", "Gibson Assembly…")).toBe(1);
    expect(matchScore("sembl", "Gibson Assembly…")).toBe(2);
    expect(matchScore("tools", "Simulate Agarose Gel…", ["Tools"])).toBe(3);
    expect(matchScore("xyz", "Gibson Assembly…")).toBeNull();
  });

  it("leaves out items that cannot be run or are hidden, and nothing for an empty query", () => {
    const items = [{ text: "Gel", path: [] }, { ...item("Gel"), hidden: true }, item("Gel")];
    expect(rankMenuItems("gel", items)).toHaveLength(1);
    expect(rankMenuItems("  ", [item("Gel")])).toHaveLength(0);
  });
});
