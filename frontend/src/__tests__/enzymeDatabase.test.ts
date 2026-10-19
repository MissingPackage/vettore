import { describe, expect, it } from "vitest";
import { cutsites } from "../enzymeLogic";
import cases from "./fixtures/cuts.json";

// Cross-engine check: OVE's cut-site engine, fed with the REBASE-derived database, must cut
// exactly where Biopython cuts (fixture written by scripts/make_enzyme_fixture.py).
describe("bundled enzyme database", () => {
  for (const c of cases as unknown as { name: string; sequence: string; circular: boolean; cuts: Record<string, number[]> }[]) {
    it(`matches Biopython cut positions on ${c.name}`, () => {
      const found = cutsites(c.sequence, c.circular);
      const ove: Record<string, number[]> = {};
      for (const [name, sites] of Object.entries(found)) {
        ove[name] = sites.map((s) => s.topSnipPosition % c.sequence.length).sort((a, b) => a - b);
      }
      expect(ove).toEqual(c.cuts);
    });
  }
});
