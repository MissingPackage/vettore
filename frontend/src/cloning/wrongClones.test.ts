import { describe, expect, it } from "vitest";
import type { CloningResult, Product } from "./types";
import { wrongClonesOf } from "./wrongClones";

const product = (sequence: string, label?: string): Product => ({
  doc: { name: `p-${sequence}`, sequence, circular: true, features: [], primers: [] },
  operation: "Restriction cloning",
  parents: [],
  details: "",
  oligos: [],
  label,
});

describe("wrongClonesOf", () => {
  it("lists the other products, then the engine's wrong clones, each molecule once", () => {
    const expected = product("AAAA", "forward");
    const result: CloningResult = {
      products: [expected, product("TTTT", "reverse")],
      warnings: [],
      oligos: [],
      screen_against: [product("TTTT", "insert in the other orientation"), product("GGGG", "parent vector")],
    };
    // The chosen product may carry the name the user typed: it is recognised by its sequence.
    const chosen = { ...expected, doc: { ...expected.doc, name: "my clone" } };
    expect(wrongClonesOf(chosen, result).map((w) => w.name)).toEqual(["reverse", "parent vector"]);
  });

  it("is empty for a simulation with one product and nothing to screen against", () => {
    const only = product("AAAA");
    expect(wrongClonesOf(only, { products: [only], warnings: [], oligos: [] })).toEqual([]);
  });
});
