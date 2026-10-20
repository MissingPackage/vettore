import { describe, expect, it } from "vitest";
import { defaultProductName } from "./productName";

describe("default product names", () => {
  it("add the method to a name made of the inputs", () => {
    expect(defaultProductName("pUC19 + GFP", "Gibson assembly", [])).toBe("pUC19 + GFP (Gibson)");
    expect(defaultProductName("pUC19 + GFP", "Restriction cloning", [])).toBe("pUC19 + GFP (restriction cloning)");
    expect(defaultProductName("pDEST + pENTR", "Gateway LR reaction", [])).toBe("pDEST + pENTR (LR)");
    // An operation without a short name is used as it is.
    expect(defaultProductName("pCR + GFP", "TOPO TA cloning", [])).toBe("pCR + GFP (TOPO TA cloning)");
  });

  it("tell apart two products of the same files made by different methods", () => {
    const gibson = defaultProductName("pUC19 + KanMX", "Gibson assembly", []);
    const restriction = defaultProductName("pUC19 + KanMX", "Restriction cloning", [gibson]);
    expect(restriction).not.toBe(gibson);
  });

  it("number a name that is already open", () => {
    expect(defaultProductName("pUC19 + GFP", "Gibson assembly", ["pUC19 + GFP (Gibson)"]))
      .toBe("pUC19 + GFP (Gibson) 2");
    expect(defaultProductName("pUC19 + GFP", "Gibson assembly", ["pUC19 + GFP (Gibson)", "pUC19 + GFP (Gibson) 2"]))
      .toBe("pUC19 + GFP (Gibson) 3");
  });

  it("keep names that already say what was done", () => {
    expect(defaultProductName("pUC19 101..900 PCR product", "PCR", [])).toBe("pUC19 101..900 PCR product");
    expect(defaultProductName("pUC19 mut", "Mutagenesis", [])).toBe("pUC19 mut");
    expect(defaultProductName("pUC19 mut", "Mutagenesis", ["pUC19 mut"])).toBe("pUC19 mut 2");
  });
});
