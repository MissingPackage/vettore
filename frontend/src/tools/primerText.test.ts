import { describe, expect, it } from "vitest";
import {
  buildRequest,
  DEFAULT_FORM,
  defaultProductRange,
  orderingText,
  primerNames,
  regionLength,
  validate,
} from "./primerText";

const form = (patch: Partial<typeof DEFAULT_FORM>) => ({ ...DEFAULT_FORM, ...patch });

describe("region", () => {
  it("counts the region across the origin of a circular sequence only", () => {
    expect(regionLength(form({ start: "101", end: "300" }), 1000, false)).toBe(200);
    expect(regionLength(form({ start: "901", end: "100" }), 1000, true)).toBe(200);
    expect(regionLength(form({ start: "901", end: "100" }), 1000, false)).toBeNull();
  });

  it("asks for a region before anything else", () => {
    expect(validate(DEFAULT_FORM, 1000, false).region).toMatch(/Pick a feature/);
    expect(validate(form({ start: "0", end: "10" }), 1000, false).region).toMatch(/between 1 and 1,000/);
    expect(validate(form({ start: "500", end: "100" }), 1000, false).region).toMatch(/linear/);
    expect(validate(form({ start: "500", end: "100" }), 1000, true).region).toBeUndefined();
  });
});

describe("constraints", () => {
  const ok = form({ start: "301", end: "500" });

  it("accepts the defaults for a standard PCR", () => {
    expect(validate(ok, 1500, false)).toEqual({});
  });

  it("checks every group next to its field", () => {
    expect(validate({ ...ok, sizeMin: "22", sizeOpt: "20" }, 1500, false).size).toMatch(/optimum/);
    expect(validate({ ...ok, sizeMax: "40" }, 1500, false).size).toMatch(/8 to 36/);
    expect(validate({ ...ok, tmMin: "65" }, 1500, false).tm).toMatch(/optimum/);
    expect(validate({ ...ok, gcMin: "70" }, 1500, false).gc).toMatch(/minimum/);
    expect(validate({ ...ok, maxTmDifference: "-1" }, 1500, false).tmDiff).toBeDefined();
  });

  it("needs room in the product for the region and both primers", () => {
    expect(defaultProductRange(200, 18)).toEqual([236, 536]);
    expect(validate({ ...ok, productMax: "230" }, 1500, false).product)
      .toBe("At least 236 bp: the 200 bp region plus two primers.");
    expect(validate({ ...ok, productMin: "400", productMax: "300" }, 1500, false).product).toMatch(/exceed/);
  });
});

it("sends 0-based positions and leaves an empty product range to the engine", () => {
  const body = buildRequest(form({ start: "301", end: "500" }),
                            { name: "p", sequence: "A".repeat(1500), circular: false }) as any;
  expect(body.start).toBe(300);
  expect(body.end).toBe(499);
  expect(body.product_min).toBeNull();
  expect(body.size_opt).toBe(20);
  expect(buildRequest(DEFAULT_FORM, { name: "p", sequence: "A".repeat(1500), circular: false })).toBeNull();
});

it("names primers for ordering without clashing with existing ones", () => {
  expect(primerNames("EGFP", [])).toEqual(["EGFP_F", "EGFP_R"]);
  expect(primerNames("my gene", ["my_gene_F"])).toEqual(["my_gene_F_2", "my_gene_R_2"]);
  expect(primerNames(" ", [])).toEqual(["primer_F", "primer_R"]);
});

it("writes name ⇥ sequence lines for ordering", () => {
  expect(orderingText([{ name: "a_F", sequence: "ACGT" }, { name: "a_R", sequence: "TTGG" }]))
    .toBe("a_F\tACGT\na_R\tTTGG");
});
