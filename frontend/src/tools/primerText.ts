// Form logic for Primers → Design PCR Primers, kept free of React so it can be unit-tested.
// The form shows positions 1-based, as biologists count; the engine takes 0-based inclusive.

export interface PrimerForm {
  start: string;          // 1-based, inclusive
  end: string;
  productMin: string;
  productMax: string;
  sizeMin: string;
  sizeOpt: string;
  sizeMax: string;
  tmMin: string;
  tmOpt: string;
  tmMax: string;
  gcMin: string;
  gcMax: string;
  maxTmDifference: string;
}

/** Standard PCR (sources in src/vettore/primer_design.py): 18–24 nt, Tm 57–63 °C, GC 40–60 %. */
export const DEFAULT_FORM: PrimerForm = {
  start: "",
  end: "",
  productMin: "",
  productMax: "",
  sizeMin: "18",
  sizeOpt: "20",
  sizeMax: "24",
  tmMin: "57",
  tmOpt: "60",
  tmMax: "63",
  gcMin: "40",
  gcMax: "60",
  maxTmDifference: "5",
};

/** How far beyond the smallest possible product the default range reaches. */
export const PRODUCT_SLACK = 300;
export const MAX_PRIMER_LENGTH = 36;

export type FormErrors = Partial<Record<"region" | "product" | "size" | "tm" | "gc" | "tmDiff", string>>;

const int = (text: string) => (/^\s*\d+\s*$/.test(text) ? Number(text) : NaN);
const num = (text: string) => (text.trim() === "" ? NaN : Number(text));
const n = (value: number) => value.toLocaleString("en-US");

/** Region length in bp, across the origin of a circular sequence; null if not set. */
export function regionLength(form: PrimerForm, length: number, circular: boolean): number | null {
  const start = int(form.start);
  const end = int(form.end);
  if (!(start >= 1 && end >= 1 && start <= length && end <= length)) return null;
  if (start <= end) return end - start + 1;
  return circular ? length - start + end + 1 : null;
}

/** The product range the engine would use when the fields are left as they are. */
export function defaultProductRange(regionBp: number, sizeMin: number): [number, number] {
  const smallest = regionBp + 2 * sizeMin;
  return [smallest, smallest + PRODUCT_SLACK];
}

/** Inline messages next to the fields; empty when the form can be sent. */
export function validate(form: PrimerForm, length: number, circular: boolean): FormErrors {
  const errors: FormErrors = {};
  const start = int(form.start);
  const end = int(form.end);
  if (form.start.trim() === "" || form.end.trim() === "") {
    errors.region = "Pick a feature or type the region's start and end. Tip: select the region in the map before opening this dialog and it is filled in for you.";
  } else if (!(start >= 1 && start <= length) || !(end >= 1 && end <= length)) {
    errors.region = `Start and end must be whole numbers between 1 and ${n(length)}.`;
  } else if (start > end && !circular) {
    errors.region = "The start must come before the end on a linear sequence.";
  }

  const [sMin, sOpt, sMax] = [int(form.sizeMin), int(form.sizeOpt), int(form.sizeMax)];
  if ([sMin, sOpt, sMax].some(Number.isNaN)) {
    errors.size = "Enter whole numbers.";
  } else if (!(sMin <= sOpt && sOpt <= sMax)) {
    errors.size = "Minimum ≤ optimum ≤ maximum.";
  } else if (sMin < 8 || sMax > MAX_PRIMER_LENGTH) {
    errors.size = `Primer3 designs primers of 8 to ${MAX_PRIMER_LENGTH} nt.`;
  }

  const [tMin, tOpt, tMax] = [num(form.tmMin), num(form.tmOpt), num(form.tmMax)];
  if ([tMin, tOpt, tMax].some(Number.isNaN)) errors.tm = "Enter temperatures in °C.";
  else if (!(tMin <= tOpt && tOpt <= tMax)) errors.tm = "Minimum ≤ optimum ≤ maximum.";

  const [gMin, gMax] = [num(form.gcMin), num(form.gcMax)];
  if ([gMin, gMax].some(Number.isNaN) || gMin < 0 || gMax > 100) errors.gc = "Enter percentages from 0 to 100.";
  else if (gMin > gMax) errors.gc = "The minimum must not exceed the maximum.";

  const diff = num(form.maxTmDifference);
  if (Number.isNaN(diff) || diff < 0) errors.tmDiff = "Enter 0 °C or more.";

  const region = regionLength(form, length, circular);
  if (!errors.region && region !== null && !errors.size) {
    const [defaultMin, defaultMax] = defaultProductRange(region, sMin);
    const pMin = form.productMin.trim() === "" ? defaultMin : int(form.productMin);
    const pMax = form.productMax.trim() === "" ? defaultMax : int(form.productMax);
    if (Number.isNaN(pMin) || Number.isNaN(pMax)) errors.product = "Enter whole numbers of base pairs.";
    else if (pMax < defaultMin) {
      errors.product = `At least ${n(defaultMin)} bp: the ${n(region)} bp region plus two primers.`;
    } else if (pMin > pMax) errors.product = "The minimum must not exceed the maximum.";
    else if (pMin > length) errors.product = `The sequence is only ${n(length)} bp long.`;
  }
  return errors;
}

/** Request body for /api/primer-design/design, or null while the form has errors. */
export function buildRequest(
  form: PrimerForm,
  template: { name: string; sequence: string; circular: boolean },
): object | null {
  if (Object.keys(validate(form, template.sequence.length, template.circular)).length) return null;
  const optional = (text: string) => (text.trim() === "" ? null : int(text));
  return {
    template,
    start: int(form.start) - 1,
    end: int(form.end) - 1,
    product_min: optional(form.productMin),
    product_max: optional(form.productMax),
    size_min: int(form.sizeMin),
    size_opt: int(form.sizeOpt),
    size_max: int(form.sizeMax),
    tm_min: num(form.tmMin),
    tm_opt: num(form.tmOpt),
    tm_max: num(form.tmMax),
    gc_min: num(form.gcMin),
    gc_max: num(form.gcMax),
    max_tm_difference: num(form.maxTmDifference),
  };
}

/** Primer names safe for ordering forms (no spaces), unique among the document's primers. */
export function primerNames(prefix: string, existing: string[]): [string, string] {
  const base = (prefix.trim() || "primer").replace(/\s+/g, "_");
  const taken = new Set(existing);
  for (let i = 1; ; i++) {
    const suffix = i === 1 ? "" : `_${i}`;
    const pair: [string, string] = [`${base}_F${suffix}`, `${base}_R${suffix}`];
    if (!taken.has(pair[0]) && !taken.has(pair[1])) return pair;
  }
}

/** "name ⇥ sequence" lines, the format order forms accept. */
export function orderingText(primers: { name: string; sequence: string }[]): string {
  return primers.map((p) => `${p.name}\t${p.sequence}`).join("\n");
}

/** A 1-based position range as written in the dialog. */
export function range(start0: number, end0: number): string {
  return `${n(start0 + 1)}–${n(end0 + 1)}`;
}
