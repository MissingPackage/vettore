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

const WHOLE_REGION = "To amplify the whole region, primers included, use Actions → Design Primers to "
  + "Amplify a Region, which places the primers at the ends of the region itself.";

/**
 * Why the region leaves no room for a primer on each side of it, with what to do; null when
 * there is room or the region is not set. The same check, in the same words, as the engine's
 * (primer_design._no_room): the dialog shows it in the result pane before asking the engine.
 */
export function roomProblem(form: PrimerForm, length: number, circular: boolean, name: string): string | null {
  const region = regionLength(form, length, circular);
  const sizeMin = int(form.sizeMin);
  if (region === null || Number.isNaN(sizeMin)) return null;
  const start = int(form.start) - 1;
  const end = int(form.end) - 1;
  if (circular) {
    if (region + 2 * sizeMin <= length) return null;
    return `The region covers ${n(region)} of the ${n(length)} bp of “${name}”, which leaves no room `
      + `for primers: each needs at least ${sizeMin} bases outside the region. Choose a region at `
      + `least ${2 * sizeMin} bp shorter than the plasmid. ${WHOLE_REGION}`;
  }
  if (start >= sizeMin && end + sizeMin <= length - 1) return null;
  return `Primers need at least ${sizeMin} bases on each side of the region, outside it, but this `
    + `region reaches ${start < sizeMin ? "the start" : "the end"} of the linear sequence “${name}”. `
    + `Choose a region from ${n(sizeMin + 1)} to ${n(length - sizeMin)}. ${WHOLE_REGION}`;
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
  // A region with no room for primers is explained in the result pane (roomProblem), not here.
  const noRoom = roomProblem(form, length, circular, "") !== null;
  if (!errors.region && region !== null && !errors.size && !noRoom) {
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
  if (roomProblem(form, template.sequence.length, template.circular, template.name)) return null;
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

/** The names under which both primers of a pair are already in the document, or null. */
export function pairAlreadyIn(
  forward: string,
  reverse: string,
  primers: { name: string; bases?: string }[],
): [string, string] | null {
  const find = (sequence: string) =>
    primers.find((p) => p.bases?.toUpperCase() === sequence.toUpperCase())?.name;
  const [f, r] = [find(forward), find(reverse)];
  return f && r ? [f, r] : null;
}

/** "name ⇥ sequence" lines, the format order forms accept. */
export function orderingText(primers: { name: string; sequence: string }[]): string {
  return primers.map((p) => `${p.name}\t${p.sequence}`).join("\n");
}

/** A 1-based position range as written in the dialog. */
export function range(start0: number, end0: number): string {
  return `${n(start0 + 1)}–${n(end0 + 1)}`;
}
