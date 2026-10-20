// AAV packaging check: find the ITRs of a transfer plasmid and measure what a capsid would
// have to package. Pure functions, so the arithmetic is covered by unit tests.
//
// Sizes, with their sources:
// - Wild-type AAV2 is 4679 nt including both ITRs (RefSeq NC_001401.2; Srivastava A, Lusby EW,
//   Berns KI, "Nucleotide sequence and organization of the adeno-associated virus 2 genome",
//   J Virol 1983;45:555-64).
// - Single-stranded vectors package up to about the wild-type size; above it the titre falls
//   and the genomes that are packaged are truncated at the 5' end (Dong JY, Fan PD, Frizzell
//   RA, "Quantitative analysis of the packaging capacity of recombinant adeno-associated
//   virus", Hum Gene Ther 1996;7:2101-12; Wu Z, Yang H, Colosi P, "Effect of genome size on
//   AAV vector packaging", Mol Ther 2010;18:80-6).
// - A self-complementary vector packages an inverted dimer of the cassette, so it holds about
//   half as much (McCarty DM et al., "Self-complementary recombinant adeno-associated virus
//   (scAAV) vectors promote efficient transduction independently of DNA synthesis",
//   Gene Ther 2001;8:1248-54).

import type { OveAnnotation, SequenceData } from "../api";

/** Genome size, ITRs included, that a single-stranded vector packages reliably. */
export const SS_LIMIT_BP = 4700;
/** Above this a single-stranded genome is packaged truncated, if at all. */
export const SS_HARD_LIMIT_BP = 5000;
/** The self-complementary vector carries the cassette twice, so the room is halved. */
export const SC_LIMIT_BP = Math.round(SS_LIMIT_BP / 2);
export const SC_HARD_LIMIT_BP = Math.round(SS_HARD_LIMIT_BP / 2);
/** Wild-type AAV2, for reference on the bar. */
export const WILD_TYPE_BP = 4679;

export type VectorKind = "ss" | "sc";

export function limitFor(kind: VectorKind): { limit: number; hardLimit: number } {
  return kind === "sc"
    ? { limit: SC_LIMIT_BP, hardLimit: SC_HARD_LIMIT_BP }
    : { limit: SS_LIMIT_BP, hardLimit: SS_HARD_LIMIT_BP };
}

// "ITR" as a word: a protein translation can contain those letters inside a peptide
// (SaCas9 has ...LVITRDEN...), and a feature called "mITRa" is not an ITR either.
const ITR_WORD = /\bITRs?\b|inverted terminal repeat/i;
// Features that mark the plasmid backbone rather than the cargo.
const BACKBONE_WORD =
  /\b(ori|origin of replication|f1 ori|pUC origin|ColE1|AmpR|KanR|CmR|SpecR|NeoR|bla|ampicillin|kanamycin|chloramphenicol|spectinomycin|resistance)\b/i;

export interface Annotation {
  name: string;
  type: string;
  start: number;          // 0-based inclusive
  end: number;            // 0-based inclusive; end < start means it crosses the origin
  notes: Record<string, string[]>;
}

export interface Itr extends Annotation {
  length: number;
  /** Where it came from: the document's own annotation, or the feature catalogue. */
  source: "annotation" | "catalogue";
}

/** OVE keeps annotations either as a list or as a map keyed by id. */
export function annotationList(value: SequenceData["features"] | undefined): OveAnnotation[] {
  if (!value) return [];
  return Array.isArray(value) ? value : Object.values(value);
}

function noteText(notes: Record<string, string[]> | undefined): string {
  if (!notes) return "";
  return Object.entries(notes)
    // A translation is protein, not annotation text: searching it gives false ITRs.
    .filter(([key]) => key !== "translation")
    .map(([key, value]) => `${key} ${Array.isArray(value) ? value.join(" ") : String(value)}`)
    .join(" ");
}

function normalise(annotation: OveAnnotation): Annotation {
  return {
    name: annotation.name ?? "",
    type: annotation.type ?? "",
    start: annotation.start,
    end: annotation.end,
    notes: (annotation.notes as Record<string, string[]>) ?? {},
  };
}

/** Length of start..end inclusive, wrapping through the origin when start > end. */
export function spanLength(start: number, end: number, sequenceLength: number): number {
  return start <= end ? end - start + 1 : sequenceLength - start + end + 1;
}

export function looksLikeItr(annotation: Annotation): boolean {
  const notes = noteText(annotation.notes);
  if (ITR_WORD.test(annotation.name) || ITR_WORD.test(notes)) return true;
  // A repeat_region flagged as an inverted repeat is an ITR in an AAV plasmid.
  return annotation.type === "repeat_region" && /\binverted\b/i.test(notes);
}

/** ITRs annotated in the document, in order of position. */
export function findItrs(
  features: SequenceData["features"] | undefined,
  sequenceLength: number,
  source: Itr["source"] = "annotation",
): Itr[] {
  return annotationList(features)
    .map(normalise)
    .filter(looksLikeItr)
    .map((annotation) => ({
      ...annotation,
      length: spanLength(annotation.start, annotation.end, sequenceLength),
      source,
    }))
    .sort((a, b) => a.start - b.start || a.end - b.end);
}

export interface Cargo {
  /** The ITR the genome starts at, and the one it ends at, in reading order. */
  from: Itr;
  to: Itr;
  start: number;           // 0-based inclusive, start of `from`
  end: number;             // 0-based inclusive, end of `to`
  /** Everything a capsid would carry: ITR to ITR, both included. */
  totalBp: number;
  /** What sits between the two ITRs: promoter, transgene, polyA … */
  insertBp: number;
  crossesOrigin: boolean;
  /** Backbone features inside this span (origin, resistance marker): it is probably not the cargo. */
  backbone: string[];
}

function backboneIn(
  features: SequenceData["features"] | undefined,
  start: number,
  end: number,
  sequenceLength: number,
): string[] {
  const inside = (position: number) =>
    start <= end
      ? position >= start && position <= end
      : position >= start || position <= end;
  return annotationList(features)
    .map(normalise)
    .filter((annotation) => {
      const text = `${annotation.name} ${annotation.type} ${noteText(annotation.notes)}`;
      const marker = annotation.type === "rep_origin" || BACKBONE_WORD.test(text);
      return marker && inside(annotation.start);
    })
    .map((annotation) => annotation.name || annotation.type);
}

function cargoFor(
  from: Itr,
  to: Itr,
  features: SequenceData["features"] | undefined,
  sequenceLength: number,
): Cargo {
  const totalBp = spanLength(from.start, to.end, sequenceLength);
  return {
    from,
    to,
    start: from.start,
    end: to.end,
    totalBp,
    insertBp: Math.max(0, totalBp - from.length - to.length),
    crossesOrigin: from.start > to.end,
    backbone: backboneIn(features, from.start, to.end, sequenceLength),
  };
}

/**
 * The two ways to read a pair of ITRs on a circular plasmid: through the cassette, or the
 * other way round through the backbone. Linear sequences have only the first.
 */
export function cargoOptions(
  pair: [Itr, Itr],
  data: Pick<SequenceData, "circular" | "features">,
  sequenceLength: number,
): Cargo[] {
  const [first, second] = pair;
  const options = [cargoFor(first, second, data.features, sequenceLength)];
  if (data.circular) options.push(cargoFor(second, first, data.features, sequenceLength));
  return options;
}

/**
 * Which side of the ITRs is the cargo. The side carrying the origin of replication and the
 * resistance marker is the plasmid backbone, so the other one is what gets packaged; without
 * such markers the span that does not cross the origin is used.
 */
export function chooseCargo(options: Cargo[]): { cargo: Cargo; reason: string } {
  if (options.length === 1) return { cargo: options[0], reason: "" };
  const [through, around] = options;
  if (through.backbone.length !== around.backbone.length) {
    const cargo = through.backbone.length < around.backbone.length ? through : around;
    const other = cargo === through ? around : through;
    return {
      cargo,
      reason: `The other side carries ${other.backbone.slice(0, 3).join(", ")}, so it is the `
        + "plasmid backbone.",
    };
  }
  const cargo = through.crossesOrigin && !around.crossesOrigin ? around : through;
  return {
    cargo,
    reason: "No origin or resistance marker is annotated, so the span that does not cross "
      + "position 1 was taken. Check the choice below.",
  };
}

export type VerdictLevel = "ok" | "tight" | "over";

export interface Verdict {
  level: VerdictLevel;
  headline: string;
  detail: string;
  /** How full the capsid is, 1 = at the limit (can be above 1). */
  load: number;
  limit: number;
}

export function verdictFor(totalBp: number, kind: VectorKind): Verdict {
  const { limit, hardLimit } = limitFor(kind);
  const load = totalBp / limit;
  const name = kind === "sc" ? "self-complementary AAV" : "single-stranded AAV";
  const spare = limit - totalBp;
  if (totalBp <= limit) {
    return {
      level: "ok",
      headline: `Fits a ${name} vector`,
      detail: `${spare.toLocaleString()} bp to spare before the ${limit.toLocaleString()} bp `
        + "limit, ITRs included.",
      load,
      limit,
    };
  }
  if (totalBp <= hardLimit) {
    return {
      level: "tight",
      headline: `Over the limit for a ${name} vector`,
      detail: `${(totalBp - limit).toLocaleString()} bp too long. Genomes of this size are `
        + "still packaged, but the titre drops and part of the preparation carries truncated "
        + "genomes. Shorten the cassette or use a smaller polyA, promoter or linker.",
      load,
      limit,
    };
  }
  return {
    level: "over",
    headline: `Too long for a ${name} vector`,
    detail: `${(totalBp - limit).toLocaleString()} bp over the `
      + `${limit.toLocaleString()} bp limit and beyond `
      + `${hardLimit.toLocaleString()} bp, where genomes are packaged truncated if at all. `
      + "Split the cassette over two vectors or remove elements.",
    load,
    limit,
  };
}

/** Position 1-based, as a biologist counts. */
export function asPosition(index: number): number {
  return index + 1;
}
