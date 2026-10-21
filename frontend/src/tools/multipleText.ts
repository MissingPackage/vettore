// Wording and choices for "All to each other" in Tools → Align Multiple Sequences (progressive
// alignment in the engine, src/vettore/msa.py), kept free of React so they can be unit-tested.
import type { Difference, DifferenceCounts } from "./alignmentText";

export interface MultipleLimits {
  max_sequences: number;
  max_length: number;
  max_total: number;
  min_length: number;
}

/** Shown until the engine answers; the engine's own values replace them. */
export const DEFAULT_MULTIPLE_LIMITS: MultipleLimits = {
  max_sequences: 50, max_length: 15000, max_total: 150000, min_length: 15,
};

export interface MultipleRow {
  index: number;
  name: string;
  length: number;
  circular: boolean;
  strand: 1 | -1;
  rotation: number;
  origin_moved_to: number | null;
  /** The sequence as aligned (turned and rotated), with gaps. */
  row: string;
  first_column: number;
  last_column: number;
  identity: number;
}

export interface MultipleResult {
  alignmentTracks: unknown[];
  sequences: MultipleRow[];
  consensus: string;
  /** Per column, 0-100; null where fewer than two sequences cover the column. */
  conservation: (number | null)[];
  columns: number;
  warnings: string[];
}

export interface Comparison {
  differences: Difference[];
  counts: DifferenceCounts;
  identity: number;
  start: number | null;
  end: number | null;
  overlap: number;
  /** Stretches between compared segments with nothing in common, in the first sequence's
   *  numbering (start and end null when it has no bases there). */
  unrelated: { start: number | null; end: number | null; a_bases: number; b_bases: number }[];
  /** Positions where one of the two has N: not counted as differences. */
  unknown: number;
}

const n = (value: number) => value.toLocaleString("en-US");

/** Why a document cannot take part, or null if it can. */
export function ineligibleForMultiple(length: number, limits: MultipleLimits): string | null {
  if (length === 0) return "empty";
  if (length < limits.min_length) return `too short (under ${limits.min_length} bp)`;
  if (length > limits.max_length) return `too long (${n(length)} bp; the limit is ${n(limits.max_length)} bp)`;
  return null;
}

/**
 * Default selection: the active document, plus the other one when only two can take part.
 * As in the star alignment, unrelated open files are not ticked for the user.
 */
export function defaultMultipleChoice(
  candidates: { editorName: string; length: number }[],
  active: string | null,
  limits: MultipleLimits,
): string[] {
  const usable = candidates.filter((c) => !ineligibleForMultiple(c.length, limits));
  if (usable.length === 2) return usable.map((c) => c.editorName);
  return usable.some((c) => c.editorName === active) && active ? [active] : [];
}

/** Inline error for the choice, a prompt while fewer than two are ticked, or null. */
export function multipleChoiceProblem(lengths: number[], limits: MultipleLimits): string | null {
  if (lengths.length > limits.max_sequences) {
    return `Too many: at most ${limits.max_sequences}. Untick ${lengths.length - limits.max_sequences}.`;
  }
  const total = lengths.reduce((a, b) => a + b, 0);
  if (total > limits.max_total) {
    return `Together these are ${n(total)} bp; at most ${n(limits.max_total)} bp can be aligned at once. `
      + "Untick some, or align the regions you need.";
  }
  return null;
}

export function multiplePrompt(chosen: number): string | null {
  if (chosen >= 2) return null;
  return "Tick at least two sequences to line up with each other: variants of a construct, capsid "
    + "genes of several serotypes, versions of one plasmid.";
}

/** "Aligning… 40 %" style line for the progress bar. */
export function progressText(message: string, progress: number): string {
  return `${message}… ${Math.round(progress * 100)} %`;
}

export interface ConservationBin {
  /** First and last alignment column of the bin, 1-based. */
  start: number;
  end: number;
  /** Average conservation of the bin's columns, 0-100; null when no column has one. */
  mean: number | null;
  /** The least conserved column's value, 0-100 (null as for ``mean``). */
  min: number | null;
}

/** The conservation of every column folded into at most ``bins`` bars. Columns covered by
 *  a single sequence (null) compare nothing and are left out of the averages. */
export function conservationBins(conservation: (number | null)[], bins: number): ConservationBin[] {
  const count = Math.max(1, Math.min(bins, conservation.length));
  const out: ConservationBin[] = [];
  for (let b = 0; b < count; b++) {
    const start = Math.floor((b * conservation.length) / count);
    const end = Math.floor(((b + 1) * conservation.length) / count);
    if (end <= start) continue;
    const slice = conservation.slice(start, end).filter((c): c is number => c !== null);
    out.push({
      start: start + 1,
      end,
      mean: slice.length ? slice.reduce((a, c) => a + c, 0) / slice.length : null,
      min: slice.length ? Math.min(...slice) : null,
    });
  }
  return out;
}

/** "Columns 1,201–1,225: on average 84 % of the sequences share the most common base
 *  (lowest 50 %)". A gap is never agreement. */
export function describeBin(bin: ConservationBin): string {
  const where = bin.start === bin.end ? `Column ${n(bin.start)}` : `Columns ${n(bin.start)}–${n(bin.end)}`;
  if (bin.mean === null) return `${where}: only one sequence covers ${bin.start === bin.end ? "it" : "them"}`;
  const lowest = bin.start === bin.end || bin.min === Math.round(bin.mean) ? "" : ` (lowest ${bin.min} %)`;
  const share = bin.start === bin.end ? "" : "on average ";
  return `${where}: ${share}${Math.round(bin.mean)} % of the sequences share the most common base${lowest}`;
}

/** Bar class by conservation: the eye goes to the variable stretches. */
export function conservationLevel(mean: number | null): "high" | "mid" | "low" | "none" {
  if (mean === null) return "none";
  if (mean >= 90) return "high";
  if (mean >= 60) return "mid";
  return "low";
}

/** Stretches of at least ``minLength`` columns below ``threshold`` %, longest first. */
export function variableRegions(conservation: (number | null)[], threshold = 60, minLength = 10, max = 5)
  : { start: number; end: number }[] {
  const regions: { start: number; end: number }[] = [];
  let start = -1;
  for (let i = 0; i <= conservation.length; i++) {
    const value = conservation[i];
    const low = i < conservation.length && value !== null && value < threshold;
    if (low && start < 0) start = i;
    if (!low && start >= 0) {
      if (i - start >= minLength) regions.push({ start: start + 1, end: i });
      start = -1;
    }
  }
  return regions.sort((a, b) => (b.end - b.start) - (a.end - a.start) || a.start - b.start).slice(0, max);
}

/** The row's note for the summary table: strand and origin, in plain words. */
export function orientationNote(row: MultipleRow): string {
  const parts = [row.strand === 1 ? "forward" : "reverse complement"];
  if (row.origin_moved_to) parts.push(`origin moved to its base ${n(row.origin_moved_to)}`);
  return parts.join(", ");
}
