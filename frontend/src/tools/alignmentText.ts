// Wording and choices for Tools → Align Multiple Sequences, kept free of React so they can be
// unit-tested. Positions come from the engine 1-based, in the reference's numbering.

export interface Difference {
  kind: "mismatch" | "insertion" | "deletion";
  /** Mismatch and deletion: first and last reference base. Insertion: the two flanking bases. */
  start: number;
  end: number;
  reference: string;
  sequence: string;
  /** All to each other only: near an end of the two, where nothing beyond holds the
   *  alignment, so the same bases could read as another mix of substitutions and indels. */
  uncertain?: boolean;
}

export interface DifferenceCounts {
  mismatches: number;
  insertions: number;
  deletions: number;
  inserted_bases: number;
  deleted_bases: number;
}

export interface AlignmentLimits {
  max_sequences: number;
  max_length: number;
  min_length: number;
}

/** Shown until the engine answers; the engine's own values replace them. */
export const DEFAULT_LIMITS: AlignmentLimits = { max_sequences: 10, max_length: 15000, min_length: 15 };

const n = (value: number) => value.toLocaleString("en-US");

function bases(text: string, max = 12): string {
  return text.length > max ? `${text.slice(0, max)}…` : text;
}

/** One difference in words, e.g. "1,204: A → G" or "between 1,001 and 1,002: 7 bp extra (GATTACA)". */
export function describeDifference(d: Difference): string {
  if (d.kind === "mismatch") return `${n(d.start)}: ${d.reference} → ${d.sequence}`;
  if (d.kind === "deletion") {
    const where = d.start === d.end ? n(d.start) : `${n(d.start)}–${n(d.end)}`;
    return `${where}: ${d.reference.length} bp missing (${bases(d.reference)})`;
  }
  return `between ${n(d.start)} and ${n(d.end)}: ${d.sequence.length} bp extra (${bases(d.sequence)})`;
}

const plural = (count: number, one: string, many: string) => `${n(count)} ${count === 1 ? one : many}`;

/** "none", or e.g. "3 mismatches · 1 insertion (7 bp) · 1 deletion (4 bp)". */
export function summarizeCounts(c: DifferenceCounts): string {
  const parts: string[] = [];
  if (c.mismatches) parts.push(plural(c.mismatches, "mismatch", "mismatches"));
  if (c.insertions) parts.push(`${plural(c.insertions, "insertion", "insertions")} (${n(c.inserted_bases)} bp)`);
  if (c.deletions) parts.push(`${plural(c.deletions, "deletion", "deletions")} (${n(c.deleted_bases)} bp)`);
  return parts.length ? parts.join(" · ") : "none";
}

/** Tab-separated, one difference per line, for pasting into a spreadsheet. */
export function differencesAsTsv(name: string, diffs: Difference[]): string {
  const notes = diffs.some((d) => d.uncertain);
  const header = "sequence\tkind\tstart\tend\treference bases\tsequence bases" + (notes ? "\tnote" : "");
  const rows = diffs.map((d) => [name, d.kind, d.start, d.end, d.reference, d.sequence,
    ...(notes ? [d.uncertain ? "near the end: uncertain" : ""] : [])].join("\t"));
  return [header, ...rows].join("\n");
}

/**
 * The name each open document goes by in the dialog and its results: its tab name, plus its tab
 * number when another open tab has the same name ("Pasted sequence (tab 3)").
 */
export function distinctNames(docs: { editorName: string; name: string }[]): Record<string, string> {
  const uses = new Map<string, number>();
  for (const d of docs) uses.set(d.name, (uses.get(d.name) ?? 0) + 1);
  return Object.fromEntries(docs.map((d, i) => [
    d.editorName,
    (uses.get(d.name) ?? 0) > 1 ? `${d.name} (tab ${i + 1})` : d.name,
  ]));
}

/** Why a document cannot take part, or null if it can. */
export function ineligible(length: number, limits: AlignmentLimits): string | null {
  if (length === 0) return "empty";
  if (length < limits.min_length) return `too short (under ${limits.min_length} bp)`;
  if (length > limits.max_length) return `too long (${n(length)} bp; the limit is ${n(limits.max_length)} bp)`;
  return null;
}

/**
 * Default selection: the other document when it is the only one that can take part, else
 * nothing. Unrelated open files only share chance matches with the reference, so the user
 * ticks what they expect to resemble it.
 */
export function defaultChoice(
  candidates: { editorName: string; length: number }[],
  reference: string | null,
  limits: AlignmentLimits,
): string[] {
  const usable = candidates.filter((c) => c.editorName !== reference && !ineligible(c.length, limits));
  return usable.length === 1 ? [usable[0].editorName] : [];
}

/** Neutral prompt while nothing is ticked yet, or null. */
export function choicePrompt(chosen: number): string | null {
  if (chosen > 0) return null;
  return "Tick the sequences to compare with the reference: versions of the same plasmid, a gene "
    + "and its source, reads of a clone. Unrelated files only share short chance matches.";
}

/** Inline error for the "Compare with" list, or null when the choice can be aligned. */
export function choiceProblem(chosen: number, limits: AlignmentLimits): string | null {
  const max = limits.max_sequences - 1;
  if (chosen > max) {
    return `Too many: at most ${max} besides the reference. Untick ${chosen - max}.`;
  }
  return null;
}

export const MODE_HELP: Record<string, string> = {
  "semi-global": "Aligned end to end; overhanging ends are not compared.",
  global: "Two circular sequences, aligned all the way round.",
  local: "Only part of this sequence resembles the reference; the best-matching region is shown.",
};
