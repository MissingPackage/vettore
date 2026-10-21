/** The plain sentences of the Sanger contig: coverage, the verdict on the consensus, differences. */
import type { SequenceData } from "./api";

/** Phred Q30, one wrong call in 1,000: one read this good confirms a difference on its own. */
export const CONFIRMING_QUALITY = 30;
/** Phred Q20: two reads this good confirm a difference together; one this good against it is enough to doubt it. */
export const SUPPORTING_QUALITY = 20;

export interface ContigDifference {
  kind: "mismatch" | "insertion" | "deletion";
  /** 1-based on the reference; for an insertion, the reference base that follows it. */
  position: number;
  reference: string;
  consensus: string;
  /** Reads that show it. */
  reads: number;
  /** The best quality among them. */
  quality: number;
  confirmed: boolean;
}

export interface ContigCall {
  read: string;
  strand: number;
  call: string;
  quality: number;
}

export interface Disagreement {
  kind: "base" | "insertion" | "no call";
  position: number;
  reference: string;
  consensus: string;
  calls: ContigCall[];
}

export interface Contig {
  reads: number;
  consensus_length: number;
  /** Runs of 0 (no read), 1 (one strand) or 2 (both strands), 1-based inclusive; start > end across the origin. */
  coverage: { start: number; end: number; strands: number }[];
  covered: { start: number; end: number }[];
  /** Covered, clearly called, nothing left to check; `reads` is the fewest reads over any base. */
  verified: { start: number; end: number; reads: number }[];
  /** Covered positions where no read makes a clear call: the reference is kept there. */
  unclear: number;
  /** Covered positions where good reads conflict, the vote is tied, or only calls under Q20 back the consensus. */
  doubtful: number;
  differences: ContigDifference[];
  disagreements: Disagreement[];
  document?: SequenceData;
}

const n = (value: number) => value.toLocaleString("en-US");

export function span(start: number, end: number): string {
  return start === end ? n(start) : `${n(start)}–${n(end)}`;
}

/** "a", "a and b", "a, b and c"; past `max` items, "a, b and 3 more stretches". */
export function listed(items: string[], max = 4): string {
  const shown = items.length > max ? [...items.slice(0, max - 1), `${n(items.length - max + 1)} more stretches`] : items;
  return shown.length <= 1 ? shown.join("") : `${shown.slice(0, -1).join(", ")} and ${shown[shown.length - 1]}`;
}

/** "223–861 covered on both strands; 862–903 on one strand only; no read over 1–222 and 904–1,125." */
export function coverageSentence(contig: Contig): string {
  const of = (strands: number) => contig.coverage.filter((r) => r.strands === strands).map((r) => span(r.start, r.end));
  const both = of(2);
  const one = of(1);
  const none = of(0);
  if (!both.length && !one.length) return "No read covers the reference.";
  if (!one.length && !none.length) return "Every base is covered on both strands.";
  const parts: string[] = [];
  if (both.length) parts.push(`${listed(both)} covered on both strands`);
  if (one.length) parts.push(both.length ? `${listed(one)} on one strand only` : `${listed(one)} covered on one strand only`);
  if (none.length) parts.push(`no read over ${listed(none)}`);
  const text = parts.join("; ");
  return `${text[0].toUpperCase()}${text.slice(1)}.`;
}

export function describeDifference(d: ContigDifference): string {
  if (d.kind === "mismatch") return `${n(d.position)}: ${d.reference} → ${d.consensus}`;
  if (d.kind === "deletion") return `${n(d.position)}: ${d.reference} missing in the clone`;
  return `${n(d.position)}: ${d.consensus} extra in the clone`;
}

/** Why a difference is trusted, or why to look at its trace. */
export function evidence(d: ContigDifference): string {
  const seen = d.reads > 1 ? `seen in ${n(d.reads)} reads, best Q${d.quality}` : `1 read at Q${d.quality}`;
  if (d.confirmed) return seen;
  if (d.reads === 1 && d.quality < CONFIRMING_QUALITY) return `${seen}, under Q${CONFIRMING_QUALITY}`;
  return `${seen}; another read shows something else, or the reads are under Q${SUPPORTING_QUALITY}`;
}

/** One sentence that answers "is my clone right?" from the consensus of every read. */
export function contigVerdict(contig: Contig, referenceName: string, coveragePercent: number): { text: string; ok: boolean } {
  const where = `${coveragePercent} % of ${referenceName}`;
  const total = contig.differences.length;
  const disputed = contig.disagreements.filter((d) => d.kind !== "no call").length;
  const settled = (disputed
    ? ` The reads disagree at ${n(disputed)} place${disputed === 1 ? "" : "s"} (listed under Contig); the consensus`
      + " takes the call with the most quality behind it, and the plasmid's base where the reads are even."
    : "")
    + (contig.doubtful
      ? ` ${n(contig.doubtful)} covered position${contig.doubtful === 1 ? " is" : "s are"} not verified: good reads`
        + ` conflict there, a difference is not confirmed, or only calls under Q${SUPPORTING_QUALITY} back it. Check ${contig.doubtful === 1 ? "it" : "them"}`
        + " in the traces."
      : "")
    + (contig.unclear
      ? ` ${n(contig.unclear)} covered position${contig.unclear === 1 ? " has" : "s have"} no clear call in any read:`
        + " the plasmid's base is kept there, unverified."
      : "");
  if (!total) {
    return {
      ok: !contig.unclear && !contig.doubtful,
      text: `No differences: the consensus of ${n(contig.reads)} reads matches ${referenceName} wherever they cover it`
        + ` (${where}).${settled}`,
    };
  }
  const confirmed = contig.differences.filter((d) => d.confirmed).length;
  const check = total - confirmed;
  let how: string;
  if (!check) how = total === 1 ? "It is confirmed" : total === 2 ? "Both are confirmed" : `All ${n(total)} are confirmed`;
  else if (!confirmed) how = total === 1 ? "It needs a look at the trace" : `${total === 2 ? "Both" : `All ${n(total)}`} need a look at the trace`;
  else how = `${n(confirmed)} ${confirmed === 1 ? "is" : "are"} confirmed, ${n(check)} need${check === 1 ? "s" : ""} a look at the trace`;
  return {
    ok: false,
    text: `${n(total)} difference${total === 1 ? "" : "s"} between the consensus of ${n(contig.reads)} reads and`
      + ` ${referenceName} (covering ${where}). ${how}: a difference counts as confirmed when two reads show it`
      + ` at Q${SUPPORTING_QUALITY} or better, or one read at Q${CONFIRMING_QUALITY} or better, and no read at`
      + ` Q${SUPPORTING_QUALITY} or better shows something else.${settled}`,
  };
}
