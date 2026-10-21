// Find Diagnostic Digest, the parts without React: which enzymes methylation blocks on the
// molecules to compare, and the protocol text copied for the bench.
import { enzymeRows, methylationWarnings } from "../enzymeLogic";

export interface PlannerLane {
  name: string;
  sizes: number[];
  uncut: boolean;
  text: string;
  bands: { size: number; apparent_size: number; rf: number; intensity: number; label: string }[];
}

export interface PlannerResult {
  enzymes: string[];
  /** The lab's buffer for each enzyme (File → Preferences → My enzymes), null when not set. */
  buffers: (string | null)[];
  /** The two enzymes go in different buffers. */
  mixed_buffers: boolean;
  same_pattern: string[][];
  sentence: string;
  clear: boolean;
  min_separation: number;
  lanes: PlannerLane[];
}

export interface PlannerResponse {
  results: PlannerResult[];
  diagnostic_count: number;
  tried: number;
  enzymes_tried: number;
  reason: string | null;
  ladder: { size: number; apparent_size: number; rf: number; intensity: number; label: string; reference?: boolean }[];
  ladder_label: string;
  percent: number;
}

/**
 * Enzymes that Dam or Dcm methylation stops at one of their sites on any of these molecules,
 * as REBASE records it (the same table as the Enzymes view). Plasmid DNA from the usual E. coli
 * strains carries both methylations, so these enzymes would leave sites uncut.
 */
export function blockedByMethylation(molecules: { sequence: string; circular: boolean }[]): string[] {
  const blocked = new Set<string>();
  for (const m of molecules) {
    if (!m.sequence) continue;
    const rows = enzymeRows("commercial", m.sequence, m.circular);
    const notes = methylationWarnings(m.sequence, m.circular, rows);
    for (const [name, list] of Object.entries(notes)) {
      if (list.some((n) => n.blocking)) blocked.add(name);
    }
  }
  return [...blocked].sort((a, b) => a.localeCompare(b));
}

/** "EcoRI (rCutSmart) + SalI (r3.1)": each enzyme with its buffer, when the lab gave one. */
export function enzymesWithBuffers(result: Pick<PlannerResult, "enzymes" | "buffers">): string {
  return result.enzymes
    .map((e, i) => (result.buffers?.[i]?.trim() ? `${e} (${result.buffers[i]!.trim()})` : e))
    .join(" + ");
}

/** What to do when the two enzymes of a double digest go in different buffers. */
export function mixedBuffersAdvice(result: PlannerResult): string {
  const [a, b] = result.enzymes;
  return `${a} and ${b} go in different buffers: digest with one, then the other, or use a buffer both are active in.`;
}

const bp = (n: number) => `${n.toLocaleString("en-US")} bp`;

function laneLine(lane: PlannerLane, circular: boolean): string {
  if (lane.uncut) {
    return `${bp(lane.sizes[0])}, not cut${circular ? " (a supercoiled plasmid runs ahead of its size)" : ""}`;
  }
  return lane.sizes.map((s) => (s < 150 ? `${bp(s)} (too small to see)` : bp(s))).join(" + ");
}

/** The digest as a short protocol: enzymes, gel, and the bands of every outcome. */
export function protocolText(
  expectedName: string,
  result: PlannerResult,
  gel: { percent: number; ladderLabel: string },
  circular: boolean[],
): string {
  const kind = result.enzymes.length === 1 ? "single digest" : "double digest";
  const lines = [
    `Diagnostic digest for “${expectedName}”`,
    `Digest: ${enzymesWithBuffers(result)} (${kind})`,
    ...(result.mixed_buffers ? [`Note: ${mixedBuffersAdvice(result)}`] : []),
    `Gel: ${gel.percent.toFixed(1)} % agarose, ladder ${gel.ladderLabel}`,
    "",
    "Expected bands:",
    ...result.lanes.map((lane, i) => `  ${lane.name}: ${laneLine(lane, circular[i] ?? true)}`),
  ];
  if (result.same_pattern.length) {
    lines.push("", `Same bands with: ${result.same_pattern.map((e) => e.join(" + ")).join(", ")}`);
  }
  lines.push("", "Sizes predicted by Vettore (Find Diagnostic Digest); fragments under 150 bp are not counted.");
  return lines.join("\n");
}
