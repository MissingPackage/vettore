// "Copy as FASTA": the selection, or the whole sequence, as a FASTA record for emails and chats.

export interface Range {
  /** 0-based, inclusive; start > end crosses the origin of a circular sequence (OVE's model). */
  start: number;
  end: number;
}

/** The selection as a range, or null when nothing is selected (only a caret, or no layer). */
export function selectionRange(
  length: number,
  layer: { start?: number; end?: number } | null | undefined,
): Range | null {
  if (!layer || layer.start === undefined || layer.end === undefined) return null;
  const { start, end } = layer;
  if (start < 0 || end < 0 || start >= length || end >= length) return null;
  return { start, end };
}

export function rangeLength(range: Range, length: number): number {
  return range.start <= range.end ? range.end - range.start + 1 : length - range.start + range.end + 1;
}

export function subsequence(sequence: string, range: Range): string {
  return range.start <= range.end
    ? sequence.slice(range.start, range.end + 1)
    : sequence.slice(range.start) + sequence.slice(0, range.end + 1);
}

/**
 * A FASTA record with 60 bases per line. The identifier is the document name without spaces;
 * for a selection it carries the 1-based positions, as in "pUC19:101-623".
 */
export function fastaRecord(name: string, sequence: string, range: Range | null): string {
  const id = (name.trim() || "sequence").replace(/\s+/g, "_");
  const bases = range ? subsequence(sequence, range) : sequence;
  const header = range ? `>${id}:${range.start + 1}-${range.end + 1}` : `>${id}`;
  const lines = [header];
  for (let i = 0; i < bases.length; i += 60) lines.push(bases.slice(i, i + 60));
  return lines.join("\n") + "\n";
}
