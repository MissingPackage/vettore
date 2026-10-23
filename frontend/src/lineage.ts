// The lineage of a document: the tree from the document back to the files it was made from.
// Built from the History, which cloning steps fill with a snapshot of every input (S6.1).
import type { SequenceData } from "./api";
import { baseName } from "./files/formats";
import { LINEAGE_DEPTH, type HistoryEntry, type HistoryParent } from "./history";
import type { Doc } from "./documents";

export interface LineageNode {
  /** Unique within the tree, so a row can be opened and collapsed. */
  id: string;
  name: string;
  /** "4,329 bp, circular", or "size unknown" when the file kept no snapshot. */
  size: string;
  /** The step this node came out of ("Restriction cloning"), or how it reached Vettore. */
  operation: string;
  /** The step that was run on this node, which is what joins it to the node above. */
  usedBy?: string;
  time: string;
  detail?: string;
  /** The input as it was when the step ran; missing when the snapshot was not kept. */
  data?: SequenceData;
  /** The file the input was open from, when it had one. */
  path?: string | null;
  /** The saved file left this snapshot out because the History did not fit. */
  dropped: boolean;
  /** How far from the document: 0 is the document itself. */
  depth: number;
  children: LineageNode[];
}

export function sizeOf(data: SequenceData | undefined): string {
  if (!data || typeof data.sequence !== "string") return "size unknown";
  const bp = data.sequence.length.toLocaleString("en-US");
  return `${bp} bp, ${data.circular ? "circular" : "linear"}`;
}

/**
 * The step a document came out of: the newest step that was run on inputs (a cloning), or
 * else the first step of all, which is how the document came into being — opened from a file,
 * created, pasted. Never the newest step of all, because saving, exporting and printing are
 * steps too: a construct that was cloned and then saved came out of the cloning.
 */
function madeBy(history: HistoryEntry[]): HistoryEntry | undefined {
  return stepsWithParents(history)[0] ?? history[0];
}

/** Every step that was run on inputs, newest first: one branch of the tree each. */
function stepsWithParents(history: HistoryEntry[]): HistoryEntry[] {
  return [...history].reverse().filter((step) => step.parents?.length);
}

function parentNode(parent: HistoryParent, step: HistoryEntry, id: string,
                    depth: number): LineageNode {
  const own = parent.history ?? [];
  const made = madeBy(own);
  return {
    id,
    name: parent.name,
    size: sizeOf(parent.data),
    // What the input itself came out of; a file saved by another program says nothing of it.
    operation: made?.action ?? "where it came from is not recorded",
    usedBy: step.action,
    time: made?.time ?? step.time,
    detail: made?.detail,
    data: parent.data,
    path: parent.path ?? null,
    dropped: !parent.data,
    depth,
    children: depth + 1 >= LINEAGE_DEPTH ? [] : childrenOf(own, id, depth + 1),
  };
}

function childrenOf(history: HistoryEntry[], prefix: string, depth: number): LineageNode[] {
  const out: LineageNode[] = [];
  stepsWithParents(history).forEach((step, i) => {
    (step.parents ?? []).forEach((parent, j) => {
      out.push(parentNode(parent, step, `${prefix}.${i}-${j}`, depth));
    });
  });
  return out;
}

/**
 * The document and everything it was made from. `data` is the document as it is now, so the
 * root shows its current size even while the snapshots below it show the inputs as they were.
 */
export function buildLineage(doc: Doc, data: SequenceData | undefined): LineageNode {
  const made = madeBy(doc.history);
  return {
    id: "root",
    name: doc.name,
    size: sizeOf(data),
    operation: made?.action ?? "Opened",
    time: made?.time ?? "",
    detail: made?.detail,
    data,
    path: doc.path,
    dropped: false,
    depth: 0,
    children: childrenOf(doc.history, "n", 1),
  };
}

/** How many steps the tree goes back, counting the document as one. */
export function lineageDepth(node: LineageNode): number {
  // Folded rather than spread: a file can carry more inputs than a call takes arguments.
  return 1 + node.children.reduce((deepest, child) => Math.max(deepest, lineageDepth(child)), 0);
}

/** The lineage as plain lines, for the clipboard (the indentation shows the tree). */
export function lineageText(node: LineageNode): string {
  const pad = "  ".repeat(node.depth);
  const where = node.path ? ` [${baseName(node.path)}]` : "";
  const used = node.usedBy ? `, used by ${node.usedBy}` : "";
  const head = `${pad}${node.name} — ${node.operation} (${node.size})${used}${where}`;
  return [head, ...node.children.map(lineageText)].join("\n");
}
