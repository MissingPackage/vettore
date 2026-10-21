import { useEffect, useMemo, useState } from "react";
import { api, type SequenceData } from "./api";
import { baseName, type Doc } from "./documents";
import { buildLineage, lineageDepth, lineageText, type LineageNode } from "./lineage";
import { useCopyText } from "./cloning/useCopyText";
import "./history.css";

/** Every path the tree mentions, so the engine is asked once which files are still there. */
function pathsIn(node: LineageNode): string[] {
  const own = node.path ? [node.path] : [];
  return [...own, ...node.children.flatMap(pathsIn)];
}

function Row({ node, existing, onReopen, onOpenFile }: {
  node: LineageNode;
  existing: Set<string>;
  onReopen: (data: SequenceData, from: string) => void;
  onOpenFile: (path: string) => void;
}) {
  const fileThere = !!node.path && existing.has(node.path);
  return (
    <li className="vt-lineage-node" data-testid="lineage-node" data-depth={node.depth}>
      <div className="vt-lineage-row">
        <span className="vt-lineage-name">{node.name}</span>
        <span className="vt-lineage-size">{node.size}</span>
        <span className="vt-lineage-step">{node.operation}</span>
        {node.usedBy && (
          <span className="vt-lineage-used">used by {node.usedBy}</span>
        )}
        {node.time && <span className="vt-lineage-time">{node.time}</span>}
        {node.depth > 0 && node.data && (
          <button
            className="bp3-button bp3-small bp3-minimal"
            data-testid="lineage-open-snapshot"
            title={`Open “${node.name}” as it was when this step ran`}
            onClick={() => onReopen(node.data as SequenceData, `the lineage of “${node.name}”`)}
          >
            Open the snapshot
          </button>
        )}
        {fileThere && (
          <button
            className="bp3-button bp3-small bp3-minimal"
            data-testid="lineage-open-file"
            title={`Open “${baseName(node.path as string)}” as it is on disk now`}
            onClick={() => onOpenFile(node.path as string)}
          >
            Open the file
          </button>
        )}
        {node.path && !fileThere && node.depth > 0 && (
          <span className="vt-lineage-gone" data-testid="lineage-file-gone">
            “{baseName(node.path)}” is no longer there
          </span>
        )}
        {node.dropped && (
          <span className="vt-lineage-gone" data-testid="lineage-no-snapshot">
            no copy was kept in the file
          </span>
        )}
      </div>
      {node.detail && <div className="vt-lineage-detail">{node.detail}</div>}
      {node.children.length > 0 && (
        <ul className="vt-lineage-children">
          {node.children.map((child) => (
            <Row key={child.id} node={child} existing={existing}
                 onReopen={onReopen} onOpenFile={onOpenFile} />
          ))}
        </ul>
      )}
    </li>
  );
}

/** The Lineage tab: the document, and the tree of what it was made from (S6.1). */
export default function LineageView({ doc, data, onReopen, onOpenFile }: {
  doc: Doc;
  /** The document as it is now, for the size of the root. */
  data: SequenceData | undefined;
  onReopen: (data: SequenceData, from: string) => void;
  onOpenFile: (path: string) => void;
}) {
  const tree = useMemo(() => buildLineage(doc, data), [doc, data]);
  const [existing, setExisting] = useState<Set<string>>(new Set());
  const { copy, fallbackDialog } = useCopyText();
  const paths = useMemo(() => [...new Set(pathsIn(tree))], [tree]);

  useEffect(() => {
    if (!paths.length) {
      setExisting(new Set());
      return;
    }
    let live = true;
    // Which files are still on disk decides "Open the file" versus "Open the snapshot".
    api.pathsExist(paths)
      .then((r) => { if (live) setExisting(new Set(r.existing)); })
      .catch(() => { if (live) setExisting(new Set()); });
    return () => { live = false; };
  }, [paths]);

  const steps = lineageDepth(tree);
  return (
    <div className="vt-history vt-lineage" data-testid="lineage-view">
      <h2>Lineage of “{doc.name}”</h2>
      <p className="vt-history-intro">
        {tree.children.length
          ? `Where this construct came from, ${steps === 1 ? "one step" : `${steps} steps`} back. `
            + "Each line is a document as it was when the step below it ran: open the snapshot "
            + "to see it exactly as it was, or the file to see it as it is now."
          : "This document was not made from other documents in Vettore, so there is nothing "
            + "above it yet. A cloning step keeps a copy of everything it was run on, and that "
            + "is what builds this tree."}
      </p>
      <ul className="vt-lineage-tree">
        <Row node={tree} existing={existing} onReopen={onReopen} onOpenFile={onOpenFile} />
      </ul>
      {tree.children.length > 0 && (
        <button
          className="bp3-button bp3-small"
          data-testid="lineage-copy"
          onClick={() => copy(lineageText(tree), "Copy the lineage")}
        >
          Copy the lineage
        </button>
      )}
      {fallbackDialog}
    </div>
  );
}
