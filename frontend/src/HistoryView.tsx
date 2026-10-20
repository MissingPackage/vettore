import type { SequenceData } from "./api";
import type { Doc } from "./documents";

/** Provenance of the document: opened, created, saved, and every cloning step with its inputs. */
export default function HistoryView({ doc, onReopen }: {
  doc: Doc;
  /** Open a snapshot of an input document in a new tab. */
  onReopen: (data: SequenceData, from: string) => void;
}) {
  return (
    <div className="vt-history">
      <h2>History of “{doc.name}”</h2>
      <table className="bp3-html-table bp3-html-table-condensed bp3-html-table-striped">
        <thead>
          <tr><th>When</th><th>Action</th><th>Details</th><th>Inputs</th></tr>
        </thead>
        <tbody>
          {[...doc.history].reverse().map((h, i) => (
            <tr key={i}>
              <td>{h.time}</td>
              <td>{h.action}</td>
              <td>{h.detail ?? ""}</td>
              <td>
                {(h.parents ?? []).map((p, j) => (
                  <button
                    key={j}
                    className="bp3-button bp3-small bp3-minimal"
                    data-testid="history-reopen"
                    title={`Open “${p.name}” as it was when this step ran`}
                    onClick={() => onReopen(p.data, `${h.action} in “${doc.name}”`)}
                  >
                    Open “{p.name}”
                  </button>
                ))}
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
