import type { Doc } from "./documents";

/** S1 history: provenance of the document (opened, created, saved). Cloning steps join in S2. */
export default function HistoryView({ doc }: { doc: Doc }) {
  return (
    <div className="vt-history">
      <h2>History of “{doc.name}”</h2>
      <table className="bp3-html-table bp3-html-table-condensed bp3-html-table-striped">
        <thead>
          <tr><th>When</th><th>Action</th><th>Details</th></tr>
        </thead>
        <tbody>
          {[...doc.history].reverse().map((h, i) => (
            <tr key={i}><td>{h.time}</td><td>{h.action}</td><td>{h.detail ?? ""}</td></tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
