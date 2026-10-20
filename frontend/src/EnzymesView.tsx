import { useMemo, useState } from "react";
import { useSelector } from "react-redux";
import { actions } from "@teselagen/ove";
import { ENZYME_SETS, enzymeRows, methylationWarnings, type EnzymeRow, type EnzymeSetId } from "./enzymeLogic";
import { store } from "./store";

type SortKey = "name" | "sites" | "location";

/** Enzymes view: enzymes of the chosen set with number of sites, cut locations and ends. */
export default function EnzymesView({ editorName, set }: { editorName: string; set: EnzymeSetId }) {
  const { sequence, circular } = useSelector(
    (state: any) => ({
      sequence: (state.VectorEditor[editorName]?.sequenceData?.sequence ?? "") as string,
      circular: !!state.VectorEditor[editorName]?.sequenceData?.circular,
    }),
    (a, b) => a.sequence === b.sequence && a.circular === b.circular,
  );
  const [query, setQuery] = useState("");
  const [sort, setSort] = useState<{ key: SortKey; desc: boolean }>({ key: "location", desc: false });
  const [selected, setSelected] = useState<string | null>(null);

  const rows = useMemo(() => enzymeRows(set, sequence, circular), [set, sequence, circular]);
  const methylation = useMemo(() => methylationWarnings(sequence, circular, rows),
                              [sequence, circular, rows]);
  const shown = useMemo(() => {
    const q = query.trim().toLowerCase();
    const filtered = q ? rows.filter((r) => r.name.toLowerCase().includes(q)) : rows;
    const value = (r: EnzymeRow) =>
      sort.key === "name" ? r.name.toLowerCase() : sort.key === "sites" ? r.cuts.length : r.cuts[0]?.location ?? 0;
    return [...filtered].sort((a, b) => {
      const va = value(a);
      const vb = value(b);
      const cmp = va < vb ? -1 : va > vb ? 1 : a.name.localeCompare(b.name);
      return sort.desc ? -cmp : cmp;
    });
  }, [rows, query, sort]);

  const header = (key: SortKey, label: string) => (
    <th
      className="vt-sortable"
      onClick={() => setSort((s) => ({ key, desc: s.key === key ? !s.desc : false }))}
    >
      {label}
      {sort.key === key ? (sort.desc ? " ▼" : " ▲") : ""}
    </th>
  );

  const select = (row: EnzymeRow, cutIndex = 0) => {
    setSelected(row.name);
    const cut = row.cuts[cutIndex];
    if (cut) store.dispatch(actions.selectionLayerUpdate({ start: cut.start, end: cut.end }, { editorName }));
  };

  const setLabel = ENZYME_SETS.find((s) => s.id === set)?.label;

  return (
    <div className="vt-table-view" data-testid="enzymes-view">
      <div className="vt-table-toolbar">
        <input
          className="bp3-input"
          placeholder="Find an enzyme by name"
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          data-testid="enzyme-search"
          title="Type a name, for example BamHI"
        />
        <span className="vt-muted">
          {rows.length === 0
            ? `No enzyme of “${setLabel}” cuts this sequence`
            : `${rows.length} enzyme${rows.length === 1 ? "" : "s"} of “${setLabel}” cut this `
              + `sequence${query.trim() && shown.length !== rows.length
                ? `, ${shown.length} shown` : ""}`}
          . Click a position to select that site.
        </span>
      </div>
      <table className="bp3-html-table bp3-html-table-condensed bp3-interactive vt-table">
        <thead>
          <tr>
            {header("name", "Enzyme")}
            <th title="5′→3′ on the top strand">Recognition site</th>
            {header("sites", "Sites")}
            {header("location", "Location(s)")}
            <th title="What the cut leaves: blunt ends, or an overhang of so many bases">
              Ends
            </th>
            <th title="Dam (GATC) or Dcm (CCWGG) methylation inside the site, in the usual
                       cloning strains">
              Methylation
            </th>
          </tr>
        </thead>
        <tbody>
          {shown.map((r) => (
            <tr key={r.name} className={r.name === selected ? "vt-row-selected" : ""} onClick={() => select(r)}>
              <td className={r.cuts.length === 1 ? "vt-bold" : ""}>{r.name}</td>
              <td className="vt-mono">{r.site}</td>
              <td>{r.cuts.length}</td>
              <td>
                {r.cuts.map((c, i) => (
                  <a
                    key={i}
                    href="#"
                    className="vt-cut"
                    onClick={(e) => {
                      e.preventDefault();
                      e.stopPropagation();
                      select(r, i);
                    }}
                  >
                    {c.location}
                  </a>
                ))}
              </td>
              <td>
                {r.ends}
                {r.overhang ? <span className="vt-mono vt-muted"> {r.overhang}</span> : null}
              </td>
              <td>
                {(methylation[r.name] ?? []).map((note) => (
                  <span
                    key={note.kind}
                    className={note.blocking ? "vt-warn" : "vt-muted"}
                    title={note.blocking
                      ? `REBASE records this enzyme as "${note.status}" by ${note.kind} methylation, and a ${note.kind} methylated base lies inside this site. Grow the plasmid in a ${note.kind === "Dam" ? "dam⁻" : "dcm⁻"} strain, or pick another enzyme.`
                      : `A ${note.kind} methylated base lies inside this site, but REBASE records no sensitivity for this enzyme. Check your supplier before relying on the digest.`}
                  >
                    {note.kind}: {note.status}{" "}
                  </span>
                ))}
              </td>
            </tr>
          ))}
        </tbody>
      </table>
      {!shown.length && (
        <p className="vt-muted" data-testid="enzymes-empty">
          {query.trim()
            ? `No enzyme in “${setLabel}” matches “${query.trim()}”. It may cut more often than `
              + "this set allows: try another set in the box at the top right, or Enzymes → "
              + "Choose Enzymes… for the full list."
            : `No enzyme of “${setLabel}” cuts this sequence. Pick a wider set in the box at `
              + "the top right (Unique & Dual Cutters, or Nonredundant Commercial)."}
        </p>
      )}
      {Object.keys(methylation).length > 0 && (
        <p className="vt-muted" data-testid="methylation-note">
          A base methylated by Dam (GATC) or Dcm (CCWGG) in the usual cloning strains lies
          inside these recognition sites, and REBASE records the enzyme as affected. Grow the
          plasmid in a dam⁻/dcm⁻ strain or choose another enzyme. Entries marked “not recorded”
          have the overlap but no sensitivity data in REBASE — check your supplier.
        </p>
      )}
    </div>
  );
}
