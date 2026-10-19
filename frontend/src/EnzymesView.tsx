import { useMemo, useState } from "react";
import { useSelector } from "react-redux";
import { actions } from "@teselagen/ove";
import { ENZYME_SETS, enzymeRows, type EnzymeRow, type EnzymeSetId } from "./enzymeLogic";
import { store } from "./store";

type SortKey = "name" | "sites" | "location";

/** Enzymes view: enzymes of the chosen set with number of sites and cut locations. */
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
          placeholder="Find enzymes"
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          data-testid="enzyme-search"
        />
        <span className="vt-muted">
          {rows.length} enzymes in “{setLabel}” cut this sequence
        </span>
      </div>
      <table className="bp3-html-table bp3-html-table-condensed bp3-interactive vt-table">
        <thead>
          <tr>
            {header("name", "Enzyme")}
            <th>Recognition site</th>
            {header("sites", "Sites")}
            {header("location", "Location(s)")}
            <th>Ends</th>
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
            </tr>
          ))}
        </tbody>
      </table>
      {!shown.length && <p className="vt-muted">No enzymes of this set cut the sequence.</p>}
    </div>
  );
}
