// Primers → Primer Library: the primers the lab already keeps, in a .csv the lab shares in a
// folder (OneDrive, Dropbox, a lab drive). Vettore reads the file from the path the user chose
// and writes it back there; the engine (src/vettore/primer_library.py) owns the shared file,
// the one binding rule and the products. This dialog lists and searches the primers, imports
// from a spreadsheet, adds and removes them, and shows which bind the open document.
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import { api, desktop, whenDesktopReady } from "../api";
import { Field } from "../cloning/inputs";
import { useCopyText } from "../cloning/useCopyText";
import { messageOf } from "../files/notices";
import {
  baseName,
  bindingText,
  cleanSequence,
  filterPrimers,
  mergeImport,
  newEntry,
  orderingText,
  toCsv,
  type BoundPrimer,
  type LibraryFile,
  type LibraryPrimer,
  type RejectedRow,
} from "./primerLibraryText";
import type { Tool, ToolProps } from "./types";
import "./primerLibrary.css";

// Blueprint 3 typings predate React 18 (no implicit `children`).
const { Dialog, Callout } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

const EMPTY: LibraryFile = { path: "", exists: false, mtime: null, primers: [], rejected: [] };
type Notice = { intent: "success" | "warning" | "danger"; text: string } | null;

function PrimerLibraryDialog(props: ToolProps) {
  const { copy, fallbackDialog } = useCopyText();
  const [file, setFile] = useState<LibraryFile>(EMPTY);
  const [primers, setPrimers] = useState<LibraryPrimer[]>([]);
  const [dirty, setDirty] = useState(false);
  const [query, setQuery] = useState("");
  const [notice, setNotice] = useState<Notice>(null);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [adding, setAdding] = useState({ name: "", sequence: "", notes: "", location: "" });
  const [importText, setImportText] = useState<string | null>(null);
  const [bound, setBound] = useState<BoundPrimer[] | null>(null);
  const [today, setToday] = useState("");
  const fileInput = useRef<HTMLInputElement>(null);

  const load = useCallback(async () => {
    setBusy(true);
    setError(null);
    try {
      const got = await api.primerLibrary.get();
      setFile(got);
      setPrimers(got.primers);
      setDirty(false);
      if (got.rejected.length) {
        setNotice({ intent: "warning", text: rejectedText(got.rejected, "in the library file") });
      }
    } catch (e) {
      setError(messageOf(e));
    } finally {
      setBusy(false);
    }
  }, []);

  useEffect(() => {
    load();
    api.primerLibrary.parse("").then((r) => setToday(r.today)).catch(() => undefined);
  }, [load]);

  // Which library primers bind the active document (the "on this sequence" panel and the map
  // toggle use the same engine rule).
  const activeData = props.active?.data();
  const activeKey = activeData ? `${props.active?.editorName}:${activeData.sequence.length}` : "";
  useEffect(() => {
    if (!activeData) {
      setBound(null);
      return;
    }
    let live = true;
    // Bindings are computed against the saved shared file, so they refresh when the document
    // changes or the file is saved or reloaded (its modification time changes).
    api.primerLibrary.bindings(activeData)
      .then((r) => live && setBound(r.primers))
      .catch(() => live && setBound(null));
    return () => { live = false; };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [activeKey, file.mtime]);

  const shown = useMemo(() => filterPrimers(primers, query), [primers, query]);

  const change = (list: LibraryPrimer[]) => {
    setPrimers(list);
    setDirty(true);
    setNotice(null);
  };

  const save = async () => {
    setBusy(true);
    setError(null);
    try {
      // The file may have changed under the user (a colleague saved it): warn rather than
      // overwrite blindly.
      const current = await api.primerLibrary.get();
      if (file.mtime && current.mtime && current.mtime !== file.mtime) {
        setNotice({ intent: "warning", text: "The library file changed on disk since you opened "
          + "it. Reload to see the other changes, or Save to keep yours." });
        setFile((f) => ({ ...f, mtime: current.mtime }));
        setBusy(false);
        return;
      }
      const saved = await api.primerLibrary.save(primers as unknown as object[]);
      setFile((f) => ({ ...f, path: saved.path, exists: true, mtime: saved.mtime }));
      setDirty(false);
      setNotice({ intent: "success", text: `Saved ${saved.count} primer${saved.count === 1 ? "" : "s"} to ${baseName(saved.path)}.` });
    } catch (e) {
      setError(messageOf(e));
    } finally {
      setBusy(false);
    }
  };

  const chooseFile = async () => {
    const d = (await whenDesktopReady(0)) ?? desktop();
    if (!d) {
      setNotice({ intent: "warning", text: "Choosing the library file needs the desktop app." });
      return;
    }
    const dir = file.path ? file.path.replace(/[\\/][^\\/]*$/, "") : "";
    const chosen = await d.save_dialog(baseName(file.path) || "lab-primers.csv", dir, "csv");
    if (!chosen) return;
    setBusy(true);
    try {
      const got = await api.primerLibrary.setPath(chosen);
      setFile(got);
      setPrimers(got.primers);
      setDirty(false);
      setNotice({ intent: "success", text: `The library is now ${baseName(got.path)}. Keep it in a shared folder so your lab shares it.` });
    } catch (e) {
      setError(messageOf(e));
    } finally {
      setBusy(false);
    }
  };

  const addPrimer = () => {
    const seq = cleanSequence(adding.sequence);
    if (!seq) {
      setNotice({ intent: "danger", text: "Enter a DNA sequence (A, C, G, T) for the primer." });
      return;
    }
    const existing = primers.find((p) => p.sequence === seq);
    if (existing) {
      setNotice({ intent: "warning", text: `That sequence is already in the library as “${existing.name}”.` });
      return;
    }
    change([...primers, newEntry(adding.name, seq, today, adding.notes.trim(), adding.location.trim())]);
    setAdding({ name: "", sequence: "", notes: "", location: "" });
  };

  const remove = (primer: LibraryPrimer) =>
    change(primers.filter((p) => p !== primer));

  const editField = (primer: LibraryPrimer, field: "notes" | "location", value: string) =>
    change(primers.map((p) => (p === primer ? { ...p, [field]: value } : p)));

  const runImport = async (text: string) => {
    setBusy(true);
    try {
      const parsed = await api.primerLibrary.parse(text);
      const incoming = parsed.primers.map((p) => ({ ...p, added: p.added || today }));
      const { merged, added, duplicates } = mergeImport(primers, incoming);
      change(merged);
      setImportText(null);
      const parts = [`Imported ${added} primer${added === 1 ? "" : "s"}`];
      if (duplicates) parts.push(`${duplicates} already in the library`);
      if (parsed.rejected.length) parts.push(rejectedText(parsed.rejected, "skipped"));
      setNotice({ intent: added ? "success" : "warning", text: parts.join("; ") + ". Save to write them to the shared file." });
    } catch (e) {
      setError(messageOf(e));
    } finally {
      setBusy(false);
    }
  };

  const importFromFile = async (f: File) => {
    try {
      await runImport(await f.text());
    } catch {
      setNotice({ intent: "danger", text: `“${f.name}” could not be read.` });
    }
  };

  const exportCsv = async () => {
    const text = toCsv(primers);
    const d = (await whenDesktopReady(0)) ?? desktop();
    if (d) {
      const chosen = await d.save_dialog("lab-primers.csv", "", "csv");
      if (!chosen) return;
      try {
        await api.saveBinary(chosen, new Blob([text], { type: "text/csv" }));
        setNotice({ intent: "success", text: `Exported ${primers.length} primers to ${baseName(chosen)}.` });
      } catch (e) {
        setError(messageOf(e));
      }
    } else {
      const url = URL.createObjectURL(new Blob([text], { type: "text/csv" }));
      const a = document.createElement("a");
      a.href = url;
      a.download = "lab-primers.csv";
      a.click();
      URL.revokeObjectURL(url);
    }
  };

  const copyOne = async (primer: LibraryPrimer) => {
    if (await copy(orderingText([primer]), primer.name)) {
      setNotice({ intent: "success", text: `Copied ${primer.name}: paste it into an order form.` });
    }
  };

  return (
    <Dialog isOpen onClose={props.onClose} title="Primer Library" style={{ width: 1180, maxWidth: "95vw" }}
            canOutsideClickClose={!dirty}>
      <div className={`${Classes.DIALOG_BODY} vt-primer-library`} data-testid="primer-library-dialog">
        <div className="vt-pl-head">
          <div>
            <p className="vt-pl-intro">
              The primers your lab already has. Keep the file in a shared folder so colleagues
              share it; Vettore reads it from there and writes it back.
            </p>
            <p className="vt-pl-file" data-testid="pl-file">
              Library file: <b>{file.path ? baseName(file.path) : "none yet"}</b>
              {file.path && <span className="vt-muted"> · {file.path}</span>}
              {!file.exists && file.path && <span className="vt-muted"> (new, written on first save)</span>}
            </p>
          </div>
          <div className="vt-pl-head-actions">
            <Button small data-testid="pl-choose-file" onClick={chooseFile}>Change file…</Button>
            <Button small data-testid="pl-reload" onClick={load} disabled={busy}
                    title="Re-read the shared file (e.g. after a colleague changed it)">Reload</Button>
          </div>
        </div>

        {error && <Callout intent="danger" data-testid="pl-error">{error}</Callout>}
        {notice && <Callout intent={notice.intent} data-testid="pl-notice">{notice.text}</Callout>}

        <div className="vt-pl-toolbar">
          <input className="bp3-input" data-testid="pl-search" placeholder="Search name or sequence…"
                 value={query} onChange={(e) => setQuery(e.target.value)} />
          <span className="vt-muted" data-testid="pl-count">
            {shown.length} of {primers.length} primer{primers.length === 1 ? "" : "s"}
          </span>
          <span className="vt-pl-spacer" />
          <Button small data-testid="pl-import" onClick={() => setImportText("")}>Import…</Button>
          <Button small data-testid="pl-export" disabled={!primers.length} onClick={exportCsv}>Export .csv</Button>
          <input ref={fileInput} type="file" accept=".csv,.tsv,.txt" style={{ display: "none" }}
                 data-testid="pl-file-input"
                 onChange={(e) => { const f = e.target.files?.[0]; if (f) importFromFile(f); e.target.value = ""; }} />
        </div>

        <div className="vt-pl-table-wrap">
          <table className="bp3-html-table bp3-html-table-condensed vt-table vt-pl-table" data-testid="pl-table">
            <thead>
              <tr>
                <th>Name</th><th>Sequence 5′→3′</th><th>nt</th><th>Tm</th>
                <th>Notes</th><th>Location</th><th>Added</th><th></th>
              </tr>
            </thead>
            <tbody>
              {shown.length === 0 && (
                <tr><td colSpan={8} className="vt-muted" data-testid="pl-empty">
                  {primers.length ? "No primer matches the search." : "No primers yet. Add one below or Import a spreadsheet."}
                </td></tr>
              )}
              {shown.map((p, i) => (
                <tr key={`${p.sequence}-${i}`} data-testid="pl-row">
                  <td>{p.name}</td>
                  <td className="vt-pl-seq">{p.sequence}</td>
                  <td>{p.length}</td>
                  <td>{p.tm === null ? "—" : `${p.tm.toFixed(1)} °C`}</td>
                  <td>
                    <input className="bp3-input bp3-small vt-pl-edit" aria-label={`Notes for ${p.name}`}
                           value={p.notes} onChange={(e) => editField(p, "notes", e.target.value)} />
                  </td>
                  <td>
                    <input className="bp3-input bp3-small vt-pl-edit" aria-label={`Location of ${p.name}`}
                           placeholder="box 3 A4" value={p.location}
                           onChange={(e) => editField(p, "location", e.target.value)} />
                  </td>
                  <td className="vt-muted">{p.added || "—"}</td>
                  <td className="vt-pl-rowactions">
                    <Button small minimal icon="duplicate" aria-label={`Copy ${p.name}`}
                            data-testid="pl-copy" onClick={() => copyOne(p)} />
                    <Button small minimal icon="trash" aria-label={`Remove ${p.name}`}
                            data-testid="pl-remove" onClick={() => remove(p)} />
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>

        <div className="vt-pl-add" data-testid="pl-add">
          <span className="vt-field-label">Add a primer</span>
          <div className="vt-pl-add-row">
            <input className="bp3-input bp3-small" data-testid="pl-add-name" placeholder="name"
                   value={adding.name} onChange={(e) => setAdding({ ...adding, name: e.target.value })} />
            <input className="bp3-input bp3-small vt-pl-add-seq" data-testid="pl-add-seq" placeholder="sequence 5′→3′"
                   value={adding.sequence} onChange={(e) => setAdding({ ...adding, sequence: e.target.value })} />
            <input className="bp3-input bp3-small" data-testid="pl-add-notes" placeholder="notes"
                   value={adding.notes} onChange={(e) => setAdding({ ...adding, notes: e.target.value })} />
            <input className="bp3-input bp3-small" data-testid="pl-add-location" placeholder="box 3 A4"
                   value={adding.location} onChange={(e) => setAdding({ ...adding, location: e.target.value })} />
            <Button small intent="primary" data-testid="pl-add-button" onClick={addPrimer}>Add</Button>
          </div>
        </div>

        {props.active && (
          <div className="vt-pl-bound" data-testid="pl-bound">
            <span className="vt-field-label">On “{props.active.name}”</span>
            {bound === null && <span className="vt-muted">Looking for library primers on this sequence…</span>}
            {bound && bound.length === 0 && (
              <span className="vt-muted" data-testid="pl-bound-none">No library primer binds this sequence.</span>
            )}
            {bound && bound.length > 0 && (
              <ul className="vt-pl-bound-list">
                {bound.map((p) => (
                  <li key={p.sequence}>
                    <b>{p.name}</b>{" "}
                    {p.bindings.map((b, i) => (
                      <span key={i} className="vt-muted">{i > 0 ? "; " : ""}{bindingText(b)}</span>
                    ))}
                  </li>
                ))}
              </ul>
            )}
          </div>
        )}
      </div>

      {importText !== null && (
        <ImportPanel text={importText} onText={setImportText} onFile={() => fileInput.current?.click()}
                     onImport={() => runImport(importText)} onCancel={() => setImportText(null)} busy={busy} />
      )}
      {fallbackDialog}

      <div className={Classes.DIALOG_FOOTER}>
        <div className="vt-pl-footer-note vt-muted">
          {dirty ? "Unsaved changes — Save writes them to the shared file." : "Up to date with the shared file."}
        </div>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button data-testid="pl-close" onClick={props.onClose}>Close</Button>
          <Button intent="primary" data-testid="pl-save" disabled={!dirty || busy} onClick={save}>Save to File</Button>
        </div>
      </div>
    </Dialog>
  );
}

function ImportPanel(props: {
  text: string; busy: boolean;
  onText: (t: string) => void; onFile: () => void; onImport: () => void; onCancel: () => void;
}) {
  return (
    <div className={Classes.DIALOG_BODY + " vt-pl-import"} data-testid="pl-import-panel">
      <Field label="Import primers"
             hint="Paste two or more columns copied from Excel (name, sequence, optional notes and location), or open a .csv/.tsv. Rows that are not DNA are listed and left out.">
        <textarea className="bp3-input bp3-fill vt-pl-import-text" rows={6} data-testid="pl-import-text"
                  placeholder={"KanR_F\tGTAAAACGACGGCCAGT\nKanR_R\tCAGGAAACAGCTATGAC"}
                  value={props.text} onChange={(e) => props.onText(e.target.value)} />
      </Field>
      <div className="vt-pl-import-actions">
        <Button small data-testid="pl-import-file" onClick={props.onFile}>Open a .csv/.tsv file…</Button>
        <span className="vt-pl-spacer" />
        <Button small data-testid="pl-import-cancel" onClick={props.onCancel}>Cancel</Button>
        <Button small intent="primary" data-testid="pl-import-run" disabled={props.busy || !props.text.trim()}
                onClick={props.onImport}>Import pasted rows</Button>
      </div>
    </div>
  );
}

function rejectedText(rejected: RejectedRow[], where: string): string {
  const names = rejected.slice(0, 3).map((r) => `“${r.name || r.sequence}”`).join(", ");
  const more = rejected.length > 3 ? ` and ${rejected.length - 3} more` : "";
  return `${rejected.length} row${rejected.length === 1 ? "" : "s"} ${where} (not DNA): ${names}${more}`;
}

export const tools: Tool[] = [
  { id: "primer-library", menu: "Primers", label: "Primer Library", Component: PrimerLibraryDialog },
];
