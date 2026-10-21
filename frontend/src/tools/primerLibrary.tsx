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
  heldText,
  rejectedText,
  saveChanges,
  today as todayDate,
} from "./primerLibraryText";
import type { Tool, ToolProps } from "./types";
import "./primerLibrary.css";

// Blueprint 3 typings predate React 18 (no implicit `children`).
const { Dialog, Callout } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

const EMPTY: LibraryFile = { path: "", exists: false, mtime: null, delimiter: ";", primers: [], rejected: [] };
const HELD_TIP = "This row of the file runs over several lines: change it in Excel.";
type Notice = { intent: "success" | "warning" | "danger"; text: string } | null;

function PrimerLibraryDialog(props: ToolProps) {
  const { copy, fallbackDialog } = useCopyText();
  const [file, setFile] = useState<LibraryFile>(EMPTY);
  const [primers, setPrimers] = useState<LibraryPrimer[]>([]);
  const [dirty, setDirty] = useState(false);
  // An action that would drop unsaved changes, waiting for the user to confirm it.
  const [pendingDiscard, setPendingDiscard] = useState<{ run: () => void; what: string } | null>(null);
  const [query, setQuery] = useState("");
  const [notice, setNotice] = useState<Notice>(null);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [adding, setAdding] = useState({ name: "", sequence: "", notes: "", location: "" });
  const [importText, setImportText] = useState<string | null>(null);
  const [bound, setBound] = useState<BoundPrimer[] | null>(null);
  const today = todayDate();
  const fileInput = useRef<HTMLInputElement>(null);

  const load = useCallback(async () => {
    setBusy(true);
    setError(null);
    try {
      const got = await api.primerLibrary.get();
      setFile(got);
      setPrimers(got.primers);
      setDirty(false);
      const warnings = [
        got.rejected.length ? rejectedText(got.rejected, "file") : "",
        got.held?.length ? heldText(got.held) : "",
      ].filter(Boolean);
      if (warnings.length) setNotice({ intent: "warning", text: warnings.join(" ") });
    } catch (e) {
      setError(messageOf(e));
    } finally {
      setBusy(false);
    }
  }, []);

  useEffect(() => {
    load();
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
    // The dialog's own list, saved or not, so a primer just added or imported shows at once.
    const t = setTimeout(() => {
      api.primerLibrary.bindings(activeData, primers as unknown as object[])
        .then((r) => live && setBound(r.primers))
        .catch(() => live && setBound(null));
    }, 250);
    return () => { live = false; clearTimeout(t); };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [activeKey, primers]);

  const shown = useMemo(() => filterPrimers(primers, query), [primers, query]);

  const change = (list: LibraryPrimer[]) => {
    setPrimers(list);
    setDirty(true);
    setNotice(null);
  };

  /** Save, merged with the file as it is now; true when it was written. */
  const save = async (): Promise<boolean> => {
    setBusy(true);
    setError(null);
    try {
      // Only this window's changes, applied by the engine to the file as it is now: what a
      // colleague added, edited or deleted since it was loaded stays as they left it.
      const changes = saveChanges(file.primers, primers, file.rejected);
      const saved = await api.primerLibrary.save(changes);
      const mine = new Set([...primers, ...file.primers].map((p) => p.sequence));
      const after = new Set(saved.primers.map((p) => p.sequence));
      const theirs = saved.primers.filter((p) => !mine.has(p.sequence)).length;
      const touched = new Set([...changes.edited, ...changes.removed].map((c) => c.sequence));
      const dropped = file.primers.filter((p) => !touched.has(p.sequence) && !after.has(p.sequence)).length;
      setFile(saved);
      setPrimers(saved.primers);
      setDirty(false);
      const plural = (n: number, one: string, many: string) => `${n} ${n === 1 ? one : many}`;
      const extra = [
        saved.skipped && `${plural(saved.skipped, "primer you added was", "primers you added were")} already in the file (someone added the same oligo meanwhile): the file's entry is kept.`,
        saved.gone && `${plural(saved.gone, "change was", "changes were")} not saved: someone had changed or removed that row of the file meanwhile.`,
        saved.held && `${plural(saved.held, "change was", "changes were")} not saved: that row of the file runs over several lines, so Vettore leaves it as it is (change it in Excel).`,
        theirs && `${plural(theirs, "primer", "primers")} someone else added meanwhile ${theirs === 1 ? "is" : "are"} now listed.`,
        dropped > 0 && `${plural(dropped, "primer", "primers")} someone else removed meanwhile ${dropped === 1 ? "is" : "are"} no longer listed.`,
      ].filter(Boolean).join(" ");
      setNotice({ intent: saved.gone || saved.skipped ? "warning" : "success",
                  text: `Saved to ${baseName(saved.path)}: ${saved.added} added, ${saved.changed} changed, ${saved.removed} removed.${extra ? " " + extra : ""}` });
      return true;
    } catch (e) {
      setError(messageOf(e));
      return false;
    } finally {
      setBusy(false);
    }
  };

  /** Run an action that re-reads the file or closes the window, asking first when there are
   * unsaved changes it would drop. */
  const guard = (what: string, run: () => void) => {
    if (busy) return;                     // never while a save is being written
    if (dirty) setPendingDiscard({ what, run });
    else run();
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

  const remove = (primer: LibraryPrimer) => {
    change(primers.filter((p) => p !== primer));
    // A second row of the file with the same oligo (a second tube) is not a primer of the list,
    // so it stays in the file and becomes the library's entry for that oligo after saving.
    const twin = file.rejected.find((r) => cleanSequence(r.sequence) === primer.sequence);
    if (twin) {
      setNotice({ intent: "warning", text: `Another row of the file has the same oligo (“${twin.name || "unnamed"}”). `
        + "It stays in the file, and after saving it is listed in place of the one you removed." });
    }
  };

  const editField = (primer: LibraryPrimer, field: "notes" | "location", value: string) =>
    change(primers.map((p) => (p === primer ? { ...p, [field]: value } : p)));

  const runImport = async (text: string) => {
    setBusy(true);
    try {
      const parsed = await api.primerLibrary.parse(text);
      const incoming = parsed.primers.map((p) => ({ ...p, added: p.added || today, fresh: true }));
      const { merged, added, duplicates } = mergeImport(primers, incoming);
      change(merged);
      setImportText(null);
      const parts = [`Imported ${added} primer${added === 1 ? "" : "s"}`];
      if (duplicates) parts.push(`${duplicates} already in the library`);
      setNotice({ intent: added && !parsed.rejected.length ? "success" : "warning",
                  text: parts.join("; ") + ". Save to write them to the shared file."
                    + (parsed.rejected.length ? " " + rejectedText(parsed.rejected, "paste") : "") });
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
    // The library file's separator (";" for Excel set to Italian) and a byte-order mark, so
    // Excel opens the export in columns with accents and µ intact.
    const text = "\ufeff" + toCsv(primers, file.delimiter);
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
    <Dialog isOpen onClose={() => guard("close the window", props.onClose)} title="Primer Library"
            style={{ width: 1180, maxWidth: "95vw" }} canOutsideClickClose={!dirty}>
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
            <Button small data-testid="pl-choose-file" disabled={busy}
                    onClick={() => guard("choose another file", chooseFile)}>Change file…</Button>
            <Button small data-testid="pl-reload" onClick={() => guard("reload the file", load)} disabled={busy}
                    title="Re-read the shared file (e.g. after a colleague changed it)">Reload</Button>
          </div>
        </div>

        {pendingDiscard && (
          <Callout intent="warning" data-testid="pl-discard">
            You have changes that are not saved to the shared file. {pendingDiscard.what[0].toUpperCase()}
            {pendingDiscard.what.slice(1)} without saving them?
            <div className="vt-pl-discard-actions">
              <Button small intent="primary" data-testid="pl-discard-save"
                      onClick={async () => { const next = pendingDiscard; setPendingDiscard(null); if (await save()) next.run(); }}>
                Save, then {pendingDiscard.what}
              </Button>
              <Button small intent="danger" data-testid="pl-discard-yes"
                      onClick={() => { const next = pendingDiscard; setPendingDiscard(null); setDirty(false); next.run(); }}>
                Discard my changes
              </Button>
              <Button small data-testid="pl-discard-no" onClick={() => setPendingDiscard(null)}>Keep editing</Button>
            </div>
          </Callout>
        )}
        {error && <Callout intent="danger" data-testid="pl-error">{error}</Callout>}
        {notice && <Callout intent={notice.intent} data-testid="pl-notice">{notice.text}</Callout>}

        <div className="vt-pl-toolbar">
          <input className="bp3-input" data-testid="pl-search" placeholder="Search name or sequence…"
                 value={query} onChange={(e) => setQuery(e.target.value)} />
          <span className="vt-muted" data-testid="pl-count">
            {shown.length} of {primers.length} primer{primers.length === 1 ? "" : "s"}
          </span>
          <span className="vt-pl-spacer" />
          <Button small data-testid="pl-import" disabled={busy} onClick={() => setImportText("")}>Import…</Button>
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
                  <td title={p.searched === false ? "Shorter than 15 nt: not looked for on sequences" : undefined}>
                    {p.length}{p.searched === false && <span className="vt-muted"> · too short to search</span>}
                    {p.degenerate && <span className="vt-muted"> · degenerate</span>}
                  </td>
                  <td>{p.tm === null ? "—" : `${p.tm.toFixed(1)} °C`}</td>
                  <td>
                    {/* A text area, not an input: a note written on several lines in Excel keeps
                        its line breaks when it is edited here. */}
                    <textarea className="bp3-input bp3-small vt-pl-edit" aria-label={`Notes for ${p.name}`}
                              rows={Math.min(3, p.notes.split("\n").length)}
                              value={p.notes} disabled={busy || p.held} title={p.held ? HELD_TIP : undefined} onChange={(e) => editField(p, "notes", e.target.value)} />
                  </td>
                  <td>
                    <textarea className="bp3-input bp3-small vt-pl-edit" aria-label={`Location of ${p.name}`}
                              rows={Math.min(3, p.location.split("\n").length)}
                              placeholder="box 3 A4" value={p.location} disabled={busy || p.held}
                              title={p.held ? HELD_TIP : undefined}
                              onChange={(e) => editField(p, "location", e.target.value)} />
                  </td>
                  <td className="vt-muted">{p.added || "—"}</td>
                  <td className="vt-pl-rowactions">
                    <Button small minimal icon="duplicate" aria-label={`Copy ${p.name}`}
                            data-testid="pl-copy" onClick={() => copyOne(p)} />
                    <Button small minimal icon="trash" aria-label={`Remove ${p.name}`}
                            data-testid="pl-remove" disabled={busy || p.held} onClick={() => remove(p)} />
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
            <Button small intent="primary" data-testid="pl-add-button" disabled={busy} onClick={addPrimer}>Add</Button>
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
                    {p.bindings.length === 0 && (
                      <span className="vt-muted" data-testid="pl-bound-many">
                        binds {p.sites} places{p.repeat ? " (a repeat)" : ""}: not drawn on the map
                      </span>
                    )}
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
          <Button data-testid="pl-close" disabled={busy} onClick={() => guard("close the window", props.onClose)}>Close</Button>
          <Button intent="primary" data-testid="pl-save" disabled={!dirty || busy} onClick={() => save()}>Save to File</Button>
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
        <Button small data-testid="pl-import-file" disabled={props.busy} onClick={props.onFile}>Open a .csv/.tsv file…</Button>
        <span className="vt-pl-spacer" />
        <Button small data-testid="pl-import-cancel" onClick={props.onCancel}>Cancel</Button>
        <Button small intent="primary" data-testid="pl-import-run" disabled={props.busy || !props.text.trim()}
                onClick={props.onImport}>Import pasted rows</Button>
      </div>
    </div>
  );
}


export const tools: Tool[] = [
  { id: "primer-library", menu: "Primers", label: "Primer Library", Component: PrimerLibraryDialog },
];
