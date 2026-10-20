import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useSelector } from "react-redux";
import { Editor } from "@teselagen/ove";
import { api, whenDesktopReady, type OpenResult, type SequenceData } from "./api";
import {
  baseName,
  canSaveInPlace,
  createDocument,
  editorState,
  isDirty,
  markSaved,
  now,
  removeDocument,
  stem,
  type Doc,
} from "./documents";
import { buildMenu } from "./menu";
import NewSequenceDialog, { type NewSequence } from "./NewSequenceDialog";
import CloningDialog from "./cloning/CloningDialog";
import GelDialog from "./GelDialog";
import SangerDialog from "./SangerDialog";
import DetectFeaturesDialog from "./DetectFeaturesDialog";
import { actionById, type CloningAction } from "./cloning/registry";
import type { CloningResult, OpenDoc, Product } from "./cloning/types";
import HistoryView from "./HistoryView";
import EditorSlot from "./EditorSlot";
import EnzymesView from "./EnzymesView";
import StatusLine from "./StatusLine";
import { exportMap } from "./mapExport";
import { applyEnzymeSet, ENZYME_SETS, type EnzymeSetId } from "./enzymeSets";
import { ENZYMES } from "./enzymeLogic";
import { showView, VIEWS, type ViewId } from "./views";

const DOCTABS_HEIGHT = 30;
const VIEWBAR_HEIGHT = 32;
const STATUS_HEIGHT = 26;

function useWindowHeight() {
  const [height, setHeight] = useState(window.innerHeight);
  useEffect(() => {
    const onResize = () => setHeight(window.innerHeight);
    window.addEventListener("resize", onResize);
    return () => window.removeEventListener("resize", onResize);
  }, []);
  return height;
}

function download(filename: string, content: string | Blob, type = "text/plain") {
  const blob = content instanceof Blob ? content : new Blob([content], { type });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

export default function App() {
  const [docs, setDocs] = useState<Doc[]>([]);
  const [active, setActive] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [newOpen, setNewOpen] = useState(false);
  const [cloning, setCloning] = useState<CloningAction | null>(null);
  const [gelOpen, setGelOpen] = useState(false);
  const [sangerOpen, setSangerOpen] = useState(false);
  const [detectOpen, setDetectOpen] = useState(false);
  const [recent, setRecent] = useState<string[]>([]);
  const fileInput = useRef<HTMLInputElement>(null);
  const height = useWindowHeight();
  const docsRef = useRef(docs);
  docsRef.current = docs;

  const activeDoc = docs.find((d) => d.editorName === active) ?? null;
  const dirtyByEditor = useSelector((state: any) => {
    const out: Record<string, boolean> = {};
    for (const d of docs) out[d.editorName] = isDirty(state.VectorEditor[d.editorName]);
    return out;
  }, (a, b) => JSON.stringify(a) === JSON.stringify(b));
  const circular = useSelector((state: any) =>
    active ? !!state.VectorEditor[active]?.sequenceData?.circular : false,
  );
  const activeSequence = useSelector((state: any) =>
    active ? (state.VectorEditor[active]?.sequenceData?.sequence as string | undefined) : undefined,
  );

  const unsavedCount = Object.values(dirtyByEditor).filter(Boolean).length;
  useEffect(() => {
    api.unsaved(unsavedCount).catch(() => undefined);
  }, [unsavedCount]);

  const updateDoc = useCallback((editorName: string, patch: Partial<Doc>) => {
    setDocs((prev) => prev.map((d) => (d.editorName === editorName ? { ...d, ...patch } : d)));
  }, []);

  const logHistory = useCallback((editorName: string, action: string, detail?: string) => {
    setDocs((prev) =>
      prev.map((d) =>
        d.editorName === editorName
          ? { ...d, history: [...d.history, { time: now(), action, detail }] }
          : d,
      ),
    );
  }, []);

  const refreshRecent = useCallback(() => {
    api.recent().then((r) => setRecent(r.paths)).catch(() => undefined);
  }, []);

  const addResult = useCallback((result: OpenResult, origin: string) => {
    const created = result.documents.map((sequenceData) =>
      createDocument(sequenceData, result.path, {
        time: now(),
        action: "Opened",
        detail: result.path ?? origin,
      }),
    );
    setDocs((prev) => [...prev, ...created]);
    setActive(created[created.length - 1]?.editorName ?? null);
  }, []);

  const report = useCallback((e: unknown) => setError(String((e as Error).message ?? e)), []);

  // ---- open ------------------------------------------------------------------------------
  const openPath = useCallback(
    async (path: string) => {
      const existing = docsRef.current.find((d) => d.path === path);
      if (existing) {
        setActive(existing.editorName);
        return;
      }
      try {
        addResult(await api.openPath(path), path);
        refreshRecent();
      } catch (e) {
        report(e);
      }
    },
    [addResult, refreshRecent, report],
  );

  const openFiles = useCallback(async () => {
    setError(null);
    const desktop = await whenDesktopReady(0);
    if (!desktop) {
      fileInput.current?.click();
      return;
    }
    for (const path of await desktop.open_dialog()) await openPath(path);
  }, [openPath]);

  const openBrowserFiles = useCallback(
    async (files: FileList | File[] | null) => {
      for (const file of Array.from(files ?? [])) {
        try {
          addResult(await api.openContent(file), file.name);
        } catch (e) {
          report(e);
        }
      }
      if (fileInput.current) fileInput.current.value = "";
    },
    [addResult, report],
  );

  useEffect(() => {
    refreshRecent();
    api
      .startupFiles()
      .then(({ results, errors }) => {
        results.forEach((r) => addResult(r, r.path ?? ""));
        if (errors.length) setError(errors.join("; "));
      })
      .catch(report);
  }, [addResult, refreshRecent, report]);

  // ---- new -------------------------------------------------------------------------------
  const createNew = useCallback((seq: NewSequence) => {
    setNewOpen(false);
    const sequenceData: SequenceData = { ...seq, features: [], primers: [] };
    const doc = createDocument(sequenceData, null, { time: now(), action: "Created", detail: `${seq.sequence.length} bp` });
    setDocs((prev) => [...prev, doc]);
    setActive(doc.editorName);
  }, []);

  // ---- cloning -------------------------------------------------------------------------
  // Memoised: the cloning forms key their effects on these objects' identity.
  const openDocs: OpenDoc[] = useMemo(() => docs.map((d) => ({
    editorName: d.editorName,
    name: d.name,
    data: () => editorState(d.editorName)?.sequenceData as SequenceData,
    selection: () => {
      const layer = editorState(d.editorName)?.selectionLayer;
      return layer && layer.start >= 0 && layer.end >= 0
        ? { start: layer.start, end: layer.end }
        : null;
    },
  })), [docs]);

  const createProduct = useCallback((product: Product, result: CloningResult) => {
    setCloning(null);
    const oligos = (result.oligos.length ? result.oligos : product.oligos)
      .map((o) => `${o.name} ${o.sequence}`)
      .join(" · ");
    const doc = createDocument(product.doc as SequenceData, null, {
      time: now(),
      action: product.operation,
      detail: `${product.details}. From: ${product.parents.join(", ")}.`
        + (oligos ? ` Primers: ${oligos}` : ""),
    });
    setDocs((prev) => [...prev, doc]);
    setActive(doc.editorName);
  }, []);

  // ---- save ------------------------------------------------------------------------------
  const saveDoc = useCallback(
    async (doc: Doc, saveAs: boolean): Promise<boolean> => {
      const sequenceData = editorState(doc.editorName)?.sequenceData as SequenceData;
      if (!sequenceData) return false;
      try {
        const desktop = await whenDesktopReady(0);
        let path = !saveAs && canSaveInPlace(doc) ? doc.path : null;
        if (!path && desktop) {
          path = await desktop.save_dialog(`${stem(doc.path ?? doc.name)}.gb`);
          if (!path) return false;
        }
        if (path) {
          const saved = await api.savePath(path, sequenceData);
          updateDoc(doc.editorName, { path: saved.path, name: stem(saved.path) });
          logHistory(doc.editorName, "Saved", saved.path);
          refreshRecent();
        } else {
          const { content } = await api.exportText(sequenceData, "genbank");
          download(`${stem(doc.name)}.gb`, content);
          logHistory(doc.editorName, "Downloaded", `${stem(doc.name)}.gb`);
        }
        markSaved(doc.editorName);
        return true;
      } catch (e) {
        report(e);
        return false;
      }
    },
    [logHistory, refreshRecent, report, updateDoc],
  );

  const exportAs = useCallback(
    async (doc: Doc, format: "genbank" | "fasta") => {
      const sequenceData = editorState(doc.editorName)?.sequenceData as SequenceData;
      try {
        const ext = format === "fasta" ? "fasta" : "gb";
        const desktop = await whenDesktopReady(0);
        if (desktop) {
          const path = await desktop.save_dialog(`${stem(doc.name)}.${ext}`);
          if (path) await api.savePath(path, sequenceData);
        } else {
          const { content } = await api.exportText(sequenceData, format);
          download(`${stem(doc.name)}.${ext}`, content);
        }
      } catch (e) {
        report(e);
      }
    },
    [report],
  );

  const exportImage = useCallback(
    async (doc: Doc, format: "svg" | "png") => {
      try {
        const blob = await exportMap(format);
        const desktop = await whenDesktopReady(0);
        const filename = `${stem(doc.name)}-map.${format}`;
        if (desktop) {
          const path = await desktop.save_dialog(filename);
          if (path) await api.saveBinary(path, blob);
        } else {
          download(filename, blob);
        }
      } catch (e) {
        report(e);
      }
    },
    [report],
  );

  // ---- close -----------------------------------------------------------------------------
  const closeDoc = useCallback(
    async (doc: Doc) => {
      if (isDirty(editorState(doc.editorName))) {
        const save = window.confirm(`Save changes to "${doc.name}" before closing?\n\nOK = Save, Cancel = more options`);
        if (save) {
          if (!(await saveDoc(doc, false))) return;
        } else if (!window.confirm(`Close "${doc.name}" and discard the changes?`)) {
          return;
        }
      }
      removeDocument(doc.editorName);
      setDocs((prev) => {
        const idx = prev.findIndex((d) => d.editorName === doc.editorName);
        const next = prev.filter((d) => d.editorName !== doc.editorName);
        setActive((cur) =>
          cur === doc.editorName ? (next[Math.min(idx, next.length - 1)]?.editorName ?? null) : cur,
        );
        return next;
      });
    },
    [saveDoc],
  );

  // ---- views -----------------------------------------------------------------------------
  const selectView = useCallback(
    (view: ViewId) => {
      if (!activeDoc) return;
      updateDoc(activeDoc.editorName, { view });
      showView(activeDoc.editorName, view, circular);
    },
    [activeDoc, circular, updateDoc],
  );

  const useEnzymeSet = useCallback(
    (set: EnzymeSetId) => {
      if (!activeDoc) return;
      updateDoc(activeDoc.editorName, { enzymeSet: set });
    },
    [activeDoc, updateDoc],
  );

  useEffect(() => {
    // Recompute the enzyme set when the sequence, topology or chosen set changes.
    if (!activeDoc || activeSequence === undefined) return;
    const t = setTimeout(() => applyEnzymeSet(activeDoc.editorName, activeDoc.enzymeSet), 200);
    return () => clearTimeout(t);
  }, [activeDoc?.editorName, activeDoc?.enzymeSet, activeSequence, circular]); // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    // Keep the map tab in sync with the topology (circular map vs linear map).
    if (activeDoc && activeDoc.view !== "history") showView(activeDoc.editorName, activeDoc.view, circular);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [active, circular]);

  // ---- keyboard shortcuts not handled by OVE --------------------------------------------
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      const mod = e.ctrlKey || e.metaKey;
      if (!mod) return;
      const key = e.key.toLowerCase();
      if (key === "o" && !e.shiftKey) {
        e.preventDefault();
        openFiles();
      } else if (key === "n" && !e.shiftKey) {
        e.preventDefault();
        setNewOpen(true);
      } else if (key === "w" && activeDoc) {
        e.preventDefault();
        closeDoc(activeDoc);
      } else if (!e.altKey && !e.shiftKey && activeDoc) {
        const view = VIEWS.find((v) => v.hotkey === key);
        if (view) {
          e.preventDefault();
          selectView(view.id);
        }
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [activeDoc, closeDoc, openFiles, selectView]);

  // ---- drag and drop ---------------------------------------------------------------------
  const onDrop = useCallback(
    (e: React.DragEvent) => {
      if (!e.dataTransfer.files.length) return;
      e.preventDefault();
      e.stopPropagation();
      openBrowserFiles(Array.from(e.dataTransfer.files));
    },
    [openBrowserFiles],
  );

  // ---- menu ------------------------------------------------------------------------------
  const handlers = useRef({
    newFile: () => setNewOpen(true),
    open: openFiles,
    openRecent: openPath,
    close: () => activeDoc && closeDoc(activeDoc),
    save: () => activeDoc && saveDoc(activeDoc, false),
    saveAs: () => activeDoc && saveDoc(activeDoc, true),
    exportGenbank: () => activeDoc && exportAs(activeDoc, "genbank"),
    exportFasta: () => activeDoc && exportAs(activeDoc, "fasta"),
    exportSvg: () => activeDoc && exportImage(activeDoc, "svg"),
    exportPng: () => activeDoc && exportImage(activeDoc, "png"),
    view: selectView,
    useEnzymeSet,
    cloning: (id: string) => setCloning(actionById(id)),
    gel: () => setGelOpen(true),
    sanger: () => setSangerOpen(true),
    detectFeatures: () => setDetectOpen(true),
    recent: [] as string[],
  });
  Object.assign(handlers.current, {
    open: openFiles,
    openRecent: openPath,
    close: () => activeDoc && closeDoc(activeDoc),
    save: () => activeDoc && saveDoc(activeDoc, false),
    saveAs: () => activeDoc && saveDoc(activeDoc, true),
    exportGenbank: () => activeDoc && exportAs(activeDoc, "genbank"),
    exportFasta: () => activeDoc && exportAs(activeDoc, "fasta"),
    exportSvg: () => activeDoc && exportImage(activeDoc, "svg"),
    exportPng: () => activeDoc && exportImage(activeDoc, "png"),
    view: selectView,
    useEnzymeSet,
    cloning: (id: string) => setCloning(actionById(id)),
    gel: () => setGelOpen(true),
    sanger: () => setSangerOpen(true),
    detectFeatures: () => setDetectOpen(true),
    recent,
  });
  // menuFilter must be stable (OVE memoizes on it); handlers are read through the ref.
  const menuFilter = useMemo(() => buildMenu(handlers), [recent]); // eslint-disable-line react-hooks/exhaustive-deps

  const onSave = useCallback(
    (opts: { isSaveAs?: boolean }) => (activeDoc ? saveDoc(activeDoc, !!opts?.isSaveAs).then(() => undefined) : undefined),
    [activeDoc, saveDoc],
  );

  const editorHeight = height - STATUS_HEIGHT;

  const docTabs = (
      <div className="vt-doctabs" style={{ height: DOCTABS_HEIGHT }}>
        {docs.map((d) => (
          <div
            key={d.editorName}
            className={`vt-doctab${d.editorName === active ? " active" : ""}`}
            onClick={() => setActive(d.editorName)}
            onAuxClick={(e) => e.button === 1 && closeDoc(d)}
            title={d.path ?? d.name}
            data-testid="doctab"
          >
            <span>{dirtyByEditor[d.editorName] ? "• " : ""}{d.name}</span>
            <button className="vt-doctab-close" title="Close (Ctrl+W)" onClick={(e) => { e.stopPropagation(); closeDoc(d); }}>×</button>
          </div>
        ))}
        {!docs.length && <span className="vt-brand">Vettore</span>}
        {error && (
          <span className="vt-error" role="alert" onClick={() => setError(null)} title="Click to dismiss">
            {error}
          </span>
        )}
      </div>
  );

  const viewBar = (
      <div className="vt-viewbar" style={{ height: VIEWBAR_HEIGHT }}>
        <div className="vt-viewtabs" role="tablist">
          {VIEWS.map((v) => (
            <button
              key={v.id}
              role="tab"
              data-testid={`view-${v.id}`}
              aria-selected={activeDoc?.view === v.id}
              className={`vt-viewtab${activeDoc?.view === v.id ? " active" : ""}`}
              disabled={!activeDoc}
              title={`${v.label} (Ctrl+${v.hotkey})`}
              onClick={() => selectView(v.id)}
            >
              {v.label}
            </button>
          ))}
        </div>
        {activeDoc && (
          <label className="vt-enzyme-set">
            Enzymes:&nbsp;
            <select
              data-testid="enzyme-set"
              value={activeDoc.enzymeSet}
              onChange={(e) => useEnzymeSet(e.target.value as EnzymeSetId)}
            >
              {ENZYME_SETS.map((set) => (
                <option key={set.id} value={set.id}>{set.label}</option>
              ))}
            </select>
          </label>
        )}
      </div>
  );

  return (
    <div className="vt-shell" onDragOver={(e) => e.preventDefault()} onDropCapture={onDrop}>
      <input
        ref={fileInput}
        type="file"
        multiple
        hidden
        data-testid="file-input"
        accept=".dna,.gb,.gbk,.gbff,.genbank,.ape,.fa,.fasta,.fna,.fas,.ab1,.abi"
        onChange={(e) => openBrowserFiles(e.target.files)}
      />
      <div className="vt-editor" style={{ height: editorHeight }}>
        {activeDoc ? (
          <>
            <EditorSlot editorName={activeDoc.editorName}>
              {docTabs}
              {viewBar}
            </EditorSlot>
            {activeDoc.view === "history" && (
              <EditorSlot editorName={activeDoc.editorName} placement="overlay">
                <HistoryView doc={activeDoc} />
              </EditorSlot>
            )}
            {activeDoc.view === "enzymes" && (
              <EditorSlot editorName={activeDoc.editorName} placement="overlay">
                <EnzymesView editorName={activeDoc.editorName} set={activeDoc.enzymeSet} />
              </EditorSlot>
            )}
            <Editor
              key={activeDoc.editorName}
              editorName={activeDoc.editorName}
              height={editorHeight}
              showMenuBar
              hideSingleImport
              allowPanelTabDraggable={false}
              menuFilter={menuFilter}
              onSave={onSave}
              onSaveAs={(opts: object) => onSave({ ...opts, isSaveAs: true })}
              alwaysAllowSave
              hideStatusBar
              additionalEnzymes={ENZYMES}
            />
          </>
        ) : (
          <>
          {docTabs}
          <Welcome recent={recent} onOpen={openFiles} onNew={() => setNewOpen(true)} onRecent={openPath} />
          </>
        )}
      </div>
      <div className="vt-status-bar" style={{ height: STATUS_HEIGHT }}>
        {activeDoc && <StatusLine editorName={activeDoc.editorName} />}
        {activeDoc && (
          <span className="vt-status-set" title="Enzyme set shown on the map">
            {ENZYME_SETS.find((e) => e.id === activeDoc.enzymeSet)?.label}
          </span>
        )}
      </div>
      {/* Remounted on every open so the form never keeps the previous file's name or topology. */}
      {newOpen && <NewSequenceDialog isOpen onClose={() => setNewOpen(false)} onCreate={createNew} />}
      {detectOpen && activeDoc && (
        <DetectFeaturesDialog
          doc={openDocs.find((d) => d.editorName === activeDoc.editorName)!}
          onClose={() => {
            setDetectOpen(false);
            logHistory(activeDoc.editorName, "Detected common features");
          }}
        />
      )}
      {sangerOpen && (
        <SangerDialog active={openDocs.find((d) => d.editorName === active) ?? null}
                      onClose={() => setSangerOpen(false)} />
      )}
      {gelOpen && (
        <GelDialog docs={openDocs} active={openDocs.find((d) => d.editorName === active) ?? null}
                   onClose={() => setGelOpen(false)} />
      )}
      {cloning && (
        <CloningDialog
          action={cloning}
          docs={openDocs}
          active={openDocs.find((d) => d.editorName === active) ?? null}
          onClose={() => setCloning(null)}
          onCreate={createProduct}
        />
      )}
    </div>
  );
}

function Welcome(props: { recent: string[]; onOpen: () => void; onNew: () => void; onRecent: (p: string) => void }) {
  return (
    <div className="vt-welcome">
      <h1>Vettore</h1>
      <p>Open a sequence file (.dna, .gb, .fasta, .ab1) or drop it here.</p>
      <div className="vt-welcome-actions">
        <button data-testid="open" onClick={props.onOpen}>Open…  (Ctrl+O)</button>
        <button data-testid="new" onClick={props.onNew}>New DNA File…  (Ctrl+N)</button>
      </div>
      {props.recent.length > 0 && (
        <div className="vt-recent">
          <h3>Recent files</h3>
          <ul>
            {props.recent.map((p) => (
              <li key={p}>
                <a href="#" onClick={(e) => { e.preventDefault(); props.onRecent(p); }} title={p}>
                  {baseName(p)}
                </a>
                <span className="vt-recent-dir">{p}</span>
              </li>
            ))}
          </ul>
        </div>
      )}
    </div>
  );
}
