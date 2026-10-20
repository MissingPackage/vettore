import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useSelector } from "react-redux";
import { Editor } from "@teselagen/ove";
import { api, EngineError, whenDesktopReady, type OpenResult, type SaveResult, type SequenceData } from "./api";
import {
  baseName,
  canSaveInPlace,
  createDocument,
  editorState,
  markSaved,
  needsSave,
  now,
  removeDocument,
  stem,
  uniqueName,
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
import Welcome from "./Welcome";
import { toolById } from "./tools";
import { exportMap } from "./mapExport";
import { applyEnzymeSet, ENZYME_SETS, type EnzymeSetId } from "./enzymeSets";
import { ENZYMES } from "./enzymeLogic";
import { showView, VIEWS, type ViewId } from "./views";
import { ACCEPT, dirName, SAVE_EXTENSION, safeFileName, type SaveKind } from "./files/formats";
import { addNotice, dismissNotice, messageOf, type Notice, type NoticeKind } from "./files/notices";
import NoticeStack from "./files/NoticeStack";
import { readClipboard, writeClipboard } from "./files/clipboard";
import { fastaRecord, rangeLength, selectionRange } from "./files/fasta";
import { resolveDrop } from "./files/drop";
import { EDITOR_TOOLS, rightClickOverrides } from "./files/editorFiles";
import PasteDialog from "./files/PasteDialog";
import CopyTextDialog from "./files/CopyTextDialog";
import SaveCopyDialog from "./files/SaveCopyDialog";
import UnsavedDialog, { type UnsavedChoice } from "./files/UnsavedDialog";
import "./files/busy.css";

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

function fromBase64(b64: string): Blob {
  return new Blob([Uint8Array.from(atob(b64), (c) => c.charCodeAt(0))], { type: "application/octet-stream" });
}

const count = (x: unknown) => (Array.isArray(x) ? x.length : x && typeof x === "object" ? Object.keys(x).length : 0);

/** Text fields keep their own copy and paste. */
function isEditable(target: EventTarget | null): boolean {
  const el = target as HTMLElement | null;
  return !!el && (el.isContentEditable || ["INPUT", "TEXTAREA", "SELECT"].includes(el.tagName));
}

function listNames(names: string[]): string {
  return names.length <= 6 ? names.join(", ") : `${names.slice(0, 5).join(", ")} and ${names.length - 5} more`;
}

export default function App() {
  const [docs, setDocs] = useState<Doc[]>([]);
  const [active, setActive] = useState<string | null>(null);
  const [notices, setNotices] = useState<Notice[]>([]);
  const [newOpen, setNewOpen] = useState(false);
  const [paste, setPaste] = useState<{ text: string; reason?: string } | null>(null);
  const [copyText, setCopyText] = useState<{ title: string; text: string } | null>(null);
  const [saveCopy, setSaveCopy] = useState<{ doc: Doc; resolve: (ok: boolean) => void } | null>(null);
  const [unsaved, setUnsaved] = useState<{ doc: Doc; resolve: (choice: UnsavedChoice) => void } | null>(null);
  const [cloning, setCloning] = useState<CloningAction | null>(null);
  const [gelOpen, setGelOpen] = useState(false);
  const [sangerOpen, setSangerOpen] = useState(false);
  const [detectOpen, setDetectOpen] = useState(false);
  const [toolId, setToolId] = useState<string | null>(null);
  const [recent, setRecent] = useState<string[]>([]);
  const fileInput = useRef<HTMLInputElement>(null);
  const height = useWindowHeight();
  const docsRef = useRef(docs);
  docsRef.current = docs;

  const activeDoc = docs.find((d) => d.editorName === active) ?? null;
  // Edited since opened or saved, or never saved at all: the tab shows a dot, closing asks.
  const unsavedByEditor = useSelector((state: any) => {
    const out: Record<string, boolean> = {};
    for (const d of docs) out[d.editorName] = needsSave(d, state.VectorEditor[d.editorName]);
    return out;
  }, (a, b) => JSON.stringify(a) === JSON.stringify(b));
  const circular = useSelector((state: any) =>
    active ? !!state.VectorEditor[active]?.sequenceData?.circular : false,
  );
  const activeSequence = useSelector((state: any) =>
    active ? (state.VectorEditor[active]?.sequenceData?.sequence as string | undefined) : undefined,
  );

  // ---- messages --------------------------------------------------------------------------
  const notify = useCallback((kind: NoticeKind, text: string) => {
    setNotices((prev) => addNotice(prev, kind, text));
  }, []);
  const dismiss = useCallback((id: number) => setNotices((prev) => dismissNotice(prev, id)), []);
  const report = useCallback((e: unknown) => notify("error", messageOf(e)), [notify]);

  // An engine error no dialog caught still reaches the message stack instead of the console.
  useEffect(() => {
    const onRejection = (e: PromiseRejectionEvent) => {
      if (!(e.reason instanceof EngineError)) return;
      e.preventDefault();
      report(e.reason);
    };
    window.addEventListener("unhandledrejection", onRejection);
    return () => window.removeEventListener("unhandledrejection", onRejection);
  }, [report]);

  const unsavedCount = Object.values(unsavedByEditor).filter(Boolean).length;
  useEffect(() => {
    // The desktop window asks before quitting while this is above zero.
    api.unsaved(unsavedCount).catch(report);
  }, [unsavedCount, report]);

  // Opening and saving large files takes a moment: show a busy pointer meanwhile.
  const [busy, setBusy] = useState(0);
  const whileBusy = useCallback(async <T,>(work: Promise<T>): Promise<T> => {
    setBusy((n) => n + 1);
    try {
      return await work;
    } finally {
      setBusy((n) => n - 1);
    }
  }, []);

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
    api.recent().then((r) => setRecent(r.paths)).catch(report);
  }, [report]);

  /** Add new tabs after the open ones and show the last. */
  const addDocs = useCallback((created: Doc[]) => {
    if (!created.length) return;
    // Updated at once, so a second document added before the next render sees these names.
    docsRef.current = [...docsRef.current, ...created];
    setDocs((prev) => [...prev, ...created]);
    setActive(created[created.length - 1].editorName);
  }, []);

  /**
   * A document that is not a file on disk gets a name no open document has ("Pasted sequence 2"),
   * so tabs and the documents lists in dialogs can tell them apart.
   */
  const distinctName = useCallback((data: SequenceData, taken: string[]): SequenceData => {
    const name = data.name ? uniqueName(data.name, taken) : data.name;
    taken.push(name);
    return name === data.name ? data : { ...data, name };
  }, []);

  /** Open the documents of a file (or of pasted text, when `origin` is null) in new tabs. */
  const addResult = useCallback((result: OpenResult, origin: string | null) => {
    const label = result.format_label;
    const fileName = result.path ? baseName(result.path) : origin ? baseName(origin) : "";
    const records = result.documents.length;
    const source = label && fileName
      ? { label, records, fileName, incomplete: !!result.warnings?.length }
      : null;
    const taken = docsRef.current.map((d) => d.name);
    const created = result.documents.map((data) => {
      const sequenceData = result.path ? data : distinctName(data, taken);
      return createDocument(
        sequenceData,
        result.path,
        origin === null
          ? { time: now(), action: "Created from clipboard", detail: `${label ?? "Text"}, ${sequenceData.sequence.length} bp` }
          : { time: now(), action: "Opened", detail: `${result.path ?? origin}${label ? ` (${label})` : ""}` },
        {
          ...(source ? { saveInPlace: !!result.path && !!result.writable, source } : {}),
          // Pasted text exists nowhere else; a file or an untouched sample does.
          neverSaved: origin === null,
        },
      );
    });
    addDocs(created);
    if (records > 1) {
      const from = origin === null ? "the clipboard" : `“${fileName}”`;
      notify("info", `${records} sequences opened from ${from}, one tab each: ${listNames(created.map((d) => d.name))}.`);
    }
    for (const warning of result.warnings ?? []) notify("warning", warning);
  }, [addDocs, distinctName, notify]);

  // ---- open ------------------------------------------------------------------------------
  const openPath = useCallback(
    async (path: string) => {
      const existing = docsRef.current.find((d) => d.path === path);
      if (existing) {
        setActive(existing.editorName);
        return;
      }
      try {
        addResult(await whileBusy(api.openPath(path)), path);
        refreshRecent();
      } catch (e) {
        report(e);
      }
    },
    [addResult, refreshRecent, report, whileBusy],
  );

  const openFiles = useCallback(async () => {
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
          addResult(await whileBusy(api.openContent(file)), file.name);
        } catch (e) {
          report(e);
        }
      }
      if (fileInput.current) fileInput.current.value = "";
    },
    [addResult, report, whileBusy],
  );

  useEffect(() => {
    refreshRecent();
    api
      .startupFiles()
      .then(({ results, errors }) => {
        results.forEach((r) => addResult(r, r.path ?? ""));
        errors.forEach((message) => notify("error", message));
      })
      .catch(report);
  }, [addResult, notify, refreshRecent, report]);

  const openSample = useCallback(
    async (id: string) => {
      try {
        addResult(await api.openSample(id), "Sample construct");
      } catch (e) {
        report(e);
      }
    },
    [addResult, report],
  );

  /** Open sequence data in a new tab (cloning products, tool results, History snapshots). */
  const openDocument = useCallback((data: SequenceData, action: string, detail?: string) => {
    const sequenceData = distinctName(data, docsRef.current.map((d) => d.name));
    addDocs([createDocument(sequenceData, null, { time: now(), action, detail }, { neverSaved: true })]);
  }, [addDocs, distinctName]);

  // ---- new -------------------------------------------------------------------------------
  const createNew = useCallback((seq: NewSequence) => {
    setNewOpen(false);
    const sequenceData = distinctName({ ...seq, features: [], primers: [] }, docsRef.current.map((d) => d.name));
    // An empty new file has nothing to lose yet; once a base is typed it counts as edited.
    const doc = createDocument(sequenceData, null, { time: now(), action: "Created", detail: `${seq.sequence.length} bp` },
      { neverSaved: seq.sequence.length > 0 });
    addDocs([doc]);
  }, [addDocs, distinctName]);

  // ---- clipboard -------------------------------------------------------------------------
  const newFromClipboard = useCallback(async () => {
    const text = await readClipboard();
    if (text === null) {
      setPaste({ text: "", reason: "Vettore could not read the clipboard directly: paste the text into the box below." });
      return;
    }
    if (!text.trim()) {
      setPaste({ text: "", reason: "The clipboard is empty. Copy a sequence first, or paste it into the box below." });
      return;
    }
    try {
      addResult(await api.openText(text), null);
    } catch {
      // Show what the clipboard held; the box explains, next to the text, what is wrong with it.
      setPaste({ text });
    }
  }, [addResult]);

  const copyAsFasta = useCallback(async () => {
    if (!activeDoc) return;
    const state = editorState(activeDoc.editorName);
    const sequence = (state?.sequenceData?.sequence as string | undefined) ?? "";
    if (!sequence) {
      notify("warning", `“${activeDoc.name}” has no sequence to copy.`);
      return;
    }
    const range = selectionRange(sequence.length, state?.selectionLayer);
    const text = fastaRecord(activeDoc.name, sequence, range);
    const bp = (range ? rangeLength(range, sequence.length) : sequence.length).toLocaleString("en-US");
    const what = range ? `bases ${range.start + 1}–${range.end + 1} (${bp} bp)` : `the whole sequence (${bp} bp)`;
    if (await writeClipboard(text)) notify("info", `Copied ${what} of “${activeDoc.name}” as FASTA.`);
    else setCopyText({ title: "Copy as FASTA", text });
  }, [activeDoc, notify]);

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
    // Snapshot the inputs as they are now, so History can reopen them after they change.
    const parents = docsRef.current
      .filter((d) => product.parents.includes(d.name))
      .map((d) => ({ name: d.name, data: editorState(d.editorName)?.sequenceData as SequenceData }))
      .filter((p) => p.data);
    const sequenceData = distinctName(product.doc as SequenceData, docsRef.current.map((d) => d.name));
    const doc = createDocument(sequenceData, null, {
      time: now(),
      action: product.operation,
      detail: `${product.details}. From: ${product.parents.join(", ")}.`
        + (oligos ? ` Primers: ${oligos}` : ""),
      parents,
    }, { neverSaved: true });
    addDocs([doc]);
  }, [addDocs, distinctName]);

  // ---- save ------------------------------------------------------------------------------
  /** Ask once whether to save a GenBank copy of a file that cannot be overwritten safely. */
  const askSaveCopy = useCallback(
    (doc: Doc) => new Promise<boolean>((resolve) => setSaveCopy({ doc, resolve })),
    [],
  );

  /** FASTA keeps only the sequence: say what the saved file does not contain. */
  const warnIfAnnotationsLost = useCallback((saved: SaveResult, data: SequenceData) => {
    const features = count(data.features);
    const primers = count(data.primers);
    if (saved.format !== "fasta" || !(features || primers)) return;
    const parts = [features && `${features} feature${features > 1 ? "s" : ""}`,
      primers && `${primers} primer${primers > 1 ? "s" : ""}`].filter(Boolean).join(" and ");
    notify("warning", `“${baseName(saved.path)}” was saved as FASTA, which keeps only the sequence: `
      + `its ${parts} are not in that file. To keep them, use File → Save As… and choose GenBank (.gb).`);
  }, [notify]);

  const saveDoc = useCallback(
    async (doc: Doc, saveAs: boolean): Promise<boolean> => {
      const sequenceData = editorState(doc.editorName)?.sequenceData as SequenceData;
      if (!sequenceData) return false;
      const current = docsRef.current.find((d) => d.editorName === doc.editorName) ?? doc;
      try {
        const desktop = await whenDesktopReady(0);
        if (!desktop) {
          // Browser: a download of a GenBank file; the original is never touched.
          const filename = `${safeFileName(current.name)}.gb`;
          const { content } = await api.exportText(sequenceData, "genbank");
          download(filename, content);
          logHistory(current.editorName, "Downloaded", filename);
          markSaved(current.editorName);
          updateDoc(current.editorName, { neverSaved: false });
          return true;
        }
        let path = !saveAs && canSaveInPlace(current) ? current.path : null;
        if (!path) {
          if (!saveAs && current.path && current.source && !current.copyNoticeShown) {
            updateDoc(current.editorName, { copyNoticeShown: true });
            if (!(await askSaveCopy(current))) return false;
          }
          path = await desktop.save_dialog(
            `${safeFileName(current.name)}${SAVE_EXTENSION.genbank}`,
            current.path ? dirName(current.path) : "",
            "genbank",
          );
          if (!path) return false;
        }
        const saved = await whileBusy(api.savePath(path, sequenceData, "genbank"));
        updateDoc(current.editorName, {
          path: saved.path,
          name: stem(saved.path),
          saveInPlace: true,
          source: { label: saved.format_label, records: 1, fileName: baseName(saved.path), incomplete: false },
          neverSaved: false,
        });
        logHistory(current.editorName, "Saved", `${saved.path} (${saved.format_label})`);
        markSaved(current.editorName);
        refreshRecent();
        notify("info", `Saved “${baseName(saved.path)}”.`);
        warnIfAnnotationsLost(saved, sequenceData);
        return true;
      } catch (e) {
        report(e);
        return false;
      }
    },
    [askSaveCopy, logHistory, notify, refreshRecent, report, updateDoc, warnIfAnnotationsLost, whileBusy],
  );

  /** Write a copy in another format; the document keeps its own file. */
  const exportAs = useCallback(
    async (doc: Doc, kind: SaveKind) => {
      const sequenceData = editorState(doc.editorName)?.sequenceData as SequenceData;
      const filename = `${safeFileName(doc.name)}${SAVE_EXTENSION[kind]}`;
      try {
        const desktop = await whenDesktopReady(0);
        if (desktop) {
          const path = await desktop.save_dialog(filename, doc.path ? dirName(doc.path) : "", kind);
          if (!path) return;
          const saved = await whileBusy(api.savePath(path, sequenceData, kind));
          logHistory(doc.editorName, "Exported", `${saved.path} (${saved.format_label})`);
          notify("info", `Exported “${baseName(saved.path)}” (${saved.format_label}).`);
        } else {
          const file = await api.exportFile(sequenceData, kind);
          download(filename, file.content_b64 ? fromBase64(file.content_b64) : file.content ?? "");
          logHistory(doc.editorName, "Downloaded", filename);
        }
      } catch (e) {
        report(e);
      }
    },
    [logHistory, notify, report, whileBusy],
  );

  const exportImage = useCallback(
    async (doc: Doc, format: "svg" | "png") => {
      try {
        const blob = await exportMap(format);
        const desktop = await whenDesktopReady(0);
        const filename = `${safeFileName(doc.name)}-map.${format}`;
        if (desktop) {
          const path = await desktop.save_dialog(filename, doc.path ? dirName(doc.path) : "", format);
          if (path) {
            const saved = await api.saveBinary(path, blob);
            notify("info", `Exported “${baseName(saved.path)}”.`);
          }
        } else {
          download(filename, blob);
        }
      } catch (e) {
        report(e);
      }
    },
    [notify, report],
  );

  // ---- close -----------------------------------------------------------------------------
  /** Save, close without saving, or keep the tab: asked when closing a tab with unsaved work. */
  const askUnsaved = useCallback(
    (doc: Doc) => new Promise<UnsavedChoice>((resolve) => setUnsaved({ doc, resolve })),
    [],
  );
  const asking = useRef(false);

  const closeDoc = useCallback(
    async (doc: Doc) => {
      if (asking.current) return;   // Ctrl+W again while the question is on screen
      const current = docsRef.current.find((d) => d.editorName === doc.editorName) ?? doc;
      if (needsSave(current, editorState(current.editorName))) {
        asking.current = true;
        let choice: UnsavedChoice;
        try {
          choice = await askUnsaved(current);
        } finally {
          asking.current = false;
        }
        if (choice === "cancel") return;
        if (choice === "save" && !(await saveDoc(current, false))) return;
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
    [askUnsaved, saveDoc],
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
      } else if (key === "v" && e.shiftKey && !isEditable(e.target)) {
        // Cancels the browser's paste too, so the open sequence is left as it is.
        e.preventDefault();
        newFromClipboard();
      } else if (key === "n" && !e.shiftKey) {
        e.preventDefault();
        setNewOpen(true);
      } else if (key === "g" && e.shiftKey && activeDoc) {
        e.preventDefault();
        setGelOpen(true);
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
  }, [activeDoc, closeDoc, newFromClipboard, openFiles, selectView]);

  // ---- paste on the start screen: Ctrl+V with no document open opens what was copied -----
  useEffect(() => {
    if (docs.length) return;
    const onPaste = (e: ClipboardEvent) => {
      if (isEditable(e.target)) return;
      const text = e.clipboardData?.getData("text/plain") ?? "";
      if (!text.trim()) return;
      e.preventDefault();
      api.openText(text).then((result) => addResult(result, null)).catch(() => setPaste({ text }));
    };
    document.addEventListener("paste", onPaste);
    return () => document.removeEventListener("paste", onPaste);
  }, [addResult, docs.length]);

  // ---- drag and drop ---------------------------------------------------------------------
  const onDrop = useCallback(
    (e: React.DragEvent) => {
      if (!e.dataTransfer.files.length) return;
      // No stopPropagation: the desktop app listens on the document for the files' full paths.
      e.preventDefault();
      resolveDrop(Array.from(e.dataTransfer.files)).then(async ({ paths, files }) => {
        for (const path of paths) await openPath(path);
        await openBrowserFiles(files);
      });
    },
    [openBrowserFiles, openPath],
  );

  // ---- menu ------------------------------------------------------------------------------
  const handlers = useRef({
    newFile: () => setNewOpen(true),
    newFromClipboard,
    open: openFiles,
    openRecent: openPath,
    close: () => activeDoc && closeDoc(activeDoc),
    save: () => activeDoc && saveDoc(activeDoc, false),
    saveAs: () => activeDoc && saveDoc(activeDoc, true),
    exportGenbank: () => activeDoc && exportAs(activeDoc, "genbank"),
    exportFasta: () => activeDoc && exportAs(activeDoc, "fasta"),
    exportDna: () => activeDoc && exportAs(activeDoc, "dna"),
    copyAsFasta,
    exportSvg: () => activeDoc && exportImage(activeDoc, "svg"),
    exportPng: () => activeDoc && exportImage(activeDoc, "png"),
    view: selectView,
    useEnzymeSet,
    cloning: (id: string) => setCloning(actionById(id)),
    gel: () => setGelOpen(true),
    sanger: () => setSangerOpen(true),
    detectFeatures: () => setDetectOpen(true),
    tool: (id: string) => setToolId(id),
    recent: [] as string[],
    hasDocument: false,
  });
  Object.assign(handlers.current, {
    newFromClipboard,
    open: openFiles,
    openRecent: openPath,
    close: () => activeDoc && closeDoc(activeDoc),
    save: () => activeDoc && saveDoc(activeDoc, false),
    saveAs: () => activeDoc && saveDoc(activeDoc, true),
    exportGenbank: () => activeDoc && exportAs(activeDoc, "genbank"),
    exportFasta: () => activeDoc && exportAs(activeDoc, "fasta"),
    exportDna: () => activeDoc && exportAs(activeDoc, "dna"),
    copyAsFasta,
    exportSvg: () => activeDoc && exportImage(activeDoc, "svg"),
    exportPng: () => activeDoc && exportImage(activeDoc, "png"),
    view: selectView,
    useEnzymeSet,
    cloning: (id: string) => setCloning(actionById(id)),
    gel: () => setGelOpen(true),
    sanger: () => setSangerOpen(true),
    detectFeatures: () => setDetectOpen(true),
    tool: (id: string) => setToolId(id),
    recent,
    hasDocument: !!activeDoc,
  });
  // menuFilter must be stable (OVE memoizes on it); handlers are read through the ref.
  const menuFilter = useMemo(() => buildMenu(handlers), [recent, !!activeDoc]); // eslint-disable-line react-hooks/exhaustive-deps
  const contextMenus = useMemo(() => rightClickOverrides(() => handlers.current.copyAsFasta()), []);

  const onSave = useCallback(
    (opts: { isSaveAs?: boolean }) => (activeDoc ? saveDoc(activeDoc, !!opts?.isSaveAs).then(() => undefined) : undefined),
    [activeDoc, saveDoc],
  );

  const editorHeight = height - STATUS_HEIGHT;

  // Each document's editor draws its own tab strip, scrolled to the start: bring the active
  // tab into sight, or a new product's tab would sit off-screen once a few files are open.
  // Runs when the tab mounts and whenever another document becomes active.
  const revealActiveTab = useCallback((tab: HTMLElement | null) => {
    const strip = tab?.parentElement;
    if (!tab || !strip) return;
    const t = tab.getBoundingClientRect();
    const s = strip.getBoundingClientRect();
    if (t.left < s.left) strip.scrollLeft -= s.left - t.left + 8;
    else if (t.right > s.right) strip.scrollLeft += t.right - s.right + 8;
  }, [active]); // eslint-disable-line react-hooks/exhaustive-deps

  const docTabs = (
      <div className="vt-doctabs" style={{ height: DOCTABS_HEIGHT }}>
        {docs.map((d) => (
          <div
            key={d.editorName}
            ref={d.editorName === active ? revealActiveTab : undefined}
            className={`vt-doctab${d.editorName === active ? " active" : ""}`}
            onClick={() => setActive(d.editorName)}
            onAuxClick={(e) => e.button === 1 && closeDoc(d)}
            title={d.path ?? (d.neverSaved ? `${d.name} (not saved to a file yet)` : d.name)}
            data-testid="doctab"
            data-unsaved={unsavedByEditor[d.editorName] ? "true" : undefined}
          >
            <span>{unsavedByEditor[d.editorName] ? "• " : ""}{d.name}</span>
            <button className="vt-doctab-close" title="Close (Ctrl+W)" onClick={(e) => { e.stopPropagation(); closeDoc(d); }}>×</button>
          </div>
        ))}
        {!docs.length && <span className="vt-brand">Vettore</span>}
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
    <div className={`vt-shell${busy ? " vt-shell-busy" : ""}`} aria-busy={busy > 0}
         onDragOver={(e) => e.preventDefault()} onDropCapture={onDrop}>
      <input
        ref={fileInput}
        type="file"
        multiple
        hidden
        data-testid="file-input"
        accept={ACCEPT}
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
                <HistoryView
                  doc={activeDoc}
                  onReopen={(data, from) => openDocument({ ...data, name: `${data.name} (from History)` }, "Reopened", from)}
                />
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
              ToolBarProps={{ toolList: EDITOR_TOOLS }}
              rightClickOverrides={contextMenus}
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
          <Welcome recent={recent} onOpen={openFiles} onNew={() => setNewOpen(true)} onRecent={openPath}
                   onPaste={newFromClipboard}
                   onSample={openSample} onError={report} />
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
      <NoticeStack notices={notices} onDismiss={dismiss} onDismissAll={() => setNotices([])} />
      {/* Remounted on every open so the form never keeps the previous file's name or topology. */}
      {newOpen && <NewSequenceDialog isOpen onClose={() => setNewOpen(false)} onCreate={createNew} />}
      {paste && (
        <PasteDialog
          initialText={paste.text}
          reason={paste.reason}
          onClose={() => setPaste(null)}
          onOpen={(result) => {
            setPaste(null);
            addResult(result, null);
          }}
        />
      )}
      {copyText && <CopyTextDialog title={copyText.title} text={copyText.text} onClose={() => setCopyText(null)} />}
      {saveCopy?.doc.source && (
        <SaveCopyDialog
          name={saveCopy.doc.name}
          source={saveCopy.doc.source}
          onCancel={() => {
            saveCopy.resolve(false);
            setSaveCopy(null);
          }}
          onConfirm={() => {
            saveCopy.resolve(true);
            setSaveCopy(null);
          }}
        />
      )}
      {unsaved && (
        <UnsavedDialog
          name={unsaved.doc.name}
          neverSaved={unsaved.doc.neverSaved}
          onChoose={(choice) => {
            unsaved.resolve(choice);
            setUnsaved(null);
          }}
        />
      )}
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
      {toolId && (() => {
        const tool = toolById(toolId);
        if (!tool) return null;
        const { Component } = tool;
        return (
          <Component
            docs={openDocs}
            active={openDocs.find((d) => d.editorName === active) ?? null}
            onClose={() => setToolId(null)}
            openDocument={openDocument}
            logHistory={logHistory}
          />
        );
      })()}
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
