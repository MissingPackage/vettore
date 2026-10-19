import { useCallback, useEffect, useRef, useState } from "react";
import { Editor, updateEditor } from "@teselagen/ove";
import { api, whenDesktopReady, type SequenceData } from "./api";
import { store } from "./store";

interface Doc {
  editorName: string;
  name: string;
  path: string | null;
}

let docCounter = 0;
const TOPBAR_HEIGHT = 34;

export default function App() {
  const [docs, setDocs] = useState<Doc[]>([]);
  const [active, setActive] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [height, setHeight] = useState(window.innerHeight - TOPBAR_HEIGHT);
  const fileInput = useRef<HTMLInputElement>(null);

  useEffect(() => {
    const onResize = () => setHeight(window.innerHeight - TOPBAR_HEIGHT);
    window.addEventListener("resize", onResize);
    return () => window.removeEventListener("resize", onResize);
  }, []);

  const addDocuments = useCallback((documents: SequenceData[], path: string | null) => {
    const created = documents.map((sequenceData) => {
      const editorName = `doc-${++docCounter}`;
      updateEditor(store, editorName, { sequenceData, sequenceDataHistory: {} });
      return { editorName, name: sequenceData.name, path };
    });
    setDocs((prev) => [...prev, ...created]);
    setActive(created[created.length - 1]?.editorName ?? null);
  }, []);

  useEffect(() => {
    api
      .startupFiles()
      .then(({ results, errors }) => {
        results.forEach((r) => addDocuments(r.documents, r.path));
        if (errors.length) setError(errors.join("; "));
      })
      .catch((e) => setError(String((e as Error).message)));
  }, [addDocuments]);

  const openFiles = useCallback(async () => {
    setError(null);
    const desktop = await whenDesktopReady(0);
    if (!desktop) {
      fileInput.current?.click();
      return;
    }
    for (const path of await desktop.open_dialog()) {
      try {
        const result = await api.openPath(path);
        addDocuments(result.documents, result.path);
      } catch (e) {
        setError(String((e as Error).message));
      }
    }
  }, [addDocuments]);

  const onBrowserFiles = useCallback(
    async (files: FileList | null) => {
      for (const file of Array.from(files ?? [])) {
        try {
          const result = await api.openContent(file);
          addDocuments(result.documents, null);
        } catch (e) {
          setError(String((e as Error).message));
        }
      }
      if (fileInput.current) fileInput.current.value = "";
    },
    [addDocuments],
  );

  return (
    <div className="vt-shell">
      <div className="vt-topbar" style={{ height: TOPBAR_HEIGHT - 9 }}>
        <button onClick={openFiles} data-testid="open">Open…</button>
        <input
          ref={fileInput}
          type="file"
          multiple
          hidden
          data-testid="file-input"
          accept=".dna,.gb,.gbk,.gbff,.ape,.fa,.fasta,.fna,.ab1"
          onChange={(e) => onBrowserFiles(e.target.files)}
        />
        <div className="vt-doctabs">
          {docs.map((d) => (
            <div
              key={d.editorName}
              className={`vt-doctab${d.editorName === active ? " active" : ""}`}
              onClick={() => setActive(d.editorName)}
              title={d.path ?? d.name}
            >
              {d.name}
            </div>
          ))}
        </div>
        {error && <span className="vt-error" role="alert">{error}</span>}
      </div>
      <div className="vt-editor">
        {active ? (
          <Editor key={active} editorName={active} height={height} showMenuBar />
        ) : (
          <div className="vt-empty">Open a sequence file (.dna, .gb, .fasta, .ab1) to start.</div>
        )}
      </div>
    </div>
  );
}
