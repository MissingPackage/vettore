import { useEffect, useState } from "react";
import { api, type Sample } from "./api";
import { baseName } from "./documents";
import { HelpDialogById } from "./tools/help";
import "./Welcome.css";

/** Shown while no document is open: what Vettore is, and the three ways in. */
export default function Welcome(props: {
  recent: string[];
  onOpen: () => void;
  onNew: () => void;
  onRecent: (path: string) => void;
  onSample: (id: string) => void;
  /** Open a sequence copied from an email or a web page. */
  onPaste?: () => void;
}) {
  const [samples, setSamples] = useState<Sample[]>([]);
  const [opening, setOpening] = useState<string | null>(null);
  // The menu bar belongs to the editor, so it is not on screen here: the help dialogs a
  // first-time user needs are opened from the welcome screen itself.
  const [help, setHelp] = useState<string | null>(null);

  useEffect(() => {
    api.samples().then((r) => setSamples(r.samples)).catch(() => undefined);
  }, []);

  const openSample = (id: string) => {
    setOpening(id);
    props.onSample(id);
  };

  return (
    <div className="vt-welcome" data-testid="welcome">
      <h1>Vettore</h1>
      <p className="vt-welcome-lead">
        Read plasmid maps, plan the cloning, order the primers and check the clone against its
        sequencing reads — on this computer, with nothing to sign in to.
      </p>

      <div className="vt-welcome-actions">
        <button className="vt-primary" data-testid="open" onClick={props.onOpen}>
          Open a file…
          <span>Ctrl+O</span>
        </button>
        <button data-testid="new" onClick={props.onNew}>
          New DNA file…
          <span>Ctrl+N</span>
        </button>
        {props.onPaste && (
          <button data-testid="paste" onClick={props.onPaste}>
            Paste a sequence…
            <span>Ctrl+Shift+V</span>
          </button>
        )}
      </div>
      <p className="vt-welcome-drop" data-testid="drop-hint">
        You can also drop a file anywhere on this window: GenBank (.gb, .gbk, .ape), .dna, EMBL,
        FASTA, .xdna, .gck, Sanger .ab1 or plain sequence text.
      </p>

      {samples.length > 0 && (
        <div className="vt-welcome-section vt-samples">
          <h3>Or try a sample construct</h3>
          <p className="vt-muted">
            Public GenBank records that come with Vettore. They open unsaved, so you can edit
            and cut them without touching anything of yours.
          </p>
          <ul>
            {samples.map((s) => (
              <li key={s.id}>
                <button
                  className="vt-link"
                  data-testid={`sample-${s.id}`}
                  disabled={opening !== null}
                  onClick={() => openSample(s.id)}
                >
                  {s.label}
                </button>
                <span className="vt-welcome-note">
                  {s.description}
                  {s.accession ? ` · ${s.accession}` : ""}
                  {opening === s.id ? " · opening…" : ""}
                </span>
              </li>
            ))}
          </ul>
        </div>
      )}

      {props.recent.length > 0 && (
        <div className="vt-welcome-section vt-recent">
          <h3>Recent files</h3>
          <ul>
            {props.recent.map((p) => (
              <li key={p}>
                <button className="vt-link" onClick={() => props.onRecent(p)} title={p}>
                  {baseName(p)}
                </button>
                <span className="vt-welcome-note vt-welcome-path">{p}</span>
              </li>
            ))}
          </ul>
        </div>
      )}

      <div className="vt-welcome-foot">
        <span className="vt-muted">First time with Vettore:</span>
        <button className="vt-link" data-testid="welcome-getting-started"
                onClick={() => setHelp("getting-started")}>
          Getting started
        </button>
        <button className="vt-link" data-testid="welcome-shortcuts"
                onClick={() => setHelp("shortcuts")}>
          Keyboard shortcuts
        </button>
        <button className="vt-link" data-testid="welcome-about" onClick={() => setHelp("about")}>
          About Vettore
        </button>
      </div>
      {help && <HelpDialogById id={help} onClose={() => setHelp(null)} />}
    </div>
  );
}
