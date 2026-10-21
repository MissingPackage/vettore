import { useState } from "react";
import type { Sample } from "./api";
import { baseName } from "./documents";
import { HelpDialogById } from "./tools/help";
import MiniMap from "./MiniMap";
import "./Welcome.css";

/** Shown while no document is open: what Vettore is, and the three ways in. */
export default function Welcome(props: {
  recent: string[];
  /** The sample constructs that come with Vettore (App lists them for the File menu too). */
  samples: Sample[];
  onOpen: () => void;
  onNew: () => void;
  onRecent: (path: string) => void;
  /** Resolves once the sample is open, or once opening it failed (the links come back). */
  onSample: (id: string) => void | Promise<void>;
  /** Open a sequence copied from an email or a web page. */
  onPaste?: () => void;
  /** File → Open from NCBI…; absent when the preference turns it off. */
  onNcbi?: () => void;
}) {
  const { samples } = props;
  const [opening, setOpening] = useState<string | null>(null);
  // The menu bar belongs to the editor, so it is not on screen here: the help dialogs a
  // first-time user needs are opened from the welcome screen itself.
  const [help, setHelp] = useState<string | null>(null);

  const openSample = (id: string) => {
    setOpening(id);
    Promise.resolve(props.onSample(id)).finally(() => setOpening(null));
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
          <ActionIcon d="M3 7h6l2 2h10v10H3z" />
          Open a file…
          <span>Ctrl+O</span>
        </button>
        <button data-testid="new" onClick={props.onNew}>
          <ActionIcon d="M6 3h9l4 4v14H6z M12 10v7 M8.5 13.5h7" />
          New DNA file…
          <span>Ctrl+N</span>
        </button>
        {props.onPaste && (
          <button data-testid="paste" onClick={props.onPaste}>
            <ActionIcon d="M8 4h8v4H8z M8 6H5v15h14V6h-3" />
            Paste a sequence…
            <span>Ctrl+Shift+V</span>
          </button>
        )}
      </div>
      <p className="vt-welcome-drop" data-testid="drop-hint">
        You can also drop a file anywhere on this window: GenBank (.gb, .gbk, .ape), .dna, EMBL,
        FASTA, .xdna, .gck, Sanger .ab1 or plain sequence text.
      </p>
      {props.onNcbi && (
        <p className="vt-welcome-drop" data-testid="welcome-ncbi-hint">
          A public GenBank record can be opened by its accession number:{" "}
          <button className="vt-link" data-testid="welcome-ncbi" onClick={props.onNcbi}>
            Open from NCBI…
          </button>{" "}
          (this uses the internet).
        </p>
      )}

      {samples.length > 0 && (
        <div className="vt-welcome-section vt-samples">
          <h3>Or try a sample construct</h3>
          <p className="vt-muted">
            Public GenBank records that come with Vettore. They open unsaved, so you can edit
            and cut them without touching anything of yours.
          </p>
          <ul className="vt-sample-cards">
            {samples.map((s) => (
              <li key={s.id}>
                <button
                  className="vt-sample-card"
                  data-testid={`sample-${s.id}`}
                  disabled={opening !== null}
                  onClick={() => openSample(s.id)}
                >
                  {s.length && s.features ? <MiniMap length={s.length} features={s.features} /> : null}
                  <span className="vt-sample-text">
                    <span className="vt-sample-name">{s.label}</span>
                    <span className="vt-welcome-note">{s.description}</span>
                    <span className="vt-sample-meta vt-num">
                      {s.length ? `${s.length.toLocaleString("en-US")} bp` : ""}
                      {s.length && s.accession ? " · " : ""}
                      {s.accession ?? ""}
                      {opening === s.id ? " · opening…" : ""}
                    </span>
                  </span>
                </button>
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

/** A line icon before each welcome action, drawn in the button's text colour. */
function ActionIcon(props: { d: string }) {
  return (
    <svg className="vt-welcome-icon" viewBox="0 0 24 24" width="18" height="18" aria-hidden="true"
         fill="none" stroke="currentColor" strokeWidth={2} strokeLinejoin="round" strokeLinecap="round">
      <path d={props.d} />
    </svg>
  );
}
