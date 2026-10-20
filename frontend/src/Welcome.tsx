import { useEffect, useState } from "react";
import { api, type Sample } from "./api";
import { baseName } from "./documents";

/** Shown while no document is open. */
export default function Welcome(props: {
  recent: string[];
  onOpen: () => void;
  onNew: () => void;
  onRecent: (path: string) => void;
  onSample: (id: string) => void;
}) {
  const [samples, setSamples] = useState<Sample[]>([]);
  useEffect(() => {
    api.samples().then((r) => setSamples(r.samples)).catch(() => undefined);
  }, []);

  return (
    <div className="vt-welcome">
      <h1>Vettore</h1>
      <p>Open a sequence file (.dna, .gb, .fasta, .ab1) or drop it here.</p>
      <div className="vt-welcome-actions">
        <button data-testid="open" onClick={props.onOpen}>Open…  (Ctrl+O)</button>
        <button data-testid="new" onClick={props.onNew}>New DNA File…  (Ctrl+N)</button>
      </div>
      {samples.length > 0 && (
        <div className="vt-samples">
          <h3>Sample constructs</h3>
          <ul>
            {samples.map((s) => (
              <li key={s.id}>
                <a href="#" data-testid={`sample-${s.id}`} onClick={(e) => { e.preventDefault(); props.onSample(s.id); }}>
                  {s.label}
                </a>
                <span className="vt-recent-dir">{s.description}</span>
              </li>
            ))}
          </ul>
        </div>
      )}
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
