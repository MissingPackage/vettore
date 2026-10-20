import { useState } from "react";
import type { SequenceData } from "./api";
import type { Doc } from "./documents";
import { useCopyText } from "./cloning/useCopyText";
import "./history.css";

/** Provenance of the document: opened, created, saved, and every cloning step with its inputs. */
export default function HistoryView({ doc, onReopen }: {
  doc: Doc;
  /** Open a snapshot of an input document in a new tab. */
  onReopen: (data: SequenceData, from: string) => void;
}) {
  const steps = [...doc.history].reverse();
  const { copy, fallbackDialog } = useCopyText();
  /** The step whose text was just put on the clipboard. */
  const [copied, setCopied] = useState<number | null>(null);
  const copyStep = async (i: number, text: string) => {
    setCopied(null);
    if (await copy(text, "Copy this step")) {
      setCopied(i);
      setTimeout(() => setCopied((current) => (current === i ? null : current)), 2500);
    }
  };
  return (
    <div className="vt-history" data-testid="history-view">
      <h2>History of “{doc.name}”</h2>
      <p className="vt-history-intro">
        What this file is made of, newest step first. A cloning step keeps a copy of everything it
        was run on, so you can open an input exactly as it was at the time — even if the file has
        changed since, or was never saved.
      </p>
      <table className="bp3-html-table bp3-html-table-condensed bp3-html-table-striped">
        <thead>
          <tr><th>When</th><th>Step</th><th>What was done</th><th>Run on</th><th /></tr>
        </thead>
        <tbody>
          {steps.map((h, i) => (
            <tr key={i}>
              <td>{h.time}</td>
              <td>{h.action}</td>
              <td className="vt-history-detail">{h.detail ?? ""}</td>
              <td>
                {h.parents?.length ? (
                  <span className="vt-history-inputs">
                    {h.parents.map((p, j) => (
                      <button
                        key={j}
                        className="bp3-button bp3-small bp3-minimal"
                        data-testid="history-reopen"
                        title={`Open “${p.name}” as it was when this step ran`}
                        onClick={() => onReopen(p.data, `${h.action} in “${doc.name}”`)}
                      >
                        Open “{p.name}”
                      </button>
                    ))}
                  </span>
                ) : (
                  <span className="vt-history-none">—</span>
                )}
              </td>
              <td>
                {h.detail && (
                  <button
                    className="bp3-button bp3-small bp3-minimal vt-history-copy"
                    data-testid="history-copy"
                    title="Copy this step (parameters and primers) to the clipboard"
                    onClick={() => copyStep(i, `${h.action}: ${h.detail}`)}
                  >
                    {copied === i ? "Copied" : "Copy"}
                  </button>
                )}
              </td>
            </tr>
          ))}
        </tbody>
      </table>
      {fallbackDialog}
    </div>
  );
}
