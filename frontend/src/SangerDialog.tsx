import { useCallback, useEffect, useRef, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import { AlignmentView, addAlignment } from "@teselagen/ove";
import { api } from "./api";
import { store } from "./store";
import type { OpenDoc } from "./cloning/types";

const { Dialog, Callout } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

const ALIGNMENT_ID = "vettore-sanger";

interface Discrepancy {
  kind: "mismatch" | "insertion" | "deletion";
  position: number;
  reference: string;
  read: string;
}

interface ReadSummary {
  name: string;
  file: string;
  strand: number;
  start: number;
  end: number;
  coverage: number;
  identity: number;
  discrepancies: Discrepancy[];
  trimmed_bases: number;
  original_length: number;
  spans_origin: boolean;
  machine: string;
}

interface SangerResult {
  alignmentTracks: unknown[];
  reads: ReadSummary[];
  warnings: string[];
  reference: { name: string; length: number; circular: boolean };
  coverage_percent: number;
  first_covered: number;
}

const TRIM_LEVELS = [
  { value: 0.01, label: "Strict (≈ Q20)" },
  { value: 0.05, label: "Normal (≈ Q13)" },
  { value: 0.2, label: "Permissive (≈ Q7)" },
];

function describe(d: Discrepancy): string {
  if (d.kind === "mismatch") return `${d.position}: ${d.reference} → ${d.read}`;
  if (d.kind === "deletion") return `${d.position}: ${d.reference || "?"} missing in the read`;
  return `${d.position}: ${d.read} extra in the read`;
}

/** Tools → Align to Reference DNA Sequence: verify a clone against the expected plasmid. */
export default function SangerDialog(props: { active: OpenDoc | null; onClose: () => void }) {
  const [files, setFiles] = useState<File[]>([]);
  const [cutoff, setCutoff] = useState(0.05);
  const [result, setResult] = useState<SangerResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const input = useRef<HTMLInputElement>(null);
  const reference = props.active?.data();

  const run = useCallback(async () => {
    if (!props.active || !files.length) return;
    setBusy(true);
    setError(null);
    try {
      const payload = await api.sangerAlign<SangerResult>({
        reference: props.active.data(),
        reads: await Promise.all(files.map(async (f) => ({
          filename: f.name,
          content_b64: await fileToBase64(f),
        }))),
        cutoff,
      });
      setResult(payload);
      if (payload.alignmentTracks.length) {
        addAlignment(store, {
          id: ALIGNMENT_ID,
          name: `${payload.reference.name} — Sanger reads`,
          alignmentType: "Sanger reads",
          alignmentTracks: payload.alignmentTracks,
          // Open Vector Editor hides chromatograms by default.
          alignmentAnnotationVisibility: { chromatogram: true },
          // Start where the reads actually are, not at base 1.
          scrollPercentageToJumpTo:
            Math.max(0, payload.first_covered - 30) / Math.max(1, payload.reference.length),
        });
      }
    } catch (e) {
      setResult(null);
      setError(String((e as Error).message));
    } finally {
      setBusy(false);
    }
  }, [props.active, files, cutoff]);

  useEffect(() => {
    if (files.length) run();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [files, cutoff]);

  return (
    <Dialog isOpen onClose={props.onClose} title="Align to Reference DNA Sequence"
            style={{ width: "92vw", maxWidth: 1500 }}>
      {/* Blueprint's Dialog does not forward unknown props, so the test id lives here. */}
      <div className={Classes.DIALOG_BODY} data-testid="sanger-dialog">
        <div className="vt-row vt-sanger-bar">
          <span>
            Reference: <b>{props.active?.name ?? "—"}</b>
            {reference ? ` (${reference.sequence.length} bp, ${reference.circular ? "circular" : "linear"})` : ""}
          </span>
          <input
            ref={input}
            type="file"
            multiple
            accept=".ab1,.abi,.fsa"
            data-testid="sanger-files"
            onChange={(e) => setFiles(Array.from(e.target.files ?? []))}
          />
          <label>
            Trim&nbsp;
            <div className="bp3-html-select bp3-minimal">
              <select data-testid="sanger-trim" value={cutoff}
                      onChange={(e) => setCutoff(Number(e.target.value))}>
                {TRIM_LEVELS.map((t) => <option key={t.value} value={t.value}>{t.label}</option>)}
              </select>
            </div>
          </label>
          {busy && <span className="vt-muted">Aligning…</span>}
        </div>

        {error && <Callout intent="danger" title="Could not align">{error}</Callout>}
        {result?.warnings.map((w, i) => (
          <Callout key={i} intent="warning" data-testid="sanger-warning">{w}</Callout>
        ))}

        {result && result.reads.length > 0 && (
          <>
            <table className="bp3-html-table bp3-html-table-condensed vt-table" data-testid="sanger-summary">
              <thead>
                <tr>
                  <th>Read</th><th>Strand</th><th>Covers</th><th>Identity</th>
                  <th>Differences</th><th>Trimmed</th>
                </tr>
              </thead>
              <tbody>
                {result.reads.map((r) => (
                  <tr key={r.file}>
                    <td title={`${r.file}${r.machine ? ` · ${r.machine}` : ""}`}>{r.name}</td>
                    <td>{r.strand === 1 ? "forward" : "reverse"}</td>
                    <td>{r.start.toLocaleString()}…{r.end.toLocaleString()}{r.spans_origin ? " (crosses the origin)" : ""}</td>
                    <td className={r.identity === 100 ? "vt-ok" : "vt-warn"}>{r.identity.toFixed(2)} %</td>
                    <td>
                      {r.discrepancies.length === 0
                        ? "none"
                        : r.discrepancies.slice(0, 6).map(describe).join("; ")
                          + (r.discrepancies.length > 6 ? ` … (${r.discrepancies.length} total)` : "")}
                    </td>
                    <td className="vt-muted">{r.trimmed_bases} of {r.original_length} nt</td>
                  </tr>
                ))}
              </tbody>
            </table>
            <p className="vt-muted">
              {result.coverage_percent} % of {result.reference.name} is covered by these reads.
              Differences are listed with reference coordinates; red marks in the view below are
              mismatches, grey marks are regions a read does not cover.
            </p>
            <div className="vt-alignment" data-testid="sanger-alignment">
              <AlignmentView id={ALIGNMENT_ID} height={420} />
            </div>
          </>
        )}

        {!result && !busy && !error && (
          <p className="vt-muted">
            Choose the .ab1 files from your sequencing provider. They are aligned to the active
            document, on either strand, after quality trimming.
          </p>
        )}
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button onClick={props.onClose}>Close</Button>
        </div>
      </div>
    </Dialog>
  );
}

async function fileToBase64(file: File): Promise<string> {
  const bytes = new Uint8Array(await file.arrayBuffer());
  let binary = "";
  for (let i = 0; i < bytes.length; i += 0x8000) {
    binary += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
  }
  return btoa(binary);
}
