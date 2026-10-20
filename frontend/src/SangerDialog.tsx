import { useCallback, useEffect, useRef, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import { AlignmentView, addAlignment } from "@teselagen/ove";
import { api } from "./api";
import { store } from "./store";
import type { OpenDoc } from "./cloning/types";

const { Dialog, Callout } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

const ALIGNMENT_ID = "vettore-sanger";
const SHOWN_DIFFERENCES = 6;
/** Phred Q20: one chance in 100 that the base call is wrong (Q = −10·log10 P). */
const GOOD_QUALITY = 20;

interface Discrepancy {
  kind: "mismatch" | "insertion" | "deletion";
  position: number;
  reference: string;
  read: string;
  /** Phred quality of the read at this difference (the lowest, for a run). */
  quality?: number;
}

interface ReadSummary {
  name: string;
  file: string;
  strand: number;
  start: number;
  end: number;
  coverage: number;
  identity: number;
  mean_quality?: number;
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

const n = (value: number) => value.toLocaleString("en-US");

function describe(d: Discrepancy): string {
  if (d.kind === "mismatch") return `${n(d.position)}: ${d.reference} → ${d.read}`;
  if (d.kind === "deletion") return `${n(d.position)}: ${d.reference || "?"} missing in the read`;
  return `${n(d.position)}: ${d.read} extra in the read`;
}

/** "about 1 chance in 16 that this base call is wrong", from the Phred definition. */
function errorChance(quality: number): string {
  const odds = Math.round(10 ** (quality / 10));
  return `about 1 chance in ${n(Math.max(1, odds))} that this base call is wrong`;
}

function Differences(props: { read: ReadSummary }) {
  const [all, setAll] = useState(false);
  const diffs = props.read.discrepancies;
  if (!diffs.length) return <span className="vt-ok">none</span>;
  const shown = all ? diffs : diffs.slice(0, SHOWN_DIFFERENCES);
  return (
    <span data-testid="sanger-differences">
      {shown.map((d, i) => {
        const weak = d.quality !== undefined && d.quality < GOOD_QUALITY;
        return (
          <span key={i}
                title={d.quality !== undefined ? `Read quality Q${d.quality}: ${errorChance(d.quality)}` : undefined}
                style={weak ? { color: "var(--vt-text-muted)" } : undefined}>
            {i > 0 && "; "}
            {describe(d)}
            {d.quality !== undefined && <span className="vt-muted"> (Q{d.quality}{weak ? ", low" : ""})</span>}
          </span>
        );
      })}
      {diffs.length > SHOWN_DIFFERENCES && (
        <Button small minimal onClick={() => setAll(!all)} data-testid="sanger-show-all">
          {all ? "show fewer" : `show all ${diffs.length}`}
        </Button>
      )}
    </span>
  );
}

/** One sentence that answers "is my clone right?" for the reads that aligned. */
function verdict(result: SangerResult): { text: string; ok: boolean } {
  const total = result.reads.reduce((sum, r) => sum + r.discrepancies.length, 0);
  const reads = `${result.reads.length} read${result.reads.length === 1 ? "" : "s"}`;
  const where = `${result.coverage_percent} % of ${result.reference.name}`;
  if (total === 0) {
    const subject = result.reads.length === 1
      ? "the read matches" : `all ${result.reads.length} reads match`;
    const cover = result.reads.length === 1 ? "it covers" : "they cover";
    return { ok: true, text: `No differences: ${subject} ${result.reference.name} wherever ${cover} it (${where}).` };
  }
  const weak = result.reads.flatMap((r) => r.discrepancies)
    .filter((d) => d.quality !== undefined && d.quality < GOOD_QUALITY).length;
  return {
    ok: false,
    text: `${n(total)} difference${total === 1 ? "" : "s"} in ${reads} (covering ${where}).`
      + (weak ? ` ${n(weak)} of them sit on low-quality bases (under Q${GOOD_QUALITY}, more than 1 chance in 100 of a`
        + " wrong call): check those peaks in the chromatogram before trusting them." : ""),
  };
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

  const summary = result && result.reads.length > 0 ? verdict(result) : null;

  return (
    <Dialog isOpen onClose={props.onClose} title="Align to Reference DNA Sequence"
            style={{ width: "92vw", maxWidth: 1500 }}>
      {/* Blueprint's Dialog does not forward unknown props, so the test id lives here. */}
      <div className={Classes.DIALOG_BODY} data-testid="sanger-dialog">
        <p style={{ marginTop: 0 }}>
          Checks a clone against the plasmid you expect: each sequencing read is trimmed of its
          low-quality ends, aligned on either strand, and every difference is listed in plasmid
          numbering.
        </p>
        {!props.active ? (
          <Callout intent="primary" data-testid="sanger-no-reference">
            Open the plasmid you expect (File → Open…) and keep its tab active, then open this
            dialog again: the reads are aligned to it.
          </Callout>
        ) : (
          <div className="vt-row vt-sanger-bar">
            <span>
              Reference: <b>{props.active.name}</b>
              {reference ? ` (${n(reference.sequence.length)} bp, ${reference.circular ? "circular" : "linear"})` : ""}
            </span>
            <label>
              Reads (.ab1)&nbsp;
              <input
                ref={input}
                type="file"
                multiple
                accept=".ab1,.abi,.fsa"
                data-testid="sanger-files"
                onChange={(e) => setFiles(Array.from(e.target.files ?? []))}
              />
            </label>
            <label title="How strictly low-quality bases are trimmed from the ends of each read">
              Trim&nbsp;
              <div className="bp3-html-select bp3-minimal">
                <select data-testid="sanger-trim" value={cutoff}
                        onChange={(e) => setCutoff(Number(e.target.value))}>
                  {TRIM_LEVELS.map((t) => <option key={t.value} value={t.value}>{t.label}</option>)}
                </select>
              </div>
            </label>
            {busy && (
              <span className="vt-muted" data-testid="sanger-busy">
                Aligning {files.length} read{files.length === 1 ? "" : "s"}…
              </span>
            )}
          </div>
        )}

        {error && <Callout intent="danger" title="Could not align">{error}</Callout>}
        {result?.warnings.map((w, i) => (
          <Callout key={i} intent="warning" data-testid="sanger-warning">{w}</Callout>
        ))}
        {result && result.reads.length === 0 && !error && (
          <Callout intent="warning" title="No read could be aligned" data-testid="sanger-none-aligned">
            The reasons are above. If a read was trimmed away, choose Permissive trimming; if it
            does not align, check that it comes from this plasmid and that the right tab is active.
          </Callout>
        )}

        {result && summary && (
          <>
            <Callout intent={summary.ok ? "success" : "warning"} data-testid="sanger-verdict">
              {summary.text}
            </Callout>
            <table className="bp3-html-table bp3-html-table-condensed vt-table" data-testid="sanger-summary">
              <thead>
                <tr>
                  <th>Read</th><th>Strand</th><th>Covers</th><th>Identity</th>
                  <th title="Mean Phred quality of the aligned bases: Q20 is 1 wrong call in 100, Q30 1 in 1,000">
                    Quality
                  </th>
                  <th>Differences (reference → read)</th><th>Trimmed</th>
                </tr>
              </thead>
              <tbody>
                {result.reads.map((r) => (
                  <tr key={r.file}>
                    <td title={`${r.file}${r.machine ? ` · ${r.machine}` : ""}`}>{r.name}</td>
                    <td>{r.strand === 1 ? "forward" : "reverse"}</td>
                    <td>{n(r.start)}…{n(r.end)}{r.spans_origin ? " (crosses the origin)" : ""}</td>
                    <td className={r.identity === 100 ? "vt-ok" : "vt-warn"}>{r.identity.toFixed(2)} %</td>
                    <td>{r.mean_quality !== undefined ? `Q${Math.round(r.mean_quality)}` : "—"}</td>
                    <td><Differences read={r} /></td>
                    <td className="vt-muted">{r.trimmed_bases} of {r.original_length} nt</td>
                  </tr>
                ))}
              </tbody>
            </table>
            <p className="vt-muted">
              Differences are listed with reference coordinates; hover one for the read quality
              behind it. In the view below, red marks are mismatches and grey marks are regions a
              read does not cover.
            </p>
            <div className="vt-alignment" data-testid="sanger-alignment">
              <AlignmentView id={ALIGNMENT_ID} height={420} />
            </div>
          </>
        )}

        {props.active && !result && !busy && !error && (
          <ol className="vt-muted" data-testid="sanger-steps">
            <li>Check that the reference above is the plasmid you expect.</li>
            <li>Choose the .ab1 files your sequencing provider sent; several at once is fine.</li>
            <li>Read the summary: every difference is listed with its position in the plasmid.</li>
          </ol>
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
