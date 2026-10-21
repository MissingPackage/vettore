// Tools → Align Multiple Sequences: open documents aligned to a reference (star alignment in
// the engine, src/vettore/alignment.py), stacked in Open Vector Editor's alignment view the
// same way as Sanger reads, with a per-sequence list of differences in reference numbering.
import { useEffect, useMemo, useRef, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import { AlignmentView, addAlignment } from "@teselagen/ove";
import { call } from "../api";
import { store } from "../store";
import { useCopyText } from "../cloning/useCopyText";
import type { Tool, ToolProps } from "./types";
import {
  choicePrompt,
  choiceProblem,
  DEFAULT_LIMITS,
  defaultChoice,
  describeDifference,
  differencesAsTsv,
  distinctNames,
  ineligible,
  MODE_HELP,
  summarizeCounts,
  type AlignmentLimits,
  type Difference,
  type DifferenceCounts,
} from "./alignmentText";
import "./alignment.css";

// Blueprint 3 typings predate React 18 (no implicit `children`).
const { Dialog, Callout } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

const ALIGNMENT_ID = "vettore-multiple-alignment";
const SHOWN_DIFFERENCES = 300;

interface SequenceSummary {
  index: number;
  name: string;
  length: number;
  strand: 1 | -1;
  origin_moved_to: number | null;
  mode: "global" | "semi-global" | "local";
  start: number;
  end: number;
  spans_origin: boolean;
  coverage: number;
  aligned_length: number;
  identity: number;
  differences: Difference[];
  counts: DifferenceCounts;
}

interface AlignResult {
  alignmentTracks: unknown[];
  reference: { name: string; length: number; circular: boolean };
  sequences: SequenceSummary[];
  warnings: string[];
  first_covered: number;
}

const bp = (value: number) => value.toLocaleString("en-US");

function AlignMultipleSequences(props: ToolProps) {
  const { docs } = props;
  const [limits, setLimits] = useState<AlignmentLimits>(DEFAULT_LIMITS);
  const lengths = useMemo(
    () => Object.fromEntries(docs.map((d) => [d.editorName, d.data()?.sequence?.length ?? 0])),
    [docs],
  );
  // Several pasted sequences are all called "Pasted sequence": the tab number tells them apart.
  const names = useMemo(() => distinctNames(docs), [docs]);
  const [reference, setReference] = useState<string | null>(
    props.active?.editorName ?? docs[0]?.editorName ?? null,
  );
  const candidates = docs.map((d) => ({ editorName: d.editorName, length: lengths[d.editorName] }));
  const [chosen, setChosen] = useState<string[]>(() => defaultChoice(candidates, reference, limits));
  const [result, setResult] = useState<AlignResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [selected, setSelected] = useState(0);
  const [copied, setCopied] = useState(false);
  const { copy: copyText, fallbackDialog } = useCopyText();
  const request = useRef(0);

  useEffect(() => {
    call<AlignmentLimits>("alignment/limits").then(setLimits).catch(() => undefined);
  }, []);

  const referenceDoc = docs.find((d) => d.editorName === reference) ?? null;
  const others = docs.filter((d) => d.editorName !== reference);
  const referenceProblem = referenceDoc ? ineligible(lengths[referenceDoc.editorName], limits) : null;
  const problem = choiceProblem(chosen.length, limits);
  const prompt = choicePrompt(chosen.length);
  const ready = !!referenceDoc && !referenceProblem && !problem && chosen.length > 0;

  const changeReference = (editorName: string) => {
    setReference(editorName);
    setChosen(defaultChoice(candidates, editorName, limits));
    setResult(null);
  };

  const toggle = (editorName: string) => {
    setChosen((prev) => prev.includes(editorName)
      ? prev.filter((x) => x !== editorName)
      : [...prev, editorName]);
    setResult(null);
  };

  const run = async () => {
    if (!referenceDoc || !ready) return;
    const id = ++request.current;
    setBusy(true);
    setError(null);
    setResult(null);
    try {
      // Tab names, not the names stored in the files, so the user recognises each track.
      const chosenDocs = others.filter((d) => chosen.includes(d.editorName));
      const payload = await call<AlignResult>("alignment/align", {
        reference: { ...referenceDoc.data(), name: names[referenceDoc.editorName] },
        sequences: chosenDocs.map((d) => ({ ...d.data(), name: names[d.editorName] })),
      });
      if (id !== request.current) return;
      setResult(payload);
      setSelected(0);
      if (payload.alignmentTracks.length) {
        addAlignment(store, {
          id: ALIGNMENT_ID,
          name: `${payload.reference.name} — alignment`,
          alignmentType: "Alignment to reference",
          alignmentTracks: payload.alignmentTracks,
          // The reference's features come along but start hidden (eye menu of the view):
          // stacked feature rows would push the compared sequences out of sight.
          alignmentAnnotationLabelVisibility: { features: true },
          // Start where the sequences are, not at base 1.
          scrollPercentageToJumpTo:
            Math.max(0, payload.first_covered - 30) / Math.max(1, payload.reference.length),
        });
      }
    } catch (e) {
      if (id === request.current) setError(String((e as Error).message));
    } finally {
      if (id === request.current) setBusy(false);
    }
  };

  const row = result?.sequences[selected] ?? null;
  const copy = async () => {
    if (!row) return;
    setCopied(false);
    if (await copyText(differencesAsTsv(row.name, row.differences), `Differences of ${row.name}`)) {
      setCopied(true);
      setTimeout(() => setCopied(false), 2500);
    }
  };

  return (
    <Dialog isOpen onClose={props.onClose} title="Align Multiple Sequences"
            style={{ width: "92vw", maxWidth: 1500 }} canOutsideClickClose={false}>
      {/* Blueprint's Dialog does not forward unknown props, so the test id lives here. */}
      <div className={`${Classes.DIALOG_BODY} vt-msa`} data-testid="msa-dialog">
        <p className="vt-msa-intro">
          Lines up open sequences against a reference, on either strand, and lists every
          difference in the reference's numbering.
        </p>

        {docs.length < 2 ? (
          <Callout intent="primary" data-testid="msa-empty">
            Open at least two sequences to compare (File → Open…). The sequence in the active tab
            is used as the reference.
          </Callout>
        ) : (
          <div className="vt-msa-setup">
            <div className="vt-msa-field">
              <span className="vt-msa-label">Reference</span>
              <div className="bp3-html-select">
                <select data-testid="msa-reference" value={reference ?? ""}
                        onChange={(e) => changeReference(e.target.value)}>
                  {docs.map((d) => (
                    <option key={d.editorName} value={d.editorName}>
                      {names[d.editorName]} ({bp(lengths[d.editorName])} bp)
                    </option>
                  ))}
                </select>
              </div>
              {referenceProblem && (
                <span className="vt-msa-invalid" data-testid="msa-reference-problem">
                  This sequence is {referenceProblem}: choose another reference.
                </span>
              )}
            </div>

            <div className="vt-msa-field">
              <span className="vt-msa-label">Compare with</span>
              <div className="vt-msa-choices" data-testid="msa-choices">
                {others.map((d) => {
                  const why = ineligible(lengths[d.editorName], limits);
                  const circular = !!d.data()?.circular;
                  return (
                    <label key={d.editorName} className={`vt-msa-choice${why ? " is-disabled" : ""}`}>
                      <input
                        type="checkbox"
                        data-testid="msa-choice"
                        disabled={!!why}
                        checked={!why && chosen.includes(d.editorName)}
                        onChange={() => toggle(d.editorName)}
                      />
                      <span>{names[d.editorName]}</span>
                      <span className="vt-msa-meta">
                        {bp(lengths[d.editorName])} bp, {circular ? "circular" : "linear"}
                      </span>
                      {why && <span className="vt-msa-invalid">{why}</span>}
                    </label>
                  );
                })}
              </div>
              {problem && (
                <span className="vt-msa-invalid" data-testid="msa-choice-problem">{problem}</span>
              )}
              {prompt && <span className="vt-msa-meta" data-testid="msa-choice-prompt">{prompt}</span>}
            </div>

            <div className="vt-msa-actions">
              <Button intent="primary" data-testid="msa-align" disabled={!ready || busy}
                      loading={busy} onClick={run}>
                Align
              </Button>
              {busy ? (
                <span className="vt-msa-busy" data-testid="msa-busy">
                  Aligning {chosen.length} sequence{chosen.length === 1 ? "" : "s"} to{" "}
                  {referenceDoc && names[referenceDoc.editorName]}… long plasmids take a few seconds each.
                </span>
              ) : (
                <span className="vt-msa-meta">
                  Up to {limits.max_sequences} sequences at a time (reference included), each up
                  to {bp(limits.max_length)} bp. Both strands are tried; circular sequences are
                  turned to start where the reference starts.
                </span>
              )}
            </div>
          </div>
        )}

        {error && (
          <Callout intent="danger" title="Could not align" data-testid="msa-error">{error}</Callout>
        )}
        {result?.warnings.map((w, i) => (
          <Callout key={i} intent="warning" data-testid="msa-warning">{w}</Callout>
        ))}
        {result && result.sequences.length === 0 && !error && (
          <Callout intent="warning">
            None of the chosen sequences could be aligned to {result.reference.name}.
          </Callout>
        )}

        {result && result.sequences.length > 0 && (
          <>
            <table className="bp3-html-table bp3-html-table-condensed bp3-interactive vt-table"
                   data-testid="msa-summary">
              <thead>
                <tr>
                  <th>Sequence</th><th>Strand</th><th>Covers (reference bases)</th>
                  <th>Coverage</th><th>Identity</th><th>Differences</th>
                </tr>
              </thead>
              <tbody>
                {result.sequences.map((s, i) => (
                  <tr key={s.index} className={i === selected ? "vt-row-selected" : ""}
                      onClick={() => setSelected(i)} title="Click to list the differences below">
                    <td>
                      {s.name}
                      {s.mode === "local" && <span className="vt-msa-note"> · partial match</span>}
                      {s.origin_moved_to && (
                        <span className="vt-msa-note"> · origin moved to its base {bp(s.origin_moved_to)}</span>
                      )}
                    </td>
                    <td>{s.strand === 1 ? "forward" : "reverse complement"}</td>
                    <td>
                      {bp(s.start)}–{bp(s.end)}
                      {s.spans_origin ? " (crosses the origin)" : ""}
                    </td>
                    <td>{s.coverage.toFixed(1)} %</td>
                    <td className={s.identity === 100 ? "vt-ok" : "vt-warn"}>{s.identity.toFixed(2)} %</td>
                    <td>{summarizeCounts(s.counts)}</td>
                  </tr>
                ))}
              </tbody>
            </table>

            {row && (
              <div className="vt-msa-details" data-testid="msa-details">
                <div className="vt-msa-details-head">
                  <b>{row.name}</b>
                  <span className="vt-msa-meta">
                    {MODE_HELP[row.mode]} {bp(row.aligned_length)} of its {bp(row.length)} bp are
                    aligned. Bases are given on the reference strand.
                  </span>
                  {row.differences.length > 0 && (
                    <Button small minimal icon={copied ? undefined : "duplicate"} className={copied ? "vt-done" : undefined}
                            data-testid="msa-copy" onClick={copy}>
                      {copied ? "Copied: paste it into a spreadsheet" : "Copy list (for a spreadsheet)"}
                    </Button>
                  )}
                </div>
                {row.differences.length === 0 ? (
                  <p className="vt-ok">Identical to {result.reference.name} over the aligned region.</p>
                ) : (
                  <ul className="vt-msa-diffs" data-testid="msa-differences">
                    {row.differences.slice(0, SHOWN_DIFFERENCES).map((d, i) => (
                      <li key={i} className={`vt-msa-diff is-${d.kind}`}>{describeDifference(d)}</li>
                    ))}
                    {row.differences.length > SHOWN_DIFFERENCES && (
                      <li className="vt-msa-meta">
                        … and {bp(row.differences.length - SHOWN_DIFFERENCES)} more; copy the list
                        to see them all.
                      </li>
                    )}
                  </ul>
                )}
              </div>
            )}

            <p className="vt-msa-meta">
              Below, the reference is on top; coloured marks on the other tracks are differences
              from it, and gaps (–) are bases a sequence lacks or regions it does not cover.
            </p>
            <div className="vt-msa-view" data-testid="msa-alignment">
              <AlignmentView id={ALIGNMENT_ID} height={420} />
            </div>
          </>
        )}
      </div>
      {fallbackDialog}
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button onClick={props.onClose}>Close</Button>
        </div>
      </div>
    </Dialog>
  );
}

export const tools: Tool[] = [
  {
    id: "align-multiple",
    menu: "Tools",
    label: "Align Multiple Sequences",
    requiresDocument: true,
    Component: AlignMultipleSequences,
  },
];
