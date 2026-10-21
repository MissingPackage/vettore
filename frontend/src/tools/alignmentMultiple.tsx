// "All to each other" in Tools → Align Multiple Sequences: a progressive alignment of the
// chosen documents (src/vettore/msa.py), run as a background job with its progress shown,
// stacked under a consensus row with a conservation bar, and the differences between any two
// sequences the user picks.
import { useEffect, useMemo, useRef, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, ProgressBar } from "@blueprintjs/core";
import { AlignmentView, addAlignment } from "@teselagen/ove";
import { call } from "../api";
import { store } from "../store";
import { useCopyText } from "../cloning/useCopyText";
import type { ToolProps } from "./types";
import { describeDifference, differencesAsTsv, summarizeCounts } from "./alignmentText";
import {
  conservationBins,
  conservationLevel,
  DEFAULT_MULTIPLE_LIMITS,
  defaultMultipleChoice,
  describeBin,
  ineligibleForMultiple,
  multipleChoiceProblem,
  multiplePrompt,
  orientationNote,
  progressText,
  variableRegions,
  type Comparison,
  type MultipleLimits,
  type MultipleResult,
} from "./multipleText";

// Blueprint 3 typings predate React 18 (no implicit `children`).
const { Callout } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

const ALIGNMENT_ID = "vettore-progressive-alignment";
const SHOWN_DIFFERENCES = 300;
const POLL_MS = 300;
const BAR_BINS = 400;

type JobState =
  | { state: "running"; progress: number; message: string }
  | { state: "done"; result: MultipleResult }
  | { state: "cancelled" };

const bp = (value: number) => value.toLocaleString("en-US");
const wait = (ms: number) => new Promise((resolve) => setTimeout(resolve, ms));

interface Props extends ToolProps {
  names: Record<string, string>;
  lengths: Record<string, number>;
}

export function AllToEachOther({ docs, active, names, lengths }: Props) {
  const [limits, setLimits] = useState<MultipleLimits>(DEFAULT_MULTIPLE_LIMITS);
  const candidates = docs.map((d) => ({ editorName: d.editorName, length: lengths[d.editorName] }));
  const [chosen, setChosen] = useState<string[]>(
    () => defaultMultipleChoice(candidates, active?.editorName ?? null, limits));
  const [progress, setProgress] = useState<{ progress: number; message: string } | null>(null);
  const [result, setResult] = useState<MultipleResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const job = useRef<string | null>(null);
  const request = useRef(0);

  useEffect(() => {
    call<MultipleLimits>("alignment/multiple/limits").then(setLimits).catch(() => undefined);
  }, []);

  // Closing the dialog stops an alignment still running.
  useEffect(() => () => { void cancelJob(); }, []);

  const cancelJob = async () => {
    const id = job.current;
    job.current = null;
    request.current++;
    if (id) await call(`alignment/multiple/job/${id}/cancel`, {}).catch(() => undefined);
  };

  const ticked = docs.filter((d) => chosen.includes(d.editorName)
    && !ineligibleForMultiple(lengths[d.editorName], limits));
  const problem = multipleChoiceProblem(ticked.map((d) => lengths[d.editorName]), limits);
  const prompt = multiplePrompt(ticked.length);
  const busy = progress !== null;
  const ready = ticked.length >= 2 && !problem;

  const toggle = (editorName: string) => {
    setChosen((prev) => prev.includes(editorName)
      ? prev.filter((x) => x !== editorName)
      : [...prev, editorName]);
    setResult(null);
  };

  const run = async () => {
    if (!ready) return;
    await cancelJob();
    const id = ++request.current;
    setError(null);
    setResult(null);
    setProgress({ progress: 0, message: "Starting" });
    try {
      // Tab names, not the names stored in the files, so the user recognises each track.
      const { job: jobId } = await call<{ job: string }>("alignment/multiple/start", {
        sequences: ticked.map((d) => ({ ...d.data(), name: names[d.editorName] })),
      });
      // Stopped, closed or started again while the request was on its way: stop it too.
      if (id !== request.current) {
        await call(`alignment/multiple/job/${jobId}/cancel`, {}).catch(() => undefined);
        return;
      }
      job.current = jobId;
      for (;;) {
        await wait(POLL_MS);
        if (id !== request.current) return;
        const state = await call<JobState>(`alignment/multiple/job/${jobId}`);
        if (id !== request.current) return;
        if (state.state === "running") {
          setProgress({ progress: state.progress, message: state.message });
          continue;
        }
        job.current = null;
        if (state.state === "done") show(state.result);
        return;
      }
    } catch (e) {
      if (id === request.current) setError(String((e as Error).message));
    } finally {
      if (id === request.current) setProgress(null);
    }
  };

  const show = (payload: MultipleResult) => {
    setResult(payload);
    if (payload.alignmentTracks.length) {
      addAlignment(store, {
        id: ALIGNMENT_ID,
        name: "Multiple alignment",
        alignmentType: "Multiple Sequence Alignment",
        alignmentTracks: payload.alignmentTracks,
      });
    }
  };

  const stop = async () => {
    await cancelJob();
    setProgress(null);
  };

  return (
    <>
      <div className="vt-msa-setup">
        <div className="vt-msa-field">
          <span className="vt-msa-label">Sequences</span>
          <div className="vt-msa-choices" data-testid="msa-all-choices">
            {docs.map((d) => {
              const why = ineligibleForMultiple(lengths[d.editorName], limits);
              const circular = !!d.data()?.circular;
              return (
                <label key={d.editorName} className={`vt-msa-choice${why ? " is-disabled" : ""}`}>
                  <input
                    type="checkbox"
                    data-testid="msa-all-choice"
                    disabled={!!why || busy}
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
          {problem && <span className="vt-msa-invalid" data-testid="msa-all-problem">{problem}</span>}
          {prompt && <span className="vt-msa-meta" data-testid="msa-all-prompt">{prompt}</span>}
        </div>

        <div className="vt-msa-actions">
          <Button intent="primary" data-testid="msa-all-align" disabled={!ready || busy}
                  onClick={run}>
            Align
          </Button>
          {progress ? (
            <div className="vt-msa-progress" data-testid="msa-progress">
              <ProgressBar value={Math.max(0.02, progress.progress)} stripes={false} />
              <span className="vt-msa-busy">{progressText(progress.message, progress.progress)}</span>
              <Button small minimal data-testid="msa-stop" onClick={stop}>Stop</Button>
            </div>
          ) : (
            <span className="vt-msa-meta">
              Up to {limits.max_sequences} sequences, each up to {bp(limits.max_length)} bp
              and {bp(limits.max_total)} bp in all. Strand and origin follow the first ticked
              circular sequence (else the first ticked); similar sequences end up next to each
              other.
            </span>
          )}
        </div>
      </div>

      {error && (
        <Callout intent="danger" title="Could not align" data-testid="msa-error">{error}</Callout>
      )}
      {result?.warnings.map((w, i) => (
        <Callout key={i} intent="warning" data-testid="msa-warning">{w}</Callout>
      ))}
      {result && result.sequences.length > 0 && <Results result={result} />}
    </>
  );
}

function Results({ result }: { result: MultipleResult }) {
  return (
    <>
      <table className="bp3-html-table bp3-html-table-condensed vt-table" data-testid="msa-all-summary">
        <thead>
          <tr>
            <th>Sequence</th><th>Length</th><th>Strand</th><th>Columns</th>
            <th title="Identical columns over the columns it covers, gaps included">Identity to the consensus</th>
          </tr>
        </thead>
        <tbody>
          {result.sequences.map((s) => (
            <tr key={s.index}>
              <td>{s.name}</td>
              <td>{bp(s.length)} bp</td>
              <td>{orientationNote(s)}</td>
              <td>{bp(s.first_column)}–{bp(s.last_column)} of {bp(result.columns)}</td>
              <td className={s.identity === 100 ? "vt-ok" : undefined}>{s.identity.toFixed(2)} %</td>
            </tr>
          ))}
        </tbody>
      </table>

      <ConservationBar conservation={result.conservation} />

      <p className="vt-msa-meta">
        Below, the consensus is on top (the most common base of each column; an ambiguity code
        where two are equally common); coloured marks on the sequences are differences from it,
        and gaps (–) are bases a sequence lacks or regions it does not cover.
      </p>
      <div className="vt-msa-view" data-testid="msa-alignment">
        <AlignmentView id={ALIGNMENT_ID} height={420} />
      </div>

      <PairDifferences rows={result.sequences} />
    </>
  );
}

function ConservationBar({ conservation }: { conservation: (number | null)[] }) {
  const bins = useMemo(() => conservationBins(conservation, BAR_BINS), [conservation]);
  const regions = useMemo(() => variableRegions(conservation), [conservation]);
  const counted = conservation.filter((c): c is number => c !== null);
  const mean = counted.reduce((a, c) => a + c, 0) / Math.max(1, counted.length);
  return (
    <div className="vt-msa-field" data-testid="msa-conservation">
      <span className="vt-msa-label">
        Conservation <span className="vt-msa-meta">
          (share of the sequences with each column's most common base, {Math.round(mean)} % on
          average; gaps never count as agreement; hover a bar for its columns)
        </span>
      </span>
      <svg className="vt-msa-conservation" viewBox={`0 0 ${bins.length} 100`}
           preserveAspectRatio="none" shapeRendering="crispEdges" role="img"
           aria-label={`Conservation along ${bp(conservation.length)} columns`}>
        {bins.map((b, i) => (
          <rect key={i} x={i} width={1} y={100 - Math.max(2, b.mean ?? 100)}
                height={Math.max(2, b.mean ?? 100)}
                className={`is-${conservationLevel(b.mean)}`}>
            <title>{describeBin(b)}</title>
          </rect>
        ))}
      </svg>
      <span className="vt-msa-meta" data-testid="msa-variable">
        {regions.length
          ? `Most variable stretches (columns): ${regions.map((r) => `${bp(r.start)}–${bp(r.end)}`).join(", ")}.`
          : "No stretch of 10 or more columns where fewer than 60 % of the sequences agree."}
      </span>
    </div>
  );
}

function PairDifferences({ rows }: { rows: MultipleResult["sequences"] }) {
  const [a, setA] = useState(0);
  const [b, setB] = useState(Math.min(1, rows.length - 1));
  const [comparison, setComparison] = useState<Comparison | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [copied, setCopied] = useState(false);
  const { copy: copyText, fallbackDialog } = useCopyText();

  useEffect(() => {
    let current = true;
    setComparison(null);
    setError(null);
    if (a === b) return;
    call<Comparison>("alignment/multiple/compare", { a: rows[a], b: rows[b] })
      .then((c) => { if (current) setComparison(c); })
      .catch((e) => { if (current) setError(String((e as Error).message)); });
    return () => { current = false; };
  }, [a, b, rows]);

  const copy = async () => {
    if (!comparison) return;
    setCopied(false);
    if (await copyText(differencesAsTsv(rows[b].name, comparison.differences),
                       `Differences of ${rows[b].name} from ${rows[a].name}`)) {
      setCopied(true);
      setTimeout(() => setCopied(false), 2500);
    }
  };

  const pick = (value: number, set: (v: number) => void, testId: string) => (
    <div className="bp3-html-select">
      <select data-testid={testId} value={value} onChange={(e) => set(Number(e.target.value))}>
        {rows.map((r, i) => <option key={r.index} value={i}>{r.name}</option>)}
      </select>
    </div>
  );

  return (
    <div className="vt-msa-details" data-testid="msa-pair">
      <div className="vt-msa-details-head">
        <b>Differences between two sequences</b>
        <span className="vt-msa-pair-pick">
          {pick(b, setB, "msa-pair-b")} compared with {pick(a, setA, "msa-pair-a")}
        </span>
        {comparison && comparison.differences.length > 0 && (
          <Button small minimal icon={copied ? undefined : "duplicate"} className={copied ? "vt-done" : undefined}
                  data-testid="msa-pair-copy" onClick={copy}>
            {copied ? "Copied: paste it into a spreadsheet" : "Copy list (for a spreadsheet)"}
          </Button>
        )}
      </div>
      {a === b && <p className="vt-msa-meta">Choose two different sequences.</p>}
      {error && <p className="vt-msa-invalid">{error}</p>}
      {comparison && (
        comparison.overlap === 0 ? (
          <p className="vt-msa-meta">These two do not overlap in the alignment.</p>
        ) : (
          <>
            <p className="vt-msa-meta" data-testid="msa-pair-summary">
              Over {rows[a].name}'s bases {bp(comparison.start ?? 0)}–{bp(comparison.end ?? 0)}:
              {" "}{comparison.identity.toFixed(2)} % identical; {summarizeCounts(comparison.counts)}.
              Positions and bases are given in {rows[a].name}'s own numbering and strand.
              {comparison.unknown > 0 && ` ${bp(comparison.unknown)} position${comparison.unknown === 1 ? "" : "s"} with N (unknown base) are not counted as differences.`}
            </p>
            {comparison.unrelated.length > 0 && (
              <Callout intent="warning" data-testid="msa-pair-unrelated">
                Not compared, because the two have nothing recognisably in common there (for
                example different inserts in the same backbone):{" "}
                {comparison.unrelated.map((u) => u.start === null
                  ? `${bp(u.b_bases)} bp of ${rows[b].name} with no counterpart`
                  : `${rows[a].name} ${bp(u.start)}–${bp(u.end ?? u.start)} against ${bp(u.b_bases)} bp of ${rows[b].name}`).join("; ")}.
              </Callout>
            )}
            {comparison.differences.length === 0 ? (
              <p className="vt-ok">Identical where they overlap.</p>
            ) : (
              <ul className="vt-msa-diffs" data-testid="msa-pair-differences">
                {comparison.differences.slice(0, SHOWN_DIFFERENCES).map((d, i) => (
                  <li key={i} className={`vt-msa-diff is-${d.kind}`}>{describeDifference(d)}</li>
                ))}
                {comparison.differences.length > SHOWN_DIFFERENCES && (
                  <li className="vt-msa-meta">
                    … and {bp(comparison.differences.length - SHOWN_DIFFERENCES)} more; copy the
                    list to see them all.
                  </li>
                )}
              </ul>
            )}
          </>
        )
      )}
      {fallbackDialog}
    </div>
  );
}
