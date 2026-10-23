import { Button } from "@blueprintjs/core";
import type { SequenceData } from "./api";
import {
  type Contig,
  type ContigCall,
  coverageSentence,
  describeDifference,
  evidence,
  span,
} from "./contigText";
import "./contig.css";

const n = (value: number) => value.toLocaleString("en-US");

const STRANDS = [
  { strands: 0, label: "no read" },
  { strands: 1, label: "one strand" },
  { strands: 2, label: "both strands" },
];

/** One bar per run of coverage, the whole reference from left to right. */
function CoverageStrip(props: { contig: Contig; length: number }) {
  const { length } = props;
  const pieces = props.contig.coverage.flatMap((run) =>
    run.start <= run.end
      ? [{ ...run }]
      : [{ ...run, end: length }, { ...run, start: 1 }]);
  return (
    <div className="vt-contig-strip-wrap">
      <svg className="vt-contig-strip" viewBox={`0 0 ${length} 1`} preserveAspectRatio="none"
           role="img" aria-label={coverageSentence(props.contig)} data-testid="sanger-coverage-strip">
        {pieces.map((run, i) => (
          <rect key={i} x={run.start - 1} width={run.end - run.start + 1} y={0} height={1}
                className={`vt-contig-cover-${run.strands}`}>
            <title>{`${span(run.start, run.end)}: ${STRANDS[run.strands].label}`}</title>
          </rect>
        ))}
      </svg>
      <div className="vt-contig-scale vt-muted">
        <span>1</span>
        <span className="vt-contig-legend">
          {STRANDS.map((s) => (
            <span key={s.strands}><i className={`vt-contig-swatch vt-contig-cover-${s.strands}`} />{s.label}</span>
          ))}
        </span>
        <span>{n(length)}</span>
      </div>
    </div>
  );
}

function call(c: ContigCall): string {
  const shown = c.call === "-" ? "gap" : c.call === "" ? "nothing" : c.call;
  return `${c.read} (${c.strand === 1 ? "fwd" : "rev"}): ${shown} Q${c.quality}`;
}

/** The Contig section of the Sanger dialog: two or more reads, one consensus over the reference. */
export default function SangerContig(props: {
  contig: Contig;
  referenceLength: number;
  onSave?: (doc: SequenceData) => void;
}) {
  const { contig } = props;
  const disputed = contig.disagreements;
  return (
    <section className="vt-contig" data-testid="sanger-contig">
      <div className="vt-row vt-contig-head">
        <h4 className="vt-contig-title">Contig: {n(contig.reads)} reads, one consensus</h4>
        {props.onSave && contig.document && (
          <Button small icon="floppy-disk" data-testid="sanger-save-consensus"
                  title="Open the consensus as a new document: every difference is annotated where it is, the clearly called and checked stretches are annotated “Sanger verified”, and the places to check are marked"
                  onClick={() => props.onSave!(contig.document!)}>
            Save consensus
          </Button>
        )}
      </div>
      <p className="vt-contig-coverage" data-testid="sanger-coverage">{coverageSentence(contig)}</p>
      <CoverageStrip contig={contig} length={props.referenceLength} />

      {contig.differences.length > 0 && (
        <table className="bp3-html-table bp3-html-table-condensed vt-table" data-testid="sanger-contig-differences">
          <thead>
            <tr><th>Difference in the consensus (reference → clone)</th><th>Evidence</th><th>Verdict</th></tr>
          </thead>
          <tbody>
            {contig.differences.map((d, i) => (
              <tr key={i}>
                <td>{describeDifference(d)}</td>
                <td>{evidence(d)}</td>
                <td className={d.confirmed ? "vt-ok" : "vt-warn"}>{d.confirmed ? "confirmed" : "check the trace"}</td>
              </tr>
            ))}
          </tbody>
        </table>
      )}

      {disputed.length > 0 && (
        <details className="vt-contig-disputed" data-testid="sanger-disagreements" open={disputed.length <= 5}>
          <summary>
            {disputed.length === 1 ? "1 place" : `${n(disputed.length)} places`} to look at in the traces: the reads
            do not agree, make no clear call, or no read runs across a join
          </summary>
          <table className="bp3-html-table bp3-html-table-condensed vt-table">
            <thead><tr><th>Position</th><th>Reference</th><th>Consensus</th><th>Each read</th></tr></thead>
            <tbody>
              {disputed.map((d, i) => (
                <tr key={i}>
                  <td>{n(d.position)}{d.kind === "insertion" || d.kind === "not read across" ? " (before)" : ""}</td>
                  <td>{d.reference || "—"}</td>
                  <td>
                    {d.kind === "no call" ? "no clear call: reference kept"
                      : d.kind === "not read across" ? "no read runs across this join: not verified"
                        : d.consensus || "—"}
                  </td>
                  <td className="vt-muted">{d.calls.map(call).join(" · ")}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </details>
      )}
    </section>
  );
}
