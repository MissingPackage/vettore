import { useEffect, useMemo, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import { api } from "./api";
import { loadPreferences } from "./preferences";
import type { OpenDoc } from "./cloning/types";
import { EnzymePicker, Field } from "./cloning/inputs";
import "./gel.css";

const { Dialog, Callout } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

interface GelBand {
  size: number;
  apparent_size: number;
  rf: number;
  intensity: number;
  label: string;
  reference?: boolean;
}

interface GelResult {
  /** `note`: why a lane is not what was asked (an enzyme that does not cut, an empty file). */
  lanes: { name: string; bands: GelBand[]; sizes: number[]; note?: string | null }[];
  ladder: GelBand[];
  ladder_label: string;
  percent: number;
}

interface LaneSetup {
  editorName: string | null;
  enzymes: string[];
}

const GEL_HEIGHT = 440;
const LANE_WIDTH = 104;
const MARGIN_LEFT = 62;
const TITLE_Y = 14;
/** Lane headers take up to three lines of 12 px above the wells. */
const LINE_HEIGHT = 12;
const HEADER_LINES = 3;
const WELL_Y = TITLE_Y + HEADER_LINES * LINE_HEIGHT + 6;
/** Characters of a 10 px header line that fit over one lane. */
const HEADER_CHARS = 17;
/** Ladder sizes closer than this (px) would overprint: only one of them is written. */
const LABEL_GAP = 10;

/** Words over at most `lines` lines of `max` characters; what does not fit ends in "…". */
function wrapWords(text: string, max: number, lines: number): string[] {
  const out: string[] = [];
  let line = "";
  for (const word of text.split(/\s+/).filter(Boolean)) {
    if (!line) line = word;
    else if (`${line} ${word}`.length <= max) line = `${line} ${word}`;
    else {
      out.push(line);
      line = word;
    }
  }
  if (line) out.push(line);
  const cut = (l: string) => (l.length > max ? `${l.slice(0, max - 1)}…` : l);
  const kept = out.slice(0, lines).map(cut);
  if (out.length > lines) kept[lines - 1] = `${kept[lines - 1].slice(0, max - 1)}…`;
  return kept;
}

/** Which ladder bands get their size written: reference bands first, then any band far
 * enough from every size already written. */
function labelledBands(ys: number[], reference: boolean[]): boolean[] {
  const keep = ys.map(() => false);
  const placed: number[] = [];
  const order = ys.map((_, i) => i).sort((a, b) => Number(reference[b]) - Number(reference[a]) || a - b);
  for (const i of order) {
    if (placed.every((y) => Math.abs(y - ys[i]) >= LABEL_GAP)) {
      keep[i] = true;
      placed.push(ys[i]);
    }
  }
  return keep;
}

/** Tools → Simulate Agarose Gel: lanes are digests of the open documents. */
export default function GelDialog(props: { docs: OpenDoc[]; active: OpenDoc | null; onClose: () => void }) {
  // The document in front is lane 1, uncut: that already draws a gel at the first click.
  const [lanes, setLanes] = useState<LaneSetup[]>([
    { editorName: props.active?.editorName ?? null, enzymes: [] },
  ]);
  // Ladder and agarose come from Preferences; the engine's defaults are the fallback.
  const preferences = loadPreferences();
  const [ladder, setLadder] = useState(preferences.gelLadder);
  const [percent, setPercent] = useState(preferences.gelPercent);
  const [ladders, setLadders] = useState<{ id: string; label: string }[]>([]);
  const [percentages, setPercentages] = useState<number[]>([0.7, 1.0, 1.2, 1.5, 2.0]);
  const [result, setResult] = useState<GelResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  useEffect(() => {
    api.gelLadders()
      .then((r) => {
        setLadders(r.ladders);
        setPercentages(r.percentages);
        if (!r.ladders.some((l) => l.id === preferences.gelLadder)) setLadder(r.default_ladder);
        if (!r.percentages.includes(preferences.gelPercent)) setPercent(r.default_percent);
      })
      .catch(() => undefined);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const body = useMemo(() => {
    const filled = lanes.filter((l) => l.editorName);
    if (!filled.length) return null;
    return {
      lanes: filled.map((l) => {
        const doc = props.docs.find((d) => d.editorName === l.editorName)!;
        return { name: doc.name, sequenceData: doc.data(), enzymes: l.enzymes };
      }),
      ladder,
      percent,
    };
  }, [lanes, ladder, percent, props.docs]);
  const bodyKey = body ? JSON.stringify(body) : "";

  useEffect(() => {
    if (!bodyKey) {
      setResult(null);
      return;
    }
    setError(null);
    const payload = JSON.parse(bodyKey);
    const t = setTimeout(() => {
      setBusy(true);
      api.gel<GelResult>(payload)
        .then(setResult)
        .catch((e) => {
          setResult(null);
          setError(String((e as Error).message));
        })
        .finally(() => setBusy(false));
    }, 250);
    return () => clearTimeout(t);
  }, [bodyKey]);

  const update = (i: number, patch: Partial<LaneSetup>) =>
    setLanes((rows) => rows.map((r, j) => (j === i ? { ...r, ...patch } : r)));

  const width = MARGIN_LEFT + ((result?.lanes.length ?? 1) + 1) * LANE_WIDTH + 14;

  return (
    <Dialog isOpen onClose={props.onClose} title="Simulate Agarose Gel" style={{ width: 1050, maxWidth: "95vw" }}>
      <div className={`${Classes.DIALOG_BODY} vt-cloning`} data-testid="gel-dialog">
        <div className="vt-cloning-form">
          <p className="vt-muted">
            Predicts the pattern of a diagnostic digest: one lane per open file, next to the
            ladder you run. Leave a lane without enzymes to load the plasmid uncut.
          </p>
          {lanes.map((lane, i) => {
            const doc = props.docs.find((d) => d.editorName === lane.editorName) ?? null;
            const data = doc?.data();
            return (
              <div key={i} className="vt-insert-row">
                <Field label={`Lane ${i + 1}`}>
                  <div className="bp3-html-select bp3-fill">
                    <select
                      data-testid={`gel-doc-${i}`}
                      value={lane.editorName ?? ""}
                      onChange={(e) => update(i, { editorName: e.target.value })}
                    >
                      <option value="" disabled>Choose an open file…</option>
                      {props.docs.map((d) => (
                        <option key={d.editorName} value={d.editorName}>{d.name}</option>
                      ))}
                    </select>
                  </div>
                </Field>
                <Field
                  label="Digest with"
                  hint={doc
                    ? "Up to three enzymes. Leave empty to load the uncut plasmid."
                    : "Choose a file first."}
                >
                  <EnzymePicker
                    value={lane.enzymes}
                    onChange={(enzymes) => update(i, { enzymes })}
                    max={3}
                    sequence={data?.sequence}
                    circular={data?.circular}
                    testId={`gel-enzymes-${i}`}
                  />
                </Field>
                {lanes.length > 1 && (
                  <button className="bp3-button bp3-small bp3-minimal"
                          onClick={() => setLanes((r) => r.filter((_, j) => j !== i))}>
                    Remove lane
                  </button>
                )}
              </div>
            );
          })}
          <button className="bp3-button bp3-small" data-testid="gel-add-lane"
                  onClick={() => setLanes((r) => [...r, { editorName: null, enzymes: [] }])}>
            + Add lane
          </button>
          <Field label="Ladder" hint="Your usual ladder: File → Preferences keeps the choice.">
            <div className="bp3-html-select bp3-fill">
              <select data-testid="gel-ladder" value={ladder} onChange={(e) => setLadder(e.target.value)}>
                {ladders.map((l) => <option key={l.id} value={l.id}>{l.label}</option>)}
              </select>
            </div>
          </Field>
          <Field label="Agarose (%)" hint="Higher percentages separate short fragments better.">
            <div className="bp3-html-select bp3-fill">
              <select data-testid="gel-percent" value={percent} onChange={(e) => setPercent(Number(e.target.value))}>
                {percentages.map((p) => <option key={p} value={p}>{p.toFixed(1)} %</option>)}
              </select>
            </div>
          </Field>
        </div>
        <div className="vt-cloning-result">
          {error && (
            <Callout intent="danger" title="No gel to draw" data-testid="gel-error">
              {error}
            </Callout>
          )}
          {!result && !error && !busy && (
            <p className="vt-muted" data-testid="gel-empty">
              Choose a file for lane 1 and the gel appears here.
            </p>
          )}
          {busy && !result && <p className="vt-muted">Running the gel…</p>}
          {result && (
            <>
              <svg width={width} height={GEL_HEIGHT + 16} className="vt-gel" data-testid="gel-image">
                <rect x={MARGIN_LEFT - 4} y={WELL_Y - 14} width={width - MARGIN_LEFT + 4}
                      height={GEL_HEIGHT - WELL_Y + 24} fill="#1b1b1b" rx={4} />
                <Lane x={MARGIN_LEFT} label={result.ladder_label} bands={result.ladder} isLadder />
                {result.lanes.map((lane, i) => (
                  <Lane key={i} x={MARGIN_LEFT + (i + 1) * LANE_WIDTH} label={lane.name}
                        bands={lane.bands} />
                ))}
              </svg>
              <p className="vt-muted">
                {result.percent.toFixed(1)} % agarose. Hover a band to read its size. Positions
                are an approximation (semi-log migration), meant for reading a digest pattern
                against the ladder rather than for measuring a real run.
              </p>
              {result.lanes.some((lane) => lane.note) && (
                <Callout intent="warning" data-testid="gel-notes">
                  {result.lanes.map((lane, i) => lane.note && (
                    <p key={i} data-testid="gel-note" style={{ margin: 0 }}>
                      <b>Lane {i + 1}:</b> {lane.note}
                    </p>
                  ))}
                </Callout>
              )}
              {result.lanes.map((lane, i) => (
                <p key={i} className="vt-muted">
                  <b>{lane.name}</b>: {lane.sizes.length
                    ? `${lane.sizes.map((s) => s.toLocaleString("en-US")).join(", ")} bp`
                    : "no bands"}
                </p>
              ))}
            </>
          )}
        </div>
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button onClick={props.onClose}>Close</Button>
        </div>
      </div>
    </Dialog>
  );
}

function Lane(props: { x: number; label: string; bands: GelBand[]; isLadder?: boolean }) {
  const width = LANE_WIDTH - 20;
  // "pUC19 + GFP (Gibson) / EcoRI + BamHI": the file on up to two lines, the enzymes in grey.
  const [first, ...rest] = props.label.split(" / ");
  const enzymes = rest.join(" / ");
  const nameLines = wrapWords(first, HEADER_CHARS, enzymes ? HEADER_LINES - 1 : HEADER_LINES);
  const lines = [
    ...nameLines.map((text) => ({ text, muted: false })),
    ...(enzymes ? wrapWords(enzymes, HEADER_CHARS, 1).map((text) => ({ text, muted: true })) : []),
  ];
  const ys = props.bands.map((band) => WELL_Y + band.rf * (GEL_HEIGHT - WELL_Y - 10));
  const labelled = props.isLadder ? labelledBands(ys, props.bands.map((b) => !!b.reference)) : [];
  return (
    <g>
      <text x={props.x + width / 2} y={TITLE_Y} fontSize={10} textAnchor="middle" fill="#374151"
            data-testid="gel-lane-header">
        <title>{props.label}</title>
        {lines.map((line, i) => (
          <tspan key={i} x={props.x + width / 2} dy={i === 0 ? 0 : LINE_HEIGHT}
                 fill={line.muted ? "#6b7280" : undefined}>
            {line.text}
          </tspan>
        ))}
      </text>
      <rect x={props.x} y={WELL_Y - 8} width={width} height={7} fill="#3f3f46" rx={1} />
      {props.bands.map((band, i) => {
        const y = ys[i];
        // Placed with a transform, not y, so a band slides to its new size when the enzymes
        // change and the eye sees which fragment split (gel.css; instant with reduced motion).
        return (
          <g key={i}>
            <rect className="vt-gel-band" x={props.x} y={0} width={width} height={band.reference ? 5 : 4} rx={1}
                  fill="#e8f1c8" style={{ transform: `translateY(${y}px)`, opacity: 0.2 + 0.8 * band.intensity }}>
              <title>{band.label}</title>
            </rect>
            {props.isLadder && labelled[i] && (
              <text x={props.x - 10} y={y + 5} fill="#9ca3af" fontSize={9} textAnchor="end"
                    data-testid="gel-ladder-label">
                {band.label}
              </text>
            )}
          </g>
        );
      })}
    </g>
  );
}
