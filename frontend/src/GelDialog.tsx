import { useEffect, useMemo, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import { api } from "./api";
import { cutsites, ENZYMES } from "./enzymeLogic";
import type { OpenDoc } from "./cloning/types";
import { EnzymePicker, Field } from "./cloning/inputs";

const { Dialog } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

interface GelBand {
  size: number;
  apparent_size: number;
  rf: number;
  intensity: number;
  label: string;
  reference?: boolean;
}

interface GelResult {
  lanes: { name: string; bands: GelBand[]; sizes: number[] }[];
  ladder: GelBand[];
  ladder_label: string;
  percent: number;
}

interface LaneSetup {
  editorName: string | null;
  enzymes: string[];
}

const GEL_HEIGHT = 430;
const LANE_WIDTH = 104;
const MARGIN_LEFT = 62;
const TITLE_Y = 14;
const WELL_Y = 40;

/** Tools → Simulate Agarose Gel: lanes are digests of the open documents. */
export default function GelDialog(props: { docs: OpenDoc[]; active: OpenDoc | null; onClose: () => void }) {
  const [lanes, setLanes] = useState<LaneSetup[]>([
    { editorName: props.active?.editorName ?? null, enzymes: [] },
  ]);
  const [ladder, setLadder] = useState("neb-1kb-plus");
  const [percent, setPercent] = useState(1.0);
  const [ladders, setLadders] = useState<{ id: string; label: string }[]>([]);
  const [percentages, setPercentages] = useState<number[]>([0.7, 1.0, 1.2, 1.5, 2.0]);
  const [result, setResult] = useState<GelResult | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    api.gelLadders()
      .then((r) => {
        setLadders(r.ladders);
        setPercentages(r.percentages);
        setLadder(r.default_ladder);
        setPercent(r.default_percent);
      })
      .catch(() => undefined);
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
      api.gel<GelResult>(payload)
        .then(setResult)
        .catch((e) => {
          setResult(null);
          setError(String((e as Error).message));
        });
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
                <Field label="Digest with" hint="Leave empty to load the uncut plasmid">
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
          <Field label="Ladder">
            <div className="bp3-html-select bp3-fill">
              <select data-testid="gel-ladder" value={ladder} onChange={(e) => setLadder(e.target.value)}>
                {ladders.map((l) => <option key={l.id} value={l.id}>{l.label}</option>)}
              </select>
            </div>
          </Field>
          <Field label="Agarose (%)">
            <div className="bp3-html-select bp3-fill">
              <select data-testid="gel-percent" value={percent} onChange={(e) => setPercent(Number(e.target.value))}>
                {percentages.map((p) => <option key={p} value={p}>{p.toFixed(1)} %</option>)}
              </select>
            </div>
          </Field>
        </div>
        <div className="vt-cloning-result">
          {error && <p className="vt-error">{error}</p>}
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
                {result.percent.toFixed(1)} % agarose. Positions are an approximation
                (semi-log migration), meant for reading a digest pattern against the ladder.
              </p>
              {result.lanes.map((lane, i) => (
                <p key={i} className="vt-muted">
                  <b>{lane.name}</b>: {lane.sizes.map((s) => s.toLocaleString()).join(", ")} bp
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
  const [first, ...rest] = props.label.split(" / ");
  return (
    <g>
      <text x={props.x + width / 2} y={TITLE_Y} fontSize={10} textAnchor="middle" fill="#374151">
        {first.length > 16 ? `${first.slice(0, 15)}…` : first}
        <tspan x={props.x + width / 2} dy={12} fill="#6b7280">
          {rest.join(" / ").slice(0, 18)}
        </tspan>
      </text>
      <rect x={props.x} y={WELL_Y - 8} width={width} height={7} fill="#3f3f46" rx={1} />
      {props.bands.map((band, i) => {
        const y = WELL_Y + band.rf * (GEL_HEIGHT - WELL_Y - 10);
        return (
          <g key={i}>
            <rect x={props.x} y={y} width={width} height={band.reference ? 5 : 4} rx={1}
                  fill="#e8f1c8" opacity={0.2 + 0.8 * band.intensity}>
              <title>{band.label}</title>
            </rect>
            {props.isLadder && (
              <text x={props.x - 10} y={y + 5} fill="#9ca3af" fontSize={9} textAnchor="end">
                {band.label}
              </text>
            )}
          </g>
        );
      })}
    </g>
  );
}
