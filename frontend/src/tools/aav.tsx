// Tools → AAV Packaging Check: how much DNA an AAV capsid would have to carry.
// The sizes and their literature sources are in aavPackaging.ts.
import { useCallback, useMemo, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import { api } from "../api";
import { Field } from "../cloning/inputs";
import {
  asPosition,
  cargoOptions,
  chooseCargo,
  findItrs,
  limitFor,
  SC_LIMIT_BP,
  SS_HARD_LIMIT_BP,
  SS_LIMIT_BP,
  WILD_TYPE_BP,
  verdictFor,
  type Cargo,
  type Itr,
  type VectorKind,
} from "./aavPackaging";
import type { Tool, ToolProps } from "./types";
import "./aav.css";

// Blueprint 3 typings predate React 18 (no implicit `children`).
const { Dialog, Callout } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

interface DetectedFeature {
  name: string;
  type: string;
  start: number;
  end: number;
  identity: number;
  source: string;
  note: string;
}

const INTENT: Record<string, string> = { ok: "success", tight: "warning", over: "danger" };

function bp(value: number): string {
  return `${value.toLocaleString()} bp`;
}

/** A position as a biologist writes it: 1-based, with thousands separated. */
function at(index: number): string {
  return asPosition(index).toLocaleString();
}

function AavPackagingDialog(props: ToolProps) {
  const doc = props.active;
  const data = doc?.data();
  const sequenceLength = data?.sequence.length ?? 0;

  const annotated = useMemo(
    () => (data ? findItrs(data.features, sequenceLength) : []),
    [data, sequenceLength],
  );
  const [detected, setDetected] = useState<Itr[] | null>(null);
  const [searching, setSearching] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [kind, setKind] = useState<VectorKind>("ss");
  const [pair, setPair] = useState<[number, number] | null>(null);
  const [sideIndex, setSideIndex] = useState(0);

  const itrs = detected ?? annotated;

  // Default pair: the first and the last ITR, which is the right answer for a plasmid with two.
  const chosenPair: [Itr, Itr] | null = useMemo(() => {
    if (itrs.length < 2) return null;
    const [i, j] = pair ?? [0, itrs.length - 1];
    const first = itrs[Math.min(i, itrs.length - 1)];
    const second = itrs[Math.min(j, itrs.length - 1)];
    return first === second ? null : [first, second];
  }, [itrs, pair]);

  const options: Cargo[] = useMemo(
    () => (chosenPair && data ? cargoOptions(chosenPair, data, sequenceLength) : []),
    [chosenPair, data, sequenceLength],
  );
  const suggestion = options.length ? chooseCargo(options) : null;
  const suggestedIndex = suggestion ? options.indexOf(suggestion.cargo) : 0;
  const cargo = options[sideIndex] ?? options[suggestedIndex] ?? null;
  const verdict = cargo ? verdictFor(cargo.totalBp, kind) : null;
  const { limit, hardLimit } = limitFor(kind);

  const searchCatalogue = useCallback(() => {
    if (!data) return;
    setSearching(true);
    setError(null);
    api
      .detectFeatures<{ features: DetectedFeature[] }>({ sequenceData: data })
      .then((payload) => {
        const found = findItrs(
          payload.features
            .filter((f) => f.type === "repeat_region")
            .map((f) => ({ name: f.name, type: f.type, start: f.start, end: f.end,
                           notes: { note: [f.note || f.name], detected_from: [f.source] } })),
          sequenceLength,
          "catalogue",
        );
        setDetected(found);
        setPair(null);
        setSideIndex(0);
        if (!found.length) {
          setError("No ITR was found in this sequence, annotated or not. If the construct does "
            + "have ITRs, they may be too divergent to recognise: annotate them by hand "
            + "(select the region, then Features → Add Feature…, with \"ITR\" in the name).");
        }
      })
      .catch((e) => setError(String((e as Error).message)))
      .finally(() => setSearching(false));
  }, [data, sequenceLength]);

  const itrLabel = (candidate: Itr) =>
    `${candidate.name || "ITR"} · ${at(candidate.start)}–${at(candidate.end)} `
    + `(${bp(candidate.length)})`;

  const sideLabel = (option: Cargo) =>
    `${at(option.start)} → ${at(option.end)}`
    + `${option.crossesOrigin ? " across the origin" : ""} · ${bp(option.totalBp)}`
    + `${option.backbone.length ? ` · contains ${option.backbone.slice(0, 2).join(", ")}` : ""}`;

  return (
    <Dialog isOpen onClose={props.onClose} title="AAV Packaging Check"
            style={{ width: 920, maxWidth: "95vw" }}>
      <div className={`${Classes.DIALOG_BODY} vt-cloning`} data-testid="aav-dialog">
        <div className="vt-cloning-form">
          <p className="vt-muted">
            Measures what an AAV capsid would have to package: everything from the start of one
            ITR to the end of the other, ITRs included, and compares it with the room a
            single-stranded or a self-complementary vector has.
          </p>

          {itrs.length >= 2 && (
            <>
              <Field
                label="ITRs"
                hint={detected
                  ? "Found with the feature catalogue; they are not annotated in the document yet."
                  : `Found in ${doc?.name ?? "the document"}'s own annotation.`}
              >
                <div className="vt-aav-itrs">
                  <div className="bp3-html-select bp3-fill">
                    <select
                      data-testid="aav-itr-from"
                      value={itrs.indexOf(chosenPair?.[0] ?? itrs[0])}
                      onChange={(e) => {
                        setPair([Number(e.target.value), pair?.[1] ?? itrs.length - 1]);
                        setSideIndex(0);
                      }}
                    >
                      {itrs.map((candidate, i) => (
                        <option key={i} value={i}>{itrLabel(candidate)}</option>
                      ))}
                    </select>
                  </div>
                  <div className="bp3-html-select bp3-fill">
                    <select
                      data-testid="aav-itr-to"
                      value={itrs.indexOf(chosenPair?.[1] ?? itrs[itrs.length - 1])}
                      onChange={(e) => {
                        setPair([pair?.[0] ?? 0, Number(e.target.value)]);
                        setSideIndex(0);
                      }}
                    >
                      {itrs.map((candidate, i) => (
                        <option key={i} value={i}>{itrLabel(candidate)}</option>
                      ))}
                    </select>
                  </div>
                </div>
              </Field>
              {itrs.length > 2 && (
                <Callout intent="warning" data-testid="aav-many-itrs">
                  {itrs.length} ITRs are annotated here. Pick the pair the genome lies between —
                  the two boxes above are the 5′ and the 3′ ITR.
                </Callout>
              )}

              {options.length > 1 && (
                <Field
                  label="Packaged side"
                  hint={suggestion?.reason || "Which way round the ITRs the genome runs."}
                >
                  <div className="bp3-html-select bp3-fill">
                    <select
                      data-testid="aav-side"
                      value={sideIndex}
                      onChange={(e) => setSideIndex(Number(e.target.value))}
                    >
                      {options.map((option, i) => (
                        <option key={i} value={i}>{sideLabel(option)}</option>
                      ))}
                    </select>
                  </div>
                </Field>
              )}

              <Field label="Vector" hint="Self-complementary vectors carry the cassette twice.">
                <div className="bp3-html-select bp3-fill">
                  <select
                    data-testid="aav-kind"
                    value={kind}
                    onChange={(e) => setKind(e.target.value as VectorKind)}
                  >
                    <option value="ss">Single-stranded AAV (ssAAV)</option>
                    <option value="sc">Self-complementary AAV (scAAV)</option>
                  </select>
                </div>
              </Field>
            </>
          )}

          {itrs.length < 2 && (
            <Callout intent="primary" data-testid="aav-no-itrs">
              <p>
                {itrs.length === 0
                  ? "No ITR is annotated in this document."
                  : "Only one ITR is annotated, and a genome needs two."}
              </p>
              <p>
                Vettore recognises a feature as an ITR when its name or its notes say “ITR” or
                “inverted terminal repeat”, or when it is a repeat region marked as inverted.
              </p>
              <Button
                intent="primary"
                data-testid="aav-detect"
                loading={searching}
                onClick={searchCatalogue}
              >
                Search for ITRs in the sequence
              </Button>
              <p className="vt-muted vt-aav-hint">
                This searches the common-feature catalogue, without changing the document. To
                keep the annotations, use Features → Detect Common Features…
              </p>
            </Callout>
          )}

          {error && (
            <Callout intent="warning" data-testid="aav-error">{error}</Callout>
          )}
        </div>

        <div className="vt-cloning-result">
          {!cargo && !searching && (
            <p className="vt-muted">The cargo size appears here once two ITRs are known.</p>
          )}
          {searching && <p className="vt-muted">Searching the sequence for ITRs…</p>}
          {cargo && verdict && (
            <div data-testid="aav-result">
              <Callout intent={INTENT[verdict.level]} title={verdict.headline}
                       data-testid="aav-verdict">
                {verdict.detail}
              </Callout>

              <div className="vt-aav-bar" data-testid="aav-bar"
                   title={`${bp(cargo.totalBp)} of the ${bp(limit)} a ${
                     kind === "sc" ? "self-complementary" : "single-stranded"} vector holds`}>
                <div className={`vt-aav-fill vt-aav-${verdict.level}`}
                     style={{ width: `${Math.min(100, 100 * cargo.totalBp / hardLimit)}%` }} />
                <span className="vt-aav-limit" style={{ left: `${100 * limit / hardLimit}%` }}
                      title={`Packaging limit: ${bp(limit)}`} />
                <span className="vt-aav-bar-text">
                  {bp(cargo.totalBp)} of {bp(limit)}
                </span>
              </div>
              <p className="vt-muted vt-aav-scale" data-testid="aav-scale">
                <span>0</span>
                <span>{bp(hardLimit)}</span>
                {/* Under the marker, on a line of its own so it never runs into the scale's end. */}
                <span className="vt-aav-limit-caption" data-testid="aav-limit-caption"
                      style={{ "--at": `${100 * limit / hardLimit}%` } as React.CSSProperties}>
                  limit {bp(limit)}
                </span>
              </p>

              <table className="bp3-html-table bp3-html-table-condensed vt-table"
                     data-testid="aav-sizes">
                <tbody>
                  <tr>
                    <td>ITR to ITR, ITRs included</td>
                    <td className="vt-bold" data-testid="aav-total">{bp(cargo.totalBp)}</td>
                    <td className="vt-muted">
                      {at(cargo.start)}–{at(cargo.end)}
                      {cargo.crossesOrigin ? ", across the origin" : ""}
                    </td>
                  </tr>
                  <tr>
                    <td>Between the ITRs (promoter, transgene, polyA…)</td>
                    <td data-testid="aav-insert">{bp(cargo.insertBp)}</td>
                    <td className="vt-muted">what you can still change</td>
                  </tr>
                  <tr>
                    <td>5′ ITR / 3′ ITR</td>
                    <td>{bp(cargo.from.length)} / {bp(cargo.to.length)}</td>
                    <td className="vt-muted">
                      {cargo.from.source === "catalogue" ? "found by the catalogue" : "annotated"}
                    </td>
                  </tr>
                  <tr>
                    <td>Wild-type AAV2, for comparison</td>
                    <td className="vt-muted">{bp(WILD_TYPE_BP)}</td>
                    <td className="vt-muted">ITRs included</td>
                  </tr>
                </tbody>
              </table>

              {cargo.backbone.length > 0 && (
                <Callout intent="warning" data-testid="aav-backbone-warning">
                  This span contains {cargo.backbone.slice(0, 3).join(", ")}, which belong to the
                  plasmid backbone rather than to the genome. If the cassette is on the other
                  side, change “Packaged side” above.
                </Callout>
              )}

              <p className="vt-muted vt-aav-sources">
                Room for a single-stranded genome is taken as {bp(SS_LIMIT_BP)}, about the size
                of wild-type AAV2; beyond {bp(SS_HARD_LIMIT_BP)} genomes are packaged truncated
                if at all (Dong 1996, Hum Gene Ther 7:2101; Wu 2010, Mol Ther 18:80). A
                self-complementary vector carries the cassette twice, so it holds about half,
                {" "}{bp(SC_LIMIT_BP)} (McCarty 2001, Gene Ther 8:1248).
              </p>
            </div>
          )}
        </div>
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button intent="primary" onClick={props.onClose}>Close</Button>
        </div>
      </div>
    </Dialog>
  );
}

export const tools: Tool[] = [
  {
    id: "aav-packaging",
    menu: "Tools",
    label: "AAV Packaging Check",
    requiresDocument: true,
    Component: AavPackagingDialog,
  },
];
