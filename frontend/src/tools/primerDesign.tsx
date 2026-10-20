// Primers → Design PCR Primers: primer3 pairs around a region of the active document
// (engine: src/vettore/primer_design.py). The chosen pair can be added to the document as
// primers, or copied in the name ⇥ sequence format order forms accept.
import { useEffect, useMemo, useRef, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import { actions } from "@teselagen/ove";
import { call, type OveAnnotation, type SequenceData } from "../api";
import { store } from "../store";
import type { OpenDoc } from "../cloning/types";
import { Field } from "../cloning/inputs";
import { useCopyText } from "../cloning/useCopyText";
import type { Tool, ToolProps } from "./types";
import {
  buildRequest,
  DEFAULT_FORM,
  defaultProductRange,
  orderingText,
  primerNames,
  range,
  regionLength,
  validate,
  type PrimerForm,
} from "./primerText";
import "./primerDesign.css";

// Blueprint 3 typings predate React 18 (no implicit `children`).
const { Dialog, Callout } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

interface PrimerInfo {
  sequence: string;
  start: number;
  end: number;
  length: number;
  tm: number;
  tm_primer3: number;
  tm_presets: { key: string; label: string; description: string; tm: number }[];
  gc: number;
}

interface Pair {
  rank: number;
  penalty: number;
  product_size: number;
  product_start: number;
  product_end: number;
  forward: PrimerInfo;
  reverse: PrimerInfo;
  structures: Record<string, number>;
  warnings: string[];
}

interface DesignResult {
  pairs: Pair[];
  candidates: number;
  product_range: [number, number];
  tm_method: string;
  dg_conditions: string;
}

const asList = <T,>(value: T[] | Record<string, T> | undefined): T[] =>
  Array.isArray(value) ? value : Object.values(value ?? {});

const bp = (value: number) => value.toLocaleString("en-US");

function tmTitle(p: PrimerInfo): string {
  return p.tm_presets.map((t) => `${t.label}: ${t.tm.toFixed(1)} °C — ${t.description}`).join("\n");
}

/** ΔG in kcal/mol; a structure with ΔG ≥ 0 does not form at 37 °C. */
function dg(value: number): string {
  return value >= 0 ? "none" : value.toFixed(1).replace("-", "−");
}

/** Start/end and name prefix for the dialog's first render, from the editor selection. */
function initialRegion(doc: OpenDoc): { start: string; end: string; prefix: string } {
  const sel = doc.selection();
  if (!sel) return { start: "", end: "", prefix: doc.name };
  const feature = asList<OveAnnotation>(doc.data()?.features)
    .find((f) => f.start === sel.start && f.end === sel.end);
  return { start: String(sel.start + 1), end: String(sel.end + 1), prefix: feature?.name ?? doc.name };
}

function DesignPcrPrimers(props: ToolProps) {
  const doc = props.active;
  if (!doc) {
    return (
      <Dialog isOpen onClose={props.onClose} title="Design PCR Primers">
        <div className={Classes.DIALOG_BODY}>
          <Callout intent="primary">Open a sequence first: primers are designed on the active tab.</Callout>
        </div>
      </Dialog>
    );
  }
  return <Designer {...props} doc={doc} />;
}

function Designer(props: ToolProps & { doc: OpenDoc }) {
  const { doc } = props;
  const data = doc.data() as SequenceData;
  const length = data.sequence.length;
  const circular = !!data.circular;
  const features = useMemo(
    () => asList<OveAnnotation>(data.features)
      .filter((f) => f.start >= 0 && f.end >= 0)
      .sort((a, b) => a.start - b.start),
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [],
  );
  const initial = useMemo(() => initialRegion(doc), []); // eslint-disable-line react-hooks/exhaustive-deps
  const [form, setForm] = useState<PrimerForm>({ ...DEFAULT_FORM, start: initial.start, end: initial.end });
  const [prefix, setPrefix] = useState(initial.prefix);
  const [result, setResult] = useState<DesignResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [chosen, setChosen] = useState(0);
  const [notice, setNotice] = useState<string | null>(null);
  const { copy: copyText, fallbackDialog } = useCopyText();
  const request = useRef(0);

  const errors = validate(form, length, circular);
  const region = regionLength(form, length, circular);
  const body = useMemo(
    () => buildRequest(form, { ...data, name: doc.name, circular } as any),
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [form],
  );
  const bodyKey = body ? JSON.stringify(body) : "";

  // Live design, as in the cloning dialogs: debounced, the latest request wins.
  useEffect(() => {
    setError(null);
    setNotice(null);
    if (!bodyKey) {
      setResult(null);
      return;
    }
    const id = ++request.current;
    const timer = setTimeout(() => {
      setBusy(true);
      call<DesignResult>("primer-design/design", JSON.parse(bodyKey))
        .then((r) => {
          if (id !== request.current) return;
          setResult(r);
          setChosen(0);
        })
        .catch((e) => {
          if (id !== request.current) return;
          setResult(null);
          setError(String((e as Error).message));
        })
        .finally(() => id === request.current && setBusy(false));
    }, 400);
    return () => clearTimeout(timer);
  }, [bodyKey]);

  const set = (patch: Partial<PrimerForm>) => setForm((f) => ({ ...f, ...patch }));

  const useSelection = () => {
    const sel = doc.selection();
    if (sel) set({ start: String(sel.start + 1), end: String(sel.end + 1) });
  };

  const useFeature = (index: string) => {
    const f = features[Number(index)];
    if (!f) return;
    set({ start: String(f.start + 1), end: String(f.end + 1) });
    setPrefix(f.name);
  };

  const pair = result?.pairs[chosen] ?? null;
  const existingNames = asList<OveAnnotation>((doc.data() as SequenceData).primers).map((p) => p.name);
  const [forwardName, reverseName] = primerNames(prefix, existingNames);

  const copy = async () => {
    if (!pair) return;
    const text = orderingText([
      { name: forwardName, sequence: pair.forward.sequence },
      { name: reverseName, sequence: pair.reverse.sequence },
    ]);
    setNotice(null);
    if (await copyText(text, "Primers to order")) {
      setNotice(`Copied ${forwardName} and ${reverseName}: paste them into the order form.`);
    }
  };

  const add = () => {
    if (!pair) return;
    const current = doc.data() as SequenceData;
    const note = (p: PrimerInfo) =>
      [`primer3 pair ${pair.rank}; product ${pair.product_size} bp; Tm ${p.tm.toFixed(1)} °C, `
        + `${result?.tm_method}`];
    const added: OveAnnotation[] = [
      { name: forwardName, type: "primer_bind", start: pair.forward.start, end: pair.forward.end,
        strand: 1, forward: true, bases: pair.forward.sequence, notes: { note: note(pair.forward) } },
      { name: reverseName, type: "primer_bind", start: pair.reverse.start, end: pair.reverse.end,
        strand: -1, forward: false, bases: pair.reverse.sequence, notes: { note: note(pair.reverse) } },
    ];
    store.dispatch(actions.updateSequenceData(
      { ...current, primers: [...asList<OveAnnotation>(current.primers), ...added] },
      { editorName: doc.editorName },
    ));
    props.logHistory(doc.editorName, "Added primers (primer3)",
      `${forwardName} ${pair.forward.sequence}, ${reverseName} ${pair.reverse.sequence}; `
      + `product ${pair.product_size} bp`);
    setNotice(`Added ${forwardName} and ${reverseName} to ${doc.name}; they are listed in the Primers view.`);
  };

  const defaults = region !== null && !errors.size
    ? defaultProductRange(region, Number(form.sizeMin)) : null;

  return (
    <Dialog isOpen onClose={props.onClose} title="Design PCR Primers"
            style={{ width: 1260, maxWidth: "95vw" }} canOutsideClickClose={false}>
      <div className={`${Classes.DIALOG_BODY} vt-pd`} data-testid="primer-design-dialog">
        <div className="vt-pd-form">
          <p className="vt-pd-intro">
            Finds the best primer pairs to amplify a region of <b>{doc.name}</b> (primer3). The
            product contains the whole region; the primers sit just outside it.
          </p>

          <Field label="Region to amplify">
            {features.length > 0 && (
              <div className="bp3-html-select bp3-fill">
                <select data-testid="pd-feature" value="" onChange={(e) => useFeature(e.target.value)}>
                  <option value="">Pick a feature…</option>
                  {features.map((f, i) => (
                    <option key={i} value={i}>{f.name} ({range(f.start, f.end)})</option>
                  ))}
                </select>
              </div>
            )}
            <div className="vt-pd-inline">
              <input className="bp3-input bp3-small" data-testid="pd-start" inputMode="numeric"
                     aria-label="Region start" placeholder="start" value={form.start}
                     onChange={(e) => set({ start: e.target.value })} />
              <span>…</span>
              <input className="bp3-input bp3-small" data-testid="pd-end" inputMode="numeric"
                     aria-label="Region end" placeholder="end" value={form.end}
                     onChange={(e) => set({ end: e.target.value })} />
              <Button small data-testid="pd-use-selection" disabled={!doc.selection()} onClick={useSelection}>
                Use selection
              </Button>
              {region !== null && <span className="vt-pd-muted">{bp(region)} bp</span>}
            </div>
            {errors.region && <span className="vt-pd-invalid" data-testid="pd-region-error">{errors.region}</span>}
          </Field>

          <Field label="Product size (bp)">
            <div className="vt-pd-inline">
              <input className="bp3-input bp3-small" data-testid="pd-product-min" inputMode="numeric"
                     aria-label="Smallest product" placeholder={defaults ? String(defaults[0]) : "min"}
                     value={form.productMin} onChange={(e) => set({ productMin: e.target.value })} />
              <span>to</span>
              <input className="bp3-input bp3-small" data-testid="pd-product-max" inputMode="numeric"
                     aria-label="Largest product" placeholder={defaults ? String(defaults[1]) : "max"}
                     value={form.productMax} onChange={(e) => set({ productMax: e.target.value })} />
            </div>
            {errors.product
              ? <span className="vt-pd-invalid" data-testid="pd-product-error">{errors.product}</span>
              : <span className="vt-field-hint">Leave empty for primers within 150 bp of the region on each side.</span>}
          </Field>

          <Field label="Primer length (nt): min, optimum, max">
            <div className="vt-pd-inline">
              {(["sizeMin", "sizeOpt", "sizeMax"] as const).map((k) => (
                <input key={k} className="bp3-input bp3-small" inputMode="numeric" aria-label={k}
                       value={form[k]} onChange={(e) => set({ [k]: e.target.value })} />
              ))}
            </div>
            {errors.size && <span className="vt-pd-invalid">{errors.size}</span>}
          </Field>

          <Field label="Tm (°C): min, optimum, max">
            <div className="vt-pd-inline">
              {(["tmMin", "tmOpt", "tmMax"] as const).map((k) => (
                <input key={k} className="bp3-input bp3-small" inputMode="decimal" aria-label={k}
                       data-testid={`pd-${k}`} value={form[k]} onChange={(e) => set({ [k]: e.target.value })} />
              ))}
            </div>
            {errors.tm
              ? <span className="vt-pd-invalid">{errors.tm}</span>
              : <span className="vt-field-hint">Standard preset: SantaLucia &amp; Hicks 2004, Na⁺ 50 mM, primer 250 nM.</span>}
          </Field>

          <div className="vt-pd-pair">
            <Field label="GC (%): min, max">
              <div className="vt-pd-inline">
                {(["gcMin", "gcMax"] as const).map((k) => (
                  <input key={k} className="bp3-input bp3-small" inputMode="decimal" aria-label={k}
                         value={form[k]} onChange={(e) => set({ [k]: e.target.value })} />
                ))}
              </div>
              {errors.gc && <span className="vt-pd-invalid">{errors.gc}</span>}
            </Field>
            <Field label="Max Tm difference (°C)">
              <input className="bp3-input bp3-small" inputMode="decimal" aria-label="Max Tm difference"
                     value={form.maxTmDifference} onChange={(e) => set({ maxTmDifference: e.target.value })} />
              {errors.tmDiff && <span className="vt-pd-invalid">{errors.tmDiff}</span>}
            </Field>
          </div>

          <Field label="Primer names" hint={`Added as ${forwardName} and ${reverseName}.`}>
            <input className="bp3-input bp3-small" data-testid="pd-prefix" value={prefix}
                   onChange={(e) => setPrefix(e.target.value)} />
          </Field>
          <p className="vt-pd-muted">
            To put the primers exactly on the region's ends (for cloning), use Actions → Design
            Primers to Amplify a Region.
          </p>
        </div>

        <div className="vt-pd-result">
          {busy && <p className="vt-pd-muted" data-testid="pd-busy">Designing…</p>}
          {!bodyKey && !busy && (
            <p className="vt-pd-muted">Choose the region on the left: the primer pairs appear here.</p>
          )}
          {error && (
            <Callout intent="danger" title="No primer pair" data-testid="pd-error">{error}</Callout>
          )}
          {result && (
            <>
              <div className="vt-pd-table-wrap">
              <table className="bp3-html-table bp3-html-table-condensed vt-table vt-pd-table" data-testid="pd-pairs">
                <thead>
                  <tr>
                    <th>Pair</th><th></th><th>Sequence 5′→3′</th><th>Position</th><th>nt</th>
                    <th>Tm</th><th>GC</th>
                    <th title="ΔG in kcal/mol at 37 °C; “none” means the structure does not form">
                      ΔG hairpin · self-dimer
                    </th>
                  </tr>
                </thead>
                {result.pairs.map((p, i) => {
                  const selected = i === chosen;
                  const cls = selected ? "vt-row-selected" : "";
                  return (
                    <tbody key={p.rank} className={`vt-pd-pair-rows${selected ? " is-selected" : ""}`}
                           onClick={() => setChosen(i)} data-testid="pd-pair">
                      <tr className={cls}>
                        <td rowSpan={p.warnings.length ? 3 : 2} className="vt-pd-rank">
                          <label>
                            <input type="radio" name="pd-pair" checked={selected} onChange={() => setChosen(i)} />{" "}
                            <b>{p.rank}</b> · {bp(p.product_size)} bp
                          </label>
                          <div className="vt-pd-muted" title="ΔG of the two primers pairing with each other, kcal/mol at 37 °C">
                            primer-dimer {dg(p.structures.hetero_dimer)}
                          </div>
                          <div className="vt-pd-muted" title="primer3 penalty: lower is closer to the optimum">
                            penalty {p.penalty.toFixed(2)}
                          </div>
                        </td>
                        <td>F</td>
                        <td className="vt-pd-seq">{p.forward.sequence}</td>
                        <td>{range(p.forward.start, p.forward.end)}</td>
                        <td>{p.forward.length}</td>
                        <td title={tmTitle(p.forward)}>{p.forward.tm.toFixed(1)} °C</td>
                        <td>{p.forward.gc.toFixed(0)} %</td>
                        <td>{dg(p.structures.forward_hairpin)} · {dg(p.structures.forward_self_dimer)}</td>
                      </tr>
                      <tr className={cls}>
                        <td>R</td>
                        <td className="vt-pd-seq">{p.reverse.sequence}</td>
                        <td>{range(p.reverse.start, p.reverse.end)}</td>
                        <td>{p.reverse.length}</td>
                        <td title={tmTitle(p.reverse)}>{p.reverse.tm.toFixed(1)} °C</td>
                        <td>{p.reverse.gc.toFixed(0)} %</td>
                        <td>{dg(p.structures.reverse_hairpin)} · {dg(p.structures.reverse_self_dimer)}</td>
                      </tr>
                      {p.warnings.length > 0 && (
                        <tr className={cls}>
                          <td colSpan={7} className="vt-pd-warnings" data-testid="pd-warning">
                            {p.warnings.map((w, j) => <div key={j}>⚠ {w}</div>)}
                          </td>
                        </tr>
                      )}
                    </tbody>
                  );
                })}
              </table>
              </div>
              <p className="vt-pd-muted">
                Best {result.pairs.length} of {result.candidates} pairs from primer3, ranked by its
                penalty, pairs without warnings first.
              </p>
              <p className="vt-pd-muted">
                Tm: {result.tm_method}; hover a Tm for the NEB Q5 and Taq values. ΔG in kcal/mol at
                37 °C, Na⁺ 50 mM; IDT advises keeping hairpins and dimers above −9.
              </p>
            </>
          )}
          {notice && <Callout intent="success" data-testid="pd-notice">{notice}</Callout>}
        </div>
      </div>
      {fallbackDialog}
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button data-testid="pd-close" onClick={props.onClose}>Close</Button>
          <Button icon="duplicate" data-testid="pd-copy" disabled={!pair || busy} onClick={copy}>
            Copy for Ordering
          </Button>
          <Button intent="primary" data-testid="pd-add" disabled={!pair || busy} onClick={add}>
            Add to Sequence
          </Button>
        </div>
      </div>
    </Dialog>
  );
}

export const tools: Tool[] = [
  {
    id: "design-pcr-primers",
    menu: "Primers",
    label: "Design PCR Primers",
    requiresDocument: true,
    Component: DesignPcrPrimers,
  },
];
