// One form per cloning action. Each form reports the JSON body of its /api/cloning endpoint.
import { useEffect, useRef, useState } from "react";
import { api } from "../api";
import { ENZYMES } from "../enzymeLogic";
import {
  GATEWAY_SITE_KINDS, defaultInserts, defaultLinear, defaultVector, docData, gatewayDefaults, openingFromSelection,
} from "./defaults";
import { FragmentEditor } from "./fragments";
import { newRow, rowBody, type FragmentRow } from "./fragmentRows";
import { DocSelect, EnzymePicker, Field, RegionInput, findDoc } from "./inputs";
import type { FormProps, OpenDoc } from "./types";

const dataOf = (d: OpenDoc | null) => (d ? d.data() : null);

// ---- Restriction and insertion cloning ---------------------------------------------------
export function RestrictionInsertForm({ docs, active, onChange }: FormProps) {
  const firstVector = defaultVector(docs, active);
  const [vector, setVector] = useState<string | null>(firstVector);
  const [insert, setInsert] = useState<string | null>(defaultInserts(docs, firstVector, 1)[0]);
  const [vEnz, setVEnz] = useState<string[]>([]);
  const [iEnz, setIEnz] = useState<string[]>([]);
  const [orientation, setOrientation] = useState<"forward" | "reverse">("forward");
  const v = findDoc(docs, vector);
  const i = findDoc(docs, insert);

  useEffect(() => {
    if (!v) return onChange(null, "Choose the vector");
    if (!vEnz.length) return onChange(null, `Choose the enzyme(s) that cut the vector “${v.name}”`);
    if (!i) return onChange(null, "Choose the file the insert is cut from");
    if (!iEnz.length) return onChange(null, `Choose the enzyme(s) that cut “${i.name}”`);
    onChange({
      vector: docData(v),
      vector_enzymes: vEnz,
      insert: docData(i),
      insert_enzymes: iEnz,
      insert_fragment_index: null,
      orientation,
    });
  }, [v, i, vEnz, iEnz, orientation]); // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <Field label="Vector">
        <DocSelect docs={docs} value={vector} onChange={setVector} testId="vector-doc" />
      </Field>
      <Field label="Cut vector with" hint="One enzyme, or two for directional cloning">
        <EnzymePicker value={vEnz} onChange={setVEnz} max={2} sequence={dataOf(v)?.sequence} circular={dataOf(v)?.circular} testId="vector-enzymes" />
      </Field>
      <Field label="Insert from">
        <DocSelect docs={docs} value={insert} onChange={setInsert} testId="insert-doc" />
      </Field>
      <Field label="Cut insert with">
        <EnzymePicker value={iEnz} onChange={setIEnz} max={2} sequence={dataOf(i)?.sequence} circular={dataOf(i)?.circular} testId="insert-enzymes" />
      </Field>
      <Field label="Orientation">
        <div className="bp3-html-select">
          <select value={orientation} onChange={(e) => setOrientation(e.target.value as "forward" | "reverse")}>
            <option value="forward">Forward</option>
            <option value="reverse">Reverse (flipped)</option>
          </select>
        </div>
      </Field>
    </>
  );
}

// ---- PCR ---------------------------------------------------------------------------------
function PrimerInput(props: {
  label: string;
  template: OpenDoc | null;
  value: { name: string; sequence: string };
  onChange: (p: { name: string; sequence: string }) => void;
  testId?: string;
  /** Inline validation shown under the primer. */
  children?: React.ReactNode;
}) {
  const primers = (() => {
    const p = props.template?.data().primers ?? [];
    return (Array.isArray(p) ? p : Object.values(p)).filter((x): x is typeof x & { bases: string } => !!x.bases);
  })();
  return (
    <Field label={props.label}>
      <div className="vt-primer-input">
        <input
          className="bp3-input bp3-small"
          placeholder="Name"
          value={props.value.name}
          onChange={(e) => props.onChange({ ...props.value, name: e.target.value })}
        />
        <input
          className="bp3-input bp3-small vt-mono"
          placeholder="5′ → 3′ sequence"
          data-testid={props.testId}
          value={props.value.sequence}
          onChange={(e) => props.onChange({ ...props.value, sequence: e.target.value.replace(/[^ACGTRYSWKMBDHVNacgtryswkmbdhvn]/g, "") })}
        />
        {primers.length > 0 && (
          <div className="bp3-html-select bp3-minimal">
            <select value="" onChange={(e) => {
              const p = primers.find((x) => x.name === e.target.value);
              if (p) props.onChange({ name: p.name, sequence: p.bases });
            }}>
              <option value="">From file…</option>
              {primers.map((p) => <option key={p.name} value={p.name}>{p.name}</option>)}
            </select>
          </div>
        )}
      </div>
      {props.children}
    </Field>
  );
}

export function PcrForm({ docs, active, onChange }: FormProps) {
  const [template, setTemplate] = useState<string | null>(active?.editorName ?? null);
  const [fwd, setFwd] = useState({ name: "Forward", sequence: "" });
  const [rev, setRev] = useState({ name: "Reverse", sequence: "" });
  const t = findDoc(docs, template);
  useEffect(() => {
    if (!t) return onChange(null, "Choose the template to amplify");
    if (fwd.sequence.length < 10 || rev.sequence.length < 10) {
      return onChange(null, "Enter both primers (at least 10 nt each), or pick them from the file");
    }
    onChange({ template: docData(t), forward: fwd, reverse: rev, mismatches: 0 });
  }, [t, fwd, rev]); // eslint-disable-line react-hooks/exhaustive-deps
  return (
    <>
      <Field label="Template">
        <DocSelect docs={docs} value={template} onChange={setTemplate} testId="template-doc" />
      </Field>
      <PrimerInput label="Forward primer" template={t} value={fwd} onChange={setFwd} />
      <PrimerInput label="Reverse primer" template={t} value={rev} onChange={setRev} />
    </>
  );
}

export function DesignPrimersForm({ docs, active, onChange }: FormProps) {
  const [template, setTemplate] = useState<string | null>(active?.editorName ?? null);
  const t = findDoc(docs, template);
  const sel = t?.selection() ?? null;
  const [start, setStart] = useState<number | null>(sel?.start ?? null);
  const [end, setEnd] = useState<number | null>(sel?.end ?? null);
  const [tm, setTm] = useState(60);
  const [fTail, setFTail] = useState("");
  const [rTail, setRTail] = useState("");
  useEffect(() => {
    if (!t) return onChange(null, "Choose the template");
    if (start === null || end === null) {
      return onChange(null, "Choose the region to amplify, or select it in the editor and press "
        + "“Use selection”");
    }
    onChange({ template: docData(t), start, end, target_tm: tm, forward_tail: fTail, reverse_tail: rTail });
  }, [t, start, end, tm, fTail, rTail]); // eslint-disable-line react-hooks/exhaustive-deps
  return (
    <>
      <Field label="Template">
        <DocSelect docs={docs} value={template} onChange={setTemplate} />
      </Field>
      <Field label="Region to amplify">
        <RegionInput doc={t} start={start} end={end} onChange={(s, e) => { setStart(s); setEnd(e); }} testId="region" features />
      </Field>
      <Field label="Target Tm of annealing parts (°C)">
        <input className="bp3-input bp3-small" type="number" value={tm} onChange={(e) => setTm(Number(e.target.value))} />
      </Field>
      <Field label="Forward primer 5′ tail" hint="e.g. a restriction site with a few extra bases">
        <input className="bp3-input bp3-small vt-mono" value={fTail} onChange={(e) => setFTail(e.target.value.toUpperCase().replace(/[^ACGTN]/g, ""))} />
      </Field>
      <Field label="Reverse primer 5′ tail">
        <input className="bp3-input bp3-small vt-mono" value={rTail} onChange={(e) => setRTail(e.target.value.toUpperCase().replace(/[^ACGTN]/g, ""))} />
      </Field>
    </>
  );
}

// ---- Gibson / NEBuilder HiFi / In-Fusion -------------------------------------------------
interface InsertRow {
  editorName: string | null;
  start: number | null;
  end: number | null;
}

export function overlapForm(method: "gibson" | "hifi" | "infusion") {
  return function OverlapAssemblyForm({ docs, active, onChange }: FormProps) {
    const firstVector = defaultVector(docs, active);
    const [vector, setVector] = useState<string | null>(firstVector);
    const [mode, setMode] = useState<"enzymes" | "pcr">("enzymes");
    const [vEnz, setVEnz] = useState<string[]>([]);
    const [vStart, setVStart] = useState<number | null>(null);
    const [vEnd, setVEnd] = useState<number | null>(null);
    const [inserts, setInserts] = useState<InsertRow[]>(() =>
      defaultInserts(docs, firstVector, 1).map((editorName) => ({ editorName, start: null, end: null })));
    const [overlap, setOverlap] = useState(method === "infusion" ? 15 : 25);
    const v = findDoc(docs, vector);

    useEffect(() => {
      if (!v) return onChange(null, "Choose the vector");
      const linearize = mode === "enzymes"
        ? vEnz.length ? { enzymes: vEnz } : null
        : vStart !== null && vEnd !== null ? { pcr: { start: vStart, end: vEnd } } : null;
      const rows = inserts.map((r) => ({ row: r, doc: findDoc(docs, r.editorName) }));
      if (!linearize) {
        return onChange(null, mode === "enzymes"
          ? `Choose the enzyme(s) that open the vector “${v.name}”`
          : "Choose the vector region to keep");
      }
      const missing = rows.findIndex((r) => !r.doc);
      if (missing >= 0) return onChange(null, `Choose the file for insert ${missing + 1}`);
      onChange({
        method,
        vector: { doc: docData(v), linearize },
        inserts: rows.map(({ row, doc }) => ({ doc: docData(doc!), start: row.start, end: row.end })),
        overlap,
        target_tm: 60,
      });
    }, [v, mode, vEnz, vStart, vEnd, inserts, overlap]); // eslint-disable-line react-hooks/exhaustive-deps

    const update = (i: number, patch: Partial<InsertRow>) =>
      setInserts((rows) => rows.map((r, j) => (j === i ? { ...r, ...patch } : r)));

    return (
      <>
        <Field label="Vector">
          <DocSelect docs={docs} value={vector} onChange={setVector} testId="vector-doc" />
        </Field>
        <Field label="Linearize vector by">
          <div className="bp3-html-select">
            <select value={mode} onChange={(e) => setMode(e.target.value as "enzymes" | "pcr")}>
              <option value="enzymes">Restriction digest</option>
              <option value="pcr">PCR (inverse PCR of the region to keep)</option>
            </select>
          </div>
        </Field>
        {mode === "enzymes" ? (
          <Field label="Enzymes">
            <EnzymePicker value={vEnz} onChange={setVEnz} max={2} sequence={dataOf(v)?.sequence} circular={dataOf(v)?.circular} testId="vector-enzymes" />
          </Field>
        ) : (
          <Field label="Vector region to keep">
            <RegionInput doc={v} start={vStart} end={vEnd} onChange={(s, e) => { setVStart(s); setVEnd(e); }} />
          </Field>
        )}
        {inserts.map((row, i) => {
          const d = findDoc(docs, row.editorName);
          return (
            <div key={i} className="vt-insert-row">
              <Field label={`Insert ${i + 1}`}>
                <DocSelect docs={docs} value={row.editorName} onChange={(e) => update(i, { editorName: e })} testId={`insert-doc-${i}`} />
              </Field>
              <Field label="Region to amplify" hint={d?.data().circular ? "Required for circular files" : "Empty = whole sequence"}>
                <RegionInput doc={d} start={row.start} end={row.end} onChange={(s, e) => update(i, { start: s, end: e })}
                             testId={`insert-region-${i}`} features />
              </Field>
              {inserts.length > 1 && (
                <button className="bp3-button bp3-small bp3-minimal" onClick={() => setInserts((r) => r.filter((_, j) => j !== i))}>
                  Remove insert
                </button>
              )}
            </div>
          );
        })}
        <button className="bp3-button bp3-small" onClick={() => setInserts((r) => [...r, { editorName: null, start: null, end: null }])}>
          + Add insert
        </button>
        <Field label="Overlap length (bp)">
          <input className="bp3-input bp3-small" type="number" min={10} max={60} value={overlap} onChange={(e) => setOverlap(Number(e.target.value))} />
        </Field>
      </>
    );
  };
}

// ---- Mutagenesis -------------------------------------------------------------------------
export function MutagenesisForm({ docs, active, onChange }: FormProps) {
  const [docName, setDocName] = useState<string | null>(active?.editorName ?? null);
  const d = findDoc(docs, docName);
  const sel = d?.selection() ?? null;
  const [kind, setKind] = useState<"substitution" | "insertion" | "deletion">("substitution");
  const [start, setStart] = useState<number | null>(sel?.start ?? null);
  const [end, setEnd] = useState<number | null>(sel?.end ?? null);
  const [newSeq, setNewSeq] = useState("");
  useEffect(() => {
    if (!d) return onChange(null, "Choose the plasmid to change");
    if (start === null || end === null) {
      return onChange(null, kind === "insertion"
        ? "Choose the base to insert after"
        : "Choose the region to change, or select it in the editor");
    }
    if (kind !== "deletion" && !newSeq) return onChange(null, "Enter the new bases");
    onChange({ doc: docData(d), kind, start, end: kind === "insertion" ? start : end, new_sequence: kind === "deletion" ? "" : newSeq });
  }, [d, kind, start, end, newSeq]); // eslint-disable-line react-hooks/exhaustive-deps
  return (
    <>
      <Field label="Plasmid">
        <DocSelect docs={docs} value={docName} onChange={setDocName} />
      </Field>
      <Field label="Change">
        <div className="bp3-html-select">
          <select value={kind} onChange={(e) => setKind(e.target.value as typeof kind)}>
            <option value="substitution">Substitution (replace the region)</option>
            <option value="insertion">Insertion (after the first position)</option>
            <option value="deletion">Deletion (remove the region)</option>
          </select>
        </div>
      </Field>
      <Field label={kind === "insertion" ? "Insert after base" : "Region"}>
        <RegionInput doc={d} start={start} end={end} onChange={(s, e) => { setStart(s); setEnd(e); }}
                     testId="region" features={kind !== "insertion"} />
      </Field>
      {kind !== "deletion" && (
        <Field label="New bases (5′ → 3′, top strand)">
          <input className="bp3-input bp3-small vt-mono" value={newSeq} onChange={(e) => setNewSeq(e.target.value.toUpperCase().replace(/[^ACGT]/g, ""))} />
        </Field>
      )}
    </>
  );
}

// ---- Golden Gate -------------------------------------------------------------------------
const GG_ENZYMES = new Set(Object.values(ENZYMES).filter((e) => e.goldenGate).map((e) => e.name));

export function GoldenGateForm({ docs, active, onChange }: FormProps) {
  const firstVector = defaultVector(docs, active);
  const [vector, setVector] = useState<string | null>(firstVector);
  const [parts, setParts] = useState<(string | null)[]>(() => defaultInserts(docs, firstVector, 1));
  const [enzyme, setEnzyme] = useState<string>("BsaI");
  const v = findDoc(docs, vector);
  useEffect(() => {
    const p = parts.map((n) => findDoc(docs, n));
    if (!v) return onChange(null, "Choose the destination vector");
    const missing = p.findIndex((x) => !x);
    if (missing >= 0) return onChange(null, `Choose the file for part ${missing + 1}`);
    onChange({ vector: docData(v), parts: p.map((x) => docData(x!)), enzyme });
  }, [v, parts, enzyme]); // eslint-disable-line react-hooks/exhaustive-deps
  return (
    <>
      <Field label="Enzyme">
        <div className="bp3-html-select">
          <select value={enzyme} onChange={(e) => setEnzyme(e.target.value)}>
            {[...GG_ENZYMES].sort().map((n) => <option key={n} value={n}>{n}</option>)}
          </select>
        </div>
      </Field>
      <Field label="Vector">
        <DocSelect docs={docs} value={vector} onChange={setVector} />
      </Field>
      {parts.map((p, i) => (
        <Field key={i} label={`Part ${i + 1}`}>
          <div className="vt-row">
            <DocSelect docs={docs} value={p} onChange={(e) => setParts((x) => x.map((y, j) => (j === i ? e : y)))} />
            {parts.length > 1 && (
              <button className="bp3-button bp3-small bp3-minimal" onClick={() => setParts((x) => x.filter((_, j) => j !== i))}>×</button>
            )}
          </div>
        </Field>
      ))}
      <button className="bp3-button bp3-small" onClick={() => setParts((x) => [...x, null])}>+ Add part</button>
    </>
  );
}

// ---- Gateway -----------------------------------------------------------------------------
const GATEWAY_LABELS = {
  BP: { first: "attB-flanked insert (PCR product or clone with attB sites)", second: "Donor vector (attP sites)" },
  LR: { first: "Entry clone (attL sites)", second: "Destination vector (attR sites)" },
};

export function gatewayForm(reaction: "BP" | "LR") {
  return function GatewayForm({ docs, active, onChange }: FormProps) {
    const labels = GATEWAY_LABELS[reaction];
    const [insertKind, vectorKind] = GATEWAY_SITE_KINDS[reaction];
    const [a, setA] = useState<string | null>(null);
    const [b, setB] = useState<string | null>(null);
    // The att sites of every open file, by editor name; null while the engine is asked, and
    // left null if it cannot answer (the user then chooses both files).
    const [sites, setSites] = useState<Record<string, string[]> | null>(null);
    const [looking, setLooking] = useState(true);
    const chosenByUser = useRef(false);
    const da = findDoc(docs, a);
    const db = findDoc(docs, b);

    // Pick the inputs by their att sites once, unless the user has already chosen.
    useEffect(() => {
      let live = true;
      api
        .cloning<{ sites: string[][] }>("gateway-sites", { docs: docs.map(docData) })
        .then((r) => {
          if (!live || !Array.isArray(r?.sites)) return;
          setSites(Object.fromEntries(docs.map((d, i) => [d.editorName, r.sites[i] ?? []])));
          if (chosenByUser.current) return;
          const [first, second] = gatewayDefaults(reaction, docs, r.sites, active);
          setA(first);
          setB(second);
        })
        .catch(() => undefined)
        .finally(() => live && setLooking(false));
      return () => {
        live = false;
      };
    }, []); // eslint-disable-line react-hooks/exhaustive-deps

    const hasKind = (kind: string) => !!sites && Object.values(sites).some((s) => s.some((x) => x.startsWith(kind)));
    // "Entry clone (attL sites)" → "entry clone (attL sites)": site names keep their capitals.
    const choose = (label: string, kind: string) => {
      const what = label.charAt(0).toLowerCase() + label.slice(1);
      return sites && !hasKind(kind)
        ? `Choose the ${what}. None of the open files has ${kind} sites: open it first (File → Open…).`
        : `Choose the ${what}`;
    };

    useEffect(() => {
      if (looking && (!da || !db)) return onChange(null, "Looking for att sites in the open files…");
      if (!da) return onChange(null, choose(labels.first, insertKind));
      if (!db) return onChange(null, choose(labels.second, vectorKind));
      if (da.editorName === db.editorName) {
        return onChange(null, "Choose two different files: one carries the insert, the other is the vector");
      }
      onChange({ reaction, inputs: [docData(da), docData(db)] });
    }, [da, db, looking, sites]); // eslint-disable-line react-hooks/exhaustive-deps

    const note = (d: OpenDoc) => (sites?.[d.editorName]?.length ? sites[d.editorName].join(", ") : null);
    return (
      <>
        <Field label={labels.first}>
          <DocSelect docs={docs} value={a} note={note} testId="insert-doc"
                     onChange={(v) => { chosenByUser.current = true; setA(v); }} />
        </Field>
        <Field label={labels.second}>
          <DocSelect docs={docs} value={b} note={note} testId="vector-doc"
                     onChange={(v) => { chosenByUser.current = true; setB(v); }} />
        </Field>
      </>
    );
  };
}

// ---- TA / GC / TOPO ----------------------------------------------------------------------
const DIRECTIONAL_TAG = "CACC";

export function taTopoForm(kinds: { value: string; label: string }[]) {
  return function TaTopoForm({ docs, active, onChange }: FormProps) {
    const firstVector = defaultVector(docs, active);
    const firstInsert = defaultLinear(docs, active, firstVector);
    const [kind, setKind] = useState(kinds[0].value);
    const [vector, setVector] = useState<string | null>(firstVector);
    const [position, setPosition] = useState<number | null>(() => openingFromSelection(findDoc(docs, firstVector)));
    const [source, setSource] = useState<"file" | "pcr">(firstInsert ? "file" : "pcr");
    const [insert, setInsert] = useState<string | null>(firstInsert);
    const [template, setTemplate] = useState<string | null>(firstInsert ?? active?.editorName ?? null);
    const [fwd, setFwd] = useState({ name: "Forward", sequence: "" });
    const [rev, setRev] = useState({ name: "Reverse", sequence: "" });
    const dv = findDoc(docs, vector);
    const di = findDoc(docs, insert);
    const dt = findDoc(docs, template);
    const vData = dataOf(dv);
    const directional = kind === "directional-TOPO";
    const missingTag = (seq: string | undefined) =>
      directional && !!seq && seq.length >= 4 && !seq.toUpperCase().startsWith(DIRECTIONAL_TAG);

    useEffect(() => {
      if (!dv) return onChange(null, "Choose the vector file");
      if (vData?.circular && position === null) {
        return onChange(null, `“${dv.name}” is circular: enter where it is opened (the cloning site)`);
      }
      const common = { kind, vector: docData(dv), position: vData?.circular ? position : null };
      if (source === "file") {
        if (!di) return onChange(null, "Choose the PCR product (a linear file), or amplify a template here");
        return onChange({ ...common, insert: docData(di) });
      }
      if (!dt) return onChange(null, "Choose the template to amplify");
      if (fwd.sequence.length < 12 || rev.sequence.length < 12) {
        return onChange(null, "Enter the forward and reverse primers (at least 12 nt each)");
      }
      onChange({ ...common, template: docData(dt), forward: fwd, reverse: rev });
    }, [kind, dv, position, source, di, dt, fwd, rev]); // eslint-disable-line react-hooks/exhaustive-deps

    const len = vData?.sequence.length ?? 0;
    return (
      <>
        <Field label="Method">
          <div className="bp3-html-select">
            <select data-testid="ta-kind" value={kind} onChange={(e) => setKind(e.target.value)}>
              {kinds.map((k) => <option key={k.value} value={k.value}>{k.label}</option>)}
            </select>
          </div>
        </Field>
        <Field
          label="Vector"
          hint="A linearised vector file, or a circular map plus where it opens"
        >
          <DocSelect
            docs={docs}
            value={vector}
            onChange={(v) => {
              setVector(v);
              setPosition(openingFromSelection(findDoc(docs, v)));
            }}
            testId="vector-doc"
          />
        </Field>
        {vData?.circular && (
          <Field label="Opened after base" hint="The insert goes between this base and the next (0 = before base 1)">
            <div className="vt-inline-row">
              <input
                className="bp3-input bp3-small"
                type="number"
                min={0}
                max={len}
                data-testid="vector-position"
                value={position ?? ""}
                onChange={(e) =>
                  setPosition(e.target.value.trim() === "" ? null : Math.max(0, Math.min(len, Math.round(Number(e.target.value)))))}
              />
              <button className="bp3-button bp3-small" disabled={!dv?.selection()} onClick={() => setPosition(openingFromSelection(dv))}>
                After the selection
              </button>
            </div>
          </Field>
        )}
        <Field label="PCR product">
          <div className="bp3-html-select bp3-fill">
            <select data-testid="insert-source" value={source} onChange={(e) => setSource(e.target.value as "file" | "pcr")}>
              <option value="file">An open file (the product as amplified)</option>
              <option value="pcr">Amplify a template with two primers</option>
            </select>
          </div>
        </Field>
        {source === "file" ? (
          <Field label="PCR product file">
            <DocSelect docs={docs} value={insert} onChange={setInsert} testId="insert-doc" />
            {missingTag(dataOf(di)?.sequence) && (
              <span className="vt-inline-error" data-testid="cacc-missing">
                This product does not start with CACC: directional TOPO needs CACC at the 5′ end of the forward primer.
              </span>
            )}
          </Field>
        ) : (
          <>
            <Field label="Template">
              <DocSelect docs={docs} value={template} onChange={setTemplate} testId="template-doc" />
            </Field>
            <PrimerInput label="Forward primer" template={dt} value={fwd} onChange={setFwd} testId="forward-primer">
              {missingTag(fwd.sequence) && (
                <span className="vt-inline-row vt-inline-error" data-testid="cacc-missing">
                  Directional TOPO needs CACC at the 5′ end of the forward primer.
                  <button
                    className="bp3-button bp3-small"
                    data-testid="add-cacc"
                    onClick={() => setFwd({ ...fwd, sequence: DIRECTIONAL_TAG + fwd.sequence })}
                  >
                    Add CACC
                  </button>
                </span>
              )}
            </PrimerInput>
            <PrimerInput label="Reverse primer" template={dt} value={rev} onChange={setRev} testId="reverse-primer">
              {directional && rev.sequence.toUpperCase().startsWith(DIRECTIONAL_TAG) && (
                <span className="vt-inline-warning">A reverse primer starting with CACC lets the insert go in backwards.</span>
              )}
            </PrimerInput>
          </>
        )}
      </>
    );
  };
}

// ---- Insert multiple fragments / Ligate fragments -----------------------------------------
export function InsertMultipleForm({ docs, active, onChange }: FormProps) {
  const [vector, setVector] = useState<FragmentRow>(() => newRow(docs, defaultVector(docs, active), "digest"));
  const [inserts, setInserts] = useState<FragmentRow[]>(() =>
    defaultInserts(docs, defaultVector(docs, active), 2).map((n) => newRow(docs, n, "digest")));

  useEffect(() => {
    const v = rowBody(docs, vector, "vector");
    if (typeof v === "string") return onChange(null, v);
    const bodies = inserts.map((r, i) => rowBody(docs, r, `insert ${i + 1}`));
    const missing = bodies.find((b): b is string => typeof b === "string");
    if (missing) return onChange(null, missing);
    onChange({ vector: v, inserts: bodies });
  }, [docs, vector, inserts]); // eslint-disable-line react-hooks/exhaustive-deps

  const setInsert = (i: number, row: FragmentRow) => setInserts((rows) => rows.map((r, j) => (j === i ? row : r)));
  return (
    <>
      <FragmentEditor docs={docs} row={vector} title="Vector" role="vector" prefer="digest" testId="vector" onChange={setVector} />
      {inserts.map((row, i) => (
        <FragmentEditor
          key={i}
          docs={docs}
          row={row}
          title={`Insert ${i + 1}`}
          role="insert"
          prefer="digest"
          testId={`insert-${i}`}
          onChange={(r) => setInsert(i, r)}
          onRemove={inserts.length > 1 ? () => setInserts((rows) => rows.filter((_, j) => j !== i)) : undefined}
        />
      ))}
      <button
        className="bp3-button bp3-small"
        data-testid="add-insert"
        onClick={() => setInserts((rows) => [...rows, newRow(docs, null, "digest")])}
      >
        + Add insert
      </button>
    </>
  );
}

export function LigateFragmentsForm({ docs, onChange }: FormProps) {
  const [rows, setRows] = useState<FragmentRow[]>(() =>
    defaultInserts(docs, null, 2).map((n) => newRow(docs, n, "pcr")));
  const [circular, setCircular] = useState(false);

  useEffect(() => {
    const bodies = rows.map((r, i) => rowBody(docs, r, `fragment ${i + 1}`));
    const missing = bodies.find((b): b is string => typeof b === "string");
    if (missing) return onChange(null, missing);
    if (bodies.length < 2 && !circular) {
      return onChange(null, "Add a second fragment, or tick “Close the product into a circle”");
    }
    onChange({ fragments: bodies, circular });
  }, [docs, rows, circular]); // eslint-disable-line react-hooks/exhaustive-deps

  const setRow = (i: number, row: FragmentRow) => setRows((all) => all.map((r, j) => (j === i ? row : r)));
  return (
    <>
      {rows.map((row, i) => (
        <FragmentEditor
          key={i}
          docs={docs}
          row={row}
          title={`Fragment ${i + 1}`}
          role="fragment"
          prefer="pcr"
          testId={`fragment-${i}`}
          onChange={(r) => setRow(i, r)}
          onRemove={rows.length > 1 ? () => setRows((all) => all.filter((_, j) => j !== i)) : undefined}
        />
      ))}
      <button
        className="bp3-button bp3-small"
        data-testid="add-fragment"
        onClick={() => setRows((all) => [...all, newRow(docs, null, "pcr")])}
      >
        + Add fragment
      </button>
      <label className="bp3-control bp3-checkbox">
        <input type="checkbox" data-testid="ligate-circular" checked={circular} onChange={(e) => setCircular(e.target.checked)} />
        <span className="bp3-control-indicator" />
        Close the product into a circle
      </label>
    </>
  );
}

// ---- Delete restriction fragment ---------------------------------------------------------
export function DeleteFragmentForm({ docs, active, onChange }: FormProps) {
  const [docName, setDocName] = useState<string | null>(defaultVector(docs, active));
  const [enz, setEnz] = useState<string[]>([]);
  const d = findDoc(docs, docName);
  useEffect(() => {
    if (!d) return onChange(null, "Choose the file");
    if (!enz.length) return onChange(null, "Choose the enzymes that flank the fragment to delete");
    onChange({ doc: docData(d), enzymes: enz, fragment_index: null });
  }, [d, enz]); // eslint-disable-line react-hooks/exhaustive-deps
  return (
    <>
      <Field label="File">
        <DocSelect docs={docs} value={docName} onChange={setDocName} />
      </Field>
      <Field label="Enzymes flanking the fragment to delete">
        <EnzymePicker value={enz} onChange={setEnz} max={2} sequence={dataOf(d)?.sequence} circular={dataOf(d)?.circular} />
      </Field>
    </>
  );
}

// ---- Anneal oligos -----------------------------------------------------------------------
export function AnnealOligosForm({ docs, active, onChange }: FormProps) {
  const [mode, setMode] = useState<"design" | "anneal">("design");
  const [vector, setVector] = useState<string | null>(defaultVector(docs, active));
  const [enzymes, setEnzymes] = useState<string[]>([]);
  const [insert, setInsert] = useState("");
  const [insertName, setInsertName] = useState("gRNA1");
  const [top, setTop] = useState("");
  const [bottom, setBottom] = useState("");
  const v = findDoc(docs, vector);

  useEffect(() => {
    if (!v) return onChange(null, "Choose the vector the oligos go into");
    if (!enzymes.length) return onChange(null, `Choose the enzyme(s) that cut “${v.name}”`);
    if (mode === "design") {
      if (insert.length < 4) return onChange(null, "Enter the insert (for example a 20 nt protospacer)");
      onChange({ mode, insert, insert_name: insertName || "insert", vector: docData(v), enzymes });
    } else {
      if (top.length < 10 || bottom.length < 10) {
        return onChange(null, "Enter both oligos (at least 10 nt each)");
      }
      onChange({ mode, top, top_name: `${insertName}_top`, bottom, bottom_name: `${insertName}_bottom`,
                 vector: docData(v), enzymes });
    }
  }, [mode, v, enzymes, insert, insertName, top, bottom]); // eslint-disable-line react-hooks/exhaustive-deps

  const dna = (value: string) => value.toUpperCase().replace(/[^ACGT]/g, "");
  return (
    <>
      <Field label="What do you have?">
        <div className="bp3-html-select">
          <select data-testid="anneal-mode" value={mode} onChange={(e) => setMode(e.target.value as "design" | "anneal")}>
            <option value="design">The insert — write the oligos for me</option>
            <option value="anneal">The two oligos already</option>
          </select>
        </div>
      </Field>
      <Field label="Vector">
        <DocSelect docs={docs} value={vector} onChange={setVector} testId="vector-doc" />
      </Field>
      <Field label="Cut vector with" hint="The overhangs are taken from the cut vector itself">
        <EnzymePicker value={enzymes} onChange={setEnzymes} max={2} sequence={v?.data().sequence}
                      circular={v?.data().circular} testId="vector-enzymes" />
      </Field>
      <Field label="Name">
        <input className="bp3-input bp3-small" value={insertName} onChange={(e) => setInsertName(e.target.value)} />
      </Field>
      {mode === "design" ? (
        <Field label="Insert (5′ → 3′, top strand)" hint="e.g. a 20 nt protospacer with its G">
          <input className="bp3-input bp3-small vt-mono" data-testid="anneal-insert"
                 value={insert} onChange={(e) => setInsert(dna(e.target.value))} />
        </Field>
      ) : (
        <>
          <Field label="Oligo 1 (top strand)">
            <input className="bp3-input bp3-small vt-mono" data-testid="anneal-top"
                   value={top} onChange={(e) => setTop(dna(e.target.value))} />
          </Field>
          <Field label="Oligo 2 (bottom strand)">
            <input className="bp3-input bp3-small vt-mono" data-testid="anneal-bottom"
                   value={bottom} onChange={(e) => setBottom(dna(e.target.value))} />
          </Field>
        </>
      )}
    </>
  );
}
