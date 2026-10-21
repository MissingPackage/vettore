// Tools → Find Diagnostic Digest: the single or double digest whose gel tells the expected clone
// apart from the usual wrong ones. The engine (vettore.digest_planner) ranks the digests; this
// dialog chooses what to compare and which enzymes to try, and leaves out the enzymes that
// methylation blocks (REBASE's table, as in the Enzymes view).
import { useEffect, useMemo, useRef, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import { api, type SequenceData } from "../api";
import { DocSelect, Field, findDoc } from "../cloning/inputs";
import type { OpenDoc } from "../cloning/types";
import { useCopyText } from "../cloning/useCopyText";
import GelDialog, { GelPicture, type LaneSetup } from "../GelDialog";
import { ENZYMES } from "../enzymeLogic";
import { messageOf } from "../files/notices";
import { usePreferences } from "../preferences";
import {
  blockedByMethylation,
  protocolText,
  type PlannerResponse,
  type PlannerResult,
} from "./digestPlannerText";
import type { Tool, ToolProps } from "./types";
import "./digestPlanner.css";

// Blueprint 3 typings predate React 18 (no implicit `children`).
const { Dialog } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

const EXPECTED = "Expected clone";
const COMMERCIAL_COUNT = Object.keys(ENZYMES).length;

/** A molecule to tell apart: a wrong clone of the cloning simulation, or an open file. */
interface WrongEntry {
  key: string;
  name: string;
  data: SequenceData;
  /** The open document it is, when it is one. */
  editorName: string | null;
}

function fromSimulation(doc: OpenDoc | null): WrongEntry[] {
  return (doc?.wrongClones?.() ?? []).map((w, i) => ({
    key: `sim-${i}`, name: w.name, data: w.data, editorName: null,
  }));
}

const size = (data: SequenceData) =>
  `${data.sequence.length.toLocaleString("en-US")} bp, ${data.circular ? "circular" : "linear"}`;

function DigestPlannerDialog(props: ToolProps) {
  const prefs = usePreferences();
  const [expectedId, setExpectedId] = useState<string | null>(props.active?.editorName ?? null);
  const expected = findDoc(props.docs, expectedId);
  const [wrong, setWrong] = useState<WrongEntry[]>(() => fromSimulation(expected));
  const [useMine, setUseMine] = useState(prefs.myEnzymes.length > 0);
  const [damFree, setDamFree] = useState(false);
  const [percent, setPercent] = useState(prefs.gelPercent);
  const [ladder, setLadder] = useState(prefs.gelLadder);
  const [ladders, setLadders] = useState<{ id: string; label: string }[]>([]);
  const [percentages, setPercentages] = useState<number[]>([0.7, 0.8, 1.0, 1.2, 1.5, 2.0]);
  const [response, setResponse] = useState<PlannerResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [chosen, setChosen] = useState(0);
  const [copied, setCopied] = useState<string | null>(null);
  const [gel, setGel] = useState<{ docs: OpenDoc[]; lanes: LaneSetup[] } | null>(null);
  const { copy, fallbackDialog } = useCopyText();
  const request = useRef(0);

  useEffect(() => {
    api.gelLadders()
      .then((r) => {
        setLadders(r.ladders);
        setPercentages(r.percentages);
      })
      .catch(() => undefined);
  }, []);

  // Another expected clone brings its own wrong clones.
  const firstExpected = useRef(expectedId);
  useEffect(() => {
    if (expectedId === firstExpected.current) return;
    firstExpected.current = expectedId;
    setWrong(fromSimulation(findDoc(props.docs, expectedId)));
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [expectedId]);

  const mine = useMine && prefs.myEnzymes.length > 0;
  const expectedData = expected?.data() ?? null;
  const molecules = useMemo(
    () => (expectedData ? [expectedData, ...wrong.map((w) => w.data)] : []),
    [expectedData, wrong],
  );
  const blocked = useMemo(() => (damFree ? [] : blockedByMethylation(molecules)), [damFree, molecules]);
  const tried = mine ? prefs.myEnzymes : null;
  const leftOut = tried ? blocked.filter((b) => tried.includes(b)) : blocked;

  const body = useMemo(() => {
    if (!expectedData || !wrong.length) return null;
    return {
      expected: { name: EXPECTED, sequenceData: expectedData },
      alternatives: wrong.map((w) => ({ name: w.name, sequenceData: w.data })),
      enzymes: tried,
      excluded: blocked,
      ladder,
      percent,
    };
  }, [expectedData, wrong, tried, blocked, ladder, percent]);
  const bodyKey = body ? JSON.stringify(body) : "";

  useEffect(() => {
    setError(null);
    if (!bodyKey) {
      setResponse(null);
      return;
    }
    const id = ++request.current;
    const t = setTimeout(() => {
      setBusy(true);
      api.diagnosticDigest(JSON.parse(bodyKey))
        .then((r) => {
          if (id !== request.current) return;
          setResponse(r);
          setChosen(0);
        })
        .catch((e) => {
          if (id !== request.current) return;
          setResponse(null);
          setError(messageOf(e));
        })
        .finally(() => id === request.current && setBusy(false));
    }, 250);
    return () => clearTimeout(t);
  }, [bodyKey]);

  const result: PlannerResult | null = response?.results[chosen] ?? null;
  useEffect(() => setCopied(null), [result]);

  const addable = props.docs.filter((d) => d.editorName !== expectedId
    && !wrong.some((w) => w.editorName === d.editorName));

  const copyProtocol = async () => {
    if (!result || !response || !expected) return;
    const circular = [expectedData?.circular ?? true, ...wrong.map((w) => w.data.circular)];
    const text = protocolText(expected.name, result, { percent: response.percent, ladderLabel: response.ladder_label }, circular);
    if (await copy(text, "Diagnostic digest")) setCopied("Copied: paste it into your notebook or protocol.");
  };

  const sendToGel = () => {
    if (!result || !expected) return;
    // Wrong clones that are not open files join the gel as lanes of their own.
    const virtual: OpenDoc[] = wrong.filter((w) => !w.editorName).map((w) => ({
      editorName: `diagnostic-${w.key}`, name: w.name, data: () => w.data, selection: () => null,
    }));
    const lanes: LaneSetup[] = [
      { editorName: expected.editorName, enzymes: result.enzymes },
      ...wrong.map((w) => ({ editorName: w.editorName ?? `diagnostic-${w.key}`, enzymes: result.enzymes })),
    ];
    setGel({ docs: [...props.docs, ...virtual], lanes });
  };

  if (gel) {
    return (
      <GelDialog docs={gel.docs} active={expected} onClose={props.onClose}
                 initial={{ lanes: gel.lanes, ladder, percent }} />
    );
  }

  return (
    <Dialog isOpen onClose={props.onClose} title="Find Diagnostic Digest" style={{ width: 1100, maxWidth: "95vw" }}
            canOutsideClickClose={false}>
      <div className={`${Classes.DIALOG_BODY} vt-cloning vt-digest-planner`} data-testid="digest-dialog">
        <div className="vt-cloning-form">
          <p className="vt-dialog-intro">
            Finds the single or double digest whose gel tells your clone apart from the usual wrong
            ones, so one gel shows which colonies are right.
          </p>
          <Field label="Expected clone" hint="The construct you want: usually the product of your cloning.">
            <DocSelect docs={props.docs} value={expectedId} onChange={setExpectedId} testId="digest-expected" />
          </Field>

          <div className="vt-field" data-testid="digest-wrong-field">
            <span className="vt-field-label">Tell apart from</span>
            {wrong.length === 0 && (
              <span className="vt-field-hint" data-testid="digest-wrong-empty">
                Nothing yet. A product of Vettore's cloning lists its wrong clones here by itself
                (insert reversed, empty vector); otherwise add the parent vector or another open file.
              </span>
            )}
            <ul className="vt-digest-wrong" data-testid="digest-wrong">
              {wrong.map((w) => (
                <li key={w.key} className="vt-card">
                  <span>
                    <b>{w.name}</b>{" "}
                    <span className="vt-muted">
                      {size(w.data)} · {w.editorName ? "open file" : "from the cloning simulation"}
                    </span>
                  </span>
                  <Button small minimal aria-label={`Remove ${w.name}`} data-testid="digest-wrong-remove"
                          onClick={() => setWrong((list) => list.filter((x) => x.key !== w.key))}>
                    Remove
                  </Button>
                </li>
              ))}
            </ul>
            {addable.length > 0 && (
              <div className="bp3-html-select bp3-fill">
                <select
                  data-testid="digest-add"
                  value=""
                  onChange={(e) => {
                    const doc = findDoc(props.docs, e.target.value);
                    if (doc) setWrong((list) => [...list, { key: `file-${doc.editorName}`, name: doc.name, data: doc.data(), editorName: doc.editorName }]);
                  }}
                >
                  <option value="" disabled>Add an open file…</option>
                  {addable.map((d) => <option key={d.editorName} value={d.editorName}>{d.name} ({size(d.data())})</option>)}
                </select>
              </div>
            )}
          </div>

          <div className="vt-field">
            <span className="vt-field-label">Enzymes to try</span>
            <label className="bp3-control bp3-radio">
              <input type="radio" name="digest-set" data-testid="digest-set-mine" checked={mine}
                     disabled={!prefs.myEnzymes.length} onChange={() => setUseMine(true)} />
              <span className="bp3-control-indicator" />
              My enzymes ({prefs.myEnzymes.length})
            </label>
            <label className="bp3-control bp3-radio">
              <input type="radio" name="digest-set" data-testid="digest-set-commercial" checked={!mine}
                     onChange={() => setUseMine(false)} />
              <span className="bp3-control-indicator" />
              All commercial enzymes ({COMMERCIAL_COUNT})
            </label>
            {!prefs.myEnzymes.length && (
              <span className="vt-field-hint" data-testid="digest-mine-hint">
                List the enzymes in your freezer in File → Preferences → My enzymes, and Vettore tries those first.
              </span>
            )}
          </div>

          <div className="vt-field">
            <span className="vt-field-label">Methylation</span>
            <label className="bp3-control bp3-checkbox">
              <input type="checkbox" data-testid="digest-dam-free" checked={damFree}
                     onChange={(e) => setDamFree(e.target.checked)} />
              <span className="bp3-control-indicator" />
              DNA from a dam−/dcm− strain
            </label>
            <span className="vt-field-hint" data-testid="digest-methylation-hint">
              {damFree
                ? "Every enzyme is tried, also those that Dam or Dcm methylation would block."
                : leftOut.length
                  ? `Left out: ${leftOut.length > 12 ? `${leftOut.slice(0, 12).join(", ")} and ${leftOut.length - 12} more` : leftOut.join(", ")}, blocked by Dam or Dcm methylation at a site on these plasmids (plasmid DNA from the usual E. coli strains).`
                  : "No enzyme is blocked by Dam or Dcm methylation on these plasmids."}
            </span>
          </div>

          <div className="vt-row vt-digest-gel-setup">
            <Field label="Agarose (%)">
              <div className="bp3-html-select bp3-fill">
                <select data-testid="digest-percent" value={percent} onChange={(e) => setPercent(Number(e.target.value))}>
                  {percentages.map((p) => <option key={p} value={p}>{p.toFixed(1)} %</option>)}
                </select>
              </div>
            </Field>
            <Field label="Ladder">
              <div className="bp3-html-select bp3-fill">
                <select data-testid="digest-ladder" value={ladder} onChange={(e) => setLadder(e.target.value)}>
                  {ladders.length === 0 && <option value={ladder}>{ladder}</option>}
                  {ladders.map((l) => <option key={l.id} value={l.id}>{l.label}</option>)}
                </select>
              </div>
            </Field>
          </div>
        </div>

        <div className="vt-cloning-result">
          {busy && <span className="vt-busy" role="status" data-testid="digest-busy">Comparing digests…</span>}
          {!expected && (
            <div className="vt-empty-state">
              <p className="vt-empty-state-title">No expected clone</p>
              <p>Open the construct you want, then choose it on the left.</p>
            </div>
          )}
          {expected && !wrong.length && (
            <div className="vt-empty-state" data-testid="digest-nothing-to-compare">
              <p className="vt-empty-state-title">Nothing to tell apart yet</p>
              <p>Add the parent vector or another wrong clone on the left. A product of Vettore's cloning brings its wrong clones with it.</p>
            </div>
          )}
          {error && <p className="vt-error" role="alert" data-testid="digest-error">{error}</p>}
          {response && !response.results.length && response.reason && (
            <div className="vt-notice vt-notice-warning" role="status" data-testid="digest-none">
              <span className="vt-notice-title">No digest tells them apart</span>
              {response.reason}
              {mine && (
                <div>
                  <Button small data-testid="digest-try-commercial" onClick={() => setUseMine(false)}>
                    Try all commercial enzymes
                  </Button>
                </div>
              )}
            </div>
          )}
          {response && response.results.length > 0 && (
            <>
              <p className="vt-muted" data-testid="digest-summary">
                {response.diagnostic_count > response.results.length
                  ? `The best ${response.results.length} of ${response.diagnostic_count} digests that tell every wrong clone apart`
                  : `${response.results.length} digest${response.results.length === 1 ? "" : "s"} tell every wrong clone apart`}
                {" "}({response.enzymes_tried} enzymes, alone and in pairs). Choose one to see its gel.
              </p>
              <ol className="vt-digest-results" data-testid="digest-results">
                {response.results.map((r, i) => (
                  <li key={r.enzymes.join("+")}>
                    <label className={`vt-digest-result${i === chosen ? " vt-digest-result-chosen" : ""}`}>
                      <input type="radio" name="digest-result" checked={i === chosen} onChange={() => setChosen(i)}
                             data-testid="digest-result" />
                      <span>
                        <b>{r.enzymes.join(" + ")}</b>
                        {!r.clear && <span className="vt-warn"> · close bands, check on the gel</span>}
                        <br />
                        <span className="vt-digest-sentence">{r.sentence}</span>
                        {r.same_pattern.length > 0 && (
                          <span className="vt-muted vt-digest-same">
                            {" "}Same bands with {r.same_pattern.slice(0, 4).map((e) => e.join(" + ")).join(", ")}
                            {r.same_pattern.length > 4 ? ` and ${r.same_pattern.length - 4} more` : ""}.
                          </span>
                        )}
                      </span>
                    </label>
                  </li>
                ))}
              </ol>
              {result && (
                <div className="vt-digest-gel">
                  <GelPicture testId="digest-gel" ladderLabel={response.ladder_label} ladder={response.ladder}
                              lanes={result.lanes.map((l) => ({ name: `${l.name} / ${result.enzymes.join(" + ")}`, bands: l.bands }))} />
                </div>
              )}
            </>
          )}
          {copied && <span className="vt-copied" role="status" data-testid="digest-copied">{copied}</span>}
        </div>
      </div>
      {fallbackDialog}
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button onClick={props.onClose}>Close</Button>
          <Button data-testid="digest-copy" disabled={!result || busy} onClick={copyProtocol}
                  title="Enzymes, gel and the expected bands of every outcome, as text">
            Copy as Protocol
          </Button>
          <Button intent="primary" data-testid="digest-send-to-gel" disabled={!result || busy} onClick={sendToGel}>
            Send to Gel Simulation
          </Button>
        </div>
      </div>
    </Dialog>
  );
}

export const tools: Tool[] = [
  { id: "diagnostic-digest", menu: "Tools", label: "Find Diagnostic Digest", requiresDocument: true, Component: DigestPlannerDialog },
];
