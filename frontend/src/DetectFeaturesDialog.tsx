import { useCallback, useEffect, useMemo, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import { actions } from "@teselagen/ove";
import { api, type OveAnnotation, type SequenceData } from "./api";
import { store } from "./store";
import { IDENTITY_CHOICES, loadPreferences } from "./preferences";
import type { OpenDoc } from "./cloning/types";
import "./DetectFeaturesDialog.css";

const { Dialog, Callout } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

interface Detected {
  name: string;
  type: string;
  start: number;
  end: number;
  strand: number;
  forward: boolean;
  identity: number;
  note: string;
  source: string;
  length: number;
  /** The catalogue entry's length, and the bases of it (1-based) that the match covers. */
  catalogueLength?: number;
  catalogueFrom?: number;
  catalogueTo?: number;
  isNew: boolean;
  /** The existing feature that already covers this match, under whatever name. */
  annotatedAs?: string | null;
  /**
   * Set when that feature is named after a different element although the match is nearly
   * exact: “Annotated as “PGK promoter”, but the sequence is the CMV enhancer …”.
   */
  labelWarning?: string | null;
}

/** "already annotated", or "annotated as “nptI”" when the map uses another name for it. */
function annotatedText(f: Detected): string {
  return f.annotatedAs && f.annotatedAs.toLowerCase() !== f.name.toLowerCase()
    ? `annotated as “${f.annotatedAs}”`
    : "already annotated";
}

/**
 * "bases 4–708 of 711" when the match covers only part of the catalogue entry, so a 100 %
 * over a shorter stretch is not read as the whole element; null when it covers all of it.
 */
function partText(f: Detected): string | null {
  const { catalogueFrom: from, catalogueTo: to, catalogueLength: total } = f;
  if (from == null || to == null || total == null || (from === 1 && to === total)) return null;
  return `bases ${from}–${to} of ${total}`;
}

/** A message as a sentence, so another sentence can follow it. */
const sentence = (text: string) => (/[.!?]$/.test(text.trim()) ? text.trim() : `${text.trim()}.`);

interface DetectResult {
  features: Detected[];
  catalogue_size: number;
  new_count: number;
  /** How many existing features carry a label that may be wrong. */
  label_warnings?: number;
}

/** Existing features whose label may be wrong, from the engine's count or from the rows. */
function doubtfulLabels(result: DetectResult): number {
  return result.label_warnings
    ?? new Set(result.features.filter((f) => f.labelWarning).map((f) => f.annotatedAs)).size;
}

/** Colours follow the convention most plasmid maps use, by feature type. */
const COLOR: Record<string, string> = {
  promoter: "#31a354",
  enhancer: "#74c476",
  CDS: "#f16913",
  repeat_region: "#9e9ac8",
  rep_origin: "#969696",
  polyA_signal: "#e6550d",
  protein_bind: "#3182bd",
  misc_feature: "#bdbdbd",
};

export default function DetectFeaturesDialog(props: { doc: OpenDoc; onClose: () => void }) {
  const [result, setResult] = useState<DetectResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [chosen, setChosen] = useState<Set<string>>(new Set());
  const [busy, setBusy] = useState(true);
  // The threshold from Preferences is the starting point; it can be loosened here for one run.
  const [minIdentity, setMinIdentity] = useState(loadPreferences().detectMinIdentity);
  const data = props.doc.data();

  const key = (f: Detected) => `${f.name}@${f.start}`;

  useEffect(() => {
    setBusy(true);
    setError(null);
    api
      .detectFeatures<DetectResult>({ sequenceData: data, min_identity: minIdentity })
      .then((payload) => {
        setResult(payload);
        // Everything not annotated yet is ticked: one click adds the useful result.
        setChosen(new Set(payload.features.filter((f) => f.isNew).map(key)));
      })
      .catch((e) => setError(String((e as Error).message)))
      .finally(() => setBusy(false));
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [minIdentity]);

  const selected = useMemo(
    () => (result?.features ?? []).filter((f) => chosen.has(key(f))),
    [result, chosen],
  );

  const toggle = (f: Detected) =>
    setChosen((previous) => {
      const next = new Set(previous);
      if (next.has(key(f))) next.delete(key(f));
      else next.add(key(f));
      return next;
    });

  const add = useCallback(() => {
    const current = props.doc.data() as SequenceData;
    const existing = Array.isArray(current.features)
      ? current.features
      : Object.values(current.features ?? {});
    const added: OveAnnotation[] = selected.map((f) => ({
      name: f.name,
      type: f.type,
      start: f.start,
      end: f.end,
      strand: f.strand,
      forward: f.forward,
      color: COLOR[f.type] ?? COLOR.misc_feature,
      notes: {
        note: [f.note || f.name],
        detected_from: [
          `${f.source} (${f.identity.toFixed(1)}% identity${partText(f) ? ` over ${partText(f)}` : ""})`,
        ],
      },
    }));
    store.dispatch(
      actions.updateSequenceData(
        { ...current, features: [...existing, ...added] },
        { editorName: props.doc.editorName },
      ),
    );
    props.onClose();
  }, [props, selected]);

  return (
    <Dialog isOpen onClose={props.onClose} title="Detect Common Features" style={{ width: 920 }}>
      <div className={Classes.DIALOG_BODY} data-testid="detect-dialog">
        {/* Wraps on a narrow screen: the identity control moves under the text, never off the edge. */}
        <div
          className="vt-row"
          style={{ alignItems: "flex-start", justifyContent: "space-between", flexWrap: "wrap",
                   gap: "var(--vt-space-2) var(--vt-space-4)" }}
        >
          <p className="vt-muted" style={{ flex: "1 1 420px", margin: 0 }}>
            Searches both strands for the elements AAV and everyday cloning work use: ITRs,
            promoters, WPRE, polyA signals, markers, origins, fluorescent proteins, tags,
            lox/FRT sites and Cas9. Tick what you want and it is added to the map.
          </p>
          <label style={{ whiteSpace: "nowrap", flex: "none" }} title="How close a match has to be to count"
                 data-testid="detect-identity-control">
            Minimum identity&nbsp;
            <div className="bp3-html-select bp3-small">
              <select
                data-testid="detect-identity"
                value={minIdentity}
                disabled={busy}
                onChange={(e) => setMinIdentity(Number(e.target.value))}
              >
                {IDENTITY_CHOICES.map((p) => <option key={p} value={p}>{p} %</option>)}
              </select>
            </div>
          </label>
        </div>
        {busy && <p className="vt-muted" data-testid="detect-busy">Searching the sequence…</p>}
        {error && (
          <Callout intent="danger" title="Could not search">
            {sentence(error)} Try again, or lower the minimum identity above.
          </Callout>
        )}
        {result && !busy && result.features.length === 0 && (
          <Callout intent="primary" data-testid="detect-empty">
            Nothing from the catalogue of {result.catalogue_size} elements was found in{" "}
            <b>{props.doc.name}</b> at {minIdentity} % identity. Lower the threshold above to
            catch variants, or annotate the region by hand with Features → Add Feature…
          </Callout>
        )}
        {result && !busy && doubtfulLabels(result) > 0 && (
          <div className="vt-notice vt-notice-warning" role="status" data-testid="detect-label-summary">
            <span className="vt-notice-title">
              {doubtfulLabels(result)} label{doubtfulLabels(result) === 1 ? "" : "s"} may be wrong
            </span>
            {doubtfulLabels(result) === 1
              ? "The sequence under it is a different element from the one the label names. "
              : "The sequence under each is a different element from the one its label names. "}
            The rows marked below say which; they are not ticked, so nothing changes unless you
            tick one to add the element under its own name.
          </div>
        )}
        {result && result.features.length > 0 && (
          <>
            <p className="vt-muted">
              {result.features.length} of {result.catalogue_size} catalogue features found in{" "}
              <b>{props.doc.name}</b>; {result.new_count} are not annotated yet. Each one comes
              from a public GenBank record, named in the last column.
            </p>
            <table className="bp3-html-table bp3-html-table-condensed vt-table" data-testid="detect-table">
              <thead>
                <tr>
                  <th>
                    <input
                      type="checkbox"
                      data-testid="detect-all"
                      checked={selected.length === result.features.length}
                      onChange={(e) =>
                        setChosen(e.target.checked
                          ? new Set(result.features.map(key))
                          : new Set())}
                    />
                  </th>
                  <th>Feature</th><th>Type</th><th>Location</th><th>Size</th>
                  <th>Identity</th><th>From</th>
                </tr>
              </thead>
              <tbody>
                {result.features.map((f) => (
                  <tr
                    key={key(f)}
                    className={f.labelWarning ? "vt-detect-doubt" : f.isNew ? "" : "vt-muted"}
                    data-testid="detect-row"
                  >
                    <td>
                      <input type="checkbox" checked={chosen.has(key(f))} onChange={() => toggle(f)}
                             aria-label={`Add ${f.name}`} />
                    </td>
                    <td>
                      <span className="vt-swatch" style={{ background: COLOR[f.type] ?? COLOR.misc_feature }} />
                      {f.name}
                      {f.labelWarning ? (
                        <span className="vt-warn vt-detect-doubt-text" data-testid="detect-label-warning">
                          {f.labelWarning}
                        </span>
                      ) : !f.isNew && (
                        <span className="vt-muted" data-testid="detect-annotated"> · {annotatedText(f)}</span>
                      )}
                    </td>
                    <td>{f.type}</td>
                    <td className="vt-num">
                      {f.start + 1}…{f.end + 1} {f.forward ? "→" : "←"}
                    </td>
                    <td className="vt-num">
                      {f.length} bp
                      {partText(f) && (
                        <span className="vt-muted" data-testid="detect-part"
                              title={`The catalogue entry is ${f.catalogueLength} bp; only these bases of it match here`}>
                          {" "}· {partText(f)}
                        </span>
                      )}
                    </td>
                    <td className={f.identity === 100 ? "vt-num vt-ok" : "vt-num"}>{f.identity.toFixed(1)} %</td>
                    <td className="vt-muted" title={f.note}>{f.source}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </>
        )}
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button onClick={props.onClose}>Cancel</Button>
          <Button intent="primary" data-testid="detect-add" disabled={!selected.length} onClick={add}>
            Add {selected.length} feature{selected.length === 1 ? "" : "s"}
          </Button>
        </div>
      </div>
    </Dialog>
  );
}
