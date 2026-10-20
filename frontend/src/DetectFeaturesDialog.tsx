import { useCallback, useEffect, useMemo, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import { actions } from "@teselagen/ove";
import { api, type OveAnnotation, type SequenceData } from "./api";
import { store } from "./store";
import type { OpenDoc } from "./cloning/types";

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
  isNew: boolean;
}

interface DetectResult {
  features: Detected[];
  catalogue_size: number;
  new_count: number;
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
  const data = props.doc.data();

  const key = (f: Detected) => `${f.name}@${f.start}`;

  useEffect(() => {
    api
      .detectFeatures<DetectResult>({ sequenceData: data })
      .then((payload) => {
        setResult(payload);
        setChosen(new Set(payload.features.filter((f) => f.isNew).map(key)));
      })
      .catch((e) => setError(String((e as Error).message)))
      .finally(() => setBusy(false));
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

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
        detected_from: [`${f.source} (${f.identity.toFixed(1)}% identity)`],
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
        {busy && <p className="vt-muted">Searching…</p>}
        {error && <Callout intent="danger" title="Could not search">{error}</Callout>}
        {result && result.features.length === 0 && (
          <Callout intent="primary">
            No feature from the catalogue was found in {props.doc.name}.
          </Callout>
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
                  <tr key={key(f)} className={f.isNew ? "" : "vt-muted"}>
                    <td>
                      <input type="checkbox" checked={chosen.has(key(f))} onChange={() => toggle(f)} />
                    </td>
                    <td>
                      <span className="vt-swatch" style={{ background: COLOR[f.type] ?? COLOR.misc_feature }} />
                      {f.name}
                      {!f.isNew && <span className="vt-muted"> · already annotated</span>}
                    </td>
                    <td>{f.type}</td>
                    <td>
                      {f.start + 1}…{f.end + 1} {f.forward ? "→" : "←"}
                    </td>
                    <td>{f.length} bp</td>
                    <td className={f.identity === 100 ? "vt-ok" : ""}>{f.identity.toFixed(1)} %</td>
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
