import { useEffect, useRef, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import { SimpleCircularOrLinearView } from "@teselagen/ove";
import { api } from "../api";
import type { CloningAction } from "./registry";
import type { CloningResult, Oligo, OpenDoc, Product } from "./types";

// Blueprint 3 typings predate React 18 (no implicit `children`).
const { Dialog, Callout } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

function oligosAsText(oligos: Oligo[]) {
  return oligos.map((o) => `${o.name}\t${o.sequence}`).join("\n");
}

export default function CloningDialog(props: {
  action: CloningAction;
  docs: OpenDoc[];
  active: OpenDoc | null;
  onClose: () => void;
  onCreate: (product: Product, result: CloningResult) => void;
}) {
  const { action } = props;
  const [body, setBody] = useState<object | null>(null);
  const [result, setResult] = useState<CloningResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [chosen, setChosen] = useState(0);
  const request = useRef(0);

  useEffect(() => {
    setError(null);
    if (!body) {
      setResult(null);
      return;
    }
    const id = ++request.current;
    const t = setTimeout(() => {
      setBusy(true);
      api
        .cloning<CloningResult>(action.endpoint, body)
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
    }, 350);
    return () => clearTimeout(t);
  }, [action.endpoint, body]);

  const product = result?.products[chosen] ?? null;
  const oligos = result ? (result.oligos.length ? result.oligos : product?.oligos ?? []) : [];
  const { Form } = action;

  return (
    <Dialog isOpen onClose={props.onClose} title={action.title} style={{ width: 1100, maxWidth: "95vw" }} canOutsideClickClose={false}>
      <div className={`${Classes.DIALOG_BODY} vt-cloning`} data-testid="cloning-dialog">
        <div className="vt-cloning-form">
          <p className="vt-muted">{action.help}</p>
          <Form docs={props.docs} active={props.active} onChange={setBody} />
        </div>
        <div className="vt-cloning-result">
          {busy && <p className="vt-muted">Simulating…</p>}
          {!body && !busy && <p className="vt-muted">Fill in the inputs on the left to simulate the product.</p>}
          {error && (
            <Callout intent="danger" title="No product" data-testid="cloning-error">
              {error}
            </Callout>
          )}
          {result && product && (
            <>
              {result.products.length > 1 && (
                <div className="vt-products">
                  {result.products.map((p, i) => (
                    <label key={i} className="bp3-control bp3-radio bp3-inline">
                      <input type="radio" checked={i === chosen} onChange={() => setChosen(i)} />
                      <span className="bp3-control-indicator" />
                      {i === 0 ? "Expected: " : `Alternative ${i}: `}
                      {p.doc.sequence.length} bp {p.doc.circular ? "circular" : "linear"}
                    </label>
                  ))}
                </div>
              )}
              <div className="vt-preview" data-testid="cloning-preview">
                <SimpleCircularOrLinearView
                  sequenceData={product.doc}
                  editorName={`preview-${action.id}`}
                  annotationVisibility={{ features: true, primers: true, cutsites: false, parts: false }}
                />
              </div>
              <p data-testid="cloning-product">
                <b>{product.doc.name}</b> — {product.doc.sequence.length} bp, {product.doc.circular ? "circular" : "linear"}
                <br />
                <span className="vt-muted">{product.details}</span>
              </p>
              {result.warnings.length > 0 && (
                <Callout intent="warning" title="Check before ordering" data-testid="cloning-warnings">
                  <ul>{result.warnings.map((w, i) => <li key={i}>{w}</li>)}</ul>
                </Callout>
              )}
              {oligos.length > 0 && (
                <div className="vt-oligos">
                  <div className="vt-row">
                    <b>Primers to order</b>
                    <Button small minimal icon="duplicate" onClick={() => navigator.clipboard?.writeText(oligosAsText(oligos))}>
                      Copy all (name ⇥ sequence)
                    </Button>
                  </div>
                  <table className="bp3-html-table bp3-html-table-condensed vt-table" data-testid="cloning-oligos">
                    <thead>
                      <tr><th>Name</th><th>Sequence 5′→3′</th><th>Length</th><th>Annealing</th><th>Tm</th></tr>
                    </thead>
                    <tbody>
                      {oligos.map((o) => (
                        <tr key={o.name}>
                          <td>{o.name}</td>
                          <td className="vt-mono vt-break">{o.sequence}</td>
                          <td>{o.sequence.length}</td>
                          <td>{o.annealing_length ?? "—"}</td>
                          <td>{o.tm != null ? `${o.tm.toFixed(1)} °C` : "—"}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}
            </>
          )}
        </div>
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button onClick={props.onClose}>Cancel</Button>
          <Button
            intent="primary"
            data-testid="cloning-create"
            disabled={!product || busy}
            onClick={() => product && result && props.onCreate(product, result)}
          >
            Create Product
          </Button>
        </div>
      </div>
    </Dialog>
  );
}
