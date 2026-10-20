import { useCallback, useEffect, useRef, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes, Spinner } from "@blueprintjs/core";
import { SimpleCircularOrLinearView } from "@teselagen/ove";
import { api } from "../api";
import { defaultProductName } from "./productName";
import type { CloningAction } from "./registry";
import type { CloningResult, Oligo, OpenDoc, Product } from "./types";
import { useCopyText } from "./useCopyText";
import "./cloning.css";

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
  const [missing, setMissing] = useState<string | null>(null);
  const [name, setName] = useState("");
  const [result, setResult] = useState<CloningResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  // Until the user changes an input, the inputs are Vettore's guesses: a failure then says what
  // to choose, in the neutral box, rather than showing a red error.
  const [touched, setTouched] = useState(false);
  const [chosen, setChosen] = useState(0);
  const [copied, setCopied] = useState<string | null>(null);
  const { copy, fallbackDialog } = useCopyText();
  const request = useRef(0);
  // Depend on the serialised body: forms rebuild it on every render.
  const bodyKey = body ? JSON.stringify(body) : "";

  useEffect(() => {
    setError(null);
    if (!bodyKey) {
      setResult(null);
      return;
    }
    const payload = JSON.parse(bodyKey);
    const id = ++request.current;
    const t = setTimeout(() => {
      setBusy(true);
      api
        .cloning<CloningResult>(action.endpoint, payload)
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
  }, [action.endpoint, bodyKey]);

  // Stable: forms report from their effects.
  const onFormChange = useCallback((next: object | null, why?: string) => {
    setBody(next);
    setMissing(next ? null : why ?? null);
  }, []);

  const product = result?.products[chosen] ?? null;
  const suggestedName = product
    ? defaultProductName(product.doc.name, product.operation, props.docs.map((d) => d.name))
    : "";
  useEffect(() => setName(suggestedName), [suggestedName]);
  const oligos = result ? (result.oligos.length ? result.oligos : product?.oligos ?? []) : [];
  const oligoText = oligosAsText(oligos);
  // A confirmation belongs to the primers it was given for.
  useEffect(() => setCopied(null), [oligoText]);
  const copyOligos = async () => {
    if (await copy(oligoText, "Primers to order")) {
      setCopied(`Copied ${oligos.length} primer${oligos.length === 1 ? "" : "s"}: paste them into the order form.`);
    }
  };
  const { Form } = action;

  /**
   * Open every product of the simulation, each in its own tab with its own History entry: a
   * clone and its by-product, or the two orientations of an insert. The product the user
   * named keeps that name, the others take the name they are listed under.
   */
  const createAll = (all: CloningResult) => {
    const taken = props.docs.map((d) => d.name);
    all.products.forEach((p, i) => {
      const chosenName = i === chosen ? name.trim() || suggestedName : "";
      const productName = chosenName || defaultProductName(p.doc.name, p.operation, taken);
      taken.push(productName);
      props.onCreate({ ...p, doc: { ...p.doc, name: productName } }, all);
    });
  };

  return (
    <Dialog isOpen onClose={props.onClose} title={action.title} style={{ width: 1100, maxWidth: "95vw" }} canOutsideClickClose={false}>
      <div className={`${Classes.DIALOG_BODY} vt-cloning`} data-testid="cloning-dialog">
        <div
          className="vt-cloning-form"
          onChangeCapture={() => setTouched(true)}
          onClickCapture={(e) => (e.target as HTMLElement).closest("button") && setTouched(true)}
        >
          <p className="vt-help" data-testid="cloning-help">{action.help}</p>
          {action.details && <p className="vt-help-more">{action.details}</p>}
          <Form docs={props.docs} active={props.active} onChange={onFormChange} />
        </div>
        <div className="vt-cloning-result">
          {busy && (
            <div className="vt-busy" data-testid="cloning-busy">
              <Spinner size={16} /> Simulating…
            </div>
          )}
          {!bodyKey && !busy && (
            <p className="vt-missing" data-testid="cloning-missing">
              {missing ?? "Fill in the inputs on the left to simulate the product."}
            </p>
          )}
          {error && touched && (
            <Callout intent="danger" title="No product" data-testid="cloning-error">
              {error}
            </Callout>
          )}
          {error && !touched && (
            <p className="vt-missing" data-testid="cloning-missing">
              Choose the inputs on the left: the open files picked for you do not give a product.
              <span className="vt-missing-why">{error}</span>
            </p>
          )}
          {result && product && (
            <>
              {result.products.length > 1 && (
                <div className="vt-product-choice" data-testid="cloning-products">
                  <b>
                    {result.products.length} possible products: choose the one to create, or
                    create them all in their own tabs
                  </b>
                  {result.products.map((p, i) => (
                    <label key={i} className="bp3-control bp3-radio">
                      <input type="radio" checked={i === chosen} onChange={() => setChosen(i)} />
                      <span className="bp3-control-indicator" />
                      {i === 0 ? "Expected: " : `Alternative ${i}: `}
                      {p.label ?? `${p.doc.sequence.length} bp ${p.doc.circular ? "circular" : "linear"}`}
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
              <label className="vt-product-name">
                <span>Name of the new file</span>
                <input
                  className="bp3-input bp3-small"
                  data-testid="cloning-name"
                  value={name}
                  onChange={(e) => setName(e.target.value)}
                />
              </label>
              {result.warnings.length > 0 && (
                <Callout intent="warning" title="Check before ordering" data-testid="cloning-warnings">
                  <ul>{result.warnings.map((w, i) => <li key={i}>{w}</li>)}</ul>
                </Callout>
              )}
              {oligos.length > 0 && (
                <div className="vt-oligos">
                  <div className="vt-row">
                    <b>Primers to order</b>
                    <Button small minimal icon="duplicate" data-testid="cloning-copy" onClick={copyOligos}>
                      Copy all (name ⇥ sequence)
                    </Button>
                    {copied && <span className="vt-copied" role="status" data-testid="cloning-copied">{copied}</span>}
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
      {fallbackDialog}
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button onClick={props.onClose}>Cancel</Button>
          {result && result.products.length > 1 && (
            <Button
              data-testid="cloning-create-all"
              disabled={busy}
              title="Open every product of this simulation in its own tab"
              onClick={() => createAll(result)}
            >
              Create All ({result.products.length})
            </Button>
          )}
          <Button
            intent="primary"
            data-testid="cloning-create"
            disabled={!product || busy}
            onClick={() => product && result
              && props.onCreate({ ...product, doc: { ...product.doc, name: name.trim() || suggestedName } }, result)}
          >
            Create Product
          </Button>
        </div>
      </div>
    </Dialog>
  );
}
