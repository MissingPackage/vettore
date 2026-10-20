import { useEffect, useRef, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes, InputGroup, Radio, TextArea } from "@blueprintjs/core";
import { api, type OpenResult } from "../api";
import { messageOf } from "./notices";
import "./PasteDialog.css";

// Blueprint 3 typings predate React 18 (no implicit `children`); container components are untyped here.
const { Dialog, FormGroup, RadioGroup } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

/** Formats whose text says nothing about topology; the dialog's choice applies to them. */
const TOPOLOGY_FROM_DIALOG = new Set(["raw", "fasta"]);

type Check =
  | { state: "empty" }
  | { state: "checking" }
  | { state: "ok"; result: OpenResult }
  | { state: "error"; message: string };

function summary(result: OpenResult): string {
  const docs = result.documents;
  const what = result.format_label ?? "sequence";
  if (docs.length > 1) return `${what}: ${docs.length} sequences (${docs.map((d) => d.name).join(", ")})`;
  const [doc] = docs;
  return `${what}: ${doc.sequence.length.toLocaleString("en-US")} bp${result.format === "raw" ? "" : `, ${doc.circular ? "circular" : "linear"}`}`;
}

/**
 * "New File from Clipboard" when the clipboard cannot be read directly, or when what it held
 * could not be opened: a box to paste into, checked by the engine as the user types.
 */
export default function PasteDialog(props: {
  initialText: string;
  /** Why the box is shown, e.g. the clipboard could not be read. */
  reason?: string;
  onClose: () => void;
  onOpen: (result: OpenResult) => void;
}) {
  const [text, setText] = useState(props.initialText);
  const [name, setName] = useState("");
  const [circular, setCircular] = useState(false);
  const [check, setCheck] = useState<Check>({ state: "empty" });
  const request = useRef(0);

  // Ask the engine what the text is, shortly after the user stops typing.
  useEffect(() => {
    const id = ++request.current;
    if (!text.trim()) {
      setCheck({ state: "empty" });
      return;
    }
    setCheck({ state: "checking" });
    const timer = setTimeout(() => {
      api
        .openText(text, name)
        .then((result) => id === request.current && setCheck({ state: "ok", result }))
        .catch((e) => id === request.current && setCheck({ state: "error", message: messageOf(e) }));
    }, 300);
    return () => clearTimeout(timer);
  }, [text, name]);

  const open = () => {
    if (check.state !== "ok") return;
    const { result } = check;
    const documents = TOPOLOGY_FROM_DIALOG.has(result.format ?? "")
      ? result.documents.map((d) => ({ ...d, circular }))
      : result.documents;
    props.onOpen({ ...result, documents });
  };

  const status =
    check.state === "ok" ? (
      <span className="vt-paste-ok" data-testid="paste-status">Recognised as {summary(check.result)}</span>
    ) : check.state === "error" ? (
      <span className="vt-paste-error" data-testid="paste-error" role="alert">{check.message}</span>
    ) : check.state === "checking" ? (
      <span data-testid="paste-status">Checking…</span>
    ) : (
      <span data-testid="paste-status">Paste with Ctrl+V.</span>
    );

  return (
    <Dialog isOpen onClose={props.onClose} title="New File from Clipboard" style={{ width: 640 }}>
      <div className={Classes.DIALOG_BODY} data-testid="paste-dialog">
        <p className="vt-paste-intro">
          Paste a DNA sequence, a FASTA record or a whole GenBank or EMBL record. Vettore recognises
          the format and opens it as a new file.
        </p>
        {props.reason && <p className="vt-paste-reason">{props.reason}</p>}
        <FormGroup label="Text" labelFor="vt-paste-text" helperText={status}
                   intent={check.state === "error" ? "danger" : "none"}>
          <TextArea
            id="vt-paste-text"
            data-testid="paste-text"
            autoFocus
            fill
            growVertically={false}
            className="vt-paste-text"
            intent={check.state === "error" ? "danger" : "none"}
            value={text}
            onChange={(e: any) => setText(e.target.value)}
          />
        </FormGroup>
        <div className="vt-paste-options">
          <FormGroup label="Name" labelFor="vt-paste-name"
                     helperText="For a plain sequence; records keep their own names.">
            <InputGroup id="vt-paste-name" data-testid="paste-name" placeholder="Pasted sequence"
                        value={name} onChange={(e: any) => setName(e.target.value)} />
          </FormGroup>
          <RadioGroup inline label="Topology" selectedValue={circular ? "circular" : "linear"}
                      onChange={(e: any) => setCircular(e.currentTarget.value === "circular")}>
            <Radio label="Linear" value="linear" />
            <Radio label="Circular" value="circular" data-testid="paste-circular" />
          </RadioGroup>
          <span className="vt-paste-hint">GenBank and EMBL records keep their own topology.</span>
        </div>
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button onClick={props.onClose}>Cancel</Button>
          <Button intent="primary" data-testid="paste-open" disabled={check.state !== "ok"}
                  loading={check.state === "checking"} onClick={open}>
            Open
          </Button>
        </div>
      </div>
    </Dialog>
  );
}
