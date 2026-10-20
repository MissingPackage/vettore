import { useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes, InputGroup, Radio, TextArea } from "@blueprintjs/core";

// Blueprint 3 typings predate React 18 (no implicit `children`); container components are untyped here.
const { Dialog, FormGroup, RadioGroup } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

export interface NewSequence {
  name: string;
  sequence: string;
  circular: boolean;
}

const IUPAC = /[^ACGTURYSWKMBDHVN]/gi;

/** "New DNA File": paste or type a sequence; whitespace, digits and FASTA headers are stripped,
 * and RNA letters are read as DNA (U as T), as when a file is opened. */
export function cleanSequence(raw: string): { sequence: string; rejected: number; rnaLetters: number } {
  const body = raw
    .split(/\r?\n/)
    .filter((line) => !line.startsWith(">"))
    .join("")
    .replace(/[\s\d]/g, "");
  const sequence = body.replace(IUPAC, "").toUpperCase();
  const rnaLetters = (sequence.match(/U/g) ?? []).length;
  return { sequence: sequence.replace(/U/g, "T"), rejected: body.length - sequence.length, rnaLetters };
}

export default function NewSequenceDialog(props: {
  isOpen: boolean;
  onClose: () => void;
  onCreate: (seq: NewSequence) => void;
}) {
  const [name, setName] = useState("New DNA");
  const [raw, setRaw] = useState("");
  const [circular, setCircular] = useState(false);
  const { sequence, rejected, rnaLetters } = cleanSequence(raw);

  return (
    <Dialog isOpen={props.isOpen} onClose={props.onClose} title="New DNA File" style={{ width: 620 }}>
      <div className={Classes.DIALOG_BODY}>
        <p className="vt-dialog-intro" data-testid="new-seq-intro">
          Starts a new file from a sequence you type or paste. FASTA headers, spaces and numbers
          are left out. The file stays in Vettore until you save it (Ctrl+S).
        </p>
        <FormGroup label="Name" labelFor="vt-new-name">
          <InputGroup id="vt-new-name" value={name} onChange={(e: any) => setName(e.target.value)} />
        </FormGroup>
        <FormGroup label="Sequence (paste DNA or FASTA)" helperText={`${sequence.length} bp${rejected ? ` · ${rejected} invalid characters ignored` : ""}${rnaLetters ? " · RNA letters: U read as T" : ""}`}>
          <TextArea
            data-testid="new-seq-text"
            fill
            growVertically={false}
            style={{ height: 180, fontFamily: "monospace" }}
            value={raw}
            onChange={(e: any) => setRaw(e.target.value)}
          />
        </FormGroup>
        <RadioGroup inline label="Topology" selectedValue={circular ? "circular" : "linear"} onChange={(e: any) => setCircular(e.currentTarget.value === "circular")}>
          <Radio label="Linear" value="linear" inputRef={undefined} data-testid="topology-linear" />
          <Radio label="Circular" value="circular" data-testid="topology-circular" />
        </RadioGroup>
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button onClick={props.onClose}>Cancel</Button>
          <Button
            intent="primary"
            data-testid="new-seq-create"
            disabled={!name.trim()}
            onClick={() => props.onCreate({ name: name.trim(), sequence, circular })}
          >
            Create
          </Button>
        </div>
      </div>
    </Dialog>
  );
}
