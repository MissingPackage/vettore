import { useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes, InputGroup } from "@blueprintjs/core";
import { api, type OpenResult } from "../api";
import { messageOf } from "./notices";

// Blueprint 3 typings predate React 18 (no implicit `children`); container components are untyped here.
const { Dialog, FormGroup } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

/**
 * File → Open from NCBI…: a public GenBank record by its accession number, opened unsaved.
 * The only place where Vettore uses the internet; the dialog says what is sent.
 */
export default function NcbiDialog(props: {
  onClose: () => void;
  onOpen: (result: OpenResult, accession: string) => void;
}) {
  const [accession, setAccession] = useState("");
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const open = async () => {
    if (!accession.trim() || busy) return;
    setBusy(true);
    setError(null);
    try {
      const result = await api.openFromNcbi(accession);
      props.onOpen(result, result.accession);
    } catch (e) {
      setError(messageOf(e));
      setBusy(false);
    }
  };

  return (
    <Dialog isOpen onClose={props.onClose} title="Open from NCBI" style={{ width: 560 }}>
      <div className={Classes.DIALOG_BODY} data-testid="ncbi-dialog">
        <p className="vt-dialog-intro">
          Downloads a GenBank record from NCBI by its accession number and opens it as a new file,
          not yet saved.
        </p>
        <FormGroup label="Accession number" labelFor="vt-ncbi-accession"
                   helperText="For example MT612431.1, M77789.2 or NC_001422.1">
          <InputGroup
            id="vt-ncbi-accession"
            data-testid="ncbi-accession"
            autoFocus
            value={accession}
            disabled={busy}
            onChange={(e: any) => {
              setAccession(e.target.value);
              setError(null);
            }}
            onKeyDown={(e: any) => e.key === "Enter" && open()}
          />
        </FormGroup>
        <p className="vt-muted" data-testid="ncbi-privacy">
          This is the only thing in Vettore that uses the internet. Only the accession number is
          sent to NCBI: no sequence and no file leaves this computer. No account is needed.
        </p>
        {busy && <p role="status" data-testid="ncbi-busy">Downloading {accession.trim().toUpperCase()} from NCBI…</p>}
        {error && <p className="vt-field-error" role="alert" data-testid="ncbi-error">{error}</p>}
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button onClick={props.onClose}>Cancel</Button>
          <Button intent="primary" data-testid="ncbi-open" disabled={!accession.trim()} loading={busy}
                  onClick={open}>
            Open
          </Button>
        </div>
      </div>
    </Dialog>
  );
}
