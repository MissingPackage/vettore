import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import type { DocSource } from "../documents";
import { copyReason } from "./saveCopy";
import "./PasteDialog.css";

// Blueprint 3 typings predate React 18 (no implicit `children`); container components are untyped here.
const { Dialog } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

/**
 * Shown once per document, on Save, when the file it came from cannot be overwritten safely:
 * explains that a GenBank copy will be saved and the original stays as it is.
 */
export default function SaveCopyDialog(props: {
  name: string;
  source: DocSource;
  onCancel: () => void;
  onConfirm: () => void;
}) {
  return (
    <Dialog isOpen onClose={props.onCancel} title="Save a GenBank copy" style={{ width: 520 }}>
      <div className={Classes.DIALOG_BODY} data-testid="save-copy-dialog">
        <p className="vt-paste-intro">{copyReason(props.name, props.source)}</p>
        <p className="vt-paste-intro">
          Vettore leaves the original file untouched and saves your work as a new GenBank file
          (.gb), the format every sequence program opens. Next, choose where to save it.
        </p>
        {props.source.label === ".dna" && (
          <p className="vt-paste-reason">
            To give a colleague a .dna file, use File → Save As… and choose the “DNA file (.dna)” type.
          </p>
        )}
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button onClick={props.onCancel}>Cancel</Button>
          <Button intent="primary" autoFocus data-testid="save-copy-confirm" onClick={props.onConfirm}>
            Choose where to save…
          </Button>
        </div>
      </div>
    </Dialog>
  );
}
