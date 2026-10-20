import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import { quitDetail, quitHeadline, type QuitDoc } from "./desktop";
import "./PasteDialog.css";
import "./QuitDialog.css";

// Blueprint 3 typings predate React 18 (no implicit `children`); container components are untyped here.
const { Dialog } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

/**
 * Asked when the window is closed while documents have unsaved changes. Save All is the
 * default (Enter); Escape and the dialog's close button keep Vettore open.
 */
export default function QuitDialog(props: {
  docs: QuitDoc[];
  /** The document being saved by Save All, with its position in the list (from 1). */
  saving: { name: string; index: number } | null;
  onSaveAll: () => void;
  onDiscard: () => void;
  onCancel: () => void;
}) {
  const { docs, saving } = props;
  const one = docs.length === 1;
  return (
    <Dialog
      isOpen
      onClose={() => !saving && props.onCancel()}
      canOutsideClickClose={false}
      isCloseButtonShown={!saving}
      title="Quit Vettore"
      style={{ width: 540 }}
    >
      <div className={Classes.DIALOG_BODY} data-testid="quit-dialog">
        <p className="vt-paste-intro"><b>{quitHeadline(docs)}</b></p>
        {!one && (
          <ul className="vt-quit-list" data-testid="quit-docs">
            {docs.map((d) => (
              <li key={d.editorName}>
                <b>{d.name}</b> <span className="vt-muted">— {quitDetail(d)}</span>
              </li>
            ))}
          </ul>
        )}
        <p className="vt-paste-intro">
          <b>{one ? "Save" : "Save All"}</b> saves {one ? "it" : "each one"} the way File → Save
          does{one ? ` (${quitDetail(docs[0])})` : ""}, then quits. <b>Don't Save</b> quits and
          discards the changes.
        </p>
        {saving && (
          <p className="vt-busy" role="status" data-testid="quit-saving">
            Saving “{saving.name}”{one ? "" : ` (${saving.index} of ${docs.length})`}…
          </p>
        )}
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button data-testid="quit-discard" disabled={!!saving} onClick={props.onDiscard}>Don't Save</Button>
          <Button data-testid="quit-cancel" disabled={!!saving} onClick={props.onCancel}>Cancel</Button>
          <Button intent="primary" autoFocus data-testid="quit-save-all" loading={!!saving} onClick={props.onSaveAll}>
            {one ? "Save" : "Save All"}
          </Button>
        </div>
      </div>
    </Dialog>
  );
}
