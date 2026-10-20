import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import "./PasteDialog.css";

// Blueprint 3 typings predate React 18 (no implicit `children`); container components are untyped here.
const { Dialog } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

export type UnsavedChoice = "save" | "discard" | "cancel";

/**
 * Asked when a tab with unsaved work is closed. Save is the default (Enter); Escape and the
 * window's close button keep the tab open.
 */
export default function UnsavedDialog(props: {
  name: string;
  /** The document was never saved to a file (a product, a paste, a new file). */
  neverSaved: boolean;
  onChoose: (choice: UnsavedChoice) => void;
}) {
  const { onChoose } = props;
  return (
    <Dialog
      isOpen
      onClose={() => onChoose("cancel")}
      canOutsideClickClose={false}
      title="Unsaved changes"
      style={{ width: 500 }}
    >
      <div className={Classes.DIALOG_BODY} data-testid="unsaved-dialog">
        <p className="vt-paste-intro">
          <b>“{props.name}” has changes that are not saved.</b>
        </p>
        <p className="vt-paste-intro">
          {props.neverSaved
            ? "It has never been saved to a file: if you close it without saving, it is gone, "
              + "together with its History."
            : "If you close it without saving, your changes are lost."}
        </p>
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button data-testid="unsaved-discard" onClick={() => onChoose("discard")}>Don't save</Button>
          <Button data-testid="unsaved-cancel" onClick={() => onChoose("cancel")}>Cancel</Button>
          <Button intent="primary" autoFocus data-testid="unsaved-save" onClick={() => onChoose("save")}>
            Save
          </Button>
        </div>
      </div>
    </Dialog>
  );
}
