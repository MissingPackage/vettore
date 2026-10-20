import { useEffect, useRef } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import "./PasteDialog.css";

// Blueprint 3 typings predate React 18 (no implicit `children`); container components are untyped here.
const { Dialog } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

/** Shown when the system does not let Vettore use the clipboard: the text, selected, to copy by hand. */
export default function CopyTextDialog(props: { title: string; text: string; onClose: () => void }) {
  const area = useRef<HTMLTextAreaElement>(null);
  // Blueprint moves the focus to the dialog itself once it has opened, so the text is selected
  // again from onOpened; before then the textarea may not even be mounted.
  const selectAll = () => {
    area.current?.focus();
    area.current?.select();
  };
  useEffect(selectAll, []);

  return (
    <Dialog isOpen onClose={props.onClose} onOpened={selectAll} title={props.title} style={{ width: 640 }}>
      <div className={Classes.DIALOG_BODY} data-testid="copy-text-dialog">
        <p className="vt-paste-intro">
          Your system did not let Vettore put this on the clipboard. The text below is already
          selected: press Ctrl+C to copy it, then paste it where you need it.
        </p>
        <textarea
          ref={area}
          readOnly
          data-testid="copy-text"
          className={`${Classes.INPUT} ${Classes.FILL} vt-paste-text`}
          style={{ height: 200 }}
          value={props.text}
          onFocus={(e) => e.currentTarget.select()}
        />
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button intent="primary" onClick={props.onClose}>Close</Button>
        </div>
      </div>
    </Dialog>
  );
}
