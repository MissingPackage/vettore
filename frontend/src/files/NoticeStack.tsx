import { useEffect } from "react";
import { INFO_TIMEOUT_MS, type Notice } from "./notices";
import "./NoticeStack.css";

const LABEL: Record<Notice["kind"], string> = { error: "Error", warning: "Warning", info: "Done" };

function NoticeItem({ notice, onDismiss }: { notice: Notice; onDismiss: (id: number) => void }) {
  const { id, kind } = notice;
  useEffect(() => {
    if (kind !== "info") return;
    const timer = setTimeout(() => onDismiss(id), INFO_TIMEOUT_MS);
    return () => clearTimeout(timer);
  }, [id, kind, onDismiss]);

  return (
    <div
      className={`vt-notice vt-notice-${kind}`}
      role={kind === "info" ? "status" : "alert"}
      data-testid="notice"
      data-kind={kind}
    >
      <span className="vt-notice-label">{LABEL[kind]}</span>
      <span className="vt-notice-text">
        {notice.text}
        {notice.count > 1 && <span className="vt-notice-count"> (×{notice.count})</span>}
      </span>
      <button
        className="vt-notice-close"
        aria-label="Dismiss this message"
        title="Dismiss"
        data-testid="notice-dismiss"
        onClick={() => onDismiss(id)}
      >
        ×
      </button>
    </div>
  );
}

/**
 * Messages stacked at the top right of the window. Each shows its whole text (wrapped, never
 * cut) until dismissed; confirmations disappear on their own after a few seconds.
 */
export default function NoticeStack(props: {
  notices: Notice[];
  onDismiss: (id: number) => void;
  onDismissAll: () => void;
}) {
  if (!props.notices.length) return null;
  return (
    <div className="vt-notices" data-testid="notices">
      {props.notices.length > 1 && (
        <button className="vt-notices-clear" data-testid="notice-dismiss-all" onClick={props.onDismissAll}>
          Dismiss all
        </button>
      )}
      {props.notices.map((n) => (
        <NoticeItem key={n.id} notice={n} onDismiss={props.onDismiss} />
      ))}
    </div>
  );
}
