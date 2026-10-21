import { useLayoutEffect, useRef, useState } from "react";
import { VIEWS, type ViewId } from "./views";

/**
 * The view tabs (Map, Sequence, Enzymes…) as a segmented control. A white pill marks the
 * open view and slides to the next one when the view changes, so the eye follows the switch.
 * Each document's editor mounts its own copy, so the pill is placed without sliding the first
 * time; only a change of view moves it (and never with reduced motion, see --vt-duration-med).
 */
export default function ViewTabs(props: { view: ViewId | null; onSelect: (view: ViewId) => void }) {
  const list = useRef<HTMLDivElement>(null);
  const [pill, setPill] = useState<{ left: number; width: number } | null>(null);
  const [placed, setPlaced] = useState(false);

  useLayoutEffect(() => {
    const active = list.current?.querySelector<HTMLElement>(".vt-viewtab.active");
    setPill(active ? { left: active.offsetLeft, width: active.offsetWidth } : null);
  }, [props.view]);

  // Turned on after the first placement has painted, so the pill does not slide in from the left.
  useLayoutEffect(() => {
    if (pill && !placed) requestAnimationFrame(() => setPlaced(true));
  }, [pill, placed]);

  return (
    <div className={`vt-viewtabs${placed ? " vt-viewtabs-placed" : ""}`} role="tablist" ref={list}>
      {pill && (
        <span
          className="vt-viewtabs-pill"
          aria-hidden="true"
          style={{ width: pill.width, transform: `translateX(${pill.left}px)` }}
        />
      )}
      {VIEWS.map((v) => (
        <button
          key={v.id}
          role="tab"
          data-testid={`view-${v.id}`}
          aria-selected={props.view === v.id}
          className={`vt-viewtab${props.view === v.id ? " active" : ""}`}
          disabled={!props.view}
          title={`${v.label} (Ctrl+${v.hotkey})`}
          onClick={() => props.onSelect(v.id)}
        >
          {v.label}
        </button>
      ))}
    </div>
  );
}
