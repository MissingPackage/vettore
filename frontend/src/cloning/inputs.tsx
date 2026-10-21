// Shared form controls for the cloning dialogs. Coordinates are shown 1-based, as biologists count, and
// converted to the engine's 0-based inclusive coordinates.
import { useMemo, useState } from "react";
import type { OveAnnotation } from "../api";
import { cutsites, ENZYMES } from "../enzymeLogic";
import type { OpenDoc } from "./types";

export function Field(props: { label: string; children: React.ReactNode; hint?: string }) {
  return (
    <label className="vt-field">
      <span className="vt-field-label">{props.label}</span>
      {props.children}
      {props.hint && <span className="vt-field-hint">{props.hint}</span>}
    </label>
  );
}

export function DocSelect(props: {
  docs: OpenDoc[];
  value: string | null;
  onChange: (editorName: string) => void;
  testId?: string;
  /** Extra words after a file's size, e.g. the att sites it carries. */
  note?: (d: OpenDoc) => string | null;
}) {
  return (
    <div className="bp3-html-select bp3-fill">
      <select
        data-testid={props.testId}
        value={props.value ?? ""}
        onChange={(e) => props.onChange(e.target.value)}
      >
        <option value="" disabled>Choose an open file…</option>
        {props.docs.map((d) => {
          const data = d.data();
          const note = props.note?.(d);
          return (
            <option key={d.editorName} value={d.editorName}>
              {d.name} ({data.sequence.length} bp, {data.circular ? "circular" : "linear"}){note ? ` — ${note}` : ""}
            </option>
          );
        })}
      </select>
    </div>
  );
}

const ENZYME_NAMES = Object.values(ENZYMES)
  .map((e) => e.name)
  .sort((a, b) => a.localeCompare(b));

/** Up to `max` enzymes, each shown with its number of sites in the given sequence. */
export function EnzymePicker(props: {
  value: string[];
  onChange: (names: string[]) => void;
  max: number;
  sequence?: string;
  circular?: boolean;
  only?: (name: string) => boolean;
  testId?: string;
}) {
  const [text, setText] = useState("");
  const counts = useMemo(() => {
    if (!props.sequence) return {} as Record<string, number>;
    const found = cutsites(props.sequence, !!props.circular);
    return Object.fromEntries(Object.entries(found).map(([n, s]) => [n.toLowerCase(), s.length]));
  }, [props.sequence, props.circular]);
  const names = props.only ? ENZYME_NAMES.filter(props.only) : ENZYME_NAMES;
  const listId = `enz-${props.testId ?? "picker"}`;

  const add = (raw: string) => {
    const match = names.find((n) => n.toLowerCase() === raw.trim().toLowerCase());
    if (!match || props.value.includes(match) || props.value.length >= props.max) return;
    props.onChange([...props.value, match]);
    setText("");
  };

  return (
    <div className="vt-enzyme-picker">
      {props.value.map((n) => {
        const c = counts[n.toLowerCase()] ?? 0;
        // Without a sequence (a list of enzymes, as in Preferences) there are no sites to count.
        const intent = !props.sequence ? "" : c === 1 ? "bp3-intent-success" : "bp3-intent-warning";
        return (
          <span key={n} className={`bp3-tag bp3-minimal ${intent}`}>
            {n}{props.sequence ? ` · ${c} site${c === 1 ? "" : "s"}` : ""}
            <button className="bp3-tag-remove" onClick={() => props.onChange(props.value.filter((x) => x !== n))} />
          </span>
        );
      })}
      {props.value.length < props.max && (
        <>
          <input
            className="bp3-input bp3-small"
            list={listId}
            placeholder="Add enzyme"
            value={text}
            data-testid={props.testId}
            onChange={(e) => {
              setText(e.target.value);
              if (names.some((n) => n.toLowerCase() === e.target.value.trim().toLowerCase())) add(e.target.value);
            }}
            onKeyDown={(e) => e.key === "Enter" && add(text)}
          />
          <datalist id={listId}>
            {names.map((n) => (
              <option key={n} value={n}>
                {counts[n.toLowerCase()] ? `${counts[n.toLowerCase()]} site(s)` : ""}
              </option>
            ))}
          </datalist>
        </>
      )}
    </div>
  );
}

/** The document's features with a position, in map order, for picking a region by name. */
export function featuresOf(doc: OpenDoc | null): OveAnnotation[] {
  const raw = doc?.data()?.features ?? [];
  return (Array.isArray(raw) ? raw : Object.values(raw))
    .filter((f) => f.start >= 0 && f.end >= 0)
    .sort((a, b) => a.start - b.start || a.end - b.end);
}

const position = (value: number) => value.toLocaleString("en-US");

/**
 * A region as 1-based inclusive start..end; reports 0-based inclusive values. With `features`,
 * the file's features are offered by name above the positions, as in Design PCR Primers.
 */
export function RegionInput(props: {
  doc: OpenDoc | null;
  start: number | null;
  end: number | null;
  onChange: (start: number | null, end: number | null) => void;
  testId?: string;
  features?: boolean;
}) {
  const len = props.doc?.data().sequence.length ?? 0;
  const num = (v: string) => (v.trim() === "" ? null : Math.max(1, Math.min(len, Number(v))) - 1);
  const features = props.features ? featuresOf(props.doc) : [];
  // The feature stays shown while the positions are its own.
  const picked = features.findIndex((f) => f.start === props.start && f.end === props.end);
  return (
    <div className="vt-region-input">
      {features.length > 0 && (
        <div className="bp3-html-select bp3-fill">
          <select
            data-testid={props.testId && `${props.testId}-feature`}
            value={picked >= 0 ? String(picked) : ""}
            onChange={(e) => {
              const f = features[Number(e.target.value)];
              if (f) props.onChange(f.start, f.end);
            }}
          >
            <option value="">Pick a feature, or type the positions below…</option>
            {features.map((f, i) => (
              <option key={i} value={i}>
                {f.name} ({position(f.start + 1)}–{position(f.end + 1)})
              </option>
            ))}
          </select>
        </div>
      )}
      <div className="vt-region">
        <input
          className="bp3-input bp3-small"
          type="number"
          min={1}
          max={len}
          data-testid={props.testId && `${props.testId}-start`}
          value={props.start === null ? "" : props.start + 1}
          onChange={(e) => props.onChange(num(e.target.value), props.end)}
        />
        <span>…</span>
        <input
          className="bp3-input bp3-small"
          type="number"
          min={1}
          max={len}
          data-testid={props.testId && `${props.testId}-end`}
          value={props.end === null ? "" : props.end + 1}
          onChange={(e) => props.onChange(props.start, num(e.target.value))}
        />
        <button
          className="bp3-button bp3-small"
          disabled={!props.doc?.selection()}
          title="The region selected in this file's map"
          onClick={() => {
            const s = props.doc?.selection();
            if (s) props.onChange(s.start, s.end);
          }}
        >
          Use selection
        </button>
        {props.start !== null && props.end !== null && (
          <span className="vt-muted">
            {position(props.start <= props.end ? props.end - props.start + 1 : len - props.start + props.end + 1)} bp
          </span>
        )}
      </div>
    </div>
  );
}

export function findDoc(docs: OpenDoc[], editorName: string | null) {
  return docs.find((d) => d.editorName === editorName) ?? null;
}
