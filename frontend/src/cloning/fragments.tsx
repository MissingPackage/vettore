// One fragment of a ligation: a piece of a restriction digest of an open file, or a region of
// a file used as a blunt PCR product. Shared by "Insert Multiple Fragments" and "Ligate
// Fragments"; reports the engine's FragmentIn body (cloning/ligation.py).
import { useEffect, useState } from "react";
import { api } from "../api";
import { docData, selectionOf } from "./defaults";
import { newRow, type FragmentMode, type FragmentRow } from "./fragmentRows";
import { DocSelect, EnzymePicker, Field, RegionInput, findDoc } from "./inputs";
import type { OpenDoc } from "./types";
import "./cloning.css";

interface DigestFragment {
  index: number;
  start: number;
  end: number;
  length: number;
  label: string;
  left_end: string;
  right_end: string;
}

export function FragmentEditor(props: {
  docs: OpenDoc[];
  row: FragmentRow;
  title: string;
  role: "vector" | "insert" | "fragment";
  /** How a newly chosen linear file is used. */
  prefer: FragmentMode;
  testId: string;
  onChange: (row: FragmentRow) => void;
  onRemove?: () => void;
}) {
  const { docs, row, testId } = props;
  const d = findDoc(docs, row.editorName);
  const data = d?.data() ?? null;
  const [digest, setDigest] = useState<DigestFragment[] | null>(null);
  const [digestError, setDigestError] = useState<string | null>(null);
  const update = (patch: Partial<FragmentRow>) => props.onChange({ ...row, ...patch });
  const sequence = data?.sequence ?? "";

  // The fragments the chosen enzymes give, so the user can see (and pick) the piece used.
  useEffect(() => {
    setDigest(null);
    setDigestError(null);
    if (!d || row.mode !== "digest" || !row.enzymes.length) return;
    let live = true;
    api
      .cloning<{ fragments: DigestFragment[] }>("digest", { doc: docData(d), enzymes: row.enzymes })
      .then((r) => live && setDigest(r.fragments))
      .catch((e) => live && setDigestError(String((e as Error).message)));
    return () => {
      live = false;
    };
  }, [row.editorName, row.mode, row.enzymes.join(","), sequence]); // eslint-disable-line react-hooks/exhaustive-deps

  const automatic = props.role === "vector"
    ? "Automatic: the largest piece (the vector backbone)"
    : "Automatic: the piece between the two sites";

  return (
    <div className="vt-fragment" data-testid={testId}>
      <div className="vt-fragment-head">
        <b>{props.title}</b>
        {props.onRemove && (
          <button className="bp3-button bp3-small bp3-minimal" data-testid={`${testId}-remove`} onClick={props.onRemove}>
            Remove
          </button>
        )}
      </div>
      <Field label="File">
        <DocSelect
          docs={docs}
          value={row.editorName}
          testId={`${testId}-doc`}
          onChange={(editorName) => {
            const fresh = newRow(docs, editorName, props.prefer);
            props.onChange({ ...fresh, enzymes: row.enzymes });
          }}
        />
      </Field>
      <Field label="Use">
        <div className="bp3-html-select bp3-fill">
          <select
            data-testid={`${testId}-mode`}
            value={row.mode}
            onChange={(e) => {
              const mode = e.target.value as FragmentMode;
              update({ mode, fragmentIndex: null, ...(mode === "pcr" ? selectionOf(d) : {}) });
            }}
          >
            <option value="digest">Cut out with restriction enzymes</option>
            <option value="pcr">As a PCR product (blunt ends)</option>
          </select>
        </div>
      </Field>
      {row.mode === "digest" ? (
        <>
          <Field label="Cut with" hint="One or two enzymes">
            <EnzymePicker
              value={row.enzymes}
              onChange={(enzymes) => update({ enzymes, fragmentIndex: null })}
              max={2}
              sequence={data?.sequence}
              circular={data?.circular}
              testId={`${testId}-enzymes`}
            />
          </Field>
          {digestError && <span className="vt-inline-error" data-testid={`${testId}-error`}>{digestError}</span>}
          {digest && (
            <Field label="Piece to use">
              <div className="bp3-html-select bp3-fill">
                <select
                  data-testid={`${testId}-fragment`}
                  value={row.fragmentIndex ?? ""}
                  onChange={(e) => update({ fragmentIndex: e.target.value === "" ? null : Number(e.target.value) })}
                >
                  <option value="">{automatic}</option>
                  {digest.map((f) => (
                    <option key={f.index} value={f.index}>
                      {f.label}, {f.start + 1}..{f.end + 1}
                    </option>
                  ))}
                </select>
              </div>
            </Field>
          )}
        </>
      ) : (
        <Field
          label="Region"
          hint={data?.circular ? "Required for a circular file" : "Leave empty to use the whole file"}
        >
          <RegionInput doc={d} start={row.start} end={row.end} onChange={(start, end) => update({ start, end })} testId={`${testId}-region`} features />
        </Field>
      )}
    </div>
  );
}
