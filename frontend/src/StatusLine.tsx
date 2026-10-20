import { useEffect, useState } from "react";
import { useSelector } from "react-redux";
import { api } from "./api";
import { messageOf } from "./files/notices";

const TM_MAX_LENGTH = 100;

/** Selection summary shown next to the view tabs: range, length, %GC and Tm for oligo-sized selections. */
export default function StatusLine({ editorName }: { editorName: string }) {
  const { sequence, circular, start, end } = useSelector((state: any) => {
    const ed = state.VectorEditor[editorName] ?? {};
    return {
      sequence: (ed.sequenceData?.sequence ?? "") as string,
      circular: !!ed.sequenceData?.circular,
      start: ed.selectionLayer?.start ?? -1,
      end: ed.selectionLayer?.end ?? -1,
    };
  }, (a, b) => a.sequence === b.sequence && a.start === b.start && a.end === b.end && a.circular === b.circular);

  const hasSelection = start >= 0 && end >= 0;
  const selected = !hasSelection
    ? ""
    : start <= end
      ? sequence.slice(start, end + 1)
      : sequence.slice(start) + sequence.slice(0, end + 1);
  const [tm, setTm] = useState<{ value: number; tooltip: string } | null>(null);
  // Why the engine gave no Tm (e.g. the selection holds N or other ambiguous bases).
  const [tmError, setTmError] = useState<string | null>(null);

  useEffect(() => {
    setTm(null);
    setTmError(null);
    if (!selected || selected.length > TM_MAX_LENGTH || selected.length < 8) return;
    const t = setTimeout(() => {
      api
        .tm(selected)
        .then((r) =>
          setTm({
            value: r.tm,
            tooltip: r.presets.map((p) => `${p.label}: ${p.tm.toFixed(1)} °C — ${p.description}`).join("\n"),
          }),
        )
        .catch((e) => {
          setTm(null);
          setTmError(messageOf(e));
        });
    }, 150);
    return () => clearTimeout(t);
  }, [selected]);

  const gc = selected
    ? (100 * (selected.match(/[GCS]/gi)?.length ?? 0)) / selected.length
    : sequence.length
      ? (100 * (sequence.match(/[GCS]/gi)?.length ?? 0)) / sequence.length
      : 0;

  return (
    <div className="vt-status" data-testid="status-line">
      {hasSelection ? (
        <>
          <span title={start > end ? "The selection runs across the origin" : "Selected range"}>
            Selected {start + 1} … {end + 1}
            {start > end ? " (across the origin)" : ""}
          </span>
          <span>{selected.length.toLocaleString()} bp</span>
          <span title="G and C as a percentage of the selection">{gc.toFixed(1)}% GC</span>
          {tm !== null && (
            <span className="vt-tm" title={tm.tooltip} data-testid="selection-tm">
              Tm {tm.value.toFixed(1)} °C
            </span>
          )}
          {tm === null && tmError && (
            <span className="vt-muted vt-tm" title={tmError} data-testid="selection-tm-error">
              Tm: not available
            </span>
          )}
          {tm === null && selected.length > TM_MAX_LENGTH && (
            <span className="vt-muted" title={`Tm is shown for selections up to ${TM_MAX_LENGTH} bases`}>
              Tm: select ≤ {TM_MAX_LENGTH} bases
            </span>
          )}
        </>
      ) : (
        <>
          <span>{sequence.length.toLocaleString()} bp · {circular ? "circular" : "linear"}</span>
          <span title="G and C as a percentage of the whole sequence">{gc.toFixed(1)}% GC</span>
          <span className="vt-muted">Select a region for its length, %GC and Tm</span>
        </>
      )}
    </div>
  );
}
