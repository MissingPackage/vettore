import { useEffect, useState } from "react";
import { useSelector } from "react-redux";
import { api } from "./api";

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

  useEffect(() => {
    setTm(null);
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
        .catch(() => setTm(null));
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
          <span>Selected {start + 1} … {end + 1}</span>
          <span>{selected.length} bp</span>
          <span>{gc.toFixed(1)}% GC</span>
          {tm !== null && (
            <span className="vt-tm" title={tm.tooltip} data-testid="selection-tm">
              Tm {tm.value.toFixed(1)} °C
            </span>
          )}
        </>
      ) : (
        <>
          <span>{sequence.length} bp · {circular ? "circular" : "linear"}</span>
          <span>{gc.toFixed(1)}% GC</span>
        </>
      )}
    </div>
  );
}
