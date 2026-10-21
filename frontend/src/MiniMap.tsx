import { getFeatureToColorMap } from "@teselagen/sequence-utils";
import type { SampleFeature } from "./api";

const COLORS: Record<string, string> = getFeatureToColorMap();
const SIZE = 56;
const C = SIZE / 2;
const R = 20;

/** Point on the ring at a fraction of the sequence, clockwise from the top like the editor's map. */
function at(fraction: number, r: number): [number, number] {
  const a = fraction * 2 * Math.PI;
  return [C + r * Math.sin(a), C - r * Math.cos(a)];
}

/**
 * A thumbnail of a plasmid: the backbone ring and its longest features as coloured arcs, in the
 * editor's feature colours. Decorative: the card next to it names the construct and its length.
 */
export default function MiniMap(props: { length: number; features: SampleFeature[] }) {
  const { length } = props;
  return (
    <svg className="vt-minimap" viewBox={`0 0 ${SIZE} ${SIZE}`} width={SIZE} height={SIZE} aria-hidden="true">
      <circle cx={C} cy={C} r={R} fill="none" stroke="var(--vt-border-strong)" strokeWidth={1} />
      {props.features.map((f, i) => {
        const span = ((f.end - f.start + length) % length) + 1;
        if (span >= length) {
          return <circle key={i} cx={C} cy={C} r={R} fill="none" stroke={COLORS[f.type] ?? "#999"} strokeWidth={5} />;
        }
        const [x0, y0] = at(f.start / length, R);
        const [x1, y1] = at((f.start + span) / length, R);
        return (
          <path
            key={i}
            d={`M${x0.toFixed(2)} ${y0.toFixed(2)} A${R} ${R} 0 ${span > length / 2 ? 1 : 0} 1 ${x1.toFixed(2)} ${y1.toFixed(2)}`}
            fill="none"
            stroke={COLORS[f.type] ?? "#999"}
            strokeWidth={5}
          />
        );
      })}
    </svg>
  );
}
