// "Show Library Primers": the lab-library primers that bind the open document, drawn on its map
// as read-only overlay primers. They are tagged with OVERLAY_NOTE so the engine never writes
// them to a saved file (src/vettore/fileio.py strips them), and removeOverlay takes them off
// again when the toggle is turned off.
import type { OveAnnotation } from "../api";
import { bindingText, OVERLAY_NOTE, type BoundPrimer } from "./primerLibraryText";

const asList = (value: OveAnnotation[] | Record<string, OveAnnotation> | undefined): OveAnnotation[] =>
  Array.isArray(value) ? value : Object.values(value ?? {});

/** True when an annotation is a library-primer overlay Vettore drew (not a real document primer). */
export function isOverlay(primer: OveAnnotation): boolean {
  return !!primer.notes && !!(primer.notes as Record<string, string[]>)[OVERLAY_NOTE];
}

/** One map primer per binding of a library primer, labelled with the primer's name and tagged
 * as an overlay so it is never saved. */
export function overlayAnnotations(bound: BoundPrimer[]): OveAnnotation[] {
  const out: OveAnnotation[] = [];
  for (const p of bound) {
    p.bindings.forEach((b, i) => {
      out.push({
        name: p.bindings.length > 1 ? `${p.name} (${i + 1})` : p.name,
        type: "primer_bind",
        start: b.start,
        end: b.end,
        strand: b.strand,
        forward: b.strand === 1,
        // The bases drawn are those that anneal (start..end); a 5′ tail does not sit on the
        // template, so the whole primer is given in the note.
        bases: p.sequence.slice(p.sequence.length - b.footprint),
        notes: {
          [OVERLAY_NOTE]: ["library"],
          note: [`Lab library primer${p.location ? ` · ${p.location}` : ""}: ${bindingText(b)}`
            + (b.tail ? `; whole primer 5′-${p.sequence}-3′` : "")],
        },
      });
    });
  }
  return out;
}

/** The document's own primers, without any library overlay. */
export function removeOverlay(primers: OveAnnotation[] | Record<string, OveAnnotation> | undefined): OveAnnotation[] {
  return asList(primers).filter((p) => !isOverlay(p));
}

/** The primers list with the library overlay refreshed: the document's own primers, then the
 * overlay for every binding. */
export function withOverlay(
  primers: OveAnnotation[] | Record<string, OveAnnotation> | undefined,
  bound: BoundPrimer[],
): OveAnnotation[] {
  return [...removeOverlay(primers), ...overlayAnnotations(bound)];
}
