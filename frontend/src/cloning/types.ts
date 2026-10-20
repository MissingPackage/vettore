import type { SequenceData } from "../api";

export interface Oligo {
  name: string;
  sequence: string;
  tm?: number | null;
  annealing_length?: number | null;
}

export interface Product {
  doc: SequenceData;
  operation: string;
  parents: string[];
  details: string;
  oligos: Oligo[];
}

export interface CloningResult {
  products: Product[];
  warnings: string[];
  oligos: Oligo[];
}

/** An open document as seen by the cloning dialogs. */
export interface OpenDoc {
  editorName: string;
  name: string;
  data: () => SequenceData;
  /** Current selection in the editor, 0-based inclusive, or null. */
  selection: () => { start: number; end: number } | null;
}

export interface FormProps {
  docs: OpenDoc[];
  active: OpenDoc | null;
  /** Report the request body (null while the form is incomplete). */
  onChange: (body: object | null) => void;
}
