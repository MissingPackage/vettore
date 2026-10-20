// Messages shown above the editor: errors, warnings and short confirmations.
// Several can be visible at once; repeating the same message counts it instead of stacking it.
import { errorText } from "./engineErrors";

export type NoticeKind = "error" | "warning" | "info";

export interface Notice {
  id: number;
  kind: NoticeKind;
  text: string;
  /** How many times the same message arrived while it was visible. */
  count: number;
}

/** Confirmations disappear on their own; errors and warnings stay until dismissed. */
export const INFO_TIMEOUT_MS = 8000;
const MAX_NOTICES = 6;

let nextId = 0;

export function addNotice(list: Notice[], kind: NoticeKind, text: string): Notice[] {
  const same = list.find((n) => n.kind === kind && n.text === text);
  if (same) {
    return [...list.filter((n) => n !== same), { ...same, id: ++nextId, count: same.count + 1 }];
  }
  const next = [...list, { id: ++nextId, kind, text, count: 1 }];
  // Keep the newest; drop the oldest confirmations first, then the oldest of any kind.
  while (next.length > MAX_NOTICES) {
    const info = next.findIndex((n) => n.kind === "info");
    next.splice(info >= 0 ? info : 0, 1);
  }
  return next;
}

export function dismissNotice(list: Notice[], id: number): Notice[] {
  return list.filter((n) => n.id !== id);
}

/** The message of anything thrown, without a JavaScript error prefix and never "[object Object]". */
export function messageOf(e: unknown): string {
  return errorText(e);
}
