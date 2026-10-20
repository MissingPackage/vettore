// Engine error responses turned into sentences for people. The engine writes its own messages
// for everything a user can cause; anything else (a malformed request from the interface, a
// crash) still has to say what happened and what to do, never "[object Object]".

const LOG_HINT = "Help → About Vettore shows where Vettore's log file is: please send it if this "
  + "happens again.";

/** "Something went wrong" for a failure Vettore cannot explain further. */
export function internalErrorText(status?: number): string {
  return `Something went wrong inside Vettore${status ? ` (error ${status})` : ""}. Your files on disk `
    + `are not affected; try again. ${LOG_HINT}`;
}

/** Where a validation error points, e.g. "sequenceData.features[0].start" (the body prefix dropped). */
function fieldOf(loc: unknown): string {
  if (!Array.isArray(loc)) return "";
  const parts = loc.filter((p, i) => !(i === 0 && (p === "body" || p === "query" || p === "path")));
  return parts.reduce<string>((out, part) => (typeof part === "number"
    ? `${out}[${part}]`
    : out ? `${out}.${String(part)}` : String(part)), "");
}

function sentence(text: string): string {
  const trimmed = text.trim();
  return /[.!?…:]$/.test(trimmed) ? trimmed : `${trimmed}.`;
}

/** A plain-text message from any value an error carries, or "" when it holds none. */
function textOf(value: unknown): string {
  if (typeof value === "string") return value.trim();
  if (Array.isArray(value)) {
    // FastAPI's validation errors: [{loc: ["body", "path"], msg: "Field required"}, ...].
    const parts = value.map((item) => {
      if (typeof item === "string") return item.trim();
      if (item && typeof item === "object") {
        const { loc, msg } = item as { loc?: unknown; msg?: unknown };
        const field = fieldOf(loc);
        const text = typeof msg === "string" ? msg.trim() : "";
        if (text) return field ? `${field}: ${text.charAt(0).toLowerCase()}${text.slice(1)}` : text;
      }
      return "";
    }).filter(Boolean);
    if (!parts.length) return "";
    const isValidation = value.some((item) => item && typeof item === "object" && "loc" in item);
    return isValidation
      ? `Vettore's engine could not use the request it received (${parts.join("; ")}). `
        + `This is a fault in Vettore, not in your data. ${LOG_HINT}`
      : parts.map(sentence).join(" ");
  }
  if (value && typeof value === "object") {
    for (const key of ["message", "detail", "msg", "error"]) {
      const text = textOf((value as Record<string, unknown>)[key]);
      if (text) return text;
    }
  }
  return "";
}

/** The message for an engine response with HTTP status `status` and JSON body `payload`. */
export function engineErrorText(status: number, payload: unknown): string {
  const detail = payload && typeof payload === "object" ? (payload as { detail?: unknown }).detail : undefined;
  return textOf(detail) || internalErrorText(status);
}

/** The message of anything thrown: an Error's message, a string, or an error-like object. */
export function errorText(e: unknown): string {
  if (e instanceof Error) return e.message.trim() || internalErrorText();
  return textOf(e) || internalErrorText();
}
