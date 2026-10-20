import { describe, expect, it } from "vitest";
import { engineErrorText, errorText, internalErrorText } from "./engineErrors";
import { messageOf } from "./notices";

describe("engine error messages", () => {
  it("passes the engine's own sentences through", () => {
    expect(engineErrorText(422, { detail: "“x.gb” is empty." })).toBe("“x.gb” is empty.");
  });

  it("turns a request the engine could not validate into a sentence (seen with open-path {path: null})", () => {
    const text = engineErrorText(422, {
      detail: [{ type: "string_type", loc: ["body", "path"], msg: "Input should be a valid string", input: null }],
    });
    expect(text).toContain("path: input should be a valid string");
    expect(text).toContain("fault in Vettore");
    expect(text).toContain("Help → About Vettore");
    expect(text).not.toContain("[object");
  });

  it("names nested fields of a validation error", () => {
    const text = engineErrorText(422, {
      detail: [
        { loc: ["body", "sequenceData", "features", 0, "start"], msg: "Field required" },
        { loc: ["body", "format"], msg: "Input should be 'genbank'" },
      ],
    });
    expect(text).toContain("sequenceData.features[0].start: field required; format: input should be 'genbank'");
  });

  it("reads a message out of an object or a list of strings", () => {
    expect(engineErrorText(400, { detail: { message: "The ladder file could not be read" } }))
      .toBe("The ladder file could not be read");
    expect(engineErrorText(400, { detail: ["First problem", "Second problem."] }))
      .toBe("First problem. Second problem.");
  });

  it("falls back to a plain sentence for a crash or an unreadable body", () => {
    // A 500 from the server arrives as the text "Internal Server Error", not JSON.
    expect(engineErrorText(500, {})).toBe(internalErrorText(500));
    expect(engineErrorText(500, { detail: { unexpected: 1 } })).toBe(internalErrorText(500));
    expect(engineErrorText(500, null)).toContain("error 500");
    expect(internalErrorText(500)).toContain("Help → About Vettore shows where Vettore's log file is");
  });

  it("never shows [object Object] for anything thrown", () => {
    for (const thrown of [{}, { foo: 1 }, [], null, undefined, 42, new Error("")]) {
      expect(messageOf(thrown)).not.toContain("[object");
      expect(messageOf(thrown).length).toBeGreaterThan(20);
    }
    expect(errorText({ message: "Read failed." })).toBe("Read failed.");
    expect(messageOf("plain")).toBe("plain");
  });
});
