// File → Preferences: the defaults new documents and dialogs start from.
// Every change is saved and applied at once (no OK button to forget).
import { useEffect, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import { api, type FileAssociations } from "../api";
import { associationHint, movedNote, overriddenNote } from "../files/associations";
import { messageOf } from "../files/notices";
import { ENZYME_SETS } from "../enzymeLogic";
import {
  CHROME_THEMES,
  defaultPreferences,
  IDENTITY_CHOICES,
  parseEnzymeList,
  savePreferences,
  usePreferences,
  type Preferences,
} from "../preferences";
import { EnzymePicker, Field } from "../cloning/inputs";
import type { Tool, ToolProps } from "./types";
import "./preferences.css";

// Blueprint 3 typings predate React 18 (no implicit `children`).
const { Dialog, Callout } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

const PERCENT_FALLBACK = [0.5, 0.7, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0];

/**
 * Open sequence files with Vettore on double-click. Shown only by the packaged application
 * (a copy run from source code has no program to register); applied at once, like the rest.
 */
function OpenWithVettore() {
  const [status, setStatus] = useState<FileAssociations | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [done, setDone] = useState<string | null>(null);

  useEffect(() => {
    api.fileAssociations().then(setStatus).catch(() => setStatus(null));
  }, []);
  if (!status?.available) return null;

  const change = async (enabled: boolean) => {
    setBusy(true);
    setError(null);
    setDone(null);
    try {
      setStatus(await api.setFileAssociations(enabled));
      setDone(enabled
        ? "Sequence files now open with Vettore when you double-click them."
        : "Sequence files open as they did before.");
    } catch (e) {
      setError(messageOf(e));
    } finally {
      setBusy(false);
    }
  };
  const moved = movedNote(status);
  const overridden = overriddenNote(status);

  return (
    <div className="vt-field" data-testid="pref-open-with-field">
      <span className="vt-field-label">Opening files</span>
      <label className="bp3-control bp3-checkbox vt-pref-check">
        <input
          type="checkbox"
          data-testid="pref-open-with"
          checked={status.enabled}
          disabled={busy}
          aria-describedby="pref-open-with-hint"
          onChange={(e) => change(e.target.checked)}
        />
        <span className="bp3-control-indicator" />
        Open sequence files with Vettore when you double-click them
      </label>
      <span className="vt-field-hint" id="pref-open-with-hint">{associationHint(status)}</span>
      {busy && <span className="vt-busy" role="status" data-testid="pref-open-with-busy">Applying…</span>}
      {moved && <span className="vt-field-hint" data-testid="pref-open-with-moved">{moved}</span>}
      {overridden && (
        <div className="vt-notice vt-notice-warning" role="status" data-testid="pref-open-with-overridden">
          {overridden}
        </div>
      )}
      {error && <span className="vt-field-error" role="alert" data-testid="pref-open-with-error">{error}</span>}
      {done && !busy && <span className="vt-pref-done" role="status" data-testid="pref-open-with-done">{done}</span>}
    </div>
  );
}

/** The enzymes in the lab's freezer, picked one by one or pasted as a list. */
function MyEnzymes(props: { value: string[]; onChange: (names: string[], what: string) => void }) {
  const [text, setText] = useState("");
  const [unknown, setUnknown] = useState<string[]>([]);
  const addList = () => {
    const { known, unknown: notFound } = parseEnzymeList(text);
    const added = known.filter((n) => !props.value.includes(n));
    setUnknown(notFound);
    if (added.length) props.onChange([...props.value, ...added], `${added.length} enzyme${added.length === 1 ? "" : "s"} added`);
    if (!notFound.length) setText("");
  };
  return (
    <div className={`vt-field${unknown.length ? " vt-field-invalid" : ""}`} data-testid="pref-my-enzymes-field">
      <span className="vt-field-label">My enzymes</span>
      <EnzymePicker
        value={props.value}
        onChange={(names) => props.onChange(names, "My enzymes saved")}
        max={400}
        testId="pref-my-enzymes"
      />
      <div className="vt-row">
        <input
          className="bp3-input bp3-fill"
          data-testid="pref-my-enzymes-list"
          aria-label="Paste a list of enzymes"
          aria-describedby="pref-my-enzymes-hint"
          placeholder="Or paste a list: EcoRI, BamHI, HindIII…"
          value={text}
          onChange={(e) => setText(e.target.value)}
          onKeyDown={(e) => e.key === "Enter" && addList()}
        />
        <Button data-testid="pref-my-enzymes-add" disabled={!text.trim()} onClick={addList}>Add</Button>
        {props.value.length > 0 && (
          <Button minimal data-testid="pref-my-enzymes-clear" onClick={() => props.onChange([], "My enzymes cleared")}>
            Clear
          </Button>
        )}
      </div>
      <span className="vt-field-hint" id="pref-my-enzymes-hint">
        The enzymes in your freezer. Find Diagnostic Digest tries these first; while the list is
        empty it tries every commercial enzyme.
      </span>
      {unknown.length > 0 && (
        <span className="vt-field-error" role="alert" data-testid="pref-my-enzymes-unknown">
          Vettore does not know {unknown.map((u) => `“${u}”`).join(", ")}. Check the spelling
          (for example HindIII, not Hind3), then add {unknown.length === 1 ? "it" : "them"} again.
        </span>
      )}
    </div>
  );
}

function PreferencesDialog(props: ToolProps) {
  const prefs = usePreferences();
  const [ladders, setLadders] = useState<{ id: string; label: string }[]>([]);
  const [percentages, setPercentages] = useState<number[]>(PERCENT_FALLBACK);
  const [loading, setLoading] = useState(true);
  const [ladderError, setLadderError] = useState<string | null>(null);
  const [saved, setSaved] = useState<string | null>(null);

  useEffect(() => {
    api.gelLadders()
      .then((r) => {
        setLadders(r.ladders);
        setPercentages(r.percentages);
      })
      .catch((e) => setLadderError(messageOf(e)))
      .finally(() => setLoading(false));
  }, []);

  const change = (patch: Partial<Preferences>, what: string) => {
    savePreferences(patch);
    setSaved(what);
  };

  return (
    <Dialog isOpen onClose={props.onClose} title="Preferences" style={{ width: 560 }}>
      <div className={`${Classes.DIALOG_BODY} vt-preferences`} data-testid="preferences-dialog">
        <p className="vt-muted">
          The defaults Vettore starts from. Changes apply straight away and are kept for the
          next time you open the app; documents already open keep what you chose for them.
        </p>

        <Field
          label="Appearance"
          hint="Dark keeps the map, the sequence and dialogs light: only menus, tabs and the status bar change."
        >
          <div className="bp3-html-select bp3-fill">
            <select
              data-testid="pref-chrome"
              value={prefs.chrome}
              onChange={(e) => change({ chrome: e.target.value as Preferences["chrome"] }, "Appearance saved")}
            >
              {CHROME_THEMES.map((t) => <option key={t.id} value={t.id}>{t.label}</option>)}
            </select>
          </div>
        </Field>

        <Field
          label="Enzyme set for documents you open"
          hint="Shown on the map and in the Enzymes view. You can still change it per document."
        >
          <div className="bp3-html-select bp3-fill">
            <select
              data-testid="pref-enzyme-set"
              value={prefs.enzymeSet}
              onChange={(e) => change({ enzymeSet: e.target.value as Preferences["enzymeSet"] },
                                      "Enzyme set saved")}
            >
              {ENZYME_SETS.map((set) => (
                <option key={set.id} value={set.id}>{set.label}</option>
              ))}
            </select>
          </div>
        </Field>

        <Field label="Gel ladder" hint="Pre-selected in Tools → Simulate Agarose Gel">
          <div className="bp3-html-select bp3-fill">
            <select
              data-testid="pref-gel-ladder"
              value={prefs.gelLadder}
              disabled={loading}
              onChange={(e) => change({ gelLadder: e.target.value }, "Ladder saved")}
            >
              {loading && <option>Loading ladders…</option>}
              {ladders.map((l) => <option key={l.id} value={l.id}>{l.label}</option>)}
            </select>
          </div>
          {ladderError && (
            <span className="vt-field-error" data-testid="pref-ladder-error">
              The list of ladders could not be read: {ladderError} Close and reopen Preferences
              to try again.
            </span>
          )}
        </Field>

        <Field label="Agarose percentage" hint="The gel you usually pour">
          <div className="bp3-html-select bp3-fill">
            <select
              data-testid="pref-gel-percent"
              value={prefs.gelPercent}
              onChange={(e) => change({ gelPercent: Number(e.target.value) }, "Agarose % saved")}
            >
              {percentages.map((p) => <option key={p} value={p}>{p.toFixed(1)} %</option>)}
            </select>
          </div>
        </Field>

        <Field
          label="Minimum identity when detecting common features"
          hint="Lower finds more, including variants; 100 % accepts exact matches only."
        >
          <div className="bp3-html-select bp3-fill">
            <select
              data-testid="pref-identity"
              value={prefs.detectMinIdentity}
              onChange={(e) => change({ detectMinIdentity: Number(e.target.value) },
                                      "Identity threshold saved")}
            >
              {IDENTITY_CHOICES.map((p) => <option key={p} value={p}>{p} %</option>)}
            </select>
          </div>
        </Field>

        <MyEnzymes value={prefs.myEnzymes} onChange={(myEnzymes, what) => change({ myEnzymes }, what)} />

        <OpenWithVettore />

        {saved && (
          <Callout intent="success" data-testid="pref-saved">{saved}.</Callout>
        )}
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button
            data-testid="pref-reset"
            onClick={() => change(defaultPreferences(), "Defaults restored")}
          >
            Restore defaults
          </Button>
          <Button intent="primary" data-testid="pref-close" onClick={props.onClose}>Close</Button>
        </div>
      </div>
    </Dialog>
  );
}

export const tools: Tool[] = [
  { id: "preferences", menu: "File", label: "Preferences", Component: PreferencesDialog },
];
