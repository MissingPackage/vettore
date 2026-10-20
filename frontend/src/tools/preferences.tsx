// File → Preferences: the defaults new documents and dialogs start from.
// Every change is saved and applied at once (no OK button to forget).
import { useEffect, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import { api } from "../api";
import { ENZYME_SETS } from "../enzymeLogic";
import {
  defaultPreferences,
  IDENTITY_CHOICES,
  savePreferences,
  usePreferences,
  type Preferences,
} from "../preferences";
import { Field } from "../cloning/inputs";
import type { Tool, ToolProps } from "./types";
import "./preferences.css";

// Blueprint 3 typings predate React 18 (no implicit `children`).
const { Dialog, Callout } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

const PERCENT_FALLBACK = [0.5, 0.7, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0];

function PreferencesDialog(props: ToolProps) {
  const prefs = usePreferences();
  const [ladders, setLadders] = useState<{ id: string; label: string }[]>([]);
  const [percentages, setPercentages] = useState<number[]>(PERCENT_FALLBACK);
  const [loading, setLoading] = useState(true);
  const [saved, setSaved] = useState<string | null>(null);

  useEffect(() => {
    api.gelLadders()
      .then((r) => {
        setLadders(r.ladders);
        setPercentages(r.percentages);
      })
      .catch(() => undefined)
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
