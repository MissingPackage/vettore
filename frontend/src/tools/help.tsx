// Help menu: keyboard shortcuts, a short guide for a first day with Vettore, and About.
import { useEffect, useState } from "react";
import * as Blueprint from "@blueprintjs/core";
import { Button, Classes } from "@blueprintjs/core";
import { call } from "../api";
import { SHORTCUT_GROUPS } from "./helpShortcuts";
import type { Tool, ToolProps } from "./types";
import "./help.css";

// Blueprint 3 typings predate React 18 (no implicit `children`).
const { Dialog, Callout } = Blueprint as unknown as Record<string, React.ComponentType<any>>;

interface About {
  version: string;
  licence: string;
  log_path: string;
  third_party: { name: string; licence: string; used_for: string }[];
}

function HelpDialog(props: {
  title: string;
  testId: string;
  lead: string;
  width?: number;
  onClose: () => void;
  children: React.ReactNode;
}) {
  return (
    <Dialog isOpen onClose={props.onClose} title={props.title} style={{ width: props.width ?? 680 }}>
      <div className={`${Classes.DIALOG_BODY} vt-help`} data-testid={props.testId}>
        <p className="vt-muted">{props.lead}</p>
        {props.children}
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button intent="primary" onClick={props.onClose}>Close</Button>
        </div>
      </div>
    </Dialog>
  );
}

function Shortcuts(props: ToolProps) {
  return (
    <HelpDialog
      title="Keyboard Shortcuts"
      testId="shortcuts-dialog"
      lead="Every shortcut below works in Vettore today. On a Mac keyboard use ⌘ instead of Ctrl."
      onClose={props.onClose}
    >
      <div className="vt-shortcut-groups">
        {SHORTCUT_GROUPS.map((group) => (
          <section key={group.title}>
            <h4>{group.title}</h4>
            <table className="bp3-html-table bp3-html-table-condensed vt-table">
              <tbody>
                {group.shortcuts.map((shortcut) => (
                  <tr key={shortcut.keys}>
                    <td className="vt-shortcut-keys"><kbd>{shortcut.keys}</kbd></td>
                    <td>{shortcut.what}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </section>
        ))}
      </div>
    </HelpDialog>
  );
}

function GettingStarted(props: ToolProps) {
  return (
    <HelpDialog
      title="Getting Started"
      testId="getting-started-dialog"
      lead="A first pass through Vettore, in the order the bench work happens."
      onClose={props.onClose}
    >
      <ol className="vt-steps">
        <li>
          <b>Open a construct.</b> File → Open… (Ctrl+O), or drag the file onto the window.
          Vettore reads GenBank (.gb, .gbk, .ape), .dna, EMBL, FASTA, .xdna, .gck, Sanger .ab1
          and plain sequence text. A sequence in an email or on a web page: copy it, then File →
          New File from Clipboard (Ctrl+Shift+V). No file to hand? File → Open Sample Construct
          has three public plasmids to try things on.
        </li>
        <li>
          <b>Read the map.</b> The Map view (Ctrl+1) shows features, primers and the cut sites of
          the enzyme set chosen at the top right. The Enzymes view (Ctrl+3) lists every enzyme of
          that set with its sites, the ends it leaves and a warning when Dam or Dcm methylation
          can block the digest. Select any stretch of DNA and the bar at the bottom gives its
          length, %GC and — up to 100 bases — its Tm; hover the Tm for the other models.
        </li>
        <li>
          <b>Annotate what arrived bare.</b> Features → Detect Common Features… searches both
          strands for ITRs, promoters, WPRE, polyA signals, markers, origins, fluorescent
          proteins, tags, lox/FRT sites and Cas9, each traceable to a public GenBank record.
        </li>
        <li>
          <b>Simulate the cloning.</b> The Actions menu holds restriction cloning (Insert
          Fragment, Insert Multiple Fragments, Delete Restriction Fragment), PCR, mutagenesis,
          annealed oligos and ligation, Gibson / NEBuilder HiFi / In-Fusion, Golden Gate,
          Gateway BP and LR reactions, and TA / GC and TOPO cloning. Fill in the left side, watch
          the product on the right, then press Create Product: it opens as a new tab, unsaved
          until you save it. The primers to order are listed underneath, and Copy all puts them
          on the clipboard ready for an order form. Warnings above the button are worth reading
          before you order anything.
        </li>
        <li>
          <b>Design primers.</b> Primers → Design PCR Primers… finds the best pairs around a
          region or a feature. To put primers exactly on a region's ends for cloning, use Actions
          → Design Primers to Amplify a Region… instead.
        </li>
        <li>
          <b>Check the clone.</b> Tools → Simulate Agarose Gel… predicts the pattern of a
          diagnostic digest against a ladder. When the sequencing comes back, Tools → Align to
          Reference DNA Sequence… aligns the .ab1 reads to the construct and lists every
          difference in plasmid coordinates, with the chromatogram under it. Tools → Align
          Multiple Sequences… lines up other open sequences against one reference, and Tools →
          AAV Packaging Check… measures the ITR-to-ITR cargo against what an AAV capsid holds.
        </li>
        <li>
          <b>Keep the result.</b> Ctrl+S saves GenBank, the format every sequence tool reads.
          File → Export also writes FASTA, .dna, and the map as SVG or PNG for a figure. A dot
          before a tab's name means unsaved work: closing the tab, or Vettore, asks first and
          offers to save it. The History view
          (Ctrl+6) keeps where a product came from, down to the parent files and primers.
        </li>
      </ol>
      <Callout intent="primary">
        Vettore does not need the internet and sends nothing anywhere: files stay on this
        computer.
      </Callout>
    </HelpDialog>
  );
}

function AboutVettore(props: ToolProps) {
  const [about, setAbout] = useState<About | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    call<About>("app/about")
      .then(setAbout)
      .catch((e) => setError(String((e as Error).message)));
  }, []);

  return (
    <HelpDialog
      title="About Vettore"
      testId="about-dialog"
      lead="An offline plasmid editor and cloning simulator for the bench."
      width={640}
      onClose={props.onClose}
    >
      {error && (
        <Callout intent="danger" title="Could not read the version">
          {error}. The rest of Vettore is unaffected; restart it if this keeps happening.
        </Callout>
      )}
      {!about && !error && <p className="vt-muted">Reading the version…</p>}
      {about && (
        <>
          <table className="bp3-html-table bp3-html-table-condensed vt-table vt-about-table">
            <tbody>
              <tr>
                <td>Version</td>
                <td data-testid="about-version">{about.version}</td>
              </tr>
              <tr>
                <td>Licence</td>
                <td>
                  {about.licence} — free software. The source and the full third-party notices
                  come with the program.
                </td>
              </tr>
              <tr>
                <td>Log file</td>
                <td>
                  <span className="vt-mono vt-break">{about.log_path}</span>
                  <br />
                  <span className="vt-muted">
                    On Windows: %USERPROFILE%\.vettore\vettore.log. Send this file with the
                    sequence that caused the problem when something goes wrong.
                  </span>
                </td>
              </tr>
            </tbody>
          </table>
          <h4>Built on</h4>
          <table className="bp3-html-table bp3-html-table-condensed vt-table">
            <thead>
              <tr><th>Component</th><th>Licence</th><th>Used for</th></tr>
            </thead>
            <tbody>
              {about.third_party.map((entry) => (
                <tr key={entry.name}>
                  <td>{entry.name}</td>
                  <td className="vt-muted">{entry.licence}</td>
                  <td className="vt-muted">{entry.used_for}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </>
      )}
    </HelpDialog>
  );
}

export const tools: Tool[] = [
  { id: "shortcuts", menu: "Help", label: "Keyboard Shortcuts", Component: Shortcuts },
  { id: "getting-started", menu: "Help", label: "Getting Started", Component: GettingStarted },
  { id: "about", menu: "Help", label: "About Vettore", Component: AboutVettore },
];

/**
 * One help dialog, by tool id, for callers outside the menu bar. The welcome screen needs
 * this: the menu bar belongs to the editor and is not on screen while no document is open.
 */
export function HelpDialogById(props: { id: string; onClose: () => void }) {
  const tool = tools.find((t) => t.id === props.id);
  if (!tool) return null;
  const { Component } = tool;
  // Help dialogs read no document and change none.
  return (
    <Component
      docs={[]}
      active={null}
      onClose={props.onClose}
      openDocument={() => undefined}
      logHistory={() => undefined}
    />
  );
}
