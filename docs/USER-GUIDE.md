# Vettore — user guide

Vettore is an offline plasmid editor and cloning simulator. It follows the way bench work flows —
open, read, plan, order, verify — and runs entirely on your own
computer: no account, no internet, no data leaves the machine.

## Install (Windows)

1. Download `Vettore-windows-x64.zip` from the repository's **Releases** page.
2. Right-click the zip → **Extract All…** and pick a folder you own (Desktop or Documents).
3. Open the extracted `Vettore` folder and double-click `Vettore.exe`.
4. Windows SmartScreen may say "Windows protected your PC" because the app is not signed:
   click **More info → Run anyway**. This happens only the first time.

Nothing is installed and no administrator rights are needed. To remove Vettore, delete the
folder.

On Linux: extract `Vettore-linux-x64.zip` and run `./Vettore`.

## Opening files

- **File → Open…** (Ctrl+O), or drag a file onto the window.
- Reads `.dna`, GenBank (`.gb`, `.gbk`, `.ape`), FASTA and Sanger `.ab1` files.
- Several files stay open as tabs. A dot (•) on a tab means unsaved changes.
- **File → Open Recent** lists the last ten files.

Your `.dna` files open with their features, colours, primers (including 5′ tails) and
topology. Vettore never writes to a `.dna` file, so the originals cannot be damaged.

## The window

| Part | What it does |
|---|---|
| Menu bar | File, Edit, View, Enzymes, Features, Primers, Actions, Tools, Help |
| Toolbar | Save, import/export, undo/redo, cut sites, features, ORFs, find, visibility |
| Document tabs | One per open file |
| View tabs | Map, Sequence, Enzymes, Features, Primers, History — also Ctrl+1 … Ctrl+6 |
| Enzyme set | The set shown on the map (top right) |
| Status bar | Length, topology and %GC; for a selection also its range, length, %GC and Tm |

## Views

- **Map** — circular or linear map with features, cut sites and primers. Click a feature to
  select it; the selection is shown in every view.
- **Sequence** — both strands with features, translations, cut sites and primers.
- **Enzymes** — every enzyme of the chosen set that cuts, with the number of sites, the cut
  positions and the kind of ends. Click a position to select that site on the map.
- **Features** / **Primers** — sortable tables; double-click a row to edit it.
- **History** — where the file came from: opened, saved, or produced by a cloning simulation
  (with the parent files, the parameters and the primers used).

## Enzymes

The enzyme list comes from REBASE (through Biopython) and contains one enzyme per specificity,
so you see `BamHI` rather than a dozen isoschizomers. Choose a set from the box at the top
right or from **Enzymes → Use Enzyme Set**:

Unique 6+ cutters (the default), Unique cutters, Unique & dual cutters, 6+ cutters,
Type IIS enzymes, Golden Gate enzymes, Nonredundant commercial.

**Enzymes → Choose Enzymes…** opens the full picker for a custom list.

## Primers and Tm

Select any stretch of DNA and the status bar shows its Tm. Because different tools use
different models, hover the Tm to see three values at once:

- **Standard** — SantaLucia & Hicks 2004, 50 mM Na⁺, 250 nM primer (the default,
  the textbook conditions for a PCR primer);
- **NEB Q5 / high-fidelity** and **NEB Taq** — the conditions NEB's Tm calculator uses.

The Tm of a primer is always the Tm of the part that anneals: 5′ tails are excluded.

## Cloning simulations (Actions menu)

Every tool works the same way: fill in the left side, watch the product appear on the right,
then press **Create Product**. The product opens as a new tab and its History records what was
done. Primers that need to be ordered are listed with their length, annealing length and Tm,
and **Copy all** puts them on the clipboard as `name ⇥ sequence`, ready to paste into an order
form.

| Tool | What it does |
|---|---|
| Restriction and Insertion Cloning → Insert Fragment… | Cuts vector and insert, ligates the insert in the chosen orientation |
| Restriction and Insertion Cloning → Delete Restriction Fragment… | Removes the fragment between two sites and re-ligates |
| PCR… | Amplifies a template with two primers (tails included) |
| Design Primers to Amplify a Region… | Picks primers for a region at a target Tm, with optional 5′ tails |
| Mutagenesis… | Substitution, insertion or deletion with back-to-back Q5-style primers |
| Gibson Assembly… / NEBuilder HiFi… / In-Fusion… | Inserts one or more fragments into a linearised vector and designs the homology arms |
| Golden Gate Assembly… | Assembles type IIS parts in order, checking the junction overhangs |

Warnings appear under the product when something deserves a check before you order anything:
an enzyme that cuts more than once, a vector that can re-ligate without insert, a
non-directional ligation, primers that became very long, or an assembly that could also give a
different product.

The biology is computed by [pydna](https://github.com/pydna-group/pydna) and
[Biopython](https://biopython.org), and every operation is covered by tests that compare the
result with a sequence worked out by hand.

## Simulating a gel (Tools → Simulate Agarose Gel…)

Add one lane per file, choose the enzymes for each lane (leave them empty to load the uncut
plasmid, drawn supercoiled), pick a ladder and the agarose percentage. Band positions follow a
semi-logarithmic migration model: good for predicting the *pattern* of a diagnostic digest and
reading sizes against the ladder, not a physical simulation of your gel rig. Hover a band to
see its size.

Not implemented yet: Gateway BP/LR, TA/GC/TOPO cloning, Sanger (.ab1) alignment to a reference,
and automatic detection of common features.

## Saving and sharing

- **Save** / **Save As…** write **GenBank** (`.gb`), which Benchling, ApE and
  Geneious all open. Features, colours, primers and topology are preserved.
- **File → Export** also writes FASTA, and the map as **SVG** (for figures) or **PNG**.

## Keyboard shortcuts

| | |
|---|---|
| Ctrl+O / Ctrl+N / Ctrl+W | Open / New / Close tab |
| Ctrl+S | Save |
| Ctrl+1 … Ctrl+6 | Map, Sequence, Enzymes, Features, Primers, History |
| Ctrl+F | Find DNA or protein, with mismatches |
| Ctrl+Z / Ctrl+Shift+Z | Undo / Redo |

## If something goes wrong

The log is at `%USERPROFILE%\.vettore\vettore.log` on Windows (`~/.vettore/vettore.log` on
Linux). Send it along with the file that caused the problem.

Vettore is an independent program built on open-source
libraries.
