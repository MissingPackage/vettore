# Vettore — user guide

Vettore is an offline plasmid editor and cloning simulator. It follows the way bench work
flows — open a construct, read its map and sites, plan the cloning, order the primers, verify
the clone — and runs entirely on your own computer: no account, no internet, no data leaves
the machine.

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

The **Methylation** column warns when a base methylated by Dam (GATC) or Dcm (CCWGG) in the
usual cloning strains falls inside a recognition site *and* REBASE records that enzyme as
blocked, impaired or variable — the reason a digest sometimes fails for no visible reason.
Enzymes that cut regardless, such as BamHI, are left alone; when REBASE has no data the column
says so rather than staying silent.

## Primers and Tm

Select any stretch of DNA and the status bar shows its Tm. Because different tools use
different models, hover the Tm to see three values at once:

- **Standard** — SantaLucia & Hicks 2004, 50 mM Na⁺, 250 nM primer: the textbook conditions
  for a PCR primer, and the default;
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
| Anneal Oligos and Insert… | Clones a short insert (gRNA, shRNA, linker, tag) from two annealed oligos; give it the insert and it writes the two oligos, taking the overhangs from the cut vector itself |

Warnings appear under the product when something deserves a check before you order anything:
an enzyme that cuts more than once, a vector that can re-ligate without insert, a
non-directional ligation, primers that became very long, or an assembly that could also give a
different product.

The biology is computed by [pydna](https://github.com/pydna-group/pydna) and
[Biopython](https://biopython.org), and every operation is covered by tests that compare the
result with a sequence worked out by hand.

## Annotating a new sequence (Features → Detect Common Features…)

Opened a plasmid that arrives bare, or built one from scratch? This searches it, on both
strands, for the elements AAV work uses every day: AAV2 ITRs (including the truncated ones
found in real vectors), CMV, CAG, hSyn, EF-1α, CaMKIIα, GFAP and U6 promoters, WPRE, bGH and
SV40 polyA, AmpR and KanR, pUC and f1 origins, EGFP, mCherry, tdTomato, loxP, lox2272, FRT,
SpCas9, SaCas9 and the FLAG and HA tags.

Every entry in the catalogue is a real subsequence of a public GenBank record, and the table
shows which one it came from along with the percent identity of the match, so you can judge a
hit yourself. Tick the ones you want and they are added to the map, coloured by type, with
their provenance kept in the feature's notes. Features already annotated are shown greyed and
left unticked.

## Verifying a clone by sequencing (Tools → Align to Reference DNA Sequence…)

Open the plasmid you expect, then choose the `.ab1` files your sequencing provider sent. Each
read is quality-trimmed (Mott's algorithm; the stringency is adjustable), aligned to the
plasmid on either strand, and stacked under the reference with its chromatogram, starting at
the region the reads cover.

The table above the view is the answer to "is my clone right?": for every read it gives the
strand, the region covered in plasmid coordinates, the identity, and every difference listed
as `position: reference → read`. Insertions and deletions are grouped into one event each, and
the last line tells you how much of the plasmid the reads cover in total. Reads that cross the
origin of a circular plasmid are handled and marked.

## Simulating a gel (Tools → Simulate Agarose Gel…)

Add one lane per file, choose the enzymes for each lane (leave them empty to load the uncut
plasmid, drawn supercoiled), pick a ladder and the agarose percentage. Band positions follow a
semi-logarithmic migration model: good for predicting the *pattern* of a diagnostic digest and
reading sizes against the ladder, not a physical simulation of your gel rig. Hover a band to
see its size.

Not implemented yet: Gateway BP/LR and TA/GC/TOPO cloning, multiple sequence alignment, and
writing `.dna` files (saving is GenBank).

## Saving and sharing

- **Save** / **Save As…** write **GenBank** (`.gb`), the exchange format every sequence tool
  reads. Features, colours, primers and topology are preserved.
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

Vettore is an independent program built on open-source libraries.
