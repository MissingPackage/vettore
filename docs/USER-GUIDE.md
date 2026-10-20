# Vettore — user guide

Vettore is an offline plasmid editor and cloning simulator. It follows the way bench work
flows — open a construct, read its map and sites, plan the cloning, order the primers, verify
the clone — and runs entirely on your own computer: no account, no internet, no data leaves
the machine.

## Contents

- [Install](#install-windows)
- [Getting started](#getting-started)
- [Opening files](#opening-files)
- [Copy and paste sequences](#copy-and-paste-sequences)
- [The window](#the-window)
- [Views](#views)
- [Enzymes](#enzymes)
- [Primers and Tm](#primers-and-tm)
- [Designing PCR primers](#designing-pcr-primers-primers--design-pcr-primers)
- [Cloning simulations](#cloning-simulations-actions-menu)
- [Annotating a new sequence](#annotating-a-new-sequence-features--detect-common-features)
- [Checking what an AAV capsid can carry](#checking-what-an-aav-capsid-can-carry-tools--aav-packaging-check)
- [Comparing several sequences](#comparing-several-sequences-tools--align-multiple-sequences)
- [Verifying a clone by sequencing](#verifying-a-clone-by-sequencing-tools--align-to-reference-dna-sequence)
- [Simulating a gel](#simulating-a-gel-tools--simulate-agarose-gel)
- [Saving and sharing](#saving-and-sharing)
- [History of a file](#history-of-a-file)
- [Preferences](#preferences-file--preferences)
- [Help menu](#help-menu)
- [Keyboard shortcuts](#keyboard-shortcuts)
- [If something goes wrong](#if-something-goes-wrong)

## Install (Windows)

1. Download `Vettore-windows-x64.zip` from the repository's **Releases** page.
2. Right-click the zip → **Extract All…** and pick a folder you own (Desktop or Documents).
3. Open the extracted `Vettore` folder and double-click `Vettore.exe`.
4. Windows SmartScreen may say "Windows protected your PC" because the app is not signed:
   click **More info → Run anyway**. This happens only the first time.

Nothing is installed and no administrator rights are needed. To remove Vettore, delete the
folder. To update, download the new zip and replace the folder: your recent files and
preferences are kept.

On Linux: extract `Vettore-linux-x64.zip` and run `./Vettore`.

## Getting started

The screen you see before opening anything is the way in: **Open a file…** (Ctrl+O),
**New DNA file…** (Ctrl+N), **Paste a sequence…** (Ctrl+Shift+V), your last ten files, and three
sample constructs that come with Vettore:

| Sample | What it is |
|---|---|
| pUC19 | The classic cloning vector. Its GenBank record carries few annotations, so it is a good construct to try **Features → Detect Common Features…** on |
| pAAV-cTNT-SaCas9-U6 | An AAV transfer plasmid with both ITRs annotated — a good one for the AAV packaging check |
| pAAV-EF1A-DIO-Gephyrin.FingR-GFP | A Cre-dependent (DIO) AAV transfer plasmid with a GFP reporter |

They are public GenBank records (the accession is shown next to each one, and the record's own
authors and references travel inside the file). A sample opens as an **unsaved** document, so
you can cut it, edit it and simulate on it; **Save** asks where to put your copy and never
touches the one inside the program. Once files are open, the samples are under **File → Open
Sample Construct**.

At the bottom of that screen, **Getting started**, **Keyboard shortcuts** and **About Vettore**
open the same help you find in the Help menu once a file is open.

## Opening files

- **File → Open…** (Ctrl+O), drag files onto the window, or **File → Open Recent**. Open…
  starts in the folder you used last.
- Vettore opens GenBank (`.gb`, `.gbk`, `.ape`), `.dna`, EMBL (`.embl`), FASTA (`.fa`,
  `.fasta`), `.xdna`, `.gck`, Sanger traces (`.ab1`) and plain sequence text (`.txt`, `.seq`:
  just the DNA letters, with or without position numbers).
- Vettore looks at what is inside a file, not only at its name. A GenBank file a colleague saved
  as `.txt`, a FASTA file called `.seq` or a file with no extension opens correctly.
- A file with several sequences, for example a FASTA file with three plasmids, opens one tab per
  sequence. A message tells you how many were opened.
- Files made by other programs often bend the rules: very long names, Windows line endings,
  accented letters, features drawn across the origin. Vettore opens them anyway. A feature that
  cannot be placed on the sequence (between two bases, beyond the end, or in another record) is
  left out, and a message lists it by name.
- Files from vector vendors keep their look. A feature colour written in a note, such as
  `color: #63b349; direction: RIGHT` or a colour name like `color: Blue`, colours that feature on
  the map. When you save, the colour is written once, in the standard colour field.
- Feature types spelled in other ways — `Promoter`, `PolyA_signal`, `Rep_origin`, `ITR`,
  `Transcription_t`, `Miscellaneous`, in any letter case, or the `regulatory` features of
  current NCBI records with their class (promoter, enhancer, terminator…) — are shown as the
  standard types, so the map colours and groups them like any other file. The type written in
  the file is kept in the feature's notes, and a `regulatory` feature is saved back in its own
  form.
- If a file cannot be opened, the message names the file, says what is wrong and what to try.
  For a plain sequence it points to the line and column of the first character that is not a
  DNA letter.
- Several files stay open as tabs. A dot (•) on a tab means unsaved changes — including a
  new product of a cloning simulation or a pasted sequence that has never been saved. Closing
  such a tab asks whether to **Save**, **Don't save** or **Cancel**.
- **Closing Vettore** with unsaved documents lists them all in one question. **Save All**
  saves each one the way **File → Save** does: it asks where to put a document that was never
  saved, and it makes a GenBank copy of a `.dna` file. Then Vettore closes. **Don't Save**
  closes and discards the changes. **Cancel** keeps everything open. If a save is cancelled or
  fails, Vettore stays open and says which document was not saved.
- **Double-click a sequence file while Vettore is open** and it opens as a new tab in the
  window you already have; no second Vettore starts. A file that is already open is brought to
  the front instead of being opened twice, and starting Vettore again just brings its window
  back. To make double-clicking open sequence files in Vettore at all, see
  [Preferences](#preferences-file--preferences).

Your `.dna` files open with their features, colours, primers (including 5′ tails) and topology.

## Copy and paste sequences

Colleagues often send sequences in emails and chats. Vettore can open them, and copy sequences
for you to send back.

- **File → New File from Clipboard** (Ctrl+Shift+V) opens what you copied as a new file: a plain
  sequence, a FASTA record or a whole GenBank or EMBL record. With no file open you can also
  just press Ctrl+V on the start screen.
- A box to paste into opens instead when Vettore is not allowed to read the clipboard, or when
  the clipboard does not hold a sequence. Paste the text there with Ctrl+V. While you type, the
  box shows what it recognised (for example "GenBank: 7,204 bp, circular") or what is wrong,
  down to the line and column. For a plain sequence you can also give a name and choose linear
  or circular.
- **Edit → Copy as FASTA** copies the selection as a FASTA record, ready for an email or a
  synthesis order. With nothing selected it copies the whole sequence. The first line carries
  the positions of the selection, for example `>pUC19:101-623`. You also find it in the
  right-click **Copy** menu. If your system does not let Vettore use the clipboard, Vettore
  shows the text already selected: press Ctrl+C to copy it.

## The window

| Part | What it does |
|---|---|
| Menu bar | File, Edit, View, Enzymes, Features, Primers, Actions, Tools, Help |
| Toolbar | Save, undo/redo, cut sites, features, ORFs, find, visibility |
| Document tabs | One per open file; the open one has a blue line along its top edge |
| View tabs | Map, Sequence, Enzymes, Features, Primers, History — also Ctrl+1 … Ctrl+6 |
| Enzyme set | The set shown on the map (top right) |
| Status bar | Length, topology and %GC; for a selection also its range, length, %GC and Tm |
| Messages | Top right; they stay until you close them (×) or click **Dismiss all**. Confirmations such as "Saved" disappear on their own |

Press Tab to move between the view tabs, buttons and fields: a blue ring shows which one is
selected, and Enter or Space uses it.

## Views

- **Map** — circular or linear map with features, cut sites and primers. Click a feature to
  select it; the selection is shown in every view.
- **Sequence** — both strands with features, translations, cut sites and primers.
- **Enzymes** — every enzyme of the chosen set that cuts, with the number of sites, the cut
  positions and the kind of ends. Click a position to select that site; switch to **Map**
  (Ctrl+1) to see it there.
- **Features** / **Primers** — sortable tables; double-click a row to edit it.
- **History** — where the file came from (see [History of a file](#history-of-a-file)).

## Enzymes

The enzyme list comes from REBASE (through Biopython) and contains one enzyme per specificity,
so you see `BamHI` rather than a dozen isoschizomers. Choose a set from the box at the top
right or from **Enzymes → Use Enzyme Set**:

Unique 6+ cutters (the default), Unique cutters, Unique & dual cutters, 6+ cutters,
Type IIS enzymes, Golden Gate enzymes, Nonredundant commercial.

**Enzymes → Choose Enzymes…** opens the full picker for a custom list. The set every new
document starts with is a preference (**File → Preferences…**).

The **Methylation** column warns when a base methylated by Dam (GATC) or Dcm (CCWGG) in the
usual cloning strains falls inside a recognition site *and* REBASE records that enzyme as
blocked, impaired or variable — the reason a digest sometimes fails for no visible reason.
Enzymes that cut regardless, such as BamHI, are left alone; when REBASE has no data the column
says so rather than staying silent.

## Primers and Tm

Select a stretch of DNA of up to 100 bases — a primer-sized region — and the status bar shows
its Tm. Because different tools use
different models, hover the Tm (it is underlined with dots) to see three values at once:

- **Standard** — SantaLucia & Hicks 2004, 50 mM Na⁺, 250 nM primer: the textbook conditions
  for a PCR primer, and the default;
- **NEB Q5 / high-fidelity** and **NEB Taq** — the conditions NEB's Tm calculator uses.

The Tm of a primer is always the Tm of the part that anneals: 5′ tails are excluded.

## Designing PCR primers (Primers → Design PCR Primers…)

To pick a region before opening the dialog, select it in the map or the sequence. Or open the dialog and pick a feature from the list, or type the start and end. Choose **Primers → Design PCR Primers…**. primer3 proposes the five best primer pairs whose product contains the whole region, with the primers just outside it. To put the primers exactly on the region's ends (for cloning), use **Actions → Design Primers to Amplify a Region** instead.

The defaults suit a standard PCR:

- primers 18–24 nt long (20 ideal);
- Tm 57–63 °C (60 ideal), the two primers within 5 °C of each other;
- GC 40–60 %;
- a product up to 150 bp longer than the region on each side.

Change any of them on the left; the pairs update as you type, and a message next to a field tells you if a value cannot work.

For each pair you see the product size, and for each primer its position, length, Tm, GC and the ΔG of its hairpin and self-dimer. The pair also shows the ΔG of the two primers pairing with each other (primer-dimer). Tm uses the Standard method (SantaLucia & Hicks 2004, Na⁺ 50 mM, primer 250 nM), the same as the status bar. Hover a Tm to see the NEB Q5 and Taq values. A warning appears when a structure is more stable than −9 kcal/mol, the limit IDT recommends, or when the 3′ end of a primer also matches another place in the plasmid. Pairs without warnings are listed first.

If no pair fits, the dialog says why (for example *Tm too low*) and what to widen.

**Add to Sequence** adds the chosen pair to the plasmid's primers, named after the feature or the file (`KanR_F`, `KanR_R`), and records it in History. **Copy for Ordering** copies `name ⇥ sequence` lines, ready to paste into an order form.

## Cloning simulations (Actions menu)

Every tool works the same way: fill in the left side, watch the product appear on the right,
then press **Create Product**. The dialog starts from the files you have open and from the
region selected in the editor, and says which input is still missing. The product opens as a
new tab and its History records what was done. Primers that need to be ordered are listed with
their length, annealing length and Tm, and **Copy all** puts them on the clipboard as
`name ⇥ sequence`, ready to paste into an order form.

| Tool | What it does |
|---|---|
| Restriction and Insertion Cloning → Insert Fragment… | Cuts vector and insert, ligates the insert in the chosen orientation |
| Restriction and Insertion Cloning → Insert Multiple Fragments… | Ligates a cut vector and two or more inserts in one reaction |
| Restriction and Insertion Cloning → Delete Restriction Fragment… | Removes the fragment between two sites and re-ligates |
| PCR… | Amplifies a template with two primers (tails included) |
| Design Primers to Amplify a Region… | Places primers exactly on a region's ends at a target Tm, with optional 5′ tails |
| Mutagenesis… | Substitution, insertion or deletion with back-to-back Q5-style primers |
| Anneal Oligos and Insert… | Clones a short insert (gRNA, shRNA, linker, tag) from two annealed oligos; give it the insert and it writes the two oligos, taking the overhangs from the cut vector itself |
| Ligate Fragments… | Joins fragments end to end into a linear molecule, or closes them into a circle |
| Gibson Assembly… / NEBuilder HiFi… / In-Fusion… | Inserts one or more fragments into a linearised vector and designs the homology arms |
| Golden Gate Assembly… | Assembles type IIS parts in order, checking the junction overhangs |
| Gateway BP Reaction… / Gateway LR Reaction… | att site recombination into an entry or expression clone; LR takes up to four entry clones (MultiSite) |
| TA or GC Cloning… | A PCR product into a T-vector or C-vector, both orientations |
| TOPO Cloning… | TA, blunt or directional TOPO; directional checks the CACC on the forward primer |

When you add an enzyme to a dialog it appears as a chip showing how many times it cuts the
chosen file: **green** ("EcoRI · 1 site") when it cuts once, which is what most cloning steps
need; **amber** ("0 sites" or "2 sites") when it does not cut or cuts more than once — check
that this is what you want. Click the × on a chip to remove the enzyme.

Warnings appear under the product when something deserves a check before you order anything:
an enzyme that cuts more than once, a vector that can re-ligate without insert, a
non-directional ligation, primers that became very long, or an assembly that could also give a
different product. If no product is possible, a red box says why and names the input to
change; the result updates as soon as you do.

The biology is computed by [pydna](https://github.com/pydna-group/pydna) and
[Biopython](https://biopython.org), and every operation is covered by tests that compare the
result with a sequence worked out by hand.

### Insert several fragments in one ligation (Actions → Restriction and Insertion Cloning → Insert Multiple Fragments…)

Use this when a construct needs more than one piece: a vector plus a promoter and a gene, say,
cut so that each junction has its own matching ends.

For the vector and for every insert you choose a file and how the piece is taken from it:

- **Cut out with restriction enzymes** — pick one or two enzymes. The pieces of the digest are
  listed with their size and position, so you can see which one will be used and choose another
  if it is not the one you meant. Left on **Automatic**, Vettore takes the largest piece for the
  vector and the piece between the two sites for an insert.
- **As a PCR product** — the whole file, or the region you type (or select in the editor). The
  ends are blunt, as a proofreading polymerase leaves them.

The preview shows the product as soon as the inputs are complete. **The order you list the
inserts in does not force the result**: the ends decide how the pieces fit, and if they only fit
another way round Vettore says so and shows that product instead. When several products are
possible, each one is listed with its order and orientation — pick the one you want before
creating it.

Before you order anything, read the warnings: a vector whose two ends match each other will
close without an insert (dephosphorylate it with rSAP or CIP), blunt ligations are less
efficient and want more insert, and PCR products carry no 5′ phosphate unless you ordered the
primers phosphorylated or treat them with T4 polynucleotide kinase.

### Join fragments into a linear molecule (Actions → Ligate Fragments…)

The same fragment chooser, without a vector: two or more pieces joined end to end. Tick **Close
the product into a circle** to circularise the result instead. This is the tool for joining two
PCR products, for building a linear construct for injection or transfection, and for
re-circularising a single cut fragment.

If the fragments cannot be joined, the message names the junction that does not fit — for
example *"Promoter" ends with 5' overhang TCGA but "CDS" starts with 5' overhang AATT* — so you
know which enzyme to change.

### TA and GC cloning (Actions → TA or GC Cloning…)

A PCR product straight into a T-vector: the polymerase leaves a single 3′ A on the product, the
vector carries a single 3′ T, and the two pair. GC cloning is the same with a 3′ G insert and a
3′ C vector.

Give Vettore the vector as DNA without its single-base overhangs — either a **linearised vector
file**, or a **circular map plus the base after which it opens** (type the position, or select
the cloning site in the editor and press *After the selection*). Vettore adds the T (or C) to
the vector and the A (or G) to the PCR product itself.

The insert is either an open file (the product as amplified) or a template plus two primers,
amplified here. Because nothing holds the insert one way round, both orientations are shown as
separate products: choose one, and screen your colonies to find it.

The warning reminds you that these overhangs come from a non-proofreading polymerase; a blunt
product from a proofreading enzyme needs an A-tailing step first.

### TOPO cloning (Actions → TOPO Cloning…)

Three methods in one dialog, with the same kind of vector file as TA cloning:

- **TOPO TA** — same ends as TA cloning, both orientations.
- **Blunt TOPO** — a blunt PCR product from a proofreading polymerase, both orientations.
- **Directional TOPO** — the forward primer starts with **CACC**, which pairs with the vector's
  GTGG overhang and holds the insert one way round. One product, in the direction of the
  vector's numbering.

For directional TOPO, Vettore checks the CACC for you. If the forward primer does not start with
it, the dialog says so next to the primer and offers an **Add CACC** button that writes the
primer you should order; if the reverse primer is the one carrying CACC, you are warned that the
insert will go in backwards.

### Gateway BP and LR reactions (Actions → Gateway BP Reaction… / Gateway LR Reaction…)

Recombination between att sites, with no enzymes to choose:

- **BP** — an attB-flanked insert (a PCR product with attB tails, or a clone carrying attB sites)
  plus a donor vector with attP sites give the **entry clone**. Under the donor you choose which
  **att pair** it carries (attP1 and attP2, attP4 and attP1r, attP2r and attP3…); that pair
  decides the sites of the entry clone you get.
- **LR** — one to four entry clones (attL sites) plus a destination vector (attR sites) give the
  **expression clone**.

**MultiSite reactions (two, three or four fragments).** Use **+ Add entry clone** to put a
second, third or fourth entry clone into the same reaction, and **Remove** to take one out. The
entry clones are listed in the order they will appear in the product. When the open files fit
together — a destination vector with attR4 and attR3, a 5′ entry clone with attL4 and attR1, a
middle clone with attL1 and attL2, a 3′ clone with attR2 and attL3 — the dialog fills them in,
in that order, as soon as it opens. Under each file you see the att sites found in it, so you
can check the chain at a glance: attL4, attR1 → attL1, attL2 → attL3, attR2.

Only reactions where **both** att sites of a junction recombined are shown; a product from a
single site is the intermediate of the reaction, not what grows on the plates. The clone you are
after is listed first and the by-product after it, each with its size, so you can create either
one — or both at once with **Create All**.

Two warnings are worth knowing. If the product still carries a **ccdB** gene it would not grow in
a normal cloning strain: that happens with a vector built from two Gateway cassettes, where the
second cassette has still to be filled (run the reaction again on the product), and with an empty
entry clone that still carries its ccdB stuffer — the warning names the file it came from. And
when a vector has more than one cassette, you are told that one clone per cassette is possible,
so you can pick the one you want.

A molecule with no att sites, or with the wrong ones for the reaction, is refused with a message
that names the file and lists the sites found in each input. When the sites are right but do not
pair up, the message says which junction is missing — "attR1 in “p5E-CAGGS” has no attL1 in the
other molecules" — so you can see which entry clone is not there.

**att sites that are not Gateway sites.** The attB and attP sites of the phiC31 integrase, used
by fly landing-site vectors such as pUASTattB, look like att sites but no Gateway clonase
recombines them. Such a file is never chosen for you, is listed as "phiC31 attB — not a Gateway
site", and if you pick it anyway the dialog says so instead of letting the reaction fail.

### Choosing between possible products

Several cloning tools can give more than one product: two orientations of an insert, several
orders of fragments, a clone and its by-product. They are listed above the preview as
*Expected* and *Alternative 1, 2…*, each described in the terms that tell them apart — the order
of the fragments, whether the insert is forward or reversed, or whether it is the clone or the
by-product. Selecting one updates the map, the description and the warnings below it. **Create
Product** creates the one selected; run the dialog again to create another.
**Create All (N)** creates every product at once, each in its own tab with its own History
entry — useful when you want to keep both the clone and its by-product, or both orientations of
an insert.

You can also change **Name of the new file** before creating it: the name is what the History of
every later step will use to refer to this file.

## Annotating a new sequence (Features → Detect Common Features…)

Opened a plasmid that arrives bare, or built one from scratch? This searches it, on both
strands, for the elements AAV work uses every day: AAV2 ITRs (including the truncated ones
found in real vectors), CMV, CAG, hSyn, EF-1α, CaMKIIα, GFAP and U6 promoters, WPRE, bGH and
SV40 polyA, AmpR and KanR, pUC and f1 origins, EGFP, mCherry, tdTomato, loxP, lox2272, FRT,
SpCas9, SaCas9 and the FLAG and HA tags.

Every entry in the catalogue is a real subsequence of a public GenBank record. The table shows
which record it came from and the percent identity of the match, so you can judge a hit
yourself. Tick the ones you want and they are added to the map, coloured by type, with their
provenance kept in the feature's notes.

Features already annotated are shown greyed and left unticked. This includes elements inside a
larger feature named after them, however long that feature is: vendors often draw one "CMV
promoter" block over both the CMV enhancer and the CMV promoter, and neither is offered again.

**Labels that may be wrong.** Sometimes most of an existing feature is a nearly exact match
(98 % identity or more) for a different element than its name says. Its rows are then marked in
amber, for example: “Annotated as “PGK promoter”, but the sequence is the CMV enhancer (100 %
identity). Check the label.” A notice at the top says how many labels to check. These rows are
not ticked, so nothing changes unless you tick one: the element is then added under its
catalogue name, and you can rename or delete the old feature in the Features view. A feature
whose name is not an element, such as "my insert", is never questioned.

When a promoter or an ITR is a few bases off, loosen the identity threshold for one run; the
default is a preference.

## Checking what an AAV capsid can carry (Tools → AAV Packaging Check…)

Open your transfer plasmid and the check answers one question: is the genome short enough to be
packaged?

- It finds the **ITRs** from the plasmid's own annotation — a feature counts as an ITR when its
  name or notes say "ITR" or "inverted terminal repeat", or when it is a repeat region marked
  as inverted.
- It reports the size **from the start of one ITR to the end of the other, ITRs included**
  (what the capsid actually carries) and, separately, what lies **between** the ITRs (promoter,
  transgene, polyA — the part you can still shorten).
- The bar compares that size with the room available: about **4,700 bp** for a single-stranded
  vector, roughly the size of wild-type AAV2; above **5,000 bp** genomes come out truncated if
  they are packaged at all. A **self-complementary** vector carries the cassette twice, so it
  holds about half, **2,350 bp**. Switch between the two with the **Vector** box. The sources
  for these numbers are printed under the result.

Three situations it handles for you:

- **The cassette crosses position 1.** On a circular plasmid the two ITRs can be read either
  way round. Vettore takes the side that does *not* carry the origin of replication and the
  resistance marker, says so, and lets you choose the other side from **Packaged side**.
- **More than two ITRs annotated.** Pick the 5′ and 3′ ITR in the two boxes.
- **No ITR annotated.** Press **Search for ITRs in the sequence**: the same catalogue used by
  Detect Common Features is searched, and the sizes are computed from what it finds, without
  changing your document. To keep those annotations, run **Features → Detect Common
  Features…**. If nothing is found, annotate the region by hand (select it, then **Features →
  Add Feature…**, with "ITR" in the name) and run the check again.

The sizes are counted from the annotated ITR ends, so they are only as exact as the annotation:
a truncated ITR annotated at its real length gives the real cargo.

## Comparing several sequences (Tools → Align Multiple Sequences…)

Open the sequences you want to compare: capsid genes of different serotypes, versions of one construct, a clone and its design. Choose **Tools → Align Multiple Sequences…**. The sequence in the active tab is the **reference**. Under **Compare with**, tick the sequences to line up against it (when only one other file is open, it is already ticked), then press **Align**.

Each sequence is lined up against the reference on whichever strand fits. A circular sequence is turned so that it starts where the reference starts, so plasmids whose base 1 sits in different places still compare correctly. A fragment that crosses the origin of a circular reference is handled too.

The table gives, for every sequence:

- **Strand**: forward, or reverse complement if it fitted the other strand.
- **Covers**: the stretch of the reference it covers.
- **Coverage**: how much of the reference that stretch is.
- **Identity**: identical bases over all aligned positions, gaps included.
- **Differences**: how many mismatches, insertions and deletions.

Click a row to list every difference in the reference's numbering, for example `1,204: A → G`, `601–604: 4 bp missing (ACGT)` or `between 1,001 and 1,002: 7 bp extra (GATTACA)`. Bases are always given on the reference strand. **Copy list** puts the differences on the clipboard as columns, ready to paste into a spreadsheet. The view below stacks every sequence under the reference and marks the differences. Reference features can be switched on from the eye menu of the view.

If a sequence resembles the reference only in part (two plasmids that share a backbone but carry different cassettes), the table says *partial match* and only the shared region is shown. A sequence that shares nothing with the reference, or only a few bases by chance, is left out, with a note saying so.

Limits: up to 10 sequences at a time (reference included), each up to 15,000 bp. A large comparison takes a few seconds per sequence. To compare part of a bigger sequence, make a new file from that region first (**Actions → New File from Selection**).

## Verifying a clone by sequencing (Tools → Align to Reference DNA Sequence…)

Open the plasmid you expect, then choose the `.ab1` files your sequencing provider sent. Several at once is fine. Each read is quality-trimmed (Mott's algorithm; the stringency is adjustable), aligned to the plasmid on either strand, and stacked under the reference with its chromatogram, starting at the region the reads cover.

A sentence above the table answers "is my clone right?". Either the reads match wherever they cover the plasmid, or it says how many differences were found. For every read the table gives:

- the strand;
- the region it covers, in plasmid coordinates;
- the identity;
- the mean read **quality** (Q20 means 1 wrong base call in 100, Q30 1 in 1,000);
- every difference, written as `position: reference → read`.

Next to each difference is the read quality at that base. Hover it to see the chance that the base call is wrong. Differences on low-quality bases (under Q20) are flagged: check those peaks in the chromatogram before trusting them. Long lists show the first few differences, and **show all** opens the rest. Insertions and deletions are grouped into one event each. Reads that cross the origin of a circular plasmid are handled and marked. If no read can be aligned, the dialog explains why and what to try (gentler trimming, or the right plasmid tab).

## Simulating a gel (Tools → Simulate Agarose Gel…)

Add one lane per file, choose the enzymes for each lane (leave them empty to load the uncut
plasmid, drawn supercoiled), pick a ladder and the agarose percentage — the dialog opens on
your usual ones, set in **File → Preferences…**. Band positions follow a semi-logarithmic
migration model: good for predicting the *pattern* of a diagnostic digest and reading sizes
against the ladder, not a physical simulation of your gel rig. Hover a band to see its size.
Shortcut: Ctrl+Shift+G.

## Saving and sharing

- **Save** (Ctrl+S) writes back to the file you opened when that is safe. That means a GenBank or FASTA file holding just this one sequence, or a file Vettore saved itself.
- Sometimes writing over the file would lose something:
  - a `.dna` file, which can also hold history and notes that Vettore does not read;
  - an EMBL, `.xdna`, `.gck` or `.ab1` file;
  - one sequence out of a file that holds several;
  - a file where some features could not be shown.

  In these cases Save explains this once, then asks where to save a **GenBank copy**. It suggests the same name with `.gb`, in the same folder. The original file is never changed.
- **Save As…** (Ctrl+Shift+S) writes:
  - **GenBank** (`.gb`), the format every sequence program reads and the one we recommend for exchanging files;
  - **FASTA** (`.fasta`);
  - a **`.dna`** file.

  The format follows the extension you type. If you type none, Vettore adds `.gb`.
- FASTA keeps only the sequence. If you save as FASTA a document that has features or primers, a message reminds you that they are not in that file.
- A `.dna` file written by Vettore contains the sequence, its topology, the features with their colours, the primers with their 5′ tails, and the description. Editing history and notes kept by other programs are not written.
- **File → Export** writes a GenBank, FASTA or `.dna` copy without changing which file the document is saved to. It also exports the map as **SVG** (for figures) or **PNG**.
- **File → Print…** (Ctrl+P) prints the map. In Vettore's window on Linux printing is not
  available yet: Vettore says so instead of printing. Export the map with
  **File → Export → Map as PNG…** and print the image from your image viewer.

## History of a file

Every file records where it came from: opened, created, saved, and every cloning step with its
parameters and the primers it produced. A cloning step also keeps a copy of everything it was run
on, so **Open "…"** shows an input exactly as it was when the step ran, even if that file has
changed since or was never saved. **Copy** puts the step — parameters and primers — on the
clipboard for your notebook.

## Preferences (File → Preferences…)

What Vettore starts from, kept between sessions:

- **Enzyme set for documents you open** — the set shown on the map and in the Enzymes view of
  every file you open next. Documents already open keep what you chose for them.
- **Gel ladder** and **agarose percentage** — what **Tools → Simulate Agarose Gel…** opens on,
  so your usual gel is one click away.
- **Minimum identity when detecting common features** — 95 % by default. Lower it to catch
  variants (a promoter a few bases off), raise it to 100 % to accept exact matches only.

- **Open sequence files with Vettore when you double-click them** — only in the Vettore you
  downloaded. Tick it and `.dna`, `.gb`, `.gbk`, `.gbff`, `.ape`, `.embl`, `.fa`, `.fasta`,
  `.xdna`, `.gck` and `.ab1` files open in Vettore when you double-click them. It applies to
  your user account only and needs no administrator rights. Untick it to undo: the files open
  as they did before. If you move the Vettore folder or download a new version, tick it again
  so files open with the copy you use now; Preferences says when they still point at another
  copy. If you once told Windows to **Always** open, say, `.dna` files with another program,
  Windows keeps that choice and Preferences says so. To change it, right-click such a file,
  choose **Open with → Choose another app**, pick Vettore and tick **Always use this app**.

Every change applies at once and is saved; **Restore defaults** puts everything back (it does
not change which program opens your files).

## Help menu

- **Keyboard Shortcuts** — every shortcut that really works, grouped by what you are doing.
- **Getting Started** — a short pass through a first construct: open a file, read the map,
  annotate it, simulate the cloning, check the clone, save it.
- **About Vettore** — the version, the licence (GPL-3.0-or-later), the components Vettore is
  built on, and where the log file is: `~/.vettore/vettore.log`
  (`%USERPROFILE%\.vettore\vettore.log` on Windows). Send that file along with the sequence
  when something goes wrong.

## Keyboard shortcuts

On a Mac keyboard, use ⌘ instead of Ctrl.

| Files and tabs | |
|---|---|
| Ctrl+O / Ctrl+N | Open a sequence file / New DNA file |
| Ctrl+Shift+V | New file from a sequence you copied |
| Ctrl+S / Ctrl+Shift+S | Save (GenBank) / Save As… |
| Ctrl+P | Print the current view |
| Ctrl+W | Close the open document |
| Ctrl+Shift+G | Simulate an agarose gel |

| Views | |
|---|---|
| Ctrl+1 … Ctrl+6 | Map, Sequence, Enzymes, Features, Primers, History |

| Finding and selecting | |
|---|---|
| Ctrl+F | Find DNA or protein, with mismatches |
| Ctrl+G | Go to a position |
| Ctrl+A / Ctrl+I | Select everything / select everything except the selection |
| Ctrl+/ | Search the menus by name |

| Editing | |
|---|---|
| Ctrl+Z / Ctrl+Shift+Z | Undo / Redo |
| Ctrl+E | Reverse complement the selection |
| Ctrl+B | Set the origin at the cursor |
| Ctrl+K | Add a feature over the selection |

| The editor's own panels | |
|---|---|
| Ctrl+Shift+D | Digest panel: pick enzymes, see the fragments. A view tab brings the map back |
| Ctrl+Shift+P | The editor's PCR panel. Vettore's own PCR is **Actions → PCR…** |

## If something goes wrong

- **A file that cannot be opened**: a red message at the top right names the file, says what
  is wrong and what to try.
- **Inside a dialog**: a red box at the top of the result area says why there is no product,
  for example a primer that does not bind the template. Change the input it names and the
  result updates at once.
- **Things to check before ordering primers** appear in an amber box under the product.
- **The window stops responding and you close it**: Vettore cannot offer to save your work,
  so after a few seconds the system asks whether to quit and discard the unsaved changes.
  Answer No to keep the window open.
- All dialogs fit a 1366×768 screen. If a form is longer than the window (for example a gel
  with several lanes), scroll inside the dialog: the Cancel and Create buttons stay at the
  bottom.

The log is at `%USERPROFILE%\.vettore\vettore.log` on Windows (`~/.vettore/vettore.log` on
Linux). Send it along with the file that caused the problem.

Vettore is an independent program built on open-source libraries.
