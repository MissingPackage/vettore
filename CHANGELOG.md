# Changelog

The download links in the README always point to the latest release.

## Unreleased — opening a record from NCBI

- File → Open from NCBI… opens a public GenBank record by its accession number (for example MT612431.1) as a new, unsaved document; the welcome screen has the same link. It is the only thing in Vettore that uses the internet, only when you press Open, and it sends nothing but the accession number; no account is needed.
- Offline, the dialog says so plainly and stays open; nothing else in Vettore changes. A mistyped, unknown or protein accession, a genome too large to open (over 2 Mb, refused before downloading) and a network that blocks secure connections each get one sentence saying what to do. An accession typed without its version opens under the version NCBI sends.
- File → Preferences → Allow Open from NCBI (on by default) hides the menu item when unticked.

## Unreleased — a lighter, clearer frame

- The frame is lighter: toolbar tools that are on are white keys with a thin blue outline instead of solid squares, document tabs lose their folder borders, and the view tabs (Map, Sequence, Enzymes…) are one segmented control whose white marker slides to the open view. Nothing moved: menus, tabs, panes and shortcuts are where they were.
- The welcome screen shows each sample construct as a card with a small map of its features, its length and its GenBank accession; Open, New and Paste have icons.
- A document that has just opened — a cloning product, a pasted sequence — slides into the tab strip and its outline glows once, and creating a product says so in a notice. Copy buttons turn green with a check when the copy worked.
- Ctrl+/ (or "Search commands" in the view bar) finds any menu command by name: typing "gel" lists Simulate Agarose Gel… first, with the menu it lives in.
- In Simulate Agarose Gel the bands slide to their new size when you change a lane's enzymes, so you see which fragment split.
- File → Preferences → Appearance can make the menus, tabs and status bar dark, or follow the system; the map, sequence and dialogs stay light.
- Movement is short (under a quarter of a second), never touches the map, sequence or tables, and is switched off when the system asks for reduced motion.

## Unreleased — the lab primer library

- Primers → Primer Library… keeps the primers your lab already has: name, sequence, notes, where they are kept ("box 3 A4") and the date added. Search by name or sequence, add a primer, or import a spreadsheet — paste the columns copied from Excel, or open a .csv/.tsv; rows that are not DNA are listed and left out. Export the library as .csv.
- The library is one .csv file. Keep it in a shared folder (OneDrive, Dropbox, a lab drive) and your colleagues share it: Vettore reads it from there and writes back only what you changed: the header, your own columns and rows it cannot read stay as they are, and primers a colleague added while you had it open are kept. An Excel workbook or a file that is not a primer list is refused, never overwritten.
- Primers → Show Library Primers draws the library primers that bind the open document on its map, on either strand, tails and up to two mismatches allowed. These are a view, not part of the file: they are never written to a saved document.
- Design PCR Primers says when a library primer pair already amplifies the region ("You may already have KanR_F and KanR_R: product 1,086 bp"), and a designed pair can be Saved to Library (a sequence already there is refused with the name it has).

## Unreleased — choosing a diagnostic digest

- Tools → Find Diagnostic Digest… finds the single or double digest that tells the expected clone apart from the wrong ones on one gel, and lists the best ten with one sentence each ("Expected clone: 4.2 kb + 1.1 kb; parent vector: 4.4 kb"). Choosing one draws its gel next to your ladder; Send to Gel Simulation opens it with one lane per outcome, and Copy as Protocol copies enzymes, gel and expected sizes.
- A product of restriction cloning remembers its wrong clones (insert in the other orientation, vector closed without insert, parent vector, and the simulation's other products), and the digest planner starts from them.
- Enzymes that Dam or Dcm methylation blocks on the plasmids are left out and named, unless the DNA comes from a dam−/dcm− strain.
- File → Preferences → My enzymes keeps the enzymes in your freezer, each with the buffer you digest it in (empty until you fill it in; a list, or two columns enzyme and buffer, can be pasted from a spreadsheet); the digest planner tries them first.
- Find Diagnostic Digest shows each enzyme's buffer next to it, ranks double digests in one buffer before those that need two, and says what to do when a pair needs two different buffers; the copied protocol carries the buffers too.

## Unreleased — aligning several sequences to each other

- Tools → Align Multiple Sequences… has a second mode, **All to each other**: a progressive alignment (guide tree from shared words, most similar sequences first) of up to 50 sequences, with no reference privileged. Strand and origin follow the first circular sequence; each sequence is turned by its closest relative, and a fragment across a plasmid's origin lines up in one piece. One that shares nothing with the others is left out with a note.
- The result shows a consensus row on top of the stacked sequences, a conservation bar with the most variable stretches listed under it, each sequence's identity to the consensus, and the differences between any two sequences you pick, in the first one's numbering, ready to copy into a spreadsheet. Unrelated stretches (two inserts in one backbone) are named as not compared rather than listed as mutations, and N is never counted as a difference. The two are realigned to each other where needed, so gaps the other sequences required never show up as an insertion and a deletion between them.
- The alignment runs in the background with a progress bar and can be stopped; 10 plasmids of 15 kb take a few seconds.

## 0.3.1 — more promoters, with their species

- Detect Common Features knows the PGK (human and mouse), UbC, H1, 7SK and SFFV promoters and the full-length CAG promoter (CMV enhancer, chicken β-actin promoter and chimeric intron), each taken from a public GenBank record and cross-checked against a second, independent one.
- Promoters carry their species in the name — "U6 promoter (human)", "CaMKIIα promoter (mouse)", "PGK promoter (human)" and "(mouse)" — and the label check uses the same sequences, so a mouse PGK labelled as human, or a CMV labelled as PGK, is flagged.
- A CMV enhancer inside a detected full-length CAG is reported as part of the CAG, not as a separate element.

## 0.3.0 — files from vendors, MultiSite Gateway, a real desktop app

**Files from vendors and labs**
- Feature colours written in a note (`color: #63b349`, or a colour name) colour the map, so vendor files look as their authors drew them; saving writes each colour once.
- Feature types spelled the vendors' way (Promoter, PolyA_signal, Rep_origin, ITR, Transcription_t, Miscellaneous, any letter case) and the `regulatory` features of current NCBI records open as standard types, so they are coloured and grouped like any other file; the original type is kept and `regulatory` features are saved back in their own form.
- Detect Common Features flags labels that may be wrong: when most of a feature is a nearly exact match for a different element than its name says — a CMV promoter labelled "PGK promoter", found in a real vendor file — the row is marked in amber and a notice counts the labels to check.
- Detect Common Features counts an element inside a larger feature named after it as already annotated (a "CMV promoter" block over the CMV enhancer and promoter is not offered again).

**Gateway**
- MultiSite Gateway: an LR reaction with up to four entry clones chained through att1–att5, filled in from the open files when they fit together; the att sites of each file are shown under its picker. Checked against published MultiSite vectors and a published three-fragment clone.
- BP shows and checks which att pair of the donor vector the reaction uses.
- When inputs cannot recombine, the message lists the att sites of each file and names the missing junction.
- phiC31 attB/attP sites (fly landing-site vectors) are recognised as not Gateway sites, and said so.
- Create All (N): every product of a simulation opens in its own tab, each with its History entry.

**Desktop**
- Closing the window with unsaved work asks in Vettore's own dialog — Save All, Don't Save or Cancel — listing the documents; a frozen window can still be closed.
- Double-clicking a sequence file while Vettore is open opens it in the window you already have, not in a second copy.
- File → Preferences can make double-clicking `.dna`, `.gb`, `.fasta`, `.ab1` and the other sequence files open them in Vettore, for your user account only and without administrator rights; untick to undo.
- On Linux, File → Print… explains how to print the exported map instead of silently doing nothing.

## 0.2.0 — Files from colleagues, more cloning, primer design and alignment

**Files from colleagues**
- Opens EMBL (`.embl`), `.xdna` and `.gck` files and plain sequence text (`.txt`, `.seq`), as
  well as `.dna`, GenBank, FASTA and `.ab1`.
- Files are recognised by what is inside them: a GenBank file saved as `.txt`, a FASTA file
  named `.seq` or a file without an extension opens correctly.
- A file with several sequences opens one tab per sequence, and a message says how many.
- GenBank files from other programs open even with very long names, no topology, Windows line
  endings, accented characters or features drawn across the origin; features that cannot be
  placed on the sequence are listed by name.
- Error messages name the file, say in plain words what is wrong and what to try; for plain
  sequences they point to the line and column of the first character that is not a DNA letter.
- Save As and Export write GenBank, FASTA or `.dna` files; the format follows the extension you
  type, and `.gb` is added when you type none.
- Save only writes over a file when nothing in it would be lost. Otherwise (a `.dna`, EMBL,
  `.xdna`, `.gck` or `.ab1` file, one sequence out of several, features that could not be
  shown) it explains this once and saves a GenBank copy next to the original, which is left
  untouched.
- Files dropped on the window open from where they are on disk, so Save writes back to them.
- File → New File from Clipboard (Ctrl+Shift+V) opens a copied sequence, FASTA or GenBank/EMBL
  record; Edit → Copy as FASTA copies the selection with its positions for emails and order
  forms.
- The Open and Save dialogs start in the folder you used last.
- Sequences written in RNA letters open as DNA (U read as T), with a note.
- Saved GenBank files carry the date of the save and the record's division, and a long name is
  shortened at a word boundary with the full name kept in the definition line.
- New cloning products, pasted sequences and new files are marked unsaved (•) until you save
  them; closing a tab with unsaved work asks once: Save, Don't save or Cancel.

**Cloning simulation**
- Insert Multiple Fragments: a cut vector and two or more inserts in one ligation; every
  product the ends allow is listed with its order and orientation.
- Ligate Fragments: join fragments end to end into a linear molecule, or close them into a
  circle.
- Ligation warnings before you order anything: a vector that can close without an insert,
  blunt ends, PCR products without 5′ phosphate, fragments that only fit in another order.
- TA and GC cloning, and TOPO cloning (TA, blunt and directional, with the CACC check and an
  Add CACC button).
- Gateway BP and LR reactions, listing the clone and the by-product, with a warning when the
  product still carries ccdB; validated against published entry, destination and expression
  clones.
- Every cloning dialog starts from the files you have open, prefills regions from the editor
  selection or a feature you pick by name, says which input is still missing, and names the
  product after the method ("pUC19 + GFP (Gibson)").
- Gateway dialogs pick the entry clone, destination, attB insert and donor by their att sites.
- PCR: a reverse primer copied from the top strand is answered with the reverse complement to
  order.
- History keeps a copy of every input of a cloning step: open it again as it was when the step
  ran, or copy the step's parameters and primers for your notebook.

**Primers and verification**
- Primers → Design PCR Primers: primer3 proposes the five best pairs around a region or feature,
  with Tm, GC, product size and hairpin/dimer ΔG warnings; add a pair to the plasmid or copy it
  for ordering.
- Tools → Align Multiple Sequences: up to 10 sequences against a reference, on either strand,
  with coverage, identity and every difference in the reference's numbering; sequences that
  share only a chance match are left out, with the reason.
- Sanger verification: a one-line verdict that counts a difference once however many reads
  show it, the read quality behind every difference, low-quality calls flagged, and a clear
  message when a read comes from another plasmid.
- Tools → AAV Packaging Check: the cargo between the ITRs and including them, against the
  single-stranded and self-complementary limits, with the sources of the numbers.

**Getting started and everyday comfort**
- Welcome screen with three sample constructs (pUC19 and two AAV transfer plasmids, public
  GenBank records), recent files, paste and drop; the samples are also under File → Open
  Sample Construct.
- Help menu: Keyboard Shortcuts, Getting Started, About Vettore — also reachable before any
  file is open.
- File → Preferences: default enzyme set, gel ladder and agarose percentage, minimum identity
  for feature detection.
- Ctrl+Shift+G opens the gel simulation; menus show their shortcuts as Ctrl+….
- Detect Common Features does not offer again an element already annotated under another name.
- The gel draws a lane whose enzyme does not cut as uncut, with a note, instead of failing.
- Every error is one plain sentence that says what to do next.
- A calmer, consistent look: one font, one accent colour, readable secondary text, visible
  keyboard focus; every dialog fits a 1366×768 screen with its buttons in view.
- New app icon.

**Fixes**
- Align to Reference DNA Sequence (Sanger) and Detect Common Features work in the packaged app:
  in the 0.1.0 download both stopped with an error because a data file was missing from the
  build. The build check now exercises both on Windows and Linux.

## 0.1.0

First public release.

**Sequences and maps**
- Opens `.dna`, GenBank (`.gb`, `.gbk`, `.ape`), FASTA and Sanger `.ab1` files, keeping features, feature colours, primers with their 5′ tails, and topology.
- Map, Sequence, Enzymes, Features, Primers and History views (Ctrl+1 … Ctrl+6); several files open as tabs, with an unsaved-changes marker and a confirmation on quit.
- Sequence and annotation editing with undo/redo, reverse complement, topology change, find (DNA or protein, with mismatches) and new files from pasted sequence or FASTA.
- Saves GenBank; exports FASTA, and maps as SVG or PNG for figures.

**Enzymes**
- REBASE commercial enzymes, one per specificity, with sets for unique 6+, unique, unique & dual, 6+, type IIS and Golden Gate cutters.
- Enzymes view with cut positions and ends per site, and a Methylation column that flags sites where Dam or Dcm methylation falls inside a recognition sequence REBASE records as blocked or impaired.

**Cloning simulation** — each with a live preview, warnings and the primers to order
- Restriction cloning: insert a fragment in either orientation or delete one; checks for directionality, vector self-ligation, blunt ends and multiple cutters. Products keep the vector's numbering and the features of both parents.
- Gibson, NEBuilder HiFi and In-Fusion assembly of one or more inserts into a vector opened by digestion or inverse PCR, with homology arms designed and each design checked by assembling it.
- Golden Gate assembly with junction-overhang checks (duplicates, palindromes, domestication).
- Annealed-oligo cloning for gRNAs, shRNAs, linkers and tags, with the oligo pair written from the cut vector's own overhangs.
- PCR with tailed primers (including amplicons across the origin), primer design to a target Tm, and site-directed mutagenesis with back-to-back primers.
- Every product opens as a new file with a History entry recording its parents, parameters and primers.

**Verification**
- Sanger reads aligned to the expected construct on either strand, quality-trimmed, stacked under it with their chromatograms, with a per-read list of every difference in plasmid coordinates.
- Simulated agarose gels with seven NEB, Thermo Fisher and Invitrogen ladders.
- Automatic annotation of AAV2 ITRs, common promoters, WPRE, polyA signals, markers, origins, fluorescent proteins, lox/FRT sites, Cas9 and epitope tags, each traceable to a public GenBank record.

**Primers**
- Tm on any selection, always shown with its model and conditions: Standard (SantaLucia & Hicks 2004), and presets matching NEB's calculator for Q5 and Taq.
