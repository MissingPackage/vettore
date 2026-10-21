# Vettore — offline plasmid editor and cloning simulator

## Overview

Vettore is a desktop application for Windows and Linux that follows the everyday molecular
cloning workflow — map / sequence / enzymes / features / primers / history views, cloning
simulations, primer tools, gel simulation, Sanger verification — with the layout and
terminology bench scientists already use, and without a licence fee. Cloud tools and
Docker-based tools are out of scope: everything runs locally.

Independence rule: functionality is built from published methods, public documentation and
open libraries only. No third-party binaries are decompiled or inspected, and no third-party
names, icons, artwork or proprietary databases are reused.

Correctness rule: every molecular operation (digest, ligation, assembly, PCR, Tm, alignment)
runs on a published, tested engine (pydna, Biopython, primer3) and is covered by tests against
known reference results. The UI never computes biology on its own.

## Architecture (decided)

| Layer | Choice | Why |
|---|---|---|
| Editor UI | React + TeselaGen Open Vector Editor (`@teselagen/ove` 0.8.42, MIT) | Circular/linear map, sequence rows, cut sites, features, primers, ORFs, digest, alignments with chromatograms already exist and are maintained |
| App shell | Own React shell around OVE: document tabs, view tabs, menus, dialogs | Familiar layout; OVE stays an embedded component |
| Biology engine | Python: pydna (BSD), Biopython (BSD-like), primer3-py (GPL-2.0+) | Peer-reviewed, tested simulation of cloning operations |
| Bridge | Local FastAPI server bound to 127.0.0.1 with a per-session random token; UI talks HTTP JSON | One code path for desktop and browser/test mode |
| Desktop window | pywebview (Edge WebView2 on Windows, Qt WebEngine on Linux); fallback: open system browser | No Electron, single Python process |
| Packaging | PyInstaller one-folder build, zipped; built by GitHub Actions on `windows-latest` and `ubuntu-latest` | Unsigned exe acceptable (personal PC) |
| Native file format | Read `.dna`, `.gb/.gbk/.ape`, `.embl`, `.fasta`, `.xdna`, `.gck`, `.ab1` and plain text, detected from the content; save GenBank (`.gb`, default), FASTA or `.dna` (written with sgffp, round trip verified with Biopython and sgffp) | Every sequence tool reads GenBank; it is the lossless open standard |
| `.dna` reading | Biopython (sequence, features, primer sites after hybridization filtering) + sgffp 0.22 (MIT) for feature colours, full primer sequences, description | Biopython alone drops colours and primer 5′ tails |
| Enzyme data | Nonredundant commercial set generated from REBASE via Biopython (`vettore.enzymes` → `frontend/src/data/enzymes.json`), fed to OVE | OVE's default list contains isoschizomers with uncharacterised cuts (seen: Atu1II GGATCC shown as blunt); OVE's cut engine is cross-checked against Biopython in tests |
| Tm | Biopython `Tm_NN` with named presets: Standard (SantaLucia & Hicks 2004, default), NEB Q5, NEB Taq; all shown in the tooltip | Values differ by several °C between tools; the method is always shown |

## Slices

Each slice ends with: tests green (unit + e2e), Windows and Linux builds produced by CI, a
release zip attached to a GitHub Release, and a short "what's new" note.

- **S0 — Spike**: prove OVE in pywebview + Python bridge + Windows build. Open a `.dna` file and see its map.
- **S1 — Daily driver**: open/view/edit/save, the six views, enzymes, features, primers, find, export images.
- **S2 — Cloning simulations**: restriction & insertion, Gibson/HiFi/In-Fusion, PCR, mutagenesis, Golden Gate, Gateway, TA/TOPO, ligation; history.
- **S3 — Verification tools**: simulate agarose gel, align Sanger reads (`.ab1`) to reference, multiple alignment, detect common features, primer design.
- **S4 — Polish**: history graph, print/PDF, preferences, keyboard-shortcut parity, file associations, AAV helpers (ITR-to-ITR cargo size).
- **S5 — Screening and lab memory**: S5.1 diagnostic digest planner (started); next S5.3 lab primer library, S5.2 Sanger contig, S5.4 human and mouse elements in the catalogue.

## User stories

### S0 — Spike

**S0.1** As the developer I want the app to start as a native window on Windows and Linux so that the stack is proven before building features.
- [ ] `vettore` (Linux) and `Vettore.exe` (Windows, built in CI) open a window titled "Vettore".
- [ ] File → Open on `tests/data/pFA-KanMX4.dna` shows a circular map with 9 features and length 3941 bp.
- [ ] The server rejects API requests without the session token (HTTP 403).
- [ ] CI uploads a Windows zip and a Linux zip as workflow artifacts.

### S1 — Daily driver

**S1.1 Open files.** As a bench scientist I want to open my existing `.dna` and GenBank files so that I keep working on them.
- [ ] Opens `.dna`, `.gb`, `.gbk`, `.ape`, `.fa/.fasta`, `.ab1` (as read-only sequence); features, primers (from `.dna` primer records), topology and name are preserved.
- [ ] Several files open as tabs; drag-and-drop of a file onto the window opens it.
- [ ] Recent files list (last 10) in File menu.
- [ ] Opening a corrupt file shows an error dialog naming the file; the app keeps running.

**S1.2 Views.** As a bench scientist I want the map, sequence, enzyme, feature and primer views I am used to, so that I find everything where I expect it.
- [ ] View tabs at the bottom: Map, Sequence, Enzymes, Features, Primers, History.
- [ ] Map: circular for circular molecules, linear otherwise; toggle available.
- [ ] Sequence: double-strand rows with features, primers, cut sites, translations.
- [ ] Enzymes: table of enzymes in the current set with number of sites and positions.
- [ ] Features / Primers: sortable tables; double-click opens the edit dialog.
- [ ] Selecting any range shows in the status bar: start..end, length, %GC, and Tm when length ≤ 60 nt.

**S1.3 Edit sequence.** As a bench scientist I want to edit the sequence and annotations so that I can design constructs by hand.
- [ ] Type/paste/delete bases at the cursor; features and primers outside the edit keep their positions; overlapping ones are adjusted as in GenBank semantics.
- [ ] Reverse-complement whole sequence or selection; change topology (circular/linear); set origin.
- [ ] Add/edit/delete feature (name, type, strand, color, location); add/edit/delete primer.
- [ ] Undo/redo (Ctrl+Z / Ctrl+Y or Ctrl+Shift+Z) for every edit.
- [ ] Unsaved changes mark the tab with "•"; closing asks to save.

**S1.4 Save and export.** As a bench scientist I want to save and share my work so that collaborators can open it in whatever tool they use.
- [ ] Save / Save As writes GenBank; re-opening the file round-trips sequence, features (type, name, strand, location, color) and primers.
- [ ] Export FASTA; copy selection as plain sequence; copy reverse complement.
- [ ] Export the map as SVG and PNG.

**S1.5 Enzymes.** As a bench scientist I want to choose which enzymes are shown so that the map shows the cutters I care about.
- [ ] Enzyme set chooser with presets: Unique cutters, Unique 6+ cutters, Unique & dual cutters, Commercially available, plus a custom list.
- [ ] Search enzyme by name; show recognition site and cut position.
- [ ] Enzyme data comes from REBASE via Biopython.

**S1.6 Find.** As a bench scientist I want to find sequences so that I can locate sites and motifs.
- [ ] Ctrl+F finds DNA (both strands) and protein (in all 6 frames), with optional mismatches (0–3); results are highlighted and navigable.
- [ ] Show ORFs with a minimum length setting (OVE's default).

### S2 — Cloning simulations

Each action opens a dialog with the inputs on the left and a live preview of the product on the right,
creates a new document for the product, and records the operation in the History view.

**S2.1 Restriction & insertion cloning.**
- [ ] Choose vector + 1–2 enzymes (or one blunt site), choose insert fragment from another document by enzymes or by PCR product; orientation selector; result matches pydna's digest + ligate result base-for-base.
- [ ] Warns when ends are incompatible, when an enzyme also cuts inside the insert, and when the site is Dam/Dcm blocked [ASSUMED: methylation warnings in S2, not S1].

**S2.2 Gibson / NEBuilder HiFi / In-Fusion.**
- [ ] Select an opened vector (linearized by enzymes or by PCR) and 1–6 inserts in order; the tool designs primers with overlaps (default 20–30 nt for Gibson/HiFi, 15 nt for In-Fusion) and shows Tm of annealing parts.
- [ ] Product sequence equals pydna's assembly of the designed PCR products; primers are added to the product's primer list.

**S2.3 PCR.**
- [ ] Choose template + forward/reverse primers (from the primer list or typed, with 5′ tails); product shown and created; primers with no binding site produce an error message.

**S2.4 Mutagenesis.** Insertion, deletion, replacement via primers (Q5 SDM-style); output primers + product.

**S2.5 Golden Gate.** Type IIS enzyme (BsaI, BsmBI/Esp3I, BbsI, SapI), parts in order, overhang validation (duplicate/palindromic overhangs flagged).

**S2.6 Gateway.** BP and LR reactions (attB/attP/attL/attR) producing entry/expression clones.

**S2.7 TA / TOPO and Ligate fragments.** Simple ligation of selected fragments; TA with A-overhang products.

**S2.8 History.** Every product stores its parents and the operation; History view lists the steps; parents can be reopened.

### S3 — Verification tools

**S3.1 Simulate agarose gel.** Lanes = digests of opened documents with chosen enzymes or uncut; ladders: NEB 1 kb, NEB 1 kb Plus, NEB 100 bp, GeneRuler 1 kb, GeneRuler 1 kb Plus, GeneRuler 100 bp Plus [ASSUMED: these are the ladders the user uses]; gel % 0.7–2.0; export PNG.

**S3.2 Align to reference (Sanger).** Import one or more `.ab1` reads, trim low-quality ends, align to the current sequence, show chromatograms, mismatches/indels highlighted, coverage summary.

**S3.3 Multiple alignment** of DNA sequences.

**S3.4 Detect common features** from an offline curated database (AAV-relevant set first: ITRs, CMV/CAG/hSyn/EF1a/CaMKIIa/GFAP, WPRE, polyA signals, AmpR/KanR, ori, f1 ori, fluorescent proteins, tags, lox/FRT, U6, Cas9).

**S3.5 Primer design** (primer3): amplify a selected region with Tm, length and product-size constraints; results ranked, added to the primer list on click.

### S4 — Polish

Print / PDF of map and sequence; preferences (Tm method, default enzyme set, colors); keyboard
shortcuts for the common actions; Windows file association for `.dna`/`.gb`;
AAV helper: ITR-to-ITR cargo size with a warning above 4.7 kb.

### S5 — Screening and lab memory

**S5.1 Diagnostic digest planner.** As a bench scientist I want Vettore to pick the restriction
digest that tells my correct clone apart from the usual wrong ones, so that one gel answers
"which colonies are right?".

- Tools → Find Diagnostic Digest… opens on the active document as the expected clone. A product
  of restriction cloning brings its wrong clones (insert in the other orientation, vector closed
  without insert, parent vector) and the simulation's other products; any open file can be
  added, any entry removed.
- Candidates: single enzymes and pairs from "My enzymes" (a preference, empty by default) or,
  while it is empty, the nonredundant commercial set. Enzymes that REBASE records as blocked by
  Dam/Dcm methylation at one of their sites on these molecules are left out, unless the DNA is
  from a dam−/dcm− strain.
- A digest is diagnostic when every wrong clone has a band that runs apart from every band of
  the expected clone (or the reverse) by at least the gel model's resolution (`gel.RESOLUTION`,
  the distance at which the gel simulator merges bands). Bands under 150 bp or above the
  resolving range of the agarose percentage do not count.
- Ranking, deterministic: well-separated digests (every wrong clone at least twice the
  resolution away, every lane cut) before digests at the limit of the gel; then fewer enzymes;
  one buffer before two (below); an expected lane of 2–5 visible bands; larger separation; enzymes most labs have. Identical
  band patterns are one result. The top 10 are listed with one sentence each.
- Buffers [DECIDED 2026-10-21: the lab's own, not a vendor table]: each enzyme of "My enzymes"
  carries the buffer the lab digests it in. Results show it next to the enzyme, pairs in one
  buffer rank before pairs in two (right after "fewer enzymes"), and a pair in two buffers is
  flagged with what to do. Enzymes without a buffer never count as different.
- A result draws its gel next to the ladder; Send to Gel Simulation and Copy as Protocol.
  When nothing works the dialog says why and offers the commercial set.
- Tested on pBR322 digests worked out by hand from the published record, sizes cross-checked
  with the gel simulator's pydna digest.

### S6.3 — Open from NCBI (online, explicit)

As a bench scientist I want to open a GenBank record by its accession number, so that I do not
have to download and save it first.

- File → Open from NCBI… (and a link on the welcome screen) asks for an accession, fetches the
  GenBank record (NCBI E-utilities `efetch`, `gbwithparts`) and opens it as an unsaved document
  whose History starts with "Opened from NCBI". Only the accession is sent; no account.
- The accession is checked before anything is sent (protein accessions refused, an accession
  picked out of pasted text); offline, an uncheckable certificate, an unknown accession, NCBI
  busy, a master record without sequence and a record over 2 Mb (read from the LOCUS line,
  before downloading) each get one sentence saying what to do. The tab takes the record's
  accession.version. Nothing else uses the
  network (`src/vettore/ncbi.py` is the only module that does).
- Preference "Allow Open from NCBI", on by default; off hides the menu item and the link.
- Not included: uploading sequences anywhere (no BLAST).

## Read first

| What | Where | Why |
|---|---|---|
| Module docstrings | `src/vettore/*.py`, `src/vettore/cloning/*.py` | Each engine states its model, its sources and its limits |
| Data builders | `scripts/build_feature_database.py`, `scripts/build_methylation_table.py` | Where the bundled NCBI and REBASE data come from |
| OVE README | `frontend/node_modules/@teselagen/ove/README.md` | Editor props, data model, alignment/chromatogram model |
| Test fixtures | `tests/data/` | `.dna` and `.ab1` samples from Biopython's test suite (BSD) |

## Non-goals

- No cloud sync, accounts, telemetry or auto-update.
- No protein structure, CRISPR guide scoring, or codon optimisation.
- No third-party names, logos, icons or proprietary databases.
- No macOS build (not requested) — the code should not block one later.

## Hard constraints

- Runs fully offline on Windows 10/11 x64 and Linux x64 (the one exception, S6.3 File → Open from NCBI…, downloads a public record only when the user asks and sends only its accession); no Docker, no admin rights to install (unzip and run).
- Unsigned executable is acceptable (personal PC).
- Repository: `MissingPackage/vettore`. All repo content in English.
- Deadline: S1 ships as soon as it is green, without waiting for S2 — the editor is useful on its own.
- Biology results must match the reference engine; any operation without tests is not shipped.

## Status (2026-10-20) — v0.3.1

| Slice | State | Missing |
|---|---|---|
| S0 Spike | done | — |
| S1 Daily driver | done (S1.1–S1.6), plus more formats (EMBL, `.xdna`, `.gck`, plain text), content-based format detection, `.dna` writing, clipboard exchange, vendor conventions (colour notes, non-standard and `regulatory` feature types) | — |
| S2 Cloning | done: S2.1–S2.8 — restriction (one or several inserts), linear and circular ligation, Gibson/HiFi/In-Fusion, PCR, mutagenesis, Golden Gate, Gateway BP/LR including MultiSite LR with up to four entry clones (checked against published MultiSite vectors and a published three-fragment clone), TA/GC/TOPO, annealed oligos; History reopens every input; Create All | — |
| S3 Verification | done: S3.1–S3.5 — gel, Sanger alignment with verdict and read quality, star alignment of several sequences, common features with label checks, primer3 pair design | progressive multiple alignment (sequences are aligned to the reference, not to each other) |
| S4 Polish | mostly: preferences, shortcut list, sample constructs, help, AAV packaging check, consistent look, in-app quit dialog, one window for double-clicked files, per-user file associations (Windows and Linux) | history graph; printing in the Linux window |

Tests: 609 pytest (plus 4 opt-in Qt window tests), 110 vitest, 130 Playwright end-to-end, plus a packaging smoke test that runs the packaged binary on both
platforms in CI (file round trips in every format, cloning, Tm, primer design, alignment,
Sanger reads, sample constructs, feature detection, file associations in dry run, and a second
start handing its file to the running copy).

## Resolved

- UI language: English, the language of the field's tools and protocols.

## Open questions

1. (resolved — see above)
2. Resolved in 0.2.0: GenBank stays the default and recommended exchange format; `.dna` is written on Save As/Export only, and Save never overwrites a `.dna` it opened. The round trip is verified with Biopython and sgffp; whether other programs open Vettore-written `.dna` files has not been checked yet.
3. [ASSUMED: S2 priority order = restriction cloning, Gibson/HiFi, PCR, mutagenesis, Golden Gate, Gateway, In-Fusion, TA/TOPO] — which methods are used most? (All but Gateway and TA/TOPO are done.)
4. (resolved — Sanger verification shipped)
5. [ASSUMED: ladders listed in S3.1] — which ladder brands are used in the lab?
6. Resolved 2026-10-21: online import from NCBI by accession, on by default as an explicit menu action (S6.3, File → Open from NCBI…); only the accession is sent, nothing is uploaded. A preference hides it.
7. [ASSUMED: delivery = the zip from the GitHub Releases page] — fine?
8. [ASSUMED: name "Vettore"] — rename freely.
9. Licence: primer3-py (GPL-2.0-or-later) is bundled in-process, so a distributed build must be GPL-3.0-or-later (Apache-2.0 components rule out GPLv2). Resolved: the repository is GPL-3.0-or-later.
