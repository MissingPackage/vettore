# Changelog

The download links in the README always point to the latest release.

## 0.6.1

- The default primer Tm preset is named for its model: **Standard (SantaLucia & Hicks 2004)**.
- Release notes now come from the release tag, so every download page lists what is new.
- README rewritten around what the program does, with the full feature list.

## 0.6.0

**New**
- Actions → Anneal Oligos and Insert: clone gRNAs, shRNAs, linkers and short tags from two annealed oligos. Give it the insert and the vector, and it writes the oligo pair from the cut vector's own overhangs — any vector, any enzyme.
- Enzymes view: a Methylation column flags sites where Dam or Dcm methylation falls inside the recognition sequence and REBASE records the enzyme as blocked or impaired; enzymes that cut regardless stay quiet.

## 0.5.0

Automatic annotation.

**New**
- Features → Detect Common Features finds AAV2 ITRs (full and truncated), CMV, CAG, hSyn, EF-1α, CaMKIIα, GFAP and U6 promoters, WPRE, bGH and SV40 polyA, AmpR and KanR, pUC and f1 origins, EGFP, mCherry, tdTomato, loxP, lox2272, FRT, SpCas9, SaCas9, FLAG and HA tags — on both strands.
- Every catalogue entry is a subsequence of a public GenBank record and the match table shows which one, with its identity, so each annotation can be checked.

## 0.4.0

Clone verification by Sanger sequencing.

**New**
- Tools → Align to Reference DNA Sequence: `.ab1` reads are quality-trimmed (adjustable stringency), aligned to the open construct on either strand, and stacked under it with their chromatograms.
- A per-read table gives strand, region covered, identity and every difference as `position: reference → read`, with insertions and deletions grouped into single events; reads across the origin of a circular plasmid are handled.

## 0.3.0

**New**
- Golden Gate assembly with type IIS enzymes: junction overhangs checked pairwise, duplicate and palindromic overhangs flagged, parts with internal sites rejected with a domestication message.
- Simulated agarose gels: lanes from digests or uncut (supercoiled) plasmids, seven ladders from NEB, Thermo Fisher and Invitrogen, 0.5–3 % agarose.

## 0.2.0

Cloning simulation.

**New**
- Restriction cloning: insert a fragment in either orientation, or delete one and re-ligate; warnings for non-directional ligation, vector self-ligation, blunt ends and enzymes cutting more than once. Products keep the vector's numbering and the features of both parents.
- Gibson, NEBuilder HiFi and In-Fusion assembly of one or more inserts into a vector opened by digestion or inverse PCR, with the homology arms designed and each design checked by assembling it.
- PCR with tailed primers (including amplicons across the origin) and primer design to a target Tm.
- Site-directed mutagenesis (substitution, insertion, deletion) with back-to-back primers.
- Every product opens as a new file, with History recording its parents, parameters and primers; the primer list copies straight into an order form.
- User guide.

## 0.1.0

First release: the everyday editor.

**New**
- Opens `.dna`, GenBank (`.gb`, `.gbk`, `.ape`), FASTA and Sanger `.ab1` files, keeping features, feature colours, primers with their 5′ tails, and topology.
- Map, Sequence, Enzymes, Features, Primers and History views (Ctrl+1 … Ctrl+6), several files open as tabs, unsaved-changes marker and confirmation on quit.
- Sequence and annotation editing with undo/redo, reverse complement, topology change, new file from pasted sequence or FASTA.
- REBASE commercial enzymes, one per specificity, with sets for unique 6+, unique, unique & dual, 6+, type IIS and Golden Gate cutters; Enzymes view with cut positions and ends.
- Primer Tm on any selection, with the model and conditions always shown.
- Save as GenBank; export FASTA, and maps as SVG or PNG.

