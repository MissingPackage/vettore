# Changelog

The download links in the README always point to the latest release.

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
