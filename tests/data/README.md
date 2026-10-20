# Test fixtures

The `*.dna` files and `3730.ab1` (ABI Sanger trace) are copied unchanged from the Biopython
test suite (`Tests/`, https://github.com/biopython/biopython),
distributed under the Biopython License Agreement / BSD 3-Clause License.

## File-exchange fixtures (`exchange-*`)

Downloaded on 2026-10-20. INSDC records (GenBank/ENA) are freely available public data.

| File | Source | Changes |
|---|---|---|
| `exchange-pBR322.embl` | ENA record J01749 (pBR322), EMBL flat file from `https://www.ebi.ac.uk/ena/browser/api/embl/J01749` | none |
| `exchange-plasmids.fasta` | NCBI nuccore L09137.2, J01749.1, X06403.1 (pUC19c, pBR322, pACYC184), E-utilities `efetch rettype=fasta` | none |
| `exchange-sample-a.xdna`, `exchange-sample-b.xdna` | Biopython `Tests/Xdna/sample-a.xdna`, `sample-b.xdna` (BSD) | none |
| `exchange-artificial.gck` | Biopython `Tests/Gck/artificial.gck` (BSD) | none |
| `exchange-pUC19-numbered.seq` | NCBI L09137.2 GenBank record | only the numbered ORIGIN block (plain-sequence text with position numbers) |
| `exchange-pUC19-raw.txt` | NCBI L09137.2 FASTA | header line removed (FASTA without header) |
| `exchange-quirks-header.gb` | NCBI J01749.1 GenBank record | hand-edited: LOCUS line replaced by a long name without topology or date; one `/note="non-ASCII test: 37 °C, 5 µl"` added to the tet CDS; final `//` removed; Windows (CRLF) line endings; encoded as cp1252 |
| `PZ629150.1-unannotated.gb` | `PZ629150.1.gb` (NCBI PZ629150.1) | every feature except `source` removed with Biopython, so Detect Common Features has everything to add; sequence and header unchanged |
| `exchange-quirks-features.gb` | NCBI J01749.1 GenBank record | hand-edited: eight `misc_feature` entries labelled `test …` added before the last feature: a region across the origin written as `join(4300..4361,1..50)` and as `4300..50`, a site between two bases (`100^101`), regions beyond (`4400..4500`) and past (`4300..4400`) the end, a location in another record (`J01636.1:1..10`), and two features with the same name and unusual qualifiers. The `test …` features mark positions only and carry no biological meaning |

`cloning-gateway-entry.gb` (GenBank LT726670.1), `cloning-gateway-destination.gb` (AY568055.1)
and `cloning-gateway-expression.gb` (LT725632.1) are public NCBI records fetched with E-utilities
and saved unchanged; they validate the Gateway BP and LR reactions (`tests/test_cloning_gateway.py`,
which cites Hilson et al. 2004, Genome Res 14:2176).

`analysis-NC_001401.2.gb` (adeno-associated virus 2, complete genome) and `analysis-AY530579.1.gb`
(AAV9 hu.14 VP1 gene) are NCBI nuccore records downloaded unchanged with E-utilities efetch
(rettype=gbwithparts) on 2026-10-20; used by `tests/test_alignment.py` and
`frontend/e2e/alignment.spec.ts`.

## Vendor-notation fixtures (`fidelity-*`)

Written for Vettore on 2026-10-20 to reproduce how vector vendors write GenBank files: the
feature colour in a note (`/note="color: #63b349; direction: RIGHT"`, also `color: Blue` and
`color: None`) and non-standard feature keys (`ITR`, `Promoter`, `PROMOTER`, `Miscellaneous`,
`Transcription_t`, `Rep_origin`, `PolyA_signal`, an unknown `sgRNA`). No vendor file was
copied: sequences come from public records, features were written by hand.

| File | Sequence | Features |
|---|---|---|
| `fidelity-notes.gb` | NCBI PZ629150.1, unchanged | at the published record's own coordinates (ITRs 1..131 and 4343..4483, cTNT promoter, Kozak, SaCas9, U6 promoter, sgRNA spacer and scaffold 4228..4323, the TTTTT Pol III terminator 4324..4328, f1 ori, AmpR promoter, AmpR, ori), rewritten in vendor notation |
| `fidelity-block-cmv.gb` | catalogue entries joined in this order: CMV enhancer + CMV promoter (PV031350.1:96..679), EGFP (U55762.1:679..1398), SV40 polyA (MH282430.1:2999..3230), ori (PZ629150.1:6431..7019), AmpR promoter + AmpR (PV031350.1:2514..3479) | one `Promoter` block “CMV promoter” over the CMV enhancer and promoter (1..584), one feature per other element |
| `fidelity-mislabelled-cmv.gb` | the same as `fidelity-block-cmv.gb` | the same, with the block labelled “PGK promoter”: a wrong label, as found in a vendor's AAV vector whose own `/official_designation` says CMV |

MultiSite Gateway fixtures, public NCBI records fetched unchanged with E-utilities on
2026-10-20 (`tests/test_cloning_multisite.py`, `tests/test_cloning_gateway.py`):
`cloning-multisite-5prime.gb` (LT726877.1, attL4/attR1), `cloning-multisite-3prime.gb`
(LT726875.1, attR2/attL3), `cloning-multisite-destination.gb` (LT726807.1, pdcDNA-R4R3-dhCMV,
attR4/attR3), `cloning-multisite-destination-r4r2.gb` (AP018631.1, R4pMpGWB495, attR4/attR2),
`cloning-multisite-entry-l1r5.gb` (PP281278.1, pGN1101CK), `cloning-multisite-entry-l5l4.gb`
(PP281279.1, pGN1102CK), `cloning-multisite-entry-l4r1.gb` (LC477681.1, pCas9_P4P1r),
`cloning-multisite-donor.gb` (HM623915.1, pA-hyg-OSCAR, attP4/attP1r),
`cloning-multisite-expression.gb` (KC800697.1, a published three-fragment product) and
`cloning-phic31-attb.gb` (EF362409.1, pUASTattB, phiC31 attB; Bischof et al. 2007).
