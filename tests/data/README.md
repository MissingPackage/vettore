# Test fixtures

The `*.dna` files and `3730.ab1` (ABI Sanger trace) are copied unchanged from the Biopython
test suite (`Tests/`, https://github.com/biopython/biopython),
distributed under the Biopython License Agreement / BSD 3-Clause License. So is
`msa-opuntia.aln` (Biopython `Tests/Clustalw/opuntia.aln`, unchanged): seven *Opuntia*
chloroplast sequences (GenBank AF191658-AF191665) and their ClustalW 1.81 alignment, the
reference alignment of `tests/test_msa.py`.

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

## Catalogue records (`features-*`)

Public NCBI nuccore records downloaded unchanged with E-utilities efetch (rettype=fasta) on
2026-10-20, used by `tests/test_features.py`. Two per element: the record a catalogue entry
was taken from, so the coordinates its `source` cites can be checked without the network, and
an independent record — another laboratory, another vector — that carries the same element, so
the entry can be checked against a copy that is not its own source.

The records the entries come from:

| File | Record | Element and the coordinates the entry cites |
|---|---|---|
| `features-OQ427862.1.fasta` | OQ427862.1, vector pLKO.1-ZsGreen1 (Tokuma, direct submission) | PGK promoter (human), 1..511 |
| `features-MH325104.1.fasta` | MH325104.1, vector pLVX.TRE3G.eGFP (Berman et al., J Am Chem Soc 2018;140:18093) | PGK promoter (mouse), 3321..3820 |
| `features-OL870463.1.fasta` | OL870463.1, vector pLenti-hACE2 (Katzman et al., Viruses 2021;14:11) | UbC promoter (human), 1096..2299 |
| `features-AY640634.1.fasta` | AY640634.1, siRNA vector pAAV9(5)-hSyn-DsRed2N1-EGFP-siRNA (Michel et al., Biochem Biophys Res Commun 2005;326:307) | H1 promoter (human), 1543..1758 |
| `features-PY046810.1.fasta` | PY046810.1, sequence 103 from patent US 12370253 B2 (Pauza et al.) | 7SK promoter (human), 1..243 |
| `features-KJ697750.1.fasta` | KJ697750.1, cloning vector pRRL-sffv-eGFP-cmv-hsCDNF (Cordero Llana, direct submission) | SFFV promoter, 6155..6637 |
| `features-ON012604.1.fasta` | ON012604.1, targeting vector AAVS1 CAG iRFP720 (Feher et al., Sci Rep 2022;12:9874) | CAG promoter (full length), 4799..6471 |

The independent records:

| File | Record | Element it cross-checks |
|---|---|---|
| `features-MT921016.1.fasta` | MT921016.1, vector pTTC_FusionRed_eGFP_3UTR (Dolskiy et al., direct submission) | PGK promoter (human), annotated there at 3754..4264 |
| `features-MT891341.1.fasta` | MT891341.1, cloning vector pT-PGK-mCherry (Frei et al., Nat Commun 2020;11:4641) | PGK promoter (mouse), annotated there at 345..848 |
| `features-MN389527.1.fasta` | MN389527.1, cloning vector pMG_DV (Stach et al., ACS Synth Biol 2019) | UbC promoter (human), annotated there at 583..1799 |
| `features-FJ687158.1.fasta` | FJ687158.1, RNAi expression vector pUC18-H1 (Mathupala et al., Neurosurgery 2004;55:1410) | H1 promoter (human), annotated there at complement(414..633) |
| `features-MQ296444.1.fasta` | MQ296444.1, sequence 123 from patent WO 2021216503 (Kmiec et al.) | 7SK promoter (human), the record's own 7SK promoter sequence |
| `features-MP040844.1.fasta` | MP040844.1, sequence 35 from patent WO 2018206791 (Duchateau et al.) | SFFV promoter, the record's own SFFV promoter sequence |
| `features-HG530138.1.fasta` | HG530138.1, cloning vector pTAL7B (Frank, Skryabin & Greber) | CAG promoter (full length), annotated there at 3272..4948 |

### The elements S5.4 added

Downloaded unchanged with E-utilities `efetch rettype=fasta` on 2026-10-21, two records per
element: the one the entry cites, and an independent one — another laboratory or another
applicant — that carries the same element. `tests/test_features.py` checks the cited
coordinates against the first and finds the element in the second at ≥ 98 % identity.

Where a row says "over the N bp the two share", the two records agree on every one of those
bases but the independent record does not carry the whole entry: a patent record's cloning
site, or a slightly differently bounded annotation, falls outside. Several rows cite a patent
sequence whose only annotation is the `source` feature's `/note`, which is the applicant's
word for what the sequence is; those entries rest on the two applicants agreeing, not on a
curated feature table.

| Element (entry) | The record it cites | The independent record |
|---|---|---|
| EFS promoter (human), `LP754738.1:1..252` | LP754738.1, sequence 14 from WO2018026976 (Editas Medicine) | MP933339.1, sequence 22 from EP3789521 (Shenzhen Istirbio), 99.6 % over the 233 bp the two share |
| CBh promoter, `PY026696.1:1..796` | PY026696.1, sequence 25 from US 12359219 B2 (Gray & Sinnett, University of Texas) | QQ167606.1, JP 2024542950-A/60 (University of Pennsylvania), 99.9 % |
| SV40 early promoter, `X85976.1:2537..2884` | X85976.1, plasmid pSG424 (Sadowski & Ptashne, Nucleic Acids Res 1989;17:7539) | M83237.1, vector RSV.5(neo) (Long et al., Hum Immunol 1991;31:229), 100 % over the 344 bp the two share |
| β-actin promoter (human), `JN964682.1:20222..20741` | JN964682.1, targeted allele Vmp1:tm1e(KOMP)Wtsi (Skarnes et al., Nature 2011;474:337) | AY952326.1, β-actin GFP expression vector (Costa et al., Nat Methods 2005;2:259), 99.8 % |
| HSV TK promoter, `X64410.1:429..598` | X64410.1, plasmid pBLCAT2 (Luckow & Schütz, Nucleic Acids Res 1987;15:5490) | AJ277960.1, cloning vector pMAR (Ciana et al., Mol Endocrinol 2001;15:1104), 100 % over the 169 bp the two share |
| RSV LTR promoter (182 bp of the 528 bp LTR), `LT727488.1:complement(3997..4178)` | LT727488.1, vector pGal4-p65 (De Schamphelaire et al., direct submission) | M83237.1, vector RSV.5(neo), whose 5′ LTR holds it, 99.5 % |
| TRE promoter (tight), `MT891336.1:554..868` | MT891336.1, vector pGLM92 (Frei et al., Nat Commun 2020;11:4641) | OQ427855.1, vector Puro-iDEST-3×FLAG-hWFS1-WT (Tokuma, direct submission), 100 % |
| Cre, `MK854762.1:215..1246` | MK854762.1, shuttle vector pBFS109 (Sakanaka et al., Sci Adv 2019;5:eaaw7696) | AB542060.1, Cre expression vector pTN75 (Nishiyama), 100 % |
| FLPe, `GU253315.1:2581..3879` | GU253315.1, vector pLV.FLPe (Gonçalves et al., PLoS One 2010;5:e10954) | KM521212.1, cloning vector pECA102 (Cameron & Collins), 99.8 % over the 1272 bp the two share |
| FLPo, `KX151728.1:1842..3140` | KX151728.1, vector DBH-FLPo (Sun & Ray, PLoS One 2016;11:e0159474) | MH493814.1, vector pDC50/pT3TS-FLPo (Grajevskaja et al., PLoS One 2018;13:e0197293), 100 % over the 1272 bp the two share |
| loxN, `OR242273.1:3825..3858` | OR242273.1, vector AGG1741 (Carabajal Paladino et al., Front Bioeng Biotechnol 2023;11:1254863) | KP686439.1, vector pC-(lox2-attB2-SA-Hsp70)₃ (Diao et al., Cell Rep 2015;10:1410), 100 % |
| F3, `LC086648.1:2498..2531` | LC086648.1, BPATrap vector (Yamanishi et al.) | HQ333206.1, cloning vector pEFFS-sigHA (Wilke et al., direct submission), 100 % |
| EMCV IRES, `Z75185.1:1118..1693` | Z75185.1, vector pIRES-P (Hobbs et al., Biochem Biophys Res Commun 1998;252:368) | OQ427855.1, vector Puro-iDEST-3×FLAG-hWFS1-WT-IRES2-mCherry (Tokuma, direct submission), 99.7 % on its minus strand |
| P2A, `LT962383.1:1899..1964` | LT962383.1, somatostatin receptor 3 vector (Paramonov, direct submission) | MG880225.1, vector pSOSV/ZsG-FL (Welch et al., PLoS Negl Trop Dis 2018;12:e0006326), 98.5 % |
| T2A, `ON417756.1:739..792` | ON417756.1, MITI donor vector (Maddalena, direct submission) | PZ285977.1, pZac2.1-GfaABC1D-Nrf2-mCherry (Faust et al., Sci Adv 2026;12:eaef6777), 100 % |
| WPRE3, `PZ036182.1:1844..2089` | PZ036182.1, clone AiP12805 (Hunker & Mich, Nat Commun, in press) | MQ017553.1, sequence 23 from WO2021081280 (Shire/Klugmann), 100 % |
| hGH polyA, `MQ026515.1:1..480` | MQ026515.1, sequence 20 from WO2021071976 (Exhaura, Trinity College Dublin) | MP470633.1, sequence 45 from WO2019243847 (Replimune), 100 % |
| rabbit β-globin polyA, `MP470629.1:1..99` | MP470629.1, sequence 41 from WO2019243847 (Replimune) | OP697988.1, vector pCAG-rRBD-xrepEBNA1 (Green et al., BMC Biotechnol 2023;23:7), 99.0 % |
| PuroR (pac), `PZ709769.1:1..600` | PZ709769.1, pac gene (de la Fuente & Simarro, direct submission) | PZ267690.1, clone UbC.WT-AcrIIA4.pLenti (Sabol et al.), 99.8 % |
| BlastR (bsd), `GU474201.1:1..399` | GU474201.1, vector pYC6/CT bsd gene (Zhang et al., direct submission) | PZ267689.1, clone U6.gRNA.pLenti (Sabol et al.), 99.5 % |
| HygR (hph), `AJ439603.2:884..1903` | AJ439603.2, chimaeric cloning vector (Benito et al., Eur J Plant Pathol 1998;104:207) | U40398.1, hygromycin cassette of a mammalian vector (Yin & Schimke, Cancer Res 1995;55:4922), 99.6 % over the 998 bp the two share |
| ZeoR (ble), `L36850.1:5896..6270` | L36850.1, vector pZEO (Drocourt et al., Nucleic Acids Res 1990;18:4009) | MW471125.1, clone pCRISPR01-FnCpf1 (Abdulrachman & Pootanakit, direct submission), 100 % |
| Myc tag, `AJ243205.1:147..179` | AJ243205.1, cloning vector pUR4546 (Harmsen et al., Mol Immunol 2000;37:579) | MF966512.1, epitope-tagged latrophilin-1 (University of Kent, direct submission), which carries it twice, 100 % |
| V5 tag, `AM700585.1:70..111` | AM700585.1, V5-tagged FcεRIγ construct (Viertlböck et al., PNAS 2007;104:11718) | AY682091.1, DEC205 construct (Gurer et al., Blood 2008;112:1231), 100 % |

### The 2A census

The P2A and T2A entries were chosen by counting, not by preference. Searching E-utilities on
2026-10-21 for records whose FEATURES table annotates one of these peptides gave the records
below; each was fetched, the annotated span extracted, and Vettore's own `detect` run with
that span as a one-entry catalogue against every other record in the list. The entries in the
catalogue are the versions found in the most records. Only the two cited records and the two
cross-checks are shipped as fixtures here — the census is repeated from the network, not
offline, which is why the accessions are written down.

**T2A** — 16 records, 7 distinct sequences; the entry (`ON417756.1:739..792`) is found in 8:

KY825159.1, MK816962.1, MW503936.1, MW503937.1, MW503938.1, MW503939.1, MZ127811.1,
ON417756.1, OQ427855.1, PQ659180.1, PQ671723.1, PQ671724.1, PQ671725.1, PV611053.1,
PV611054.1, PZ285977.1

**P2A** — 15 records, 13 distinct sequences; the entry (`LT962383.1:1899..1964`) is found in 6:

KX151729.1, LT962381.1, LT962382.1, LT962383.1, MG880225.1, MK816962.1, MT559572.1,
OK413188.1, OK413189.1, OP243710.1, OR961518.1, PV863446.1, PV863447.1, PZ036169.1,
PZ259196.1

The runner-up T2A version, which the review proposed as the entry (`MK816962.1:2103..2156`),
is found in 2 of the 16.

## Sanger contig fixtures (`contig-*.ab1`)

A forward and a reverse Sanger read of the same template (a COI barcode of the earthworm
*Allolobophora chlorotica*, sample ACHLO006-09), copied unchanged on 2026-10-21 from the
sangeranalyseR package (Chao, Lanfear et al.; https://github.com/roblanf/sangeranalyseR,
`inst/extdata/Allolobophora_chlorotica/ACHLO/`, commit 4c5ebc9). The repository's LICENSE file
is MIT; its package DESCRIPTION says GPL-2.

| File | Original name | SHA-256 |
|---|---|---|
| `contig-ACHLO006-09_1_F.ab1` | `Achl_ACHLO006-09_1_F.ab1` | `fa16320fddcaf0a0d9003551c188192c0bcfb3efc93716017c8443bb28054646` |
| `contig-ACHLO006-09_2_R.ab1` | `Achl_ACHLO006-09_2_R.ab1` | `e885871a776c70292efa45df9634b84812784d95191b174d8a5ef403e19d05c5` |

`tests/test_contig.py` builds the reference from the forward read's trimmed calls between
filler, with its N calls replaced by the reverse read's; the known differences are then planted
in that reference.
