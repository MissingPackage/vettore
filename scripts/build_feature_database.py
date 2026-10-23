"""Build the common-feature database from public NCBI records.

Every entry is a real subsequence of a public GenBank record, taken at the coordinates that
record's own FEATURES table gives for the element, so each annotation Vettore makes can be
traced back to a record and a coordinate range. Nothing here is typed from memory.

Run (needs network):
    python -I scripts/build_feature_database.py > src/vettore/data/common_features.json

NCBI asks callers to identify themselves and to stay under three requests per second without
an API key; both are honoured below. Set NCBI_EMAIL to add a contact address.
"""

from __future__ import annotations

import json
import os
import sys
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TOOL = "vettore"
PAUSE = 0.4


@dataclass(frozen=True)
class Entry:
    name: str
    type: str
    accession: str
    start: int           # 1-based inclusive, as printed in the record's FEATURES table
    end: int
    note: str = ""
    strand: int = 1


#: Only elements whose coordinates were checked against the record's FEATURES table. A note
#: that says "over the N bp the two share" means the independent record carries the element
#: with a few bases less at one end: the identity is over the bases both records hold, and
#: tests/test_features.py also requires that to be at least nine tenths of the entry.
ENTRIES: list[Entry] = [
    # The two ends of the AAV genome carry the same inverted repeat, so one entry covers both.
    Entry("AAV2 ITR", "repeat_region", "NC_001401.2", 1, 145,
          "wild-type AAV2 inverted terminal repeat"),
    Entry("AAV2 ITR (truncated)", "repeat_region", "PZ629150.1", 1, 131,
          "pSub201-type ITR, 14 nt shorter than wild type"),
    Entry("CMV enhancer", "enhancer", "PV031350.1", 96, 475, "human cytomegalovirus enhancer"),
    Entry("CMV promoter", "promoter", "PV031350.1", 476, 679,
          "human cytomegalovirus immediate early promoter"),
    # Full-length CAG (Niwa H, Yamamura K, Miyazaki J, Gene 1991;108:193-199): the whole CMV
    # enhancer (the catalogue entry matches its first 380 bp exactly), the chicken β-actin
    # promoter and the chimeric chicken β-actin/rabbit β-globin intron, ending at the splice
    # acceptor. ON012604.1 (Feher et al., Sci Rep 2022;12:9874) and HG530138.1 (pTAL7B) give
    # the same boundaries; KM279353.1 and LC169511.1 annotate 26 bp more at one end and
    # ~30 bp more at the other, so the shared span is the one taken here.
    Entry("CAG promoter", "promoter", "ON012604.1", 4799, 6471,
          "CMV enhancer, chicken β-actin promoter and the chimeric chicken β-actin/rabbit "
          "β-globin intron; cross-check HG530138.1 (99.2 %)"),
    # The shorter CAG of an AAV vector: it starts 131 bp inside the CMV enhancer and its
    # intron is ~650 bp shorter than the full-length one, so the two are separate entries
    # rather than variants of one. ON785707.1 gives the range in the promoter's own note.
    Entry("CAG promoter (truncated)", "promoter", "ON785707.1", 189, 1125,
          "shortened CAG: the last 249 bp of the CMV enhancer, the chicken β-actin promoter "
          "and a shortened intron"),
    Entry("hSyn promoter", "promoter", "MH282430.1", 449, 920,
          "human synapsin I promoter, neuron-specific"),
    # Human PGK1 promoter, the 511 bp fragment lentiviral vectors use (pLKO.1 and its
    # derivatives). MT921016.1, from another laboratory and another backbone, differs at one
    # base; PX675714.1 and PY500167.1 carry it unchanged.
    Entry("PGK promoter (human)", "promoter", "OQ427862.1", 1, 511,
          "human phosphoglycerate kinase 1 promoter; cross-check MT921016.1 (99.8 %)"),
    # Mouse Pgk-1 promoter, the 500 bp fragment of PGK-neo and PGK-puro selection cassettes.
    # It is not the human promoter: the two share only ~70 % identity, hence two entries.
    # MH325104.1 (Berman et al., J Am Chem
    # Soc 2018;140:18093), LC744045.1 and MT891341.1 (Frei et al., Nat Commun 2020;11:4641)
    # annotate exactly the same 500 bp.
    Entry("PGK promoter (mouse)", "promoter", "MH325104.1", 3321, 3820,
          "mouse Pgk-1 promoter of PGK-neo/PGK-puro cassettes; "
          "cross-check MT891341.1 (100 %)"),
    # Human ubiquitin C promoter with its first exon and intron, as in UbC-driven lentiviral
    # vectors. OL870463.1 (Katzman et al., Viruses 2021;14:11) is the shared span: MN389527.1,
    # MK873089.1 and MW987535.2 annotate the same region with a few extra bases at the start
    # and, in MW987535.2, a polylinker at the end.
    Entry("UbC promoter (human)", "promoter", "OL870463.1", 1096, 2299,
          "human ubiquitin C promoter with its first intron; "
          "cross-check MN389527.1 (99.2 %)"),
    # Human H1 RNA promoter, the 216 bp of the pSUPER-type shRNA cassette. AY640634.1 is
    # Michel et al., Biochem Biophys Res Commun 2005;326:307-312; the same 216 bp is annotated
    # in the pSUPER vectors AY640627.1 and HQ416701.1, and FJ687158.1, from another
    # laboratory, annotates -230 to -9 of the same promoter.
    Entry("H1 promoter (human)", "promoter", "AY640634.1", 1543, 1758,
          "human H1 RNA promoter, RNA polymerase III; cross-check FJ687158.1 (99.5 %)"),
    # Human 7SK RNA promoter: the 243 bp directly upstream of the 7SK transcription start.
    # Checked against the human locus as well — the fragment aligns to the 243 bp in front of
    # RN7SK on NC_000006.12 (98.8 %), where NR_001445.2 starts at the next base. PY046810.1
    # and MQ296444.1 come from unrelated applicants.
    Entry("7SK promoter (human)", "promoter", "PY046810.1", 1, 243,
          "human 7SK RNA promoter, RNA polymerase III; cross-check MQ296444.1 (99.2 %)"),
    # Spleen focus-forming virus U3 promoter from the viral LTR, the 483 bp lentiviral
    # vectors use. KJ697750.1 is the pRRL-sffv vector series; MP040844.1 carries the same
    # fragment with four extra bases at each end.
    Entry("SFFV promoter", "promoter", "KJ697750.1", 6155, 6637,
          "spleen focus-forming virus U3 (LTR) promoter; cross-check MP040844.1 (100 %)"),
    Entry("EF-1α promoter (human)", "promoter", "MT612431.1", 174, 1357,
          "human elongation factor 1α promoter with intron A"),
    Entry("CaMKIIα promoter (mouse)", "promoter", "MW139900.1", 576, 1866,
          "mouse CaMKIIα promoter, excitatory neurons"),
    Entry("GFAP promoter (gfaABC1D)", "promoter", "PZ398262.1", 3384, 4070,
          "compact GFAP promoter, astrocytes"),
    Entry("U6 promoter (human)", "promoter", "PZ629150.1", 3986, 4226,
          "human U6 snRNA promoter, RNA polymerase III"),
    Entry("WPRE", "misc_feature", "MH282430.1", 2377, 2965,
          "woodchuck hepatitis virus post-transcriptional regulatory element"),
    Entry("bGH polyA", "polyA_signal", "PV031350.1", 1551, 1775,
          "bovine growth hormone polyadenylation signal"),
    Entry("SV40 polyA", "polyA_signal", "MH282430.1", 2999, 3230,
          "SV40 polyadenylation signal"),
    Entry("AmpR", "CDS", "PV031350.1", 2619, 3479, "β-lactamase (TEM-1), ampicillin resistance"),
    Entry("AmpR promoter", "promoter", "PV031350.1", 2514, 2618, ""),
    Entry("KanR", "CDS", "PX994933.1", 3733, 4548,
          "aph(3′)-Ia from Tn903, kanamycin resistance"),
    Entry("NeoR/KanR", "CDS", "U55762.1", 2629, 3423,
          "nptII from Tn5, neomycin and kanamycin resistance"),
    Entry("ori", "rep_origin", "PZ629150.1", 6431, 7019, "pUC/ColE1 high-copy origin"),
    Entry("f1 ori", "rep_origin", "PZ629150.1", 4558, 5013, "f1 phage origin"),
    Entry("EGFP", "CDS", "U55762.1", 679, 1398, "enhanced green fluorescent protein"),
    Entry("mCherry", "CDS", "AY678264.1", 1, 711, "red fluorescent protein"),
    Entry("tdTomato", "CDS", "AY678269.1", 1, 1431, "tandem dimer red fluorescent protein"),
    Entry("loxP", "protein_bind", "MT612431.1", 1539, 1572, "Cre recombinase site"),
    Entry("lox2272", "protein_bind", "MT612431.1", 1455, 1488,
          "Cre site with the GGATACTT spacer"),
    Entry("FRT", "misc_feature", "PP645326.1", 1350, 1383, "Flp recombinase site, minimal"),
    Entry("SpCas9", "CDS", "KX151730.1", 808, 5079,
          "Streptococcus pyogenes Cas9, human codon-optimised"),
    Entry("SaCas9", "CDS", "PZ629150.1", 581, 3925,
          "Staphylococcus aureus Cas9 with NLS and 3×HA"),
    Entry("3×FLAG", "misc_feature", "PV031350.1", 956, 1027, "epitope tag"),
    Entry("3×HA", "misc_feature", "PZ629150.1", 3842, 3922, "epitope tag"),

    # ---- Ubiquitous promoters: what a construct carries in the cell-line pretests --------
    # The intron-less core of the human EF-1α promoter, the one AAV and lentiviral cassettes
    # take when the 1.2 kb intron-bearing promoter is too long for the payload. LP754738.1
    # (WO2018026976, Editas) is the record's whole sequence, annotated "EFS promoter";
    # MP933339.1 (EP3789521, an unrelated applicant) carries the same core with 23 bp more at
    # one end, and MP864627.1 has it in reverse orientation.
    Entry("EFS promoter (human)", "promoter", "LP754738.1", 1, 252,
          "intron-less core of the human EF-1α promoter, as the patent record gives it: its "
          "first 5 bases are a cloning site that neither cross-check carries; "
          "cross-check MP933339.1 (99.6 % over the 233 bp the two share)"),
    # CBh (Gray SJ et al., Hum Gene Ther 2011;22:1143; PY026696.1 is Gray,S. and Sinnett,S.): a CMV enhancer fragment, the chicken
    # β-actin promoter and a shortened hybrid intron — the CAG of AAV vectors, made shorter
    # again. PY026696.1 (University of Texas) annotates the 796 bp as one misc_feature;
    # QQ167606.1 (University of Pennsylvania) carries the same promoter in an AAV cassette.
    Entry("CBh promoter", "promoter", "PY026696.1", 1, 796,
          "chicken β-actin hybrid promoter, a compact CAG; cross-check QQ167606.1 (99.9 %)"),
    # The SV40 early promoter with its 21 bp repeats, the 72 bp enhancer repeats and the
    # origin: the block that drives the resistance marker of most mammalian vectors.
    # X85976.1 (pSG424) annotates it; M83237.1 (RSV.5(neo), another laboratory) carries the
    # same 348 bp as "SV40 early promoter and origin of replication".
    Entry("SV40 early promoter", "promoter", "X85976.1", 2537, 2884,
          "SV40 early promoter with the 72 bp enhancer repeats and the origin; "
          "cross-check M83237.1 (100 % over the 344 bp the two share)"),
    # The human ACTB promoter, the 520 bp fragment of β-actin-driven cassettes. Not the
    # chicken β-actin promoter of CAG and CBh, which is a different species' promoter and a
    # separate entry. JN964682.1 annotates the fragment; AY952326.1, from another laboratory,
    # holds it inside a 2.6 kb β-actin fragment.
    Entry("β-actin promoter (human)", "promoter", "JN964682.1", 20222, 20741,
          "human ACTB promoter; cross-check AY952326.1 (99.8 %)"),
    # Herpes simplex virus thymidine kinase promoter, the weak constitutive promoter of
    # selection cassettes. X64410.1 (pBLCAT2) and AJ277960.1 (pMAR) are unrelated vectors
    # and annotate the same 170 bp.
    Entry("HSV TK promoter", "promoter", "X64410.1", 429, 598,
          "herpes simplex virus thymidine kinase promoter; "
          "cross-check AJ277960.1 (100 % over the 169 bp the two share)"),
    # The promoter part of the Rous sarcoma virus LTR, annotated on the minus strand in
    # LT727488.1 ("RSV-LTR", regulatory_class promoter), so the entry is the reverse complement
    # of that span. It is not the whole LTR: M83237.1 annotates that as 528 bp (4581..5108) with
    # these 182 bp inside it. Vectors carry this fragment, which is why it is the entry; which
    # parts of U3, R and U5 fall inside it is not annotated in either record, so it is not said.
    Entry("RSV LTR promoter", "promoter", "LT727488.1", 3997, 4178,
          "promoter fragment of the Rous sarcoma virus LTR, 182 bp of the 528 bp LTR; "
          "cross-check M83237.1 (99.5 %)", strand=-1),
    # The tetracycline-responsive promoter: seven tet operators in front of a minimal CMV
    # promoter, the "tight" version with no residual activity. MT891336.1 and OQ427855.1 come
    # from different laboratories and give the same 315 bp.
    Entry("TRE promoter (tight)", "promoter", "MT891336.1", 554, 868,
          "7× tet operator and a minimal CMV promoter, doxycycline-inducible; this is the "
          "\"tight\" version both records annotate, not the 271 bp classic TRE; "
          "cross-check OQ427855.1 (100 %)"),

    # ---- Recombinases and their sites ---------------------------------------------------
    # Cre from bacteriophage P1. Public records carry two unrelated codon versions of this
    # gene (the P1 sequence and codon-optimised ones that share only ~80 % of their bases);
    # this is the P1 version, which five unrelated vectors annotate identically.
    Entry("Cre", "CDS", "MK854762.1", 215, 1246,
          "Cre recombinase from bacteriophage P1; cross-check AB542060.1 (100 %)"),
    # Flp comes in two versions vectors use side by side: FLPe, the thermostable yeast gene,
    # and FLPo, codon-optimised for mammalian cells. Both spans are the nls-FLP fusion the two
    # records annotate — an SV40 nuclear localisation signal in front of the 423-residue
    # recombinase — so the two encode the same protein in two codon versions. They share 73 %
    # of their bases, far too little for one entry to find the other, so each has its own.
    Entry("FLPe", "CDS", "GU253315.1", 2581, 3879,
          "thermostable Flp recombinase (FLPe), yeast gene, behind an SV40 NLS; "
          "cross-check KM521212.1 (99.8 % over the 1272 bp the two share)"),
    Entry("FLPo", "CDS", "KX151728.1", 1842, 3140,
          "Flp recombinase codon-optimised for mammalian cells (FLPo), behind an SV40 NLS; "
          "cross-check MH493814.1 (100 % over the 1272 bp the two share)"),
    # loxN: a Cre site whose spacer differs from loxP, so Cre recombines it with itself and
    # not with loxP — the second site of a Brainbow-type cassette.
    Entry("loxN", "protein_bind", "OR242273.1", 3825, 3858,
          "Cre site with the GTATACCT spacer; cross-check KP686439.1 (100 %)"),
    # F3: the Flp site variant used opposite FRT, as lox2272 is used opposite loxP.
    Entry("F3", "misc_feature", "LC086648.1", 2498, 2531,
          "Flp site variant 3 (F3), recombined with itself and not with FRT; "
          "cross-check HQ333206.1 (100 %)"),

    # ---- Expression parts ---------------------------------------------------------------
    # EMCV IRES, the 576 bp of the bicistronic vectors. Z75185.1 (pIRES-P) annotates it as
    # "EMC IRES sequence"; OQ427855.1, an unrelated laboratory's vector, carries it on its
    # minus strand annotated "IRES2". Cross-checked against that vector rather than against
    # the patent sequence PG196517.1, which is the element alone with nothing around it.
    Entry("EMCV IRES", "misc_feature", "Z75185.1", 1118, 1693,
          "encephalomyocarditis virus internal ribosome entry site; "
          "cross-check OQ427855.1 (99.7 %)"),
    # P2A and T2A: the self-cleaving peptides of multi-gene cassettes. Both are short and
    # written in many codon versions, so no single entry finds them all. The two below were
    # picked by census rather than by preference: of the public records that annotate one of
    # these peptides, 16 annotate T2A in 7 distinct sequences and 15 annotate P2A in 13, and
    # these are the versions detection finds in the most of those records — T2A in 8 of the
    # 16, P2A in 6 of the 15. tests/data/README.md lists every accession in that census, so
    # the count can be repeated. A construct whose 2A is written another way is not found,
    # and the user guide says so.
    Entry("P2A", "misc_feature", "LT962383.1", 1899, 1964,
          "porcine teschovirus-1 2A self-cleaving peptide, with its GSG linker; "
          "cross-check MG880225.1 (98.5 %)"),
    Entry("T2A", "misc_feature", "ON417756.1", 739, 792,
          "Thosea asigna virus 2A self-cleaving peptide; cross-check PZ285977.1 (100 %)"),
    # WPRE3, the 246 bp WPRE: the shortened element of AAV cassettes, where the full 589 bp
    # WPRE does not fit. PZ036182.1 annotates it; MQ017553.1 (WO2021081280) is an unrelated
    # applicant's copy.
    Entry("WPRE3", "misc_feature", "PZ036182.1", 1844, 2089,
          "shortened woodchuck hepatitis virus post-transcriptional regulatory element; "
          "cross-check MQ017553.1 (100 %)"),
    Entry("hGH polyA", "polyA_signal", "MQ026515.1", 1, 480,
          "human growth hormone polyadenylation signal; cross-check MP470633.1 (100 %)"),
    Entry("rabbit β-globin polyA", "polyA_signal", "MP470629.1", 1, 99,
          "rabbit β-globin polyadenylation signal (Gil A, Proudfoot NJ, Cell 1987;49:399); "
          "cross-check OP697988.1 (99.0 %)"),

    # ---- Selection markers and tags -----------------------------------------------------
    Entry("PuroR (pac)", "CDS", "PZ709769.1", 1, 600,
          "pac from Streptomyces alboniger, puromycin resistance; "
          "cross-check PZ267690.1 (99.8 %)"),
    Entry("BlastR (bsd)", "CDS", "GU474201.1", 1, 399,
          "bsd from Aspergillus terreus, blasticidin S resistance; "
          "cross-check PZ267689.1 (99.5 %)"),
    Entry("HygR (hph)", "CDS", "AJ439603.2", 884, 1903,
          "hph from Escherichia coli, hygromycin B resistance; "
          "cross-check U40398.1 (99.6 % over the 998 bp the two share)"),
    Entry("ZeoR (ble)", "CDS", "L36850.1", 5896, 6270,
          "Sh ble from Streptoalloteichus hindustanus, zeocin and phleomycin resistance; "
          "cross-check MW471125.1 (100 %)"),
    # The two tags, like the 2A peptides, exist in more than one codon version; these are the
    # versions two unrelated vectors share.
    Entry("Myc tag", "misc_feature", "AJ243205.1", 147, 179,
          "c-Myc epitope tag; cross-check MF966512.1 (100 %)"),
    Entry("V5 tag", "misc_feature", "AM700585.1", 70, 111,
          "V5 epitope tag; cross-check AY682091.1 (100 %)"),
]

COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def fetch(accession: str) -> str:
    params = {"db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text",
              "tool": TOOL}
    if email := os.environ.get("NCBI_EMAIL"):
        params["email"] = email
    url = f"{EUTILS}?{urllib.parse.urlencode(params)}"
    # E-utilities answers 429 when several callers share an address, and the whole run is
    # wasted if one record fails, so a refusal is waited out rather than raised.
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                text = response.read().decode()
            break
        except urllib.error.URLError as error:
            if attempt == 4:
                raise SystemExit(f"NCBI did not answer for {accession}: {error}") from error
            print(f"{accession}: {error}; waiting {2 ** attempt} s", file=sys.stderr)
            time.sleep(2 ** attempt)
    lines = [line.strip() for line in text.splitlines() if line and not line.startswith(">")]
    sequence = "".join(lines).upper()
    if not sequence:
        raise SystemExit(f"NCBI returned no sequence for {accession}")
    return sequence


def main() -> None:
    cache: dict[str, str] = {}
    features = []
    for entry in ENTRIES:
        if entry.accession not in cache:
            cache[entry.accession] = fetch(entry.accession)
            time.sleep(PAUSE)
        record = cache[entry.accession]
        if entry.end > len(record):
            raise SystemExit(f"{entry.accession} is {len(record)} bp, shorter than "
                             f"{entry.name} at {entry.start}..{entry.end}")
        sequence = record[entry.start - 1 : entry.end]
        if set(sequence) - set("ACGTRYSWKMBDHVN"):
            raise SystemExit(f"{entry.name} contains unexpected characters")
        # An element a record annotates on its minus strand (the RSV LTR promoter below) is
        # stored the way it reads in a vector, and its source says so in GenBank's own
        # notation, so that the claim stays checkable against the record.
        span = f"{entry.start}..{entry.end}"
        if entry.strand == -1:
            sequence = sequence.translate(COMPLEMENT)[::-1]
            span = f"complement({span})"
        features.append({
            "name": entry.name,
            "type": entry.type,
            "sequence": sequence,
            "length": len(sequence),
            "note": entry.note,
            "source": f"{entry.accession}:{span}",
        })
    payload = {
        "source": "NCBI GenBank/RefSeq records fetched with E-utilities; coordinates from "
                  "each record's FEATURES table",
        "features": sorted(features, key=lambda f: (-f["length"], f["name"])),
    }
    json.dump(payload, sys.stdout, indent=1)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
