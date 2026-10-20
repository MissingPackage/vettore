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


#: Only elements whose coordinates were checked against the record's FEATURES table.
ENTRIES: list[Entry] = [
    # The two ends of the AAV genome carry the same inverted repeat, so one entry covers both.
    Entry("AAV2 ITR", "repeat_region", "NC_001401.2", 1, 145,
          "wild-type AAV2 inverted terminal repeat"),
    Entry("AAV2 ITR (truncated)", "repeat_region", "PZ629150.1", 1, 131,
          "pSub201-type ITR, 14 nt shorter than wild type"),
    Entry("CMV enhancer", "enhancer", "PV031350.1", 96, 475, "human cytomegalovirus enhancer"),
    Entry("CMV promoter", "promoter", "PV031350.1", 476, 679,
          "human cytomegalovirus immediate early promoter"),
    Entry("CAG promoter", "promoter", "ON785707.1", 189, 1125,
          "CMV enhancer, chicken β-actin promoter and rabbit β-globin intron"),
    Entry("hSyn promoter", "promoter", "MH282430.1", 449, 920,
          "human synapsin I promoter, neuron-specific"),
    Entry("EF-1α promoter", "promoter", "MT612431.1", 174, 1357,
          "human elongation factor 1α promoter with intron A"),
    Entry("CaMKIIα promoter", "promoter", "MW139900.1", 576, 1866,
          "mouse CaMKIIα promoter, excitatory neurons"),
    Entry("GFAP promoter (gfaABC1D)", "promoter", "PZ398262.1", 3384, 4070,
          "compact GFAP promoter, astrocytes"),
    Entry("U6 promoter", "promoter", "PZ629150.1", 3986, 4226,
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
]


def fetch(accession: str) -> str:
    params = {"db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text",
              "tool": TOOL}
    if email := os.environ.get("NCBI_EMAIL"):
        params["email"] = email
    url = f"{EUTILS}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=60) as response:
        text = response.read().decode()
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
        features.append({
            "name": entry.name,
            "type": entry.type,
            "sequence": sequence,
            "length": len(sequence),
            "note": entry.note,
            "source": f"{entry.accession}:{entry.start}..{entry.end}",
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
