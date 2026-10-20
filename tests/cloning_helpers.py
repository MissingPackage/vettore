"""Helpers for cloning tests: site-free random DNA and rotation-invariant comparison."""

import random
import re

from Bio import Restriction
from Bio.Seq import Seq
from pydna.dseqrecord import Dseqrecord

COMMON = [Restriction.EcoRI, Restriction.BamHI, Restriction.BglII, Restriction.KpnI,
          Restriction.XbaI, Restriction.PstI, Restriction.HindIII, Restriction.NotI,
          Restriction.XhoI, Restriction.SalI, Restriction.NcoI, Restriction.BsaI,
          Restriction.BsmBI, Restriction.BbsI, Restriction.SapI, Restriction.SmaI,
          Restriction.EcoRV, Restriction.NdeI]


IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "S": "[CG]",
         "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]",
         "V": "[ACG]", "N": "[ACGT]"}


def _site_patterns(enzymes) -> list:
    patterns = []
    for e in enzymes:
        for site in {e.site, str(Seq(e.site).reverse_complement())}:
            patterns.append(re.compile("".join(IUPAC[c] for c in site) + "$"))
    return patterns


def clean_dna(n: int, seed: int, avoid=COMMON) -> str:
    """Random DNA of length n with no site of the given enzymes on either strand, built base by
    base so that no recognition site is ever completed."""
    rng = random.Random(seed)
    patterns = _site_patterns(avoid)
    window = max(len(e.site) for e in avoid)
    out: list[str] = []
    while len(out) < n:
        for base in rng.sample("ACGT", 4):
            tail = "".join(out[-(window - 1):]) + base
            if not any(p.search(tail) for p in patterns):
                out.append(base)
                break
        else:  # dead end (practically impossible): backtrack one base
            out.pop()
    s = "".join(out)
    assert not any(e.search(Seq(s), linear=True) for e in avoid)
    return s


def rc(s: str) -> str:
    return str(Seq(s).reverse_complement())


def same_circular(a: str, b: str) -> bool:
    return Dseqrecord(a, circular=True).seguid() == Dseqrecord(b, circular=True).seguid()


def doc(name: str, sequence: str, circular: bool, features=(), primers=()) -> dict:
    return {"name": name, "sequence": sequence, "circular": circular,
            "features": list(features), "primers": list(primers)}


def feature(name: str, start: int, end: int, strand: int = 1, type_: str = "CDS") -> dict:
    """OVE feature: 0-based, end-inclusive."""
    return {"name": name, "type": type_, "start": start, "end": end, "strand": strand,
            "forward": strand == 1}
