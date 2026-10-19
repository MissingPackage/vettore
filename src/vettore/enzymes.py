"""Restriction-enzyme database for the editor, generated from Biopython's REBASE data.

Only commercially available enzymes with a characterised, single cleavage position are kept,
and only one enzyme per specificity (recognition site + cut positions) — the equivalent of a
"nonredundant commercial" set. NEB and Thermo Fisher names are preferred, so users see BamHI
rather than one of its many isoschizomers.

The output uses Open Vector Editor's enzyme model:
- ``site``: IUPAC recognition sequence (lower case)
- ``forwardRegex`` / ``reverseRegex``: regexes for the site on the top / bottom strand
- ``topSnipOffset`` / ``bottomSnipOffset``: cut positions relative to the first base of the
  site, in top-strand coordinates (EcoRI G^AATTC -> 1 and 5)

Regenerate the bundled copy with ``python -m vettore.enzymes > frontend/src/data/enzymes.json``
(a test fails if the committed copy is stale).

REBASE data: Roberts RJ et al., REBASE, http://rebase.neb.com (redistributed free of charge,
with citation, per REBASE terms).
"""

from __future__ import annotations

import json
import sys
from typing import Any

import Bio
from Bio import Restriction
from Bio.Restriction.Restriction_Dictionary import suppliers as SUPPLIERS
from Bio.Seq import Seq

# Supplier preference when choosing the representative of a group of isoschizomers.
SUPPLIER_RANK = ["N", "B", "R", "K", "M", "S", "X", "V", "I", "J", "O", "C", "E", "Q", "Y"]

IUPAC = {
    "a": "a", "c": "c", "g": "g", "t": "t",
    "r": "[ag]", "y": "[ct]", "s": "[cg]", "w": "[at]", "k": "[gt]", "m": "[ac]",
    "b": "[cgt]", "d": "[agt]", "h": "[act]", "v": "[acg]", "n": "[acgt]",
}

GOLDEN_GATE = {"BsaI", "BsmBI", "Esp3I", "BbsI", "BpiI", "SapI", "BspQI", "PaqCI", "BtgZI"}

# Modification-dependent enzymes (cut only methylated / hydroxymethylated DNA); their degenerate
# REBASE sites (e.g. "C" for AbaSI) are meaningless on plasmid maps.
MODIFICATION_DEPENDENT = {"AbaSI", "FspEI", "LpnPI", "MspJI"}

# Names used in most protocols win over their isoschizomers (e.g. BsmBI over Esp3I).
PREFERRED_NAMES = {
    "AatII", "Acc65I", "AflII", "AgeI", "ApaI", "AscI", "AvrII", "BamHI", "BbsI", "BglII",
    "BsaI", "BsmBI", "BsrGI", "BstBI", "ClaI", "EcoRI", "EcoRV", "FseI", "HindIII", "HpaI",
    "KpnI", "MfeI", "MluI", "NcoI", "NdeI", "NheI", "NotI", "NsiI", "PacI", "PmeI", "PstI",
    "SacI", "SacII", "SalI", "SapI", "SbfI", "SfiI", "SmaI", "SpeI", "SphI", "SwaI", "XbaI",
    "XhoI", "XmaI", "ZraI",
}


def site_regex(site: str) -> str:
    return "".join(IUPAC[c] for c in site.lower())


def _rank(enzyme: Any) -> tuple:
    best_supplier = min((SUPPLIER_RANK.index(s) for s in enzyme.suppl if s in SUPPLIER_RANK),
                        default=len(SUPPLIER_RANK))
    return (str(enzyme) not in PREFERRED_NAMES, best_supplier, -len(enzyme.suppl), str(enzyme))


def _offsets(enzyme: Any) -> tuple[str, int, int] | None:
    """Site and OVE cut offsets, re-oriented on the reverse strand if the cut is upstream."""
    site = enzyme.site
    top = enzyme.fst5
    bottom = enzyme.size + enzyme.fst3
    if top < 0 or bottom < 0:
        site = str(Seq(site).reverse_complement())
        top, bottom = enzyme.size - bottom, enzyme.size - top
    if top < 0 or bottom < 0:
        return None
    return site, top, bottom


def build_database() -> dict[str, dict[str, Any]]:
    groups: dict[tuple, list[Any]] = {}
    for enzyme in Restriction.CommOnly:
        if enzyme.fst5 is None or enzyme.scd5 is not None or enzyme.is_unknown():
            continue  # uncharacterised or two-sided cutters are not supported yet
        if str(enzyme) in MODIFICATION_DEPENDENT:
            continue
        oriented = _offsets(enzyme)
        if oriented is None:
            continue
        site, top, bottom = oriented
        # Same specificity regardless of which strand the site is written on; only
        # orientations with downstream (non-negative) cut offsets are valid for the editor.
        rc = str(Seq(site).reverse_complement())
        candidates = [(site, top, bottom)]
        if len(site) - bottom >= 0 and len(site) - top >= 0:
            candidates.append((rc, len(site) - bottom, len(site) - top))
        key = min(candidates)
        groups.setdefault(key, []).append(enzyme)

    database: dict[str, dict[str, Any]] = {}
    for (site, top, bottom), members in groups.items():
        members.sort(key=_rank)
        chosen = members[0]
        site_l = site.lower()
        rc = str(Seq(site).reverse_complement()).lower()
        size = len(site)
        name = str(chosen)
        isoschizomers = sorted(str(m) for m in members[1:])
        database[name.lower()] = {
            "name": name,
            "site": site_l,
            "forwardRegex": site_regex(site_l),
            "reverseRegex": site_regex(rc),
            "topSnipOffset": top,
            "bottomSnipOffset": bottom,
            "isType2S": top > size or bottom > size,
            "suppliers": sorted(SUPPLIERS[s][0] for s in chosen.suppl),
            "isoschizomers": isoschizomers,
            "goldenGate": bool(GOLDEN_GATE & ({name} | set(isoschizomers))),
        }
    return dict(sorted(database.items()))


def main() -> None:
    payload = {
        "source": f"REBASE via Biopython {Bio.__version__} (Bio.Restriction.CommOnly)",
        "enzymes": build_database(),
    }
    json.dump(payload, sys.stdout, indent=1, ensure_ascii=False)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
