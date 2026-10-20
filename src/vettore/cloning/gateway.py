"""Gateway BP and LR reactions, single pair and MultiSite. Routes under /api/cloning.

Site-specific recombination between att sites, run by pydna (``pydna.assembly2.gateway_assembly``
with the att consensus sequences of ``pydna.gateway``, which pydna takes from the published
site definitions; Vettore adds no att sequence of its own):

- **BP**: attB x attP -> attL + attR. An attB-flanked insert (a PCR product with attB tails, or
  a clone carrying attB sites) plus a donor vector with attP sites give the **entry clone**
  (attL) and an attR by-product.
- **LR**: attL x attR -> attB + attP. One to four entry clones plus a destination vector with
  attR sites give the **expression clone** (attB) and an attP by-product. With more than one
  entry clone the reaction is a MultiSite reaction: the fragments line up through the numbered
  sites (att1 to att5, each attL pairing with the attR of the same number), so an entry clone
  with attL4 and attR1 joins a destination vector's attR4 on one side and an attL1 entry clone
  on the other.

Each att site has its recombination core on one side: attL and attB carry it at their 3' end,
attR and attP on the 5' side of their phage arm. The core side is the side that survives in the
product, so it tells the insert of an entry clone from the backbone of a vector whatever the
orientation of the site — including the reversed sites of the MultiSite systems (the attR1 of a
5' entry clone, attB5r, attL4r and so on), which read the same way round as their unreversed
partner but face the other way.

Only reactions where **both** att sites of a junction recombine are reported: a product from a
single site is the co-integrate intermediate of the reaction, not what grows on the selection
plates.

Validated end to end against public GenBank records that carry annotated att sites
(``tests/test_cloning_gateway.py``, ``tests/test_cloning_multisite.py``).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import pairwise
from typing import Any, Literal

from Bio.Seq import Seq
from fastapi import APIRouter
from pydantic import BaseModel
from pydna.assembly2 import gateway_assembly
from pydna.dseqrecord import Dseqrecord
from pydna.gateway import find_gateway_sites

from .base import CloningError, CloningResult, Product, SequenceDoc, from_dseqrecord, to_dseqrecord

router = APIRouter()

Reaction = Literal["BP", "LR"]

OPERATIONS: dict[str, str] = {"BP": "Gateway BP reaction", "LR": "Gateway LR reaction"}
WANTED: dict[str, str] = {"BP": "Entry clone", "LR": "Expression clone"}
BY_PRODUCT: dict[str, str] = {"BP": "attR by-product", "LR": "attP by-product"}
INPUT_LABELS: dict[str, tuple[str, str]] = {
    "BP": ("attB-flanked insert", "donor vector (attP sites)"),
    "LR": ("entry clone (attL sites)", "destination vector (attR sites)"),
}
EXPECTED_SITES: dict[str, tuple[str, str]] = {"BP": ("attB", "attP"), "LR": ("attL", "attR")}
# ccdB kills the cells that keep the unreacted vector; a product that still carries it would
# not grow, so finding it in the product means the reaction was read the wrong way round.
COUNTER_SELECTION = ("ccdb", "ccdb gene", "ccdb protein")

MAX_SITE = 300      # the longest att consensus is attP, about 230 bp
MAX_ENTRY_CLONES = 4    # the largest MultiSite reaction in use takes four entry clones

# The core of an attL or attB site is at its 3' end; an attR carries it on the 5' side of its
# phage arm, an attP in the middle, with an arm on each side. The side of the core is the side
# of the molecule that ends up in the product.
CORE_AT_END = ("attL", "attB")

# phiC31 integrase uses att sites of its own, which no Gateway clonase recognises. The 34 bp
# attB and the attP core below are the published phiC31 sites (Groth et al. 2000, PNAS
# 97:5995), read from two public records that annotate them: the attB of pUASTattB (GenBank
# EF362409.1, Bischof et al. 2007, PNAS 104:3312) and the attB and attP of GenBank PX994930.1.
PHIC31_SITES = {
    "attB": "GTGCCAGGGCGTGCCCTTGGGCTCCCCGGGCGCG",
    "attP": "GTAGTGCCCCAACTGGGGTAACCTTTGAGTTCTCTCAGTTGGGGG",
}


@dataclass(frozen=True)
class Site:
    """An att site found in a molecule: ``attL1`` is kind ``attL``, number 1."""

    name: str
    start: int
    end: int
    length: int
    strand: int

    @property
    def kind(self) -> str:
        return self.name[:4]

    @property
    def number(self) -> int:
        return int(self.name[4:])


def _bounds(location) -> tuple[int, int, int]:
    """(start, end, length) of a site along the top strand; sites may wrap the origin, and
    pydna reports those as a two-part location."""
    parts = list(location.parts)
    length = sum(len(part) for part in parts)
    if len(parts) == 1:
        return int(parts[0].start), int(parts[0].end), length
    ordered = sorted(parts, key=lambda part: int(part.start))
    return int(ordered[-1].start), int(ordered[0].end), length   # wraps: (last part, first part)


def _matches(record: Dseqrecord) -> list[Site]:
    """Every site-sized att match in a molecule, in order along the top strand."""
    found = []
    for name, locations in find_gateway_sites(record, False).items():
        for location in locations:
            start, end, length = _bounds(location)
            if length <= MAX_SITE:
                found.append(Site(name, start, end, length, location.strand or 1))
    return sorted(found, key=lambda site: (site.start, -site.length))


def sites(record: Dseqrecord) -> list[Site]:
    """The att sites of a molecule as the reaction sees them. A match inside a longer one is
    left out: the attB core inside every attL, attR and attP, and the attL and attR halves
    inside an attP. So an entry clone reads attL1, attL2 and a donor vector attP1, attP2."""
    found = _matches(record)
    n = len(record.seq)

    def inside(site: Site, other: Site) -> bool:
        return (site.start - other.start) % n + site.length <= other.length

    return [site for site in found
            if not any(other.length > site.length and inside(site, other) for other in found)]


def site_names(record: Dseqrecord) -> list[str]:
    """The names of the att sites of a molecule, sorted (``["attL1", "attL2"]``)."""
    return sorted({site.name for site in sites(record)})


def _reacting(record: Dseqrecord, reaction: str) -> list[Site]:
    """The sites of a molecule that take part in this reaction: attL and attR in an LR
    reaction (a MultiSite entry clone carries one of each), attB and attP in a BP reaction."""
    return [site for site in sites(record) if site.kind in EXPECTED_SITES[reaction]]


def _faces_after(site: Site) -> bool:
    """Does the part of the molecule that survives the reaction lie after this site?"""
    return (site.kind in CORE_AT_END) == (site.strand != -1)


def _article(reaction: str) -> str:
    """"An LR reaction" (read "el-ar"), "a BP reaction"."""
    return "An" if reaction == "LR" else "A"


def _site_text(record: Dseqrecord, doc_name: str) -> str:
    found = site_names(record)
    if not found:
        return f"“{doc_name}”: no att site found"
    return f"“{doc_name}”: " + ", ".join(found)


def _all_site_text(docs: Sequence[SequenceDoc], records: Sequence[Dseqrecord]) -> str:
    return "; ".join(_site_text(record, doc.name) for doc, record in zip(docs, records))


def _kept_region(record: Dseqrecord, reaction: str) -> str | None:
    """What this molecule puts into the product: the piece its att sites both face.

    Each site faces the DNA on the side of its recombination core, which is the side the
    reaction keeps. In an entry clone or an attB-flanked insert both sites face the insert; in
    a destination or donor vector they face the backbone (the ccdB cassette between them is
    replaced). The longest such piece is the one that tells molecules apart."""
    found = _reacting(record, reaction)
    if len(found) < 2:
        return None
    sequence = str(record.seq).upper()
    n = len(sequence)
    doubled = sequence * 2
    pairs = list(pairwise(found))
    if record.circular:
        pairs.append((found[-1], found[0]))
    kept = []
    for site, following in pairs:
        length = (following.start - site.end) % n if record.circular \
            else following.start - site.end
        if not 0 < length <= n - site.length - following.length:   # overlapping sites
            continue
        if _faces_after(site) and not _faces_after(following):
            kept.append(doubled[site.end: site.end + length])
    if not record.circular:        # the ends of a linear molecule are kept by an outward site
        if not _faces_after(found[0]):
            kept.append(sequence[: found[0].start])
        if _faces_after(found[-1]):
            kept.append(sequence[found[-1].end:])
    return max(kept, key=len) if kept else None


def _probe(region: str | None) -> str | None:
    """A 60 bp piece from the middle of a region, or None when the region is too short to tell
    molecules apart."""
    if not region or len(region) < 20:
        return None
    middle = max(0, (len(region) - 60) // 2)
    return region[middle: middle + 60]


def _carries(product: Dseqrecord, region: str | None) -> bool:
    """Is that piece of a parent still in the product (either strand)?"""
    probe = _probe(region)
    if probe is None:
        return False
    doubled = (str(product.seq).upper() * 2)
    return probe in doubled or str(Seq(probe).reverse_complement()).upper() in doubled


def _has_counter_selection(record: Dseqrecord) -> bool:
    for feature in record.features:
        for values in feature.qualifiers.values():
            for value in (values if isinstance(values, (list, tuple)) else [values]):
                if str(value).strip().lower() in COUNTER_SELECTION:
                    return True
    return False


def phic31_sites(record: Dseqrecord) -> list[str]:
    """The phiC31 att sites of a molecule, by name. These are the sites of the phiC31
    integrase used in fly and mammalian landing-site transgenesis: they look like "att" sites
    but no Gateway clonase recombines them."""
    sequence = str(record.seq).upper()
    window = sequence + sequence[:200] if record.circular else sequence
    found = []
    for name, core in PHIC31_SITES.items():
        reverse = str(Seq(core).reverse_complement()).upper()
        if core in window or reverse in window:
            found.append(name)
    return found


def _missing_site_message(reaction: str, docs: Sequence[SequenceDoc],
                          records: Sequence[Dseqrecord], index: int, kind: str | None) -> str:
    """Why an input cannot take part, and what the reaction needs instead."""
    want_left, want_right = EXPECTED_SITES[reaction]
    label = _input_label(reaction, index, len(docs))
    doc, record = docs[index], records[index]
    if len(docs) == 2:
        needs = (f"{_article(reaction)} {reaction} reaction needs {want_left} sites in the "
                 f"first molecule and {want_right} sites in the second.")
    else:
        needs = (f"{_article(reaction)} {reaction} reaction needs {want_left} sites in the "
                 f"entry clones and {want_right} sites in the destination vector.")
    if kind is None:
        phic31 = phic31_sites(record)
        extra = ""
        if phic31:
            extra = (f" It carries a phiC31 {' and '.join(phic31)} site, which the phiC31 "
                     "integrase recombines; Gateway clonase does not use it.")
        return f"The {label} “{doc.name}” has no Gateway att site.{extra} {needs}"
    other = "LR" if reaction == "BP" else "BP"
    found = [site.name for site in _reacting(record, other)]
    hint = ""
    if found and not _reacting(record, reaction):
        hint = (f"It carries {', '.join(sorted(set(found)))}, which recombine in "
                f"{_article(other).lower()} {other} reaction. ")
    return (f"The {label} “{doc.name}” has no {kind} site. {hint}{needs} "
            f"Sites found — {_all_site_text(docs, records)}")


def _input_label(reaction: str, index: int, count: int) -> str:
    left, right = INPUT_LABELS[reaction]
    if index == count - 1:
        return right
    if reaction == "LR" and count > 2:
        return f"entry clone {index + 1} (attL or attR sites)"
    return left


def _pairing_message(reaction: str, docs: Sequence[SequenceDoc],
                     records: Sequence[Dseqrecord]) -> str:
    """Why the molecules do not recombine, naming the sites found in each one."""
    want_left, want_right = EXPECTED_SITES[reaction]
    numbers = sorted({site.number for record in records for site in _reacting(record, reaction)})
    unpaired = []
    for number in numbers:
        left = [doc.name for doc, record in zip(docs, records)
                for site in _reacting(record, reaction)
                if site.number == number and site.kind == want_left]
        right = [doc.name for doc, record in zip(docs, records)
                 for site in _reacting(record, reaction)
                 if site.number == number and site.kind == want_right]
        if left and not right:
            unpaired.append(f"{want_left}{number} in “{left[0]}” has no {want_right}{number} "
                            "in the other molecules")
        elif right and not left:
            unpaired.append(f"{want_right}{number} in “{right[0]}” has no {want_left}{number} "
                            "in the other molecules")
    reason = (" ".join(f"{sentence}." for sentence in unpaired) if unpaired
              else (f"{_article(reaction)} {reaction} reaction pairs {want_left}1 with "
                    f"{want_right}1, {want_left}2 with {want_right}2, and so on up to site 5."))
    molecules = "two molecules" if len(docs) == 2 else f"{len(docs)} molecules"
    return (f"No {reaction} product: the att sites of the {molecules} do not pair. "
            f"Sites found — {_all_site_text(docs, records)}. {reason}")


def _check_inputs(reaction: str, docs: Sequence[SequenceDoc],
                  records: Sequence[Dseqrecord]) -> None:
    """Every molecule must carry the sites its place in the reaction needs."""
    want_left, want_right = EXPECTED_SITES[reaction]
    for index, record in enumerate(records):
        wanted = want_right if index == len(records) - 1 else want_left
        if not sites(record):
            raise CloningError(_missing_site_message(reaction, docs, records, index, None))
        # attB is a part of attL and of attR, so only the vector-side site tells the two
        # reactions apart; a molecule with no att site at all is reported on its own.
        if not any(site.kind == wanted for site in _reacting(record, reaction)):
            raise CloningError(_missing_site_message(reaction, docs, records, index, wanted))


def _check_donor_pair(docs: Sequence[SequenceDoc], records: Sequence[Dseqrecord],
                      donor_sites: Sequence[int]) -> None:
    """The att pair the user picked for the donor vector must be the pair it carries."""
    donor, record = docs[-1], records[-1]
    found = sorted({site.number for site in _reacting(record, "BP") if site.kind == "attP"})
    wanted = sorted(set(donor_sites))
    if found != wanted:
        raise CloningError(
            f"The donor vector “{donor.name}” does not carry the att pair you chose "
            f"(attP{' and attP'.join(str(number) for number in wanted)}). "
            f"Sites found — {_site_text(record, donor.name)}. Choose the pair the donor "
            "carries, or another donor vector.")


def gateway(reaction: Reaction, *docs: SequenceDoc,
            donor_sites: Sequence[int] | None = None) -> CloningResult:
    """Run a BP or an LR reaction.

    ``docs`` are the molecules in reaction order: for LR one to four entry clones followed by
    the destination vector, for BP the attB-flanked insert followed by the donor vector.
    ``donor_sites`` is the att pair chosen for a BP donor (``[4, 1]`` for a donor with attP4
    and attP1), checked against the sites the donor carries.
    """
    if len(docs) < 2:
        raise CloningError("Choose both molecules for the reaction")
    if reaction == "BP" and len(docs) > 2:
        raise CloningError(
            "A BP reaction takes one attB-flanked molecule and one donor vector. Build the "
            "entry clones one at a time, then combine them in an LR reaction.")
    if len(docs) > MAX_ENTRY_CLONES + 1:
        raise CloningError(
            f"An LR reaction takes at most {MAX_ENTRY_CLONES} entry clones and one "
            f"destination vector ({len(docs)} molecules chosen).")

    records = [to_dseqrecord(doc) for doc in docs]
    vector_doc = docs[-1]
    _check_inputs(reaction, docs, records)
    if donor_sites and reaction == "BP":
        _check_donor_pair(docs, records, donor_sites)

    try:
        assembled = gateway_assembly(list(records), reaction, multi_site_only=True)
    except ValueError as exc:
        raise CloningError(_pairing_message(reaction, docs, records)) from exc

    # The clone the user is after carries what every molecule puts into the product: the
    # insert of each entry clone *and* the backbone of the vector. Both tests are needed: a
    # hairpin vector carries the insert twice, so the insert alone does not tell them apart.
    regions = [_kept_region(record, reaction) for record in records]
    distinct: list[Dseqrecord] = []
    seen: set[str] = set()
    for product in assembled:      # pydna lists a molecule once per orientation it found it in
        key = product.seguid()
        if key not in seen:
            seen.add(key)
            distinct.append(product)
    # A region too short to probe (a linear vector cut close to its sites) cannot tell the
    # products apart, so the other molecules decide alone.
    carries = {id(product): all(_probe(region) is None or _carries(product, region)
                                for region in regions) for product in distinct}
    ranked = sorted(distinct, key=lambda p: (not carries[id(p)], not p.circular, -len(p)))
    if not ranked:
        raise CloningError(f"No {reaction} product uses both att sites")

    warnings: list[str] = []
    carriers = [p for p in ranked if carries[id(p)]]
    if not carriers:
        missing = next((doc.name for doc, region in zip(docs, regions)
                        if _probe(region) is not None
                        and not any(_carries(p, region) for p in ranked)), docs[0].name)
        warnings.append("No product carries the sequence between the att sites of "
                        f"“{missing}”: check that the molecules are the right way "
                        "round (first the inserts or entry clones, then the vector)")
    if len(carriers) > 1:
        warnings.append(f"{len(carriers)} products carry the insert, one per Gateway cassette "
                        f"of “{vector_doc.name}”: choose the one you want, or run the reaction "
                        "again on the product to fill the other cassette")
    if carriers and _has_counter_selection(carriers[0]):
        from_entry = next((doc.name for doc, record in zip(docs[:-1], records[:-1])
                           if _has_counter_selection(record)), None)
        if from_entry and not _has_counter_selection(records[-1]):
            warnings.append(f"The product still carries the ccdB counter-selection gene of "
                            f"“{from_entry}”, so it would not grow in a standard cloning "
                            "strain: clone your insert into that entry clone first")
        else:
            warnings.append("The product still carries a ccdB counter-selection gene, so it "
                            "would not grow in a standard cloning strain: a vector with a "
                            "second Gateway cassette needs that cassette replaced as well")
    warnings.append("Only products where both att sites recombined are shown; a product from a "
                    "single site is the co-integrate intermediate of the reaction")

    return CloningResult(products=_products(reaction, docs, records, ranked, carries),
                         warnings=warnings)


def _products(reaction: str, docs: Sequence[SequenceDoc], records: Sequence[Dseqrecord],
              ranked: Sequence[Dseqrecord], carries: dict[int, bool]) -> list[Product]:
    left, right = INPUT_LABELS[reaction]
    operation = OPERATIONS[reaction]
    entries, vector_doc = list(docs[:-1]), docs[-1]
    entry_names = " + ".join(doc.name for doc in entries)
    base_name = f"{vector_doc.name} + {entry_names}"
    entry_sites = ", ".join(", ".join(site_names(record)) for record in records[:-1])
    products: list[Product] = []
    seen = {True: 0, False: 0}
    for record in ranked:
        wanted = carries[id(record)]
        seen[wanted] += 1
        role = WANTED[reaction] if wanted else BY_PRODUCT[reaction]
        name = base_name if wanted else f"{base_name} {role}"
        if seen[wanted] > 1:
            name = f"{name}, product {seen[wanted]}"
        topology = "circular" if record.circular else "linear"
        inputs = (f"{entry_names} ({left})" if len(entries) == 1
                  else f"{len(entries)} entry clones ({entry_sites})")
        details = (f"{operation}: {inputs} x {vector_doc.name} ({right}) — "
                   f"{role.lower()}, {len(record)} bp")
        products.append(Product(
            doc=from_dseqrecord(record, name), operation=operation,
            parents=[doc.name for doc in docs], details=details,
            label=f"{role} — {len(record)} bp {topology}"))
    return products


class GatewayRequest(BaseModel):
    reaction: Reaction = "LR"
    inputs: list[SequenceDoc]
    donor_sites: list[int] | None = None     # BP: the att pair chosen for the donor vector


@router.post("/gateway")
def gateway_route(body: GatewayRequest) -> CloningResult:
    if len(body.inputs) < 2:
        raise CloningError("Choose both molecules for the reaction")
    return gateway(body.reaction, *body.inputs, donor_sites=body.donor_sites)


class GatewaySitesRequest(BaseModel):
    docs: list[dict[str, Any]]


def _record(doc: dict[str, Any]) -> Dseqrecord | None:
    sequence = doc.get("sequence")
    if not isinstance(sequence, str) or not sequence.strip():
        return None
    letters = "".join(base for base in sequence.upper() if "A" <= base <= "Z")
    return Dseqrecord(letters, circular=bool(doc.get("circular")))


def _doc_sites(doc: dict[str, Any]) -> list[str]:
    record = _record(doc)
    return site_names(record) if record is not None else []


def _doc_other_sites(doc: dict[str, Any]) -> list[str]:
    """Recombination sites of other systems, so the dialog can say why a molecule that looks
    like an attB clone is not offered: ``["phiC31 attB"]``."""
    record = _record(doc)
    return [f"phiC31 {name}" for name in phic31_sites(record)] if record is not None else []


@router.post("/gateway-sites")
def gateway_sites_route(body: GatewaySitesRequest) -> dict[str, list[list[str]]]:
    """The att sites of each document, so the Gateway dialogs can pick sensible inputs: entry
    clones (attL) and a destination vector (attR) for LR, an attB insert and a donor (attP)
    for BP. ``other`` lists att sites of other recombination systems found in the same
    document (phiC31), which Gateway cannot use."""
    return {"sites": [_doc_sites(doc) for doc in body.docs],
            "other": [_doc_other_sites(doc) for doc in body.docs]}
