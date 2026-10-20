"""Gateway BP and LR reactions. Routes under /api/cloning.

Site-specific recombination between att sites, run by pydna (``pydna.assembly2.gateway_assembly``
with the att consensus sequences of ``pydna.gateway``, which pydna takes from the published
site definitions; Vettore adds no att sequence of its own):

- **BP**: attB x attP -> attL + attR. An attB-flanked insert (a PCR product with attB tails, or
  a clone carrying attB sites) plus a donor vector with attP sites give the **entry clone**
  (attL) and an attR by-product.
- **LR**: attL x attR -> attB + attP. An entry clone plus a destination vector with attR sites
  give the **expression clone** (attB) and an attP by-product.

Only reactions where **both** att sites recombine are reported: a product from a single site is
the co-integrate intermediate of the reaction, not what grows on the selection plates.

Validated end to end against three public GenBank records that carry annotated att sites
(``tests/test_cloning_gateway.py``): the LR product of a published entry clone and a published
destination vector reproduces the published expression clone's attB cassette base for base, and
the BP reaction back from it regenerates the published entry clone exactly.
"""

from __future__ import annotations

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


def _bounds(location) -> tuple[int, int, int]:
    """(start, end, length) of a site along the top strand; sites may wrap the origin, and
    pydna reports those as a two-part location."""
    parts = list(location.parts)
    length = sum(len(part) for part in parts)
    if len(parts) == 1:
        return int(parts[0].start), int(parts[0].end), length
    ordered = sorted(parts, key=lambda part: int(part.start))
    return int(ordered[-1].start), int(ordered[0].end), length   # wraps: (last part, first part)


def _sites(record: Dseqrecord) -> dict[str, int]:
    """How many att sites of each kind are in the molecule (site-sized matches only)."""
    found: dict[str, int] = {}
    for name, locations in find_gateway_sites(record, False).items():
        count = sum(1 for loc in locations if _bounds(loc)[2] <= MAX_SITE)
        if count:
            found[name] = count
    return found


def _article(reaction: str) -> str:
    """"An LR reaction" (read "el-ar"), "a BP reaction"."""
    return "An" if reaction == "LR" else "A"


def site_names(record: Dseqrecord) -> list[str]:
    """The att sites of a molecule, by name, as the reaction sees them. attB is the core of
    every attL, attR and attP site, so an attB match inside a longer site is not listed: an
    entry clone reads attL1, attL2, not attB1, attB2, attL1, attL2."""
    found = [(name, *_bounds(loc)) for name, locations in find_gateway_sites(record, False).items()
             for loc in locations if _bounds(loc)[2] <= MAX_SITE]
    n = len(record.seq)

    def inside(start: int, length: int, other_start: int, other_length: int) -> bool:
        return (start - other_start) % n + length <= other_length

    names = {name for name, start, _, length in found
             if not (name.startswith("attB") and any(
                 other_length > length and inside(start, length, other_start, other_length)
                 for _, other_start, _, other_length in found))}
    return sorted(names)


def _site_text(record: Dseqrecord, doc_name: str) -> str:
    found = site_names(record)
    if not found:
        return f"“{doc_name}”: no att site found"
    return f"“{doc_name}”: " + ", ".join(found)


def _located(record: Dseqrecord, kind: str) -> list[tuple[int, int, int, int]]:
    """(start, end, length, strand) of every site of that kind, in order along the top
    strand."""
    found = [(*_bounds(loc), loc.strand or 1)
             for name, locations in find_gateway_sites(record, False).items()
             if name.startswith(kind)
             for loc in locations if _bounds(loc)[2] <= MAX_SITE]
    return sorted(found)


def _arcs(record: Dseqrecord, kind: str) -> tuple[list[str], list[str]]:
    """The pieces between neighbouring att sites, as (cassettes, outside), top strand.

    Every att site faces the DNA it brings into the reaction: a site on the top strand faces
    the DNA after it, a site on the bottom strand the DNA before it (attX1 on the top strand
    and attX2 on the bottom strand in a cassette read the usual way round). A piece between two
    sites that face each other is a cassette (the insert of an entry or expression clone, the
    ccdB cassette of a donor or destination vector); a piece between two sites that face away
    is the vector side. A linear molecule has no piece through its ends between the sites, but
    its two ends are on the vector side."""
    located = _located(record, kind)
    if len(located) < 2:
        return [], []
    sequence = str(record.seq).upper()
    n = len(sequence)
    doubled = sequence * 2
    pairs = list(pairwise(located))
    if record.circular:
        pairs.append((located[-1], located[0]))
    cassettes, outside = [], []
    for (_, end, size, strand), (start, _, next_size, next_strand) in pairs:
        length = (start - end) % n if record.circular else start - end
        if not 0 < length <= n - size - next_size:     # overlapping sites have nothing between
            continue
        piece = doubled[end: end + length]
        if strand == 1 and next_strand == -1:
            cassettes.append(piece)
        elif strand == -1 and next_strand == 1:
            outside.append(piece)
    if not record.circular:
        outside += [sequence[: located[0][0]], sequence[max(site[1] for site in located):]]
    return cassettes, outside


def _insert_region(record: Dseqrecord, kind: str) -> str | None:
    """What the reaction moves across: the piece between two att sites that face each other.
    The site orientation decides it, not the length: an insert may be longer than the
    backbone."""
    cassettes, _ = _arcs(record, kind)
    return max(cassettes, key=len) if cassettes else None


def _backbone_region(record: Dseqrecord, kind: str) -> str | None:
    """The vector side: the longest piece between two att sites that face away from each other
    (origin, marker), or an end of a linear molecule."""
    _, outside = _arcs(record, kind)
    return max(outside, key=len) if outside else None


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


def gateway(reaction: Reaction, first_doc: SequenceDoc, second_doc: SequenceDoc) -> CloningResult:
    first = to_dseqrecord(first_doc)
    second = to_dseqrecord(second_doc)
    left, right = INPUT_LABELS[reaction]
    want_left, want_right = EXPECTED_SITES[reaction]

    # attB is a part of attL and of attR, so only the vector-side site tells the two
    # reactions apart; a molecule with no att site at all is reported on its own.
    for record, doc, label, kind in ((first, first_doc, left, want_left),
                                     (second, second_doc, right, want_right)):
        if not _sites(record):
            raise CloningError(
                f"The {label} “{doc.name}” has no Gateway att site. {_article(reaction)} {reaction} reaction needs "
                f"{want_left} sites in the first molecule and {want_right} sites in the second.")
        if not any(name.startswith(kind) for name in _sites(record)):
            raise CloningError(
                f"The {label} “{doc.name}” has no {kind} site. {_article(reaction)} {reaction} reaction needs "
                f"{want_left} sites in the first molecule and {want_right} sites in the second. "
                f"Sites found — {_site_text(first, first_doc.name)}; "
                f"{_site_text(second, second_doc.name)}")

    try:
        assembled = gateway_assembly([first, second], reaction, multi_site_only=True)
    except ValueError as exc:
        raise CloningError(
            f"No {reaction} product: the att sites of the two molecules do not pair. "
            f"Sites found — {_site_text(first, first_doc.name)}; "
            f"{_site_text(second, second_doc.name)}. {_article(reaction)} {reaction} reaction pairs "
            f"{want_left}1 with {want_right}1 and {want_left}2 with {want_right}2.") from exc

    # The clone the user is after carries the insert of the first molecule *and* the backbone
    # of the second; the by-product is the other half of the swap. Both tests are needed: a
    # hairpin vector carries the insert twice, so the insert alone does not tell them apart.
    region = _insert_region(first, want_left)
    backbone = _backbone_region(second, want_right)
    distinct: list[Dseqrecord] = []
    seen: set[str] = set()
    for product in assembled:      # pydna lists a molecule once per orientation it found it in
        key = product.seguid()
        if key not in seen:
            seen.add(key)
            distinct.append(product)
    # A backbone too short to probe (a linear vector cut close to its sites) cannot tell the
    # products apart, so the insert decides alone.
    carries = {id(p): _carries(p, region)
               and (_probe(backbone) is None or _carries(p, backbone)) for p in distinct}
    ranked = sorted(distinct, key=lambda p: (not carries[id(p)], not p.circular, -len(p)))
    if not ranked:
        raise CloningError(f"No {reaction} product uses both att sites")

    warnings: list[str] = []
    carriers = [p for p in ranked if carries[id(p)]]
    if not carriers:
        warnings.append("No product carries the sequence between the att sites of "
                        f"“{first_doc.name}”: check that the two molecules are the right way "
                        "round (first the insert or entry clone, then the vector)")
    if len(carriers) > 1:
        warnings.append(f"{len(carriers)} products carry the insert, one per Gateway cassette "
                        f"of “{second_doc.name}”: choose the one you want, or run the reaction "
                        "again on the product to fill the other cassette")
    if carriers and _has_counter_selection(carriers[0]):
        warnings.append("The product still carries a ccdB counter-selection gene, so it would "
                        "not grow in a standard cloning strain: a vector with a second Gateway "
                        "cassette needs that cassette replaced as well")
    warnings.append("Only products where both att sites recombined are shown; a product from a "
                    "single site is the co-integrate intermediate of the reaction")

    operation = OPERATIONS[reaction]
    base_name = f"{second_doc.name} + {first_doc.name}"
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
        details = (f"{operation}: {first_doc.name} ({left}) x {second_doc.name} ({right}) — "
                   f"{role.lower()}, {len(record)} bp")
        products.append(Product(
            doc=from_dseqrecord(record, name), operation=operation,
            parents=[first_doc.name, second_doc.name], details=details,
            label=f"{role} — {len(record)} bp {topology}"))
    return CloningResult(products=products, warnings=warnings)


class GatewayRequest(BaseModel):
    reaction: Reaction = "LR"
    inputs: list[SequenceDoc]


@router.post("/gateway")
def gateway_route(body: GatewayRequest) -> CloningResult:
    if len(body.inputs) != 2:
        raise CloningError("Choose both molecules for the reaction")
    return gateway(body.reaction, body.inputs[0], body.inputs[1])


class GatewaySitesRequest(BaseModel):
    docs: list[dict[str, Any]]


def _doc_sites(doc: dict[str, Any]) -> list[str]:
    sequence = doc.get("sequence")
    if not isinstance(sequence, str) or not sequence.strip():
        return []
    letters = "".join(base for base in sequence.upper() if "A" <= base <= "Z")
    return site_names(Dseqrecord(letters, circular=bool(doc.get("circular"))))


@router.post("/gateway-sites")
def gateway_sites_route(body: GatewaySitesRequest) -> dict[str, list[list[str]]]:
    """The att sites of each document, so the Gateway dialogs can pick sensible inputs: an
    entry clone (attL) and a destination vector (attR) for LR, an attB insert and a donor
    (attP) for BP."""
    return {"sites": [_doc_sites(doc) for doc in body.docs]}
