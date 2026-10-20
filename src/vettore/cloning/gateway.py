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

from typing import Literal

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


def _site_text(record: Dseqrecord, doc_name: str) -> str:
    found = _sites(record)
    if not found:
        return f"“{doc_name}”: no att site found"
    return f"“{doc_name}”: " + ", ".join(sorted(found))


def _arcs(record: Dseqrecord, kind: str) -> list[str]:
    """The pieces the molecule's att sites cut it into, shortest first."""
    located = [_bounds(loc)[:2] for name, locations in
               find_gateway_sites(record, False).items() if name.startswith(kind)
               for loc in locations if _bounds(loc)[2] <= MAX_SITE]
    if len(located) < 2:
        return []
    sequence = str(record.seq).upper()
    n = len(sequence)
    doubled = sequence * 2
    pieces = []
    for site in located:
        for other in located:
            if site is other:
                continue          # the whole molecule read past its own site is not a piece
            length = (other[0] - site[1]) % n
            if 0 < length < n:
                pieces.append(doubled[site[1]: site[1] + length])
    return sorted(pieces, key=len)


def _insert_region(record: Dseqrecord, kind: str) -> str | None:
    """What the reaction moves across: the shorter piece between the att sites — the insert is
    smaller than the backbone in every entry and expression clone in normal use."""
    arcs = _arcs(record, kind)
    return arcs[0] if arcs else None


def _backbone_region(record: Dseqrecord, kind: str) -> str | None:
    """The vector side: the longest piece outside the att cassette (origin, marker)."""
    arcs = _arcs(record, kind)
    return arcs[-1] if arcs else None


def _carries(product: Dseqrecord, region: str | None) -> bool:
    """Is that piece of a parent still in the product (either strand)?"""
    if not region or len(region) < 20:
        return False
    middle = max(0, (len(region) - 60) // 2)
    probe = region[middle: middle + 60]
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
                f"The {label} “{doc.name}” has no Gateway att site. A {reaction} reaction needs "
                f"{want_left} sites in the first molecule and {want_right} sites in the second.")
        if not any(name.startswith(kind) for name in _sites(record)):
            raise CloningError(
                f"The {label} “{doc.name}” has no {kind} site. A {reaction} reaction needs "
                f"{want_left} sites in the first molecule and {want_right} sites in the second. "
                f"Sites found — {_site_text(first, first_doc.name)}; "
                f"{_site_text(second, second_doc.name)}")

    try:
        assembled = gateway_assembly([first, second], reaction, multi_site_only=True)
    except ValueError as exc:
        raise CloningError(
            f"No {reaction} product: the att sites of the two molecules do not pair. "
            f"Sites found — {_site_text(first, first_doc.name)}; "
            f"{_site_text(second, second_doc.name)}. A {reaction} reaction pairs "
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
    carries = {id(p): _carries(p, region) and _carries(p, backbone) for p in distinct}
    ranked = sorted(distinct, key=lambda p: (not carries[id(p)], not p.circular, -len(p)))
    if not ranked:
        raise CloningError(f"No {reaction} product uses both att sites")

    warnings: list[str] = []
    carriers = [p for p in ranked if carries[id(p)]]
    if region and not carriers:
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
