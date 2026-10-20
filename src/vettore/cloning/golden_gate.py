"""Golden Gate assembly with type IIS enzymes.

Each part carries two outward-facing sites of the chosen enzyme, so cutting releases a
fragment whose ends are the 4 nt (usually) fusion overhangs. Parts are ligated in the order
given; the junction overhangs must match pairwise, and a set of overhangs that repeats itself
(or is palindromic) gives wrong clones on the bench — both are reported.

The product is also checked against pydna's own Golden Gate assembly.
"""

from __future__ import annotations

from Bio.Seq import Seq
from fastapi import APIRouter
from pydantic import BaseModel
from pydna.assembly2 import golden_gate_assembly

from .base import (
    CloningError,
    CloningResult,
    Product,
    SequenceDoc,
    from_dseqrecord,
    plural,
    to_dseqrecord,
)
from .restriction import (
    Fragment,
    IncompatibleEnds,
    circularise,
    cut_counts,
    enzyme_objects,
    fragments,
    ligate,
)

router = APIRouter()


def _overhang(end: tuple[str, str]) -> str:
    return end[1].upper()


def _rc(bases: str) -> str:
    return str(Seq(bases).reverse_complement()).upper()


def _is_palindromic(overhang: str) -> bool:
    return bool(overhang) and overhang == str(Seq(overhang).reverse_complement()).upper()


def _pick_fragment(frags: list[Fragment], enzyme_name: str, what: str, name: str) -> Fragment:
    """The piece between the two outward-facing sites: both ends cut by the enzyme."""
    flanked = [f for f in frags if f.left_enzyme == enzyme_name and f.right_enzyme == enzyme_name]
    if not flanked:
        raise CloningError(f"The {what} “{name}” has no fragment flanked by two {enzyme_name} "
                           "sites; add outward-facing sites or domesticate the part")
    if what == "vector":
        return max(flanked, key=lambda f: f.length)
    return max(flanked, key=lambda f: f.length)


def golden_gate(vector_doc: SequenceDoc, part_docs: list[SequenceDoc],
                enzyme_name: str = "BsaI") -> CloningResult:
    if not part_docs:
        raise CloningError("Add at least one part")
    [enzyme] = enzyme_objects([enzyme_name])
    name = str(enzyme)
    warnings: list[str] = []

    vector = to_dseqrecord(vector_doc)
    vector_frags = fragments(vector, [enzyme])
    backbone = _pick_fragment(vector_frags, name, "vector", vector_doc.name)
    dropout = [f for f in vector_frags if f is not backbone]
    if dropout:
        sizes = ", ".join(f"{f.length} bp" for f in dropout)
        warnings.append(f"The vector releases a {sizes} dropout fragment; gel-purify the "
                        "backbone or rely on the dropout marker for screening")

    pieces: list[Fragment] = []
    for doc in part_docs:
        record = to_dseqrecord(doc)
        count = cut_counts(record, [enzyme])[name]
        if count == 0:
            raise CloningError(f"{name} does not cut the part “{doc.name}”")
        if count > 2:
            raise CloningError(
                f"{name} cuts the part “{doc.name}” {count} times: it has internal sites and "
                "must be domesticated (remove the internal sites) before assembly")
        pieces.append(_pick_fragment(fragments(record, [enzyme]), name, "part", doc.name))

    # --- junction check ---------------------------------------------------------------------
    order = [backbone, *pieces]
    labels = [vector_doc.name, *[d.name for d in part_docs]]
    overhangs: list[str] = []
    for index, piece in enumerate(order):
        left = piece
        right = order[(index + 1) % len(order)]
        # pydna reports a 3' (right) end on the opposite strand, so compare its complement.
        donor, acceptor = _rc(_overhang(left.right_end)), _overhang(right.left_end)
        if donor != acceptor:
            raise CloningError(
                f"The junction between {labels[index]} and {labels[(index + 1) % len(labels)]} "
                f"does not match: overhang {donor or 'blunt'} against {acceptor or 'blunt'}")
        if not donor:
            raise CloningError(f"The junction between {labels[index]} and "
                               f"{labels[(index + 1) % len(labels)]} is blunt; Golden Gate needs "
                               "sticky overhangs")
        overhangs.append(donor)
        del piece

    duplicates = {o for o in overhangs if overhangs.count(o) > 1}
    if duplicates:
        warnings.append("Duplicate junction overhangs (" + ", ".join(sorted(duplicates)) +
                        "): parts can assemble in the wrong order")
    palindromic = {o for o in overhangs if _is_palindromic(o)}
    if palindromic:
        warnings.append("Palindromic overhangs (" + ", ".join(sorted(palindromic)) +
                        "): parts can insert in either orientation")

    # --- assemble -----------------------------------------------------------------------------
    try:
        linear = order[0].record
        for piece in order[1:]:
            linear = ligate(linear, piece.record)
        product = circularise(linear)
    except IncompatibleEnds as exc:
        raise CloningError(f"The parts do not ligate in this order: {exc}") from None

    cross_check = golden_gate_assembly([f.record for f in order], [enzyme], allow_blunt=False,
                                       circular_only=True)
    if not any(a.seguid() == product.seguid() for a in cross_check):
        warnings.append("pydna's Golden Gate assembly does not reproduce this product; check "
                        "the parts before using it")
    elif len(cross_check) > 1:
        warnings.append(f"{len(cross_check)} different circular products are possible with "
                        "these overhangs")

    product_name = f"{vector_doc.name} + {' + '.join(d.name for d in part_docs)}"
    details = (f"Golden Gate assembly of {plural(len(part_docs), 'part')} into {vector_doc.name} with "
               f"{name}; junction overhangs {', '.join(overhangs)}")
    return CloningResult(
        products=[Product(doc=from_dseqrecord(product, product_name),
                          operation="Golden Gate assembly",
                          parents=[vector_doc.name, *[d.name for d in part_docs]],
                          details=details)],
        warnings=warnings)


class GoldenGateRequest(BaseModel):
    vector: SequenceDoc
    parts: list[SequenceDoc]
    enzyme: str = "BsaI"


@router.post("/golden-gate")
def golden_gate_route(body: GoldenGateRequest) -> CloningResult:
    return golden_gate(body.vector, body.parts, body.enzyme)
