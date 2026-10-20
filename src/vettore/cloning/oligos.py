"""Annealed oligos: the route used to clone gRNAs, shRNAs, linkers and short tags.

Two oligos are annealed into a short duplex with sticky ends and ligated into a vector cut so
that its ends match. Vettore can work in either direction:

- **anneal**: you have the two oligos — it shows the duplex and, if a vector is given, the
  product;
- **design**: you have the insert (a 20 nt protospacer, say) and the cut vector — it writes the
  two oligos to order, taking the overhangs from the vector itself rather than from a protocol,
  so it works with any vector and any type IIS enzyme.
"""

from __future__ import annotations

from typing import Literal

from Bio.Seq import Seq
from fastapi import APIRouter
from pydantic import BaseModel
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord

from ..tm import primer_tm
from .base import (
    CloningError,
    CloningResult,
    Oligo,
    Product,
    SequenceDoc,
    from_dseqrecord,
    to_dseqrecord,
)
from .pcr import clean_bases
from .restriction import Fragment, IncompatibleEnds, circularise, enzyme_objects, fragments, ligate

router = APIRouter()

MIN_OLIGO = 10
PHOSPHATE_REMINDER = ("Annealed oligos carry no 5′ phosphate: phosphorylate them with T4 PNK "
                      "(or order them phosphorylated), otherwise the ligation leaves nicks")


def _rc(bases: str) -> str:
    return str(Seq(bases).reverse_complement()).upper()


def _overhang(end: tuple[str, str]) -> str:
    return end[1].upper()


def _backbone(vector_doc: SequenceDoc, enzyme_names: list[str]) -> Fragment:
    record = to_dseqrecord(vector_doc)
    if not record.circular:
        raise CloningError(f"The vector “{vector_doc.name}” must be a circular plasmid")
    enzymes = enzyme_objects(enzyme_names)
    pieces = fragments(record, enzymes)
    backbone = max(pieces, key=lambda f: f.length)
    if len(pieces) > 2:
        raise CloningError(f"{' and '.join(str(e) for e in enzymes)} cut the vector "
                           f"{len(pieces)} times; the vector must be cut once or twice")
    return backbone


def duplex(top: str, bottom: str) -> Dseq:
    try:
        sequence = Dseq(top, bottom)
    except (ValueError, TypeError, IndexError) as exc:  # pydna's errors for strands that
        # share no complementary stretch
        raise CloningError("The two oligos do not anneal: no complementary stretch was found. "
                           "Check that the second oligo is the reverse complement of the "
                           f"first, overhangs aside ({exc})") from None
    paired = min(len(top), len(bottom)) - abs(sequence.ovhg)
    if paired < MIN_OLIGO:
        raise CloningError("The two oligos do not anneal: their complementary stretch is only "
                           f"{max(0, paired)} bp. Check that the second oligo is the reverse "
                           "complement of the first, overhangs aside")
    return sequence


def describe_ends(sequence: Dseq) -> str:
    left, right = sequence.five_prime_end(), sequence.three_prime_end()
    def text(end: tuple[str, str]) -> str:
        return "blunt" if end[0] == "blunt" else f"{end[0]} {end[1].upper()}"
    return f"{text(left)} / {text(right)}"


def anneal(top_name: str, top: str, bottom_name: str, bottom: str,
           vector_doc: SequenceDoc | None = None,
           enzyme_names: list[str] | None = None) -> CloningResult:
    top_sequence = clean_bases(top, f"Oligo “{top_name}”")
    bottom_sequence = clean_bases(bottom, f"Oligo “{bottom_name}”")
    sequence = duplex(top_sequence, bottom_sequence)
    warnings = [PHOSPHATE_REMINDER]
    oligos = [
        Oligo(name=top_name, sequence=top_sequence, tm=_safe_tm(top_sequence),
              annealing_length=len(top_sequence)),
        Oligo(name=bottom_name, sequence=bottom_sequence, tm=_safe_tm(bottom_sequence),
              annealing_length=len(bottom_sequence)),
    ]

    if vector_doc is None:
        insert = Dseqrecord(sequence, circular=False)
        product = Product(
            doc=from_dseqrecord(insert, f"{top_name} + {bottom_name}"),
            operation="Anneal oligos", parents=[],
            details=f"Annealed {len(top_sequence)} nt and {len(bottom_sequence)} nt oligos into "
                    f"a {len(sequence)} bp duplex with {describe_ends(sequence)} ends",
            oligos=oligos)
        return CloningResult(products=[product], warnings=warnings, oligos=oligos)

    backbone = _backbone(vector_doc, enzyme_names or [])
    insert = Dseqrecord(sequence, circular=False)
    try:
        product_record = circularise(ligate(backbone.record, insert))
    except IncompatibleEnds:
        needed_left = _rc(_overhang(backbone.right_end))
        needed_right = _rc(_overhang(backbone.left_end))
        raise CloningError(
            f"The duplex ends ({describe_ends(sequence)}) do not match the cut vector. "
            f"The insert needs a 5′ {needed_left or 'blunt'} overhang on the top strand and a "
            f"5′ {needed_right or 'blunt'} overhang on the bottom strand.") from None

    name = f"{vector_doc.name} + {top_name}"
    details = (f"Ligated the annealed {len(sequence)} bp duplex into {vector_doc.name} cut with "
               f"{backbone.right_enzyme}" +
               (f" and {backbone.left_enzyme}" if backbone.left_enzyme != backbone.right_enzyme
                else ""))
    return CloningResult(
        products=[Product(doc=from_dseqrecord(product_record, name), operation="Anneal oligos",
                          parents=[vector_doc.name], details=details, oligos=oligos)],
        warnings=warnings, oligos=oligos)


def design(insert: str, vector_doc: SequenceDoc, enzyme_names: list[str],
           name: str = "oligo") -> CloningResult:
    """Write the two oligos that put ``insert`` into the vector's cut ends."""
    core = clean_bases(insert, "The insert")
    if len(core) < 4:
        raise CloningError("The insert is too short")
    backbone = _backbone(vector_doc, enzyme_names)
    left = _rc(_overhang(backbone.right_end))     # 5′ overhang of the duplex, top strand
    right = _rc(_overhang(backbone.left_end))     # 5′ overhang of the duplex, bottom strand
    if not left and not right:
        raise CloningError("The vector has blunt ends, so annealed oligos cannot be oriented; "
                           "use a vector cut to leave overhangs")

    top = left + core
    bottom = right + _rc(core)
    result = anneal(f"{name}_top", top, f"{name}_bottom", bottom, vector_doc, enzyme_names)
    result.products[0].details = (
        f"Designed oligos to insert {len(core)} bp into {vector_doc.name} cut with "
        f"{' and '.join(enzyme_names)}; overhangs {left or 'blunt'} and {right or 'blunt'}")
    result.warnings.insert(0, "Anneal the two oligos (95 °C, then cool slowly) before ligating")
    return result


def _safe_tm(sequence: str) -> float | None:
    try:
        return primer_tm(sequence)
    except ValueError:
        return None


class AnnealRequest(BaseModel):
    mode: Literal["anneal", "design"] = "anneal"
    top: str = ""
    top_name: str = "Oligo 1"
    bottom: str = ""
    bottom_name: str = "Oligo 2"
    insert: str = ""
    insert_name: str = "insert"
    vector: SequenceDoc | None = None
    enzymes: list[str] = []


@router.post("/anneal-oligos")
def anneal_route(body: AnnealRequest) -> CloningResult:
    if body.mode == "design":
        if body.vector is None:
            raise CloningError("Choose the vector the oligos go into")
        return design(body.insert, body.vector, body.enzymes, body.insert_name)
    return anneal(body.top_name, body.top, body.bottom_name, body.bottom, body.vector,
                  body.enzymes)
