"""Overlap-directed assembly: Gibson, NEBuilder HiFi and In-Fusion.

All three kits join fragments that share terminal homology, so the *design* is the same and
only the recommended homology length differs (Gibson/HiFi 20–30 bp, In-Fusion 15 bp). The
enzymes chew back opposite strands (Gibson/HiFi: 5' exonuclease, In-Fusion: 3' exonuclease),
but the repaired product is the same double-stranded molecule, so one model covers all three.

Homology is placed on one primer per junction: each insert's forward primer carries the last
``overlap`` bases of the part before it, and the last insert's reverse primer carries the
reverse complement of the first ``overlap`` bases of the vector. A vector linearised by
restriction is used at full length (both strands repaired), which is what the kits produce.
"""

from __future__ import annotations

from typing import Literal

from Bio.Seq import Seq
from fastapi import APIRouter
from pydantic import BaseModel
from pydna.assembly2 import gibson_assembly
from pydna.dseqrecord import Dseqrecord

from .base import (
    CloningError,
    CloningResult,
    Oligo,
    Product,
    SequenceDoc,
    from_dseqrecord,
    to_dseqrecord,
)
from .pcr import design_amplicon
from .regions import check_region, features_in_region, region_sequence
from .restriction import enzyme_objects, fragments

router = APIRouter()

Method = Literal["gibson", "hifi", "infusion"]
DEFAULT_OVERLAP: dict[str, int] = {"gibson": 25, "hifi": 25, "infusion": 15}
METHOD_NAMES: dict[str, str] = {"gibson": "Gibson assembly",
                                "hifi": "NEBuilder HiFi assembly",
                                "infusion": "In-Fusion cloning"}
MIN_OVERLAP = 10
MAX_PRIMER_LENGTH = 60


class Linearize(BaseModel):
    enzymes: list[str] | None = None
    pcr: dict[str, int] | None = None


class VectorIn(BaseModel):
    doc: SequenceDoc
    linearize: Linearize


class InsertIn(BaseModel):
    doc: SequenceDoc
    start: int | None = None
    end: int | None = None


class Part:
    """A piece of the final molecule, with the features it carries."""

    def __init__(self, name: str, sequence: str, source: Dseqrecord, start: int, end: int,
                 amplified: bool):
        self.name = name
        self.sequence = sequence
        self.source = source
        self.start = start
        self.end = end
        self.amplified = amplified


def _linear_vector(vector: VectorIn) -> tuple[Part, str]:
    record = to_dseqrecord(vector.doc)
    if vector.linearize.enzymes:
        if not record.circular:
            raise CloningError(f"“{vector.doc.name}” is already linear; choose PCR instead")
        enzymes = enzyme_objects(vector.linearize.enzymes)
        frags = fragments(record, enzymes)
        backbone = max(frags, key=lambda f: f.length)
        piece = backbone.record
        sequence = str(piece.seq).upper()
        part = Part(vector.doc.name, sequence, piece, 0, len(piece) - 1, amplified=False)
        how = f"cut with {' and '.join(str(e) for e in enzymes)}"
    elif vector.linearize.pcr:
        start = vector.linearize.pcr.get("start")
        end = vector.linearize.pcr.get("end")
        if start is None or end is None:
            raise CloningError("Choose the vector region to amplify")
        check_region(record, start, end, "vector")
        sequence = region_sequence(str(record.seq).upper(), start, end)
        part = Part(vector.doc.name, sequence, record, start, end, amplified=True)
        how = f"amplified from {start + 1} to {end + 1}"
    else:
        raise CloningError("Choose how to linearise the vector")
    return part, how


def _insert_part(insert: InsertIn) -> Part:
    record = to_dseqrecord(insert.doc)
    sequence = str(record.seq).upper()
    start, end = insert.start, insert.end
    if start is None or end is None:
        if record.circular:
            raise CloningError(f"“{insert.doc.name}” is circular: choose the region to amplify")
        start, end = 0, len(sequence) - 1
    check_region(record, start, end, "insert")
    return Part(insert.doc.name, region_sequence(sequence, start, end), record, start, end,
                amplified=True)


def overlap_assembly(method: Method, vector: VectorIn, inserts: list[InsertIn],
                     overlap: int | None = None, target_tm: float = 60.0) -> CloningResult:
    if not inserts:
        raise CloningError("Add at least one insert")
    homology = overlap or DEFAULT_OVERLAP[method]
    if homology < MIN_OVERLAP:
        raise CloningError(f"The overlap must be at least {MIN_OVERLAP} bp")

    vector_part, how = _linear_vector(vector)
    insert_parts = [_insert_part(i) for i in inserts]
    parts = [vector_part, *insert_parts]
    for part in parts:
        if len(part.sequence) < 2 * homology:
            raise CloningError(f"“{part.name}” is only {len(part.sequence)} bp, too short for "
                               f"{homology} bp homology arms")

    # --- design the PCR products ------------------------------------------------------------
    designed = []
    warnings: list[str] = []
    for index, part in enumerate(insert_parts):
        upstream = parts[index]  # vector for the first insert, previous insert otherwise
        f_tail = upstream.sequence[-homology:]
        r_tail = (str(Seq(vector_part.sequence[:homology]).reverse_complement())
                  if index == len(insert_parts) - 1 else "")
        amplicon = design_amplicon(part.sequence, target_tm, f_tail, r_tail, name=part.name)
        designed.append(amplicon)
        for oligo in (amplicon.forward, amplicon.reverse):
            if len(oligo.sequence) > MAX_PRIMER_LENGTH:
                warnings.append(f"Primer {oligo.name} is {len(oligo.sequence)} nt; consider a "
                                "shorter homology arm")

    vector_amplicon = None
    if vector_part.amplified:
        vector_amplicon = design_amplicon(vector_part.sequence, target_tm, "", "",
                                          name=vector_part.name)

    # --- expected product -------------------------------------------------------------------
    expected = "".join(part.sequence for part in parts)
    product = Dseqrecord(expected, circular=True)
    offset = 0
    for part in parts:
        product.features.extend(features_in_region(part.source, part.start, part.end, offset))
        offset += len(part.sequence)

    # --- check the design with pydna --------------------------------------------------------
    pieces = [Dseqrecord(vector_part.sequence, circular=False)] + [
        Dseqrecord(a.sequence, circular=False) for a in designed]
    try:
        assembled = gibson_assembly(pieces, limit=homology, circular_only=True)
    except Exception as exc:  # pydna raises various errors on impossible assemblies
        raise CloningError(f"The fragments could not be assembled: {exc}") from exc
    matching = [a for a in assembled
                if a.seguid() == Dseqrecord(expected, circular=True).seguid()]
    if not matching:
        raise CloningError(
            "The designed fragments do not assemble into the expected plasmid; check the "
            "regions and the homology length (repeated sequences at the junctions can also "
            "prevent a unique assembly)")
    if len(assembled) > len(matching):
        warnings.append(f"pydna also finds {len(assembled) - len(matching)} alternative "
                        "assembly product(s); repeated sequences may give wrong clones")

    oligos: list[Oligo] = []
    if vector_amplicon:
        oligos += [vector_amplicon.forward, vector_amplicon.reverse]
    for amplicon in designed:
        oligos += [amplicon.forward, amplicon.reverse]

    name = f"{vector_part.name} + {' + '.join(p.name for p in insert_parts)}"
    details = (f"{METHOD_NAMES[method]} of {len(insert_parts)} insert(s) into "
               f"{vector_part.name} ({how}); {homology} bp homology arms")
    return CloningResult(
        products=[Product(doc=from_dseqrecord(product, name), operation=METHOD_NAMES[method],
                          parents=[p.name for p in parts], details=details, oligos=oligos)],
        warnings=warnings, oligos=oligos)


class OverlapRequest(BaseModel):
    method: Method = "gibson"
    vector: VectorIn
    inserts: list[InsertIn]
    overlap: int | None = None
    target_tm: float = 60.0


@router.post("/overlap-assembly")
def overlap_route(body: OverlapRequest) -> CloningResult:
    return overlap_assembly(body.method, body.vector, body.inserts, body.overlap,
                            body.target_tm)


from .mutagenesis import router as _mutagenesis_router
from .oligos import router as _oligos_router
from .pcr import router as _pcr_router

router.include_router(_pcr_router)
router.include_router(_mutagenesis_router)
router.include_router(_oligos_router)
