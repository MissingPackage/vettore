"""Restriction cloning: digests, "Insert Fragment" and "Delete Restriction Fragment".

Biology runs on pydna (``Dseqrecord.cut``, sticky-end aware ``+`` and ``looped``). Coordinates
returned to the UI follow OVE: 0-based, end-inclusive, along the top strand, with
``start > end`` for a fragment that spans the origin of a circular molecule. A fragment's
``length`` is measured between the two top-strand cuts (overhangs are reported separately).

Products keep the vector's numbering when the vector origin lies in the backbone, so a
construct keeps the coordinates people already know for its parent vector.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from Bio import Restriction
from fastapi import APIRouter
from pydantic import BaseModel
from pydna.dseqrecord import Dseqrecord

from ..enzymes import build_database
from .base import (
    CloningError,
    CloningResult,
    Product,
    SequenceDoc,
    from_dseqrecord,
    plural,
    times,
    to_dseqrecord,
)

router = APIRouter()

_DATABASE: dict | None = None


def enzyme_objects(names: list[str]) -> list:
    """Bio.Restriction enzymes for names offered by the UI (case-insensitive, deduplicated)."""
    global _DATABASE
    if _DATABASE is None:
        _DATABASE = build_database()
    if not names:
        raise CloningError("Choose at least one enzyme")
    out = []
    for raw in names:
        entry = _DATABASE.get(raw.strip().lower())
        if entry is None:
            raise CloningError(f"Unknown or unsupported enzyme: {raw}")
        enzyme = getattr(Restriction, entry["name"])
        if enzyme not in out:
            out.append(enzyme)
    return out


def _end_text(end: tuple[str, str]) -> str:
    kind, bases = end
    if kind == "blunt":
        return "blunt"
    return f"{kind} overhang {bases.upper()}"


@dataclass
class Fragment:
    index: int
    record: Dseqrecord
    start: int                 # top-strand coordinate of the first base, 0-based
    end: int                   # inclusive
    length: int                # between the two top-strand cuts
    left_enzyme: str | None    # None = end of a linear molecule
    right_enzyme: str | None

    @property
    def left_end(self) -> tuple[str, str]:
        return self.record.seq.five_prime_end()

    @property
    def right_end(self) -> tuple[str, str]:
        return self.record.seq.three_prime_end()

    def ends_cut_by(self, enzyme: str) -> int:
        """How many of this fragment's two ends that enzyme made (0, 1 or 2)."""
        return (self.left_enzyme == enzyme) + (self.right_enzyme == enzyme)

    def label(self) -> str:
        left = self.left_enzyme or "end"
        right = self.right_enzyme or "end"
        return f"{left}–{right} ({self.length} bp)"

    def describe(self) -> dict:
        return {
            "index": self.index,
            "start": self.start,
            "end": self.end,
            "length": self.length,
            "left_enzyme": self.left_enzyme,
            "right_enzyme": self.right_enzyme,
            "left_end": "end of molecule" if self.left_enzyme is None else _end_text(self.left_end),
            "right_end": "end of molecule" if self.right_enzyme is None else _end_text(self.right_end),
            "label": self.label(),
        }


def cut_counts(record: Dseqrecord, enzymes: list) -> dict[str, int]:
    return {str(e): len(record.seq.get_cutsites(e)) for e in enzymes}


def fragments(record: Dseqrecord, enzymes: list) -> list[Fragment]:
    """Complete digest of ``record`` with ``enzymes``, fragments in pydna's order."""
    cuts = sorted(record.seq.get_cutsites(*enzymes), key=lambda c: c[0][0])
    if not cuts:
        names = ", ".join(str(e) for e in enzymes)
        raise CloningError(f"{names} {'does' if len(enzymes) == 1 else 'do'} not cut "
                           f"“{record.name}”")
    pieces = record.cut(*enzymes)
    n = len(record)
    out: list[Fragment] = []
    if record.circular:
        # pydna orders circular fragments from the first cut: fragment j spans cut j -> cut j+1.
        for j, piece in enumerate(pieces):
            (pos_a, _), enz_a = cuts[j]
            (pos_b, _), enz_b = cuts[(j + 1) % len(cuts)]
            length = (pos_b - pos_a) % n or n
            out.append(Fragment(j, piece, pos_a % n, (pos_b - 1) % n, length, str(enz_a),
                                str(enz_b)))
    else:
        bounds = [(0, None)] + [(c[0][0], str(c[1])) for c in cuts] + [(n, None)]
        for j, piece in enumerate(pieces):
            (pos_a, enz_a), (pos_b, enz_b) = bounds[j], bounds[j + 1]
            out.append(Fragment(j, piece, pos_a, pos_b - 1, pos_b - pos_a, enz_a, enz_b))
    if len(out) != len(pieces):
        raise CloningError("Unexpected digest result")  # defensive: pydna and cut list disagree
    return out


def _rotate_to_origin(product: Dseqrecord, vector: Dseqrecord) -> Dseqrecord:
    """Rotate a circular product so that numbering starts at the vector's origin, when the
    vector's first bases survive in the product, so familiar coordinates are preserved."""
    probe_len = min(30, len(vector))
    probe = str(vector.seq)[:probe_len].upper()
    doubled = (str(product.seq) * 2).upper()
    pos = doubled.find(probe)
    if pos < 0 or pos >= len(product):
        return product
    if doubled.find(probe, pos + 1) not in (-1, pos + len(product)):
        return product  # ambiguous, leave pydna's numbering
    return product.shifted(pos) if pos else product


def _pick_backbone(frags: list[Fragment], enzymes: list) -> Fragment:
    names = {str(e) for e in enzymes}
    candidates = [f for f in frags if {f.left_enzyme, f.right_enzyme} <= names]
    if not candidates:
        raise CloningError("No vector fragment is flanked by the chosen enzymes")
    return max(candidates, key=lambda f: f.length)


def _pick_insert(frags: list[Fragment], enzymes: list, circular: bool,
                 index: int | None) -> Fragment:
    if index is not None:
        if not 0 <= index < len(frags):
            raise CloningError(f"Insert fragment {index + 1} does not exist "
                               f"(the digest gives {plural(len(frags), 'fragment')})")
        return frags[index]
    names = {str(e) for e in enzymes}
    flanked = [f for f in frags if f.left_enzyme and f.right_enzyme
               and {f.left_enzyme, f.right_enzyme} == names]
    if not flanked and len(names) == 1:
        flanked = [f for f in frags if f.left_enzyme and f.right_enzyme]
    if not flanked:
        raise CloningError("No insert fragment is flanked by both chosen enzymes")
    # From a donor plasmid the insert is usually the smaller piece; from a PCR product the
    # fragment between the two sites is the only internal one.
    return min(flanked, key=lambda f: f.length) if circular else max(flanked, key=lambda f: f.length)


def ligate(left: Dseqrecord, right: Dseqrecord) -> Dseqrecord:
    try:
        return left + right
    except TypeError as exc:
        raise IncompatibleEnds(str(exc)) from exc


class IncompatibleEnds(Exception):
    pass


def circularise(linear: Dseqrecord) -> Dseqrecord:
    try:
        return linear.looped()
    except TypeError as exc:
        raise IncompatibleEnds(str(exc)) from exc


def _join(backbone: Dseqrecord, insert: Dseqrecord) -> Dseqrecord:
    return circularise(ligate(backbone, insert))


def can_join(backbone: Dseqrecord, insert: Dseqrecord) -> bool:
    try:
        _join(backbone, insert)
        return True
    except IncompatibleEnds:
        return False


# ---- operations ----------------------------------------------------------------------------

def digest(doc: SequenceDoc, enzyme_names: list[str]) -> list[dict]:
    record = to_dseqrecord(doc)
    return [f.describe() for f in fragments(record, enzyme_objects(enzyme_names))]


def insert_fragment(vector_doc: SequenceDoc, vector_enzymes: list[str], insert_doc: SequenceDoc,
                    insert_enzymes: list[str], insert_fragment_index: int | None = None,
                    orientation: Literal["forward", "reverse"] = "forward") -> CloningResult:
    vector = to_dseqrecord(vector_doc)
    insert = to_dseqrecord(insert_doc)
    if not vector.circular:
        raise CloningError(f"The vector “{vector_doc.name}” must be a circular plasmid")
    v_enz = enzyme_objects(vector_enzymes)
    i_enz = enzyme_objects(insert_enzymes)
    if len(v_enz) > 2 or len(i_enz) > 2:
        raise CloningError("Use one or two enzymes for each molecule")

    warnings: list[str] = []
    for name, count in cut_counts(vector, v_enz).items():
        if count == 0:
            raise CloningError(f"{name} does not cut the vector “{vector_doc.name}”")
        if count > 1:
            warnings.append(f"{name} cuts the vector {count} times; the backbone used is the "
                            "largest fragment between the chosen sites")
    for name, count in cut_counts(insert, i_enz).items():
        if count == 0:
            raise CloningError(f"{name} does not cut the insert source “{insert_doc.name}”")

    v_frags = fragments(vector, v_enz)
    backbone = _pick_backbone(v_frags, v_enz)
    i_frags = fragments(insert, i_enz)
    chosen = _pick_insert(i_frags, i_enz, insert.circular, insert_fragment_index)
    for name, count in cut_counts(insert, i_enz).items():
        # One cut per end of the insert is expected: a PCR product with an EcoRI site at each
        # end is cut twice by EcoRI. More cuts mean a site inside the source as well.
        if count > chosen.ends_cut_by(name):
            warnings.append(f"{name} cuts “{insert_doc.name}” {times(count)}; check that the "
                            f"insert fragment ({chosen.label()}) is the one you want")

    piece = chosen.record if orientation == "forward" else chosen.record.reverse_complement()
    other = chosen.record.reverse_complement() if orientation == "forward" else chosen.record
    try:
        product = _join(backbone.record, piece)
    except IncompatibleEnds:
        hint = (" The ends are compatible in the other orientation." if can_join(backbone.record, other)
                else "")
        raise CloningError(
            f"The ends do not match: vector backbone ends are {_end_text(backbone.right_end)} / "
            f"{_end_text(backbone.left_end)}, insert ends ({orientation}) are "
            f"{_end_text(piece.seq.five_prime_end())} / {_end_text(piece.seq.three_prime_end())}."
            + hint) from None

    if can_join(backbone.record, other):
        warnings.append("Non-directional cloning: the insert can ligate in both orientations; "
                        "screen colonies for the orientation you want")
    try:
        circularise(backbone.record)
        warnings.append("The vector ends are compatible with each other: the backbone can "
                        "re-ligate without insert (dephosphorylate the vector)")
    except IncompatibleEnds:
        pass
    if "blunt" in (backbone.left_end[0], backbone.right_end[0]):
        warnings.append("Blunt-end ligation is less efficient than sticky-end ligation")

    product = _rotate_to_origin(product, vector)
    name = f"{vector_doc.name} + {insert_doc.name}"
    flip = " (reverse orientation)" if orientation == "reverse" else ""
    details = (f"Inserted the {chosen.length} bp {chosen.left_enzyme}–{chosen.right_enzyme} "
               f"fragment of {insert_doc.name}{flip} into the {backbone.right_enzyme}–"
               f"{backbone.left_enzyme} sites of {vector_doc.name}")
    return CloningResult(
        products=[Product(doc=from_dseqrecord(product, name), operation="Restriction cloning",
                          parents=[vector_doc.name, insert_doc.name], details=details)],
        warnings=warnings,
    )


def delete_fragment(doc: SequenceDoc, enzyme_names: list[str],
                    fragment_index: int | None = None) -> CloningResult:
    record = to_dseqrecord(doc)
    enzymes = enzyme_objects(enzyme_names)
    frags = fragments(record, enzymes)
    if record.circular:
        if len(frags) != 2:
            raise CloningError(f"The enzymes cut {times(len(frags))}; choose enzymes that cut "
                               "exactly twice in total")
        removed = frags[fragment_index] if fragment_index is not None else min(frags, key=lambda f: f.length)
        kept = next(f for f in frags if f is not removed)
        try:
            product = circularise(kept.record)
        except IncompatibleEnds:
            raise CloningError(
                f"Incompatible ends ({_end_text(kept.right_end)} / {_end_text(kept.left_end)}) — "
                "blunting is not supported yet") from None
        product = _rotate_to_origin(product, record)
    else:
        internal = [f for f in frags if f.left_enzyme and f.right_enzyme]
        if not internal:
            raise CloningError("No internal fragment between two sites to delete")
        removed = frags[fragment_index] if fragment_index is not None else min(internal, key=lambda f: f.length)
        if removed.index in (0, len(frags) - 1):
            raise CloningError("Choose a fragment between two sites, not at an end")
        left = frags[removed.index - 1].record
        right = frags[removed.index + 1].record
        if len(frags) > 3:
            raise CloningError("The enzymes cut more than twice; choose enzymes that flank only "
                               "the fragment to delete")
        try:
            product = ligate(left, right)
        except IncompatibleEnds:
            raise CloningError("Incompatible ends — blunting is not supported yet") from None
    details = f"Deleted the {removed.label()} fragment of {doc.name}"
    return CloningResult(products=[Product(
        doc=from_dseqrecord(product, f"{doc.name} Δ{removed.left_enzyme}–{removed.right_enzyme}"),
        operation="Delete restriction fragment", parents=[doc.name], details=details)])


# ---- HTTP ----------------------------------------------------------------------------------

class DigestRequest(BaseModel):
    doc: SequenceDoc
    enzymes: list[str]


class InsertRequest(BaseModel):
    vector: SequenceDoc
    vector_enzymes: list[str]
    insert: SequenceDoc
    insert_enzymes: list[str]
    insert_fragment_index: int | None = None
    orientation: Literal["forward", "reverse"] = "forward"


class DeleteRequest(BaseModel):
    doc: SequenceDoc
    enzymes: list[str]
    fragment_index: int | None = None


@router.post("/digest")
def digest_route(body: DigestRequest) -> dict:
    return {"fragments": digest(body.doc, body.enzymes)}


@router.post("/restriction-insert")
def insert_route(body: InsertRequest) -> CloningResult:
    return insert_fragment(body.vector, body.vector_enzymes, body.insert, body.insert_enzymes,
                           body.insert_fragment_index, body.orientation)


@router.post("/delete-fragment")
def delete_route(body: DeleteRequest) -> CloningResult:
    return delete_fragment(body.doc, body.enzymes, body.fragment_index)
