"""Ligation of several fragments in one reaction. Routes under /api/cloning.

Two actions share one engine:

- **Insert Multiple Fragments**: a vector opened with one or two enzymes plus two or more
  inserts, each cut with its own enzymes (or used as a blunt PCR product), ligated into a
  plasmid;
- **Ligate Fragments**: two or more fragments joined into a linear molecule (or closed into a
  circle on request).

A fragment is either one piece of a restriction digest of an open document (its ends are what
the enzymes leave) or a region of a document taken as a PCR product made with a proofreading
polymerase (blunt ends, no 5′ phosphate).

Every product T4 DNA ligase can make from *all* the fragments, each used once, is found by
pydna's ``assembly2`` machinery with the same algorithm as ``pydna.assembly2.ligation_assembly``:
sticky ends join only when their overhangs pair completely, blunt ends join any blunt end. The
order and orientation of each fragment in a product follow from the ends alone; the order typed
in the form only decides which product is listed first. Molecules that use a subset of the
fragments (a vector closing on itself, insert concatemers) are not products, but the ones that
matter on the bench are reported as warnings.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise

from fastapi import APIRouter
from pydantic import BaseModel
from pydna.assembly2 import (
    Assembly,
    assemble,
    blunt_overlap,
    combine_algorithms,
    sticky_end_sub_strings,
)
from pydna.dseqrecord import Dseqrecord

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
from .regions import check_region, features_in_region, region_length, region_sequence
from .restriction import (
    Fragment,
    IncompatibleEnds,
    _end_text,
    _pick_backbone,
    _pick_insert,
    _rotate_to_origin,
    circularise,
    cut_counts,
    enzyme_objects,
    fragments,
)

router = APIRouter()

MAX_PRODUCTS = 50   # pydna refuses to enumerate more; past this the reaction is a random mix

PCR_PHOSPHATE_NOTE = (
    "PCR products carry no 5′ phosphate unless the primers were ordered phosphorylated: treat "
    "them with T4 polynucleotide kinase before ligation, or the joints next to them will not "
    "seal (essential if the vector is dephosphorylated)")
BLUNT_NOTE = ("Blunt ends ligate less efficiently than sticky ends: use more insert (3:1 to "
              "5:1 insert:vector) and a ligase buffer with PEG, or allow a longer ligation")


class FragmentIn(BaseModel):
    """One fragment of the reaction.

    With ``enzymes``: a piece of the digest of ``doc`` (``fragment_index`` picks it, counted
    from the digest list; ``None`` picks it automatically). Without enzymes: the region
    ``start``..``end`` (0-based, inclusive) of ``doc`` as a blunt PCR product, or the whole
    document if it is linear and no region is given.
    """

    doc: SequenceDoc
    enzymes: list[str] = []
    fragment_index: int | None = None
    start: int | None = None
    end: int | None = None


@dataclass
class Piece:
    """A fragment ready for ligation."""

    name: str              # document name (a History parent)
    label: str             # what the fragment is, for messages
    record: Dseqrecord
    from_pcr: bool
    digest: Fragment | None = None

    @property
    def left(self) -> str:
        return _end_text(self.record.seq.five_prime_end())

    @property
    def right(self) -> str:
        return _end_text(self.record.seq.three_prime_end())


def _digest_piece(item: FragmentIn, role: str, warnings: list[str]) -> Piece:
    record = to_dseqrecord(item.doc)
    enzymes = enzyme_objects(item.enzymes)
    if len(enzymes) > 2:
        raise CloningError(f"Use one or two enzymes for “{item.doc.name}”")
    counts = cut_counts(record, enzymes)
    for name, count in counts.items():
        if count == 0:
            raise CloningError(f"{name} does not cut “{item.doc.name}”; choose another enzyme "
                               "or use the fragment as a PCR product")
    frags = fragments(record, enzymes)
    if item.fragment_index is not None:
        if not 0 <= item.fragment_index < len(frags):
            raise CloningError(f"Fragment {item.fragment_index + 1} of “{item.doc.name}” does "
                               f"not exist (the digest gives {plural(len(frags), 'fragment')})")
        chosen = frags[item.fragment_index]
    elif role == "vector":
        chosen = _pick_backbone(frags, enzymes)
    else:
        chosen = _pick_insert(frags, enzymes, record.circular, None)
    for name, count in counts.items():
        # One cut per end of the fragment is expected (EcoRI at both ends of a PCR product cuts
        # twice); more cuts mean a site inside the molecule as well.
        if count > chosen.ends_cut_by(name) and item.fragment_index is None:
            warnings.append(f"{name} cuts “{item.doc.name}” {times(count)}; check that the "
                            f"fragment used ({chosen.label()}) is the one you want")
    label = f"{item.doc.name} {chosen.label()}"
    return Piece(item.doc.name, label, chosen.record, from_pcr=False, digest=chosen)


def _pcr_piece(item: FragmentIn) -> Piece:
    record = to_dseqrecord(item.doc)
    start, end = item.start, item.end
    if start is None or end is None:
        if record.circular:
            raise CloningError(f"“{item.doc.name}” is circular: choose enzymes to cut it, or "
                               "the region to use as a PCR product")
        start, end = 0, len(record) - 1
    check_region(record, start, end, "fragment")
    sequence = region_sequence(str(record.seq), start, end)
    piece = Dseqrecord(sequence, circular=False)
    piece.features.extend(features_in_region(record, start, end, 0))
    whole = start == 0 and end == len(record) - 1 and not record.circular
    where = "" if whole else f" {start + 1}..{end + 1}"
    label = (f"{item.doc.name}{where} (blunt PCR product, "
             f"{region_length(len(record), start, end)} bp)")
    return Piece(item.doc.name, label, piece, from_pcr=True)


def make_piece(item: FragmentIn, role: str, warnings: list[str]) -> Piece:
    if item.enzymes:
        return _digest_piece(item, role, warnings)
    return _pcr_piece(item)


# ---- the ligation itself -------------------------------------------------------------------

def _ligase(seqx: Dseqrecord, seqy: Dseqrecord, _limit):
    """The pairing rule of ``pydna.assembly2.ligation_assembly(allow_blunt=True)``."""
    return combine_algorithms(lambda x, y, _l: sticky_end_sub_strings(x, y, False),
                              blunt_overlap)(seqx, seqy, _limit)


@dataclass
class Ligated:
    record: Dseqrecord
    order: list[int]       # 1-based fragment numbers, negative = reverse complement

    @property
    def flipped(self) -> int:
        return sum(1 for n in self.order if n < 0)


def ligation_products(records: list[Dseqrecord], circular: bool) -> list[Ligated]:
    """Every distinct molecule that uses each fragment exactly once."""
    try:
        if len(records) == 1:
            if not circular:
                return []
            try:
                return [Ligated(circularise(records[0]), [1])]
            except IncompatibleEnds:
                return []
        asm = Assembly(records, None, _ligase, use_fragment_order=False, use_all_fragments=True)
        if circular:
            found = asm.get_circular_assemblies(max_assemblies=MAX_PRODUCTS)
            orders = [[edge[0] for edge in a] for a in found]
        else:
            found = asm.get_linear_assemblies(max_assemblies=MAX_PRODUCTS)
            orders = [[a[0][0]] + [edge[1] for edge in a] for a in found]
    except ValueError as exc:
        if "Too many" not in str(exc):  # pydna's limits on assemblies and on graph paths
            raise
        raise CloningError(
            f"These fragments can join in more than {MAX_PRODUCTS} different ways (their ends "
            "are blunt or identical), so a single reaction would give a random mix. Ligate "
            "them in steps, or use enzymes that leave different ends") from None
    out: list[Ligated] = []
    seen: set[str] = set()
    for a, order in zip(found, orders):
        record = assemble(records, a)
        key = record.seguid()
        if key not in seen:
            seen.add(key)
            out.append(Ligated(record, order))
    return out


def _expected_order(n: int, circular: bool) -> list[list[int]]:
    """Orders that mean "as listed, all forward" (a linear molecule reads either way)."""
    forward = list(range(1, n + 1))
    return [forward] if circular else [forward, [-k for k in reversed(forward)]]


def _sort(products: list[Ligated], n: int, circular: bool) -> list[Ligated]:
    expected = _expected_order(n, circular)

    def key(p: Ligated):
        flipped = p.flipped if circular else min(p.flipped, n - p.flipped)
        return (p.order not in expected, flipped, [abs(k) for k in p.order])

    return sorted(products, key=key)


def _order_text(order: list[int], pieces: list[Piece]) -> str:
    def one(k: int) -> str:
        name = pieces[abs(k) - 1].name
        return f"{name} (reversed)" if k < 0 else name
    return " → ".join(one(k) for k in order)


def _junction_report(pieces: list[Piece], circular: bool) -> str:
    """Why the fragments do not join in the order listed."""
    lines = []
    pairs = list(pairwise(pieces))
    if circular:
        pairs.append((pieces[-1], pieces[0]))
    for left, right in pairs:
        ok = bool(_ligase(left.record, right.record, None))
        if not ok:
            lines.append(f"“{left.name}” ends with {left.right} but “{right.name}” starts with "
                         f"{right.left}")
    if not lines:
        return ""
    return " Ends that do not match in the order listed: " + "; ".join(lines) + "."


def _describe_pieces(pieces: list[Piece]) -> str:
    return "; ".join(f"{p.label}: {p.left} / {p.right}" for p in pieces)


def _common_warnings(pieces: list[Piece], products: list[Ligated], circular: bool,
                     warnings: list[str]) -> None:
    if len(products) > 1:
        warnings.append(f"{len(products)} different products can form from these fragments "
                        "(different orders or orientations): screen colonies by colony PCR or a "
                        "diagnostic digest for the one you want")
    if any(p.from_pcr for p in pieces):
        warnings.append(PCR_PHOSPHATE_NOTE)
    ends = [e for p in pieces for e in (p.left, p.right)]
    if circular:
        blunt = "blunt" in ends
    else:
        blunt = "blunt" in ends[1:-1]  # outer ends of a linear product are not joined
    if blunt:
        warnings.append(BLUNT_NOTE)


def _dedupe(names: list[str]) -> list[str]:
    out: list[str] = []
    for n in names:
        if n not in out:
            out.append(n)
    return out


def _products(chosen: list[Ligated], pieces: list[Piece], base_name: str, operation: str,
              details: str, rotate_to: Dseqrecord | None) -> list[Product]:
    out = []
    parents = _dedupe([p.name for p in pieces])
    for index, item in enumerate(chosen):
        record = item.record
        if rotate_to is not None and record.circular:
            record = _rotate_to_origin(record, rotate_to)
        name = base_name if index == 0 else f"{base_name}, product {index + 1}"
        topology = "circular" if record.circular else "linear"
        label = f"{_order_text(item.order, pieces)} — {len(record)} bp {topology}"
        text = f"{details}. Order: {_order_text(item.order, pieces)}"
        if not record.circular:
            left = _end_text(record.seq.five_prime_end())
            right = _end_text(record.seq.three_prime_end())
            text += f". Ends: {left} / {right}"
        out.append(Product(doc=from_dseqrecord(record, name), operation=operation,
                           parents=parents, details=text, label=label))
    return out


# ---- operations ----------------------------------------------------------------------------

def insert_multiple(vector: FragmentIn, inserts: list[FragmentIn]) -> CloningResult:
    if not inserts:
        raise CloningError("Add at least one insert")
    warnings: list[str] = []
    vector_record = to_dseqrecord(vector.doc)
    if vector_record.circular and not vector.enzymes:
        raise CloningError(f"Choose the enzyme or enzymes that open the vector "
                           f"“{vector.doc.name}”")
    backbone = make_piece(vector, "vector", warnings)
    pieces = [backbone] + [make_piece(item, "insert", warnings) for item in inserts]

    products = _sort(ligation_products([p.record for p in pieces], circular=True),
                     len(pieces), circular=True)
    if not products:
        raise CloningError("These fragments cannot ligate into a plasmid." +
                           _junction_report(pieces, True) +
                           f" Fragment ends: {_describe_pieces(pieces)}.")
    if products[0].order != list(range(1, len(pieces) + 1)):
        warnings.insert(0, "The fragments do not ligate in the order listed with every insert "
                           "forward; the products below are the ones their ends allow")
    _common_warnings(pieces, products, True, warnings)
    try:
        circularise(backbone.record)
        if backbone.from_pcr:
            which = "the insert" if len(inserts) == 1 else "the inserts"
            warnings.append("The vector's two blunt ends could close on each other: phosphorylate "
                            f"only {which}, never the vector PCR product")
        else:
            warnings.append(f"The two ends of the cut vector are compatible ({backbone.left} / "
                            f"{backbone.right}): it can close without an insert. Dephosphorylate "
                            "the vector (rSAP or CIP) to reduce empty-vector colonies")
    except IncompatibleEnds:
        pass

    names = [p.name for p in pieces]
    base_name = " + ".join(names)
    how = (f"cut with {' and '.join(vector.enzymes)}" if vector.enzymes
           else "used as a linear fragment")
    details = (f"Ligated {plural(len(inserts), 'insert')} into {vector.doc.name} ({how}): "
               + ", ".join(p.label for p in pieces[1:]))
    return CloningResult(
        products=_products(products, pieces, base_name, "Insert multiple fragments", details,
                           vector_record),
        warnings=warnings)


def ligate_fragments(items: list[FragmentIn], circular: bool = False) -> CloningResult:
    if not items:
        raise CloningError("Add the fragments to ligate")
    if len(items) < 2 and not circular:
        raise CloningError("Add a second fragment, or close the fragment into a circle")
    warnings: list[str] = []
    pieces = [make_piece(item, "insert", warnings) for item in items]
    shape = "a circle" if circular else "a linear molecule"
    products = _sort(ligation_products([p.record for p in pieces], circular), len(pieces),
                     circular)
    if not products:
        raise CloningError(f"These fragments cannot ligate into {shape}." +
                           _junction_report(pieces, circular) +
                           f" Fragment ends: {_describe_pieces(pieces)}.")
    expected = _expected_order(len(pieces), circular)
    if products[0].order not in expected:
        warnings.insert(0, "The fragments do not ligate in the order listed; the products "
                           "below are the ones their ends allow")
    _common_warnings(pieces, products, circular, warnings)
    if not circular and len(pieces) > 1:
        closed = ligation_products([p.record for p in pieces], circular=True)
        if closed:
            warnings.append("The outer ends of the product are compatible with each other: "
                            "the same fragments can also close into a circle")

    base_name = " + ".join(p.name for p in pieces)
    details = (f"Ligated {plural(len(pieces), 'fragment')} into {shape}: "
               + ", ".join(p.label for p in pieces))
    return CloningResult(
        products=_products(products, pieces, base_name, "Ligate fragments", details,
                           pieces[0].record if circular else None),
        warnings=warnings)


# ---- HTTP ----------------------------------------------------------------------------------

class InsertMultipleRequest(BaseModel):
    vector: FragmentIn
    inserts: list[FragmentIn]


class LigateRequest(BaseModel):
    fragments: list[FragmentIn]
    circular: bool = False


@router.post("/insert-multiple")
def insert_multiple_route(body: InsertMultipleRequest) -> CloningResult:
    return insert_multiple(body.vector, body.inserts)


@router.post("/ligate")
def ligate_route(body: LigateRequest) -> CloningResult:
    return ligate_fragments(body.fragments, body.circular)
