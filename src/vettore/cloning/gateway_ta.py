"""TA, GC and TOPO cloning of PCR products. Routes under /api/cloning.

The PCR product is either an open linear document (the amplicon as amplified, primers and
tails included) or is made here from a template and two primers. The vector is an open
document holding the vector's double-stranded sequence only: a linearised vector as a linear
document, or a circular map plus the position where it is opened. Vettore adds the
single-stranded ends each method relies on, then ligates with pydna (the rule of
``pydna.assembly2.ligation_assembly``, see :mod:`.ligation`), so the orientations follow from
the ends:

- **TA**: Taq-type polymerases add one non-templated 3′ A to each strand of the PCR product
  (Clark 1988, Nucleic Acids Res 16:9677); T-vectors carry one 3′ T on each end (Promega
  pGEM-T Vector Systems manual TM042: pGEM-5Zf(+) cut with EcoRV, then a 3′ T added to both
  ends). The insert ligates in both orientations.
- **GC**: vector with a single 3′ C on each end, insert with a single 3′ G (Lucigen GC cloning,
  US patent 7,723,103 B2). Both orientations.
- **TOPO TA**: a T-vector charged with vaccinia topoisomerase I (Shuman 1994, J Biol Chem
  269:32678; Invitrogen TOPO TA Cloning manual): the same ends as TA cloning.
- **Blunt TOPO**: blunt vector, blunt PCR product from a proofreading polymerase; both
  orientations (Invitrogen Zero Blunt TOPO manual).
- **Directional TOPO**: one vector end is blunt, the other has the 5′ overhang GTGG; the
  forward primer starts with CACC, the overhang anneals to it and fixes the orientation, and a
  reverse primer starting with CACC would let the insert go in backwards (Invitrogen pENTR
  Directional TOPO Cloning Kits manual). The CACC stays in the clone, written once.

Products are numbered like the vector: a circular vector keeps its own numbering, a linearised
vector starts at its first base.
"""

from __future__ import annotations

from typing import Literal

from Bio.Seq import Seq
from fastapi import APIRouter
from pydantic import BaseModel
from pydna.dseq import Dseq
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
from .ligation import ligation_products
from .pcr import PrimerIn, run_pcr
from .regions import features_in_region

router = APIRouter()

Kind = Literal["TA", "GC", "TOPO-TA", "blunt-TOPO", "directional-TOPO"]

OPERATIONS: dict[str, str] = {
    "TA": "TA cloning",
    "GC": "GC cloning",
    "TOPO-TA": "TOPO TA cloning",
    "blunt-TOPO": "Blunt TOPO cloning",
    "directional-TOPO": "Directional TOPO cloning",
}

# Single 3′ bases: (vector, insert). Sources in the module docstring.
SINGLE_BASE_ENDS: dict[str, tuple[str, str]] = {
    "TA": ("T", "A"),
    "TOPO-TA": ("T", "A"),
    "GC": ("C", "G"),
}
DIRECTIONAL_TAG = "CACC"      # 5′ end of the forward primer (pENTR Directional TOPO manual)

BOTH_ORIENTATIONS = ("The insert ligates in either orientation, so both products are shown: "
                     "check colonies by colony PCR or sequencing for the one you want")
POLYMERASE_NOTES: dict[str, str] = {
    "TA": ("Amplify with Taq or another non-proofreading polymerase, which leaves the 3′ A "
           "overhangs; products of a proofreading polymerase are blunt and need an A-tailing "
           "step with Taq and dATP first (see the vector manual)"),
    "GC": ("The insert needs a single 3′ G overhang: amplify with the polymerase the kit "
           "specifies, or G-tail a blunt product as the kit manual describes"),
    "blunt-TOPO": ("Amplify with a proofreading polymerase: its blunt products clone into blunt "
                   "TOPO vectors, while Taq products with 3′ A overhangs do not clone well"),
    "directional-TOPO": ("Amplify with a proofreading polymerase (blunt ends); a reverse primer "
                         "that also starts with CACC would let the insert go in backwards"),
}
POLYMERASE_NOTES["TOPO-TA"] = POLYMERASE_NOTES["TA"]


def _rc(bases: str) -> str:
    return str(Seq(bases).reverse_complement()).upper()


class Opened:
    """The vector as a double-stranded linear body, with the map it came from."""

    def __init__(self, doc: SequenceDoc, position: int | None):
        record = to_dseqrecord(doc)
        n = len(record)
        sequence = str(record.seq).upper()
        if record.circular:
            if position is None:
                raise CloningError(f"“{doc.name}” is circular: enter the position where the "
                                   "vector is opened (the cloning site), or open the "
                                   "linearised vector file instead")
            if not 0 <= position <= n:
                raise CloningError(f"The opening position must be between 0 and {n}")
            position %= n
            self.body = sequence[position:] + sequence[:position]
            self.start, self.end = position, (position - 1) % n
        else:
            self.body = sequence
            self.start, self.end = 0, n - 1
            position = 0
        self.name = doc.name
        self.position = position
        self.record = record
        self.circular = record.circular

    def features(self, offset: int):
        return features_in_region(self.record, self.start, self.end, offset)


class Amplicon:
    """The PCR product, as amplified."""

    def __init__(self, record: Dseqrecord, name: str, parent: str, oligos: list[Oligo],
                 forward: str | None, reverse: str | None):
        self.record = record
        self.sequence = str(record.seq).upper()
        self.name = name
        self.parent = parent
        self.oligos = oligos
        self.forward = forward
        self.reverse = reverse

    def features(self, offset: int):
        return features_in_region(self.record, 0, len(self.record) - 1, offset)


def _amplicon(insert: SequenceDoc | None, template: SequenceDoc | None,
              forward: PrimerIn | None, reverse: PrimerIn | None) -> Amplicon:
    if template is not None:
        if forward is None or reverse is None:
            raise CloningError("Enter both primers to amplify the template")
        pcr = run_pcr(template, forward, reverse)
        product = pcr.products[0]
        record = to_dseqrecord(SequenceDoc(**product.doc))
        return Amplicon(record, product.doc["name"], template.name, pcr.oligos,
                        pcr.oligos[0].sequence, pcr.oligos[1].sequence)
    if insert is None:
        raise CloningError("Choose the PCR product, or a template and two primers")
    record = to_dseqrecord(insert)
    if record.circular:
        raise CloningError(f"“{insert.name}” is circular: choose the linear PCR product, or "
                           "amplify it here from a template and two primers")
    return Amplicon(record, insert.name, insert.name, [], None, None)


def _single_base(body: str, base: str, features) -> Dseqrecord:
    """``body`` with one unpaired base on each 3′ end (crick reads 5′→3′)."""
    record = Dseqrecord(Dseq(body + base, _rc(body) + base, ovhg=1))
    record.features = list(features)
    return record


def _ends(kind: str, vector: Opened, amplicon: Amplicon) -> tuple[Dseqrecord,
                                                                 list[tuple[Dseqrecord, bool]]]:
    """The vector and the insert molecule(s) with the ends the method relies on; each insert
    comes with ``True`` when it was built from the reverse complement of the PCR product."""
    if kind in SINGLE_BASE_ENDS:
        v_base, i_base = SINGLE_BASE_ENDS[kind]
        return (_single_base(vector.body, v_base, vector.features(1)),
                [(_single_base(amplicon.sequence, i_base, amplicon.features(1)), False)])
    if kind == "blunt-TOPO":
        v = Dseqrecord(vector.body)
        v.features = vector.features(0)
        i = Dseqrecord(amplicon.sequence)
        i.features = amplicon.features(0)
        return v, [(i, False)]
    # Directional TOPO: the vector's right end has the 5′ overhang complementary to CACC; the
    # overhang displaces the first four bases of the strand opposite CACC in the PCR product.
    v = Dseqrecord(Dseq(vector.body, _rc(vector.body + DIRECTIONAL_TAG), ovhg=0))
    v.features = vector.features(0)
    inserts = []
    for source, from_rc in ((amplicon.record, False),
                            (amplicon.record.reverse_complement(), True)):
        sequence = str(source.seq).upper()
        if sequence.startswith(DIRECTIONAL_TAG):
            i = Dseqrecord(Dseq(sequence, _rc(sequence[len(DIRECTIONAL_TAG):]),
                                ovhg=-len(DIRECTIONAL_TAG)))
            i.features = features_in_region(source, 0, len(source) - 1, 0)
            inserts.append((i, from_rc))
    return v, inserts


def _missing_cacc(amplicon: Amplicon) -> CloningError:
    if amplicon.forward is not None:
        fix = (f" Order the forward primer as {DIRECTIONAL_TAG}{amplicon.forward} and run the "
               "PCR again.")
    else:
        fix = (" Add CACC to the 5′ end of the forward primer and amplify again.")
    return CloningError(f"The PCR product “{amplicon.name}” does not start with "
                        f"{DIRECTIONAL_TAG}: directional TOPO cloning needs these four bases "
                        "at the 5′ end of the forward primer, where they pair with the vector's "
                        "GTGG overhang." + fix)


def _rotate(product: Dseqrecord, vector: Opened) -> Dseqrecord:
    """Number the product like the vector map it came from."""
    doubled = str(product.seq).upper() * 2
    at = doubled.find(vector.body)
    if at < 0 or at >= len(product):
        return product
    shift = (at + (len(vector.body) - vector.position) % len(vector.body)) % len(product)
    return product.shifted(shift) if shift else product


def ta_topo(kind: Kind, vector_doc: SequenceDoc, position: int | None = None,
            insert: SequenceDoc | None = None, template: SequenceDoc | None = None,
            forward: PrimerIn | None = None, reverse: PrimerIn | None = None) -> CloningResult:
    vector = Opened(vector_doc, position)
    amplicon = _amplicon(insert, template, forward, reverse)
    if kind == "directional-TOPO" and amplicon.forward is not None \
            and not amplicon.forward.upper().startswith(DIRECTIONAL_TAG):
        raise _missing_cacc(amplicon)
    v, inserts = _ends(kind, vector, amplicon)
    if not inserts:
        raise _missing_cacc(amplicon)

    found: list[tuple[Dseqrecord, bool]] = []    # (product, insert reversed)
    for insert_record, from_rc in inserts:
        for item in ligation_products([v, insert_record], circular=True):
            found.append((item.record, from_rc != (item.order[1] < 0)))
    if not found:  # defensive: the ends built above are compatible by construction
        raise CloningError("The vector and the insert ends do not match")
    found.sort(key=lambda f: f[1])                # forward first

    warnings: list[str] = []
    if kind == "directional-TOPO":
        if len(inserts) > 1:
            warnings.append("Both ends of the PCR product start with CACC, so the insert can go "
                            "in either orientation: remove CACC from the reverse primer")
        elif inserts[0][1]:
            warnings.append("CACC is at the 5′ end of the reverse primer, not the forward one: "
                            "the insert goes in reversed")
    if len(found) > 1:
        warnings.append(BOTH_ORIENTATIONS)
    warnings.append(POLYMERASE_NOTES[kind])

    operation = OPERATIONS[kind]
    base_name = f"{vector.name} + {amplicon.name}"
    where = (f"opened after base {vector.position}" if vector.circular
             else "linearised vector")
    products: list[Product] = []
    for index, (record, flipped) in enumerate(found):
        record = _rotate(record, vector)
        orientation = "reversed" if flipped else "forward"
        name = base_name if index == 0 else f"{base_name}, product {index + 1}"
        details = (f"{operation} of {amplicon.name} ({len(amplicon.sequence)} bp) into "
                   f"{vector.name} ({where}); insert {orientation}")
        products.append(Product(
            doc=from_dseqrecord(record, name), operation=operation,
            parents=list(dict.fromkeys([vector.name, amplicon.parent])), details=details,
            oligos=amplicon.oligos,
            label=f"Insert {orientation} — {len(record)} bp circular"))
    return CloningResult(products=products, warnings=warnings, oligos=amplicon.oligos)


class TaTopoRequest(BaseModel):
    kind: Kind = "TA"
    vector: SequenceDoc
    position: int | None = None      # circular vector: bases before the opening point
    insert: SequenceDoc | None = None
    template: SequenceDoc | None = None
    forward: PrimerIn | None = None
    reverse: PrimerIn | None = None


@router.post("/ta-topo")
def ta_topo_route(body: TaTopoRequest) -> CloningResult:
    return ta_topo(body.kind, body.vector, body.position, body.insert, body.template,
                   body.forward, body.reverse)
