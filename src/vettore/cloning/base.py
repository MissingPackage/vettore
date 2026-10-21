"""Shared models and conversions for cloning simulations.

Contract for operation modules (``restriction.py``, ``overlap.py``, ``golden_gate.py``, ...):

- Inputs arrive as OVE ``sequenceData`` dicts (the same shape the editor uses) inside pydantic
  request models; convert them with :func:`to_dseqrecord`.
- Run the biology with pydna (``pydna.assembly2``; ``pydna.assembly`` is deprecated).
- Return a :class:`CloningResult`. Each product is converted back with :func:`from_dseqrecord`
  and carries provenance (operation, parent names, human-readable parameters) for History.
- Raise :class:`CloningError` with a message a biologist understands for any user mistake
  (incompatible ends, primer does not bind, no product); the API turns it into HTTP 422.
- Never call the network (no ``pydna.genbank``/``tm_neb``).
"""

from __future__ import annotations

import re
from typing import Any

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pydantic import BaseModel, Field
from pydna.dseqrecord import Dseqrecord

from ..convert import ove_to_record, record_to_ove


class CloningError(ValueError):
    """A user-facing error: the message is shown in the dialog as is."""


class SequenceDoc(BaseModel):
    """An OVE ``sequenceData`` document (only the fields the engine needs are typed)."""

    model_config = {"extra": "allow"}

    name: str = "Untitled"
    sequence: str
    circular: bool = False
    features: list[dict[str, Any]] | dict[str, dict[str, Any]] = Field(default_factory=list)
    primers: list[dict[str, Any]] | dict[str, dict[str, Any]] = Field(default_factory=list)


class Oligo(BaseModel):
    """A primer: full 5'->3' sequence (tail + annealing part)."""

    name: str
    sequence: str
    tm: float | None = None          # Tm of the annealing part, Vettore's default preset
    annealing_length: int | None = None


class Product(BaseModel):
    doc: dict[str, Any]              # OVE sequenceData of the product
    operation: str                   # e.g. "Restriction cloning", "Gibson assembly"
    parents: list[str]               # names of the input documents
    details: str                     # human-readable parameters for History
    oligos: list[Oligo] = Field(default_factory=list)
    label: str | None = None         # tells alternative products apart in the product list


class CloningResult(BaseModel):
    products: list[Product]
    warnings: list[str] = Field(default_factory=list)
    oligos: list[Oligo] = Field(default_factory=list)   # primers to order, if any were designed
    # Wrong clones a colony screen must tell apart from the product (labelled: "insert in the
    # other orientation", "empty vector, re-ligated", "parent vector"); the diagnostic digest
    # planner starts from them.
    screen_against: list[Product] = Field(default_factory=list)


def plural(count: int, noun: str, many: str | None = None) -> str:
    """'1 insert', '2 inserts' (or the plural given)."""
    return f"{count} {noun if count == 1 else many or noun + 's'}"


def times(count: int) -> str:
    """'once', '2 times'."""
    return "once" if count == 1 else f"{count} times"


_NOT_DNA = re.compile("[^ACGTURYSWKMBDHVN]", re.IGNORECASE)
_RNA_TO_DNA = str.maketrans("Uu", "Tt")


def dna_sequence(name: str, sequence: str) -> str:
    """The sequence as the engine works on it: an RNA letter U is read as T (so enzymes cut
    and primers bind), and a character that is not a DNA letter is a mistake the user can fix."""
    if bad := _NOT_DNA.search(sequence):
        shown = "a space" if bad.group().isspace() else f"“{bad.group()}”"
        raise CloningError(
            f"“{name}” has {shown} at position {bad.start() + 1}, which is not a DNA letter "
            "(A, C, G, T or an IUPAC ambiguity code such as N). Correct the sequence and try "
            "again.")
    return sequence.translate(_RNA_TO_DNA)


def to_dseqrecord(doc: SequenceDoc | dict[str, Any]) -> Dseqrecord:
    data = doc.model_dump() if isinstance(doc, SequenceDoc) else doc
    name = str(data.get("name") or "Untitled")
    sequence = data.get("sequence")
    if not sequence or not isinstance(sequence, str):
        raise CloningError(f"“{name}” has no sequence")
    data = {**data, "sequence": dna_sequence(name, sequence)}
    record: SeqRecord = ove_to_record(data)
    molecule = Dseqrecord(record, circular=bool(data.get("circular")))
    molecule.name = name           # messages name the document as the user sees it
    return molecule


def from_dseqrecord(record: Dseqrecord, name: str) -> dict[str, Any]:
    """Convert a pydna product to OVE sequenceData (features and primer_bind sites kept)."""
    seq_record = SeqRecord(
        Seq(str(record.seq)),
        id=name[:16], name=name[:16], description=name,
        features=list(record.features),
        annotations={"molecule_type": "DNA",
                     "topology": "circular" if record.circular else "linear"},
    )
    doc = record_to_ove(seq_record, fallback_name=name)
    doc["name"] = name
    return doc


def dedupe(records: list[Dseqrecord]) -> list[Dseqrecord]:
    """Drop identical molecules (pydna returns every chemically possible product)."""
    seen: set[str] = set()
    out = []
    for r in records:
        key = r.seguid()
        if key not in seen:
            seen.add(key)
            out.append(r)
    return out
