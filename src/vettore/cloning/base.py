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


def to_dseqrecord(doc: SequenceDoc | dict[str, Any]) -> Dseqrecord:
    data = doc.model_dump() if isinstance(doc, SequenceDoc) else doc
    if not data.get("sequence"):
        raise CloningError(f"“{data.get('name') or 'Untitled'}” has no sequence")
    record: SeqRecord = ove_to_record(data)
    return Dseqrecord(record, circular=bool(data.get("circular")))


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
