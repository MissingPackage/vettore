"""Reading and writing sequence files."""

from __future__ import annotations

import io
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .convert import ove_to_record, record_to_ove

FORMATS = {
    ".dna": "snapgene",
    ".gb": "genbank",
    ".gbk": "genbank",
    ".gbff": "genbank",
    ".genbank": "genbank",
    ".ape": "genbank",
    ".fa": "fasta",
    ".fasta": "fasta",
    ".fna": "fasta",
    ".fas": "fasta",
    ".ab1": "abi",
    ".abi": "abi",
}
BINARY_FORMATS = {"snapgene", "abi"}


class UnsupportedFile(ValueError):
    pass


def detect_format(filename: str) -> str:
    fmt = FORMATS.get(Path(filename).suffix.lower())
    if fmt is None:
        raise UnsupportedFile(f"Unsupported file type: {Path(filename).suffix or filename}")
    return fmt


def read_records(filename: str, content: bytes) -> list[SeqRecord]:
    fmt = detect_format(filename)
    handle = io.BytesIO(content) if fmt in BINARY_FORMATS else io.StringIO(
        content.decode("utf-8", errors="replace")
    )
    records = list(SeqIO.parse(handle, fmt))
    if not records:
        raise UnsupportedFile(f"No sequence found in {Path(filename).name}")
    return records


def read_documents(filename: str, content: bytes) -> list[dict]:
    """Parse a file into OVE ``sequenceData`` documents (one per record)."""
    stem = Path(filename).stem
    records = read_records(filename, content)
    docs = []
    for i, record in enumerate(records):
        fallback = stem if len(records) == 1 else f"{stem} ({i + 1})"
        doc = record_to_ove(record, fallback_name=fallback)
        if detect_format(filename) == "snapgene":
            # .dna files carry no usable locus name; the file name is what users see.
            doc["name"] = fallback
        docs.append(doc)
    return docs


def write_genbank(data: dict) -> str:
    out = io.StringIO()
    SeqIO.write(ove_to_record(data), out, "genbank")
    return out.getvalue()


def write_fasta(data: dict) -> str:
    out = io.StringIO()
    SeqIO.write(ove_to_record(data), out, "fasta")
    return out.getvalue()
