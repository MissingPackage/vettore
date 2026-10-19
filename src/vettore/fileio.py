"""Reading and writing sequence files."""

from __future__ import annotations

import html
import io
import logging
import re
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .convert import ove_to_record, record_to_ove

log = logging.getLogger(__name__)

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
            _enrich_dna_file(doc, content)
        docs.append(doc)
    return docs


_NEAR_WHITE = re.compile(r"^#?(f[a-f0-9]){3}$", re.IGNORECASE)


def _enrich_dna_file(doc: dict, content: bytes) -> None:
    """Add what Biopython does not read from .dna files: feature colours, full primer
    sequences (with 5' tails) and the description. Locations stay Biopython's, which already
    applies the file's primer hybridization parameters."""
    try:
        from sgffp import SgffReader

        sg = SgffReader.from_bytes(content)
    except Exception:  # enrichment is best effort; the Biopython parse already succeeded
        log.exception("sgffp could not read the .dna file; colours and primer tails skipped")
        return

    colors = {}
    for feature in sg.features:
        color = feature.color or (feature.segments[0].color if feature.segments else None)
        if feature.segments and color and not _NEAR_WHITE.match(color):
            colors[(feature.name, feature.segments[0].start)] = color
    for ann in doc["features"]:
        start = ann["locations"][0]["start"] if ann.get("locations") else ann["start"]
        if color := colors.get((ann["name"], start)):
            ann["color"] = color

    sequences = {p.name: p.sequence for p in sg.primers if p.sequence}
    for ann in doc["primers"]:
        if bases := sequences.get(ann["name"]):
            ann["bases"] = bases.upper()
            ann.setdefault("notes", {})["primer_sequence"] = [bases.upper()]

    description = sg.notes.get("Description") if hasattr(sg.notes, "get") else None
    if description and not doc.get("description"):
        doc["description"] = html.unescape(re.sub(r"<[^>]+>", " ", str(description))).strip()


def write_genbank(data: dict) -> str:
    out = io.StringIO()
    SeqIO.write(ove_to_record(data), out, "genbank")
    return out.getvalue()


def write_fasta(data: dict) -> str:
    out = io.StringIO()
    SeqIO.write(ove_to_record(data), out, "fasta")
    return out.getvalue()
