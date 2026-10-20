"""Conversion between Biopython SeqRecords and the Open Vector Editor (OVE) JSON model.

Coordinate conventions:
- Biopython: 0-based, end-exclusive (``SimpleLocation(start, end)``).
- OVE: 0-based, end-inclusive (``start``/``end``). A feature that crosses the origin of a
  circular sequence has ``start > end``.
"""

from __future__ import annotations

import uuid
from typing import Any

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

PRIMER_TYPES = {"primer_bind", "primer"}
LOCUS_LENGTH = 16            # characters of a GenBank LOCUS name that every reader accepts
# GenBank divisions (NCBI GenBank release notes, section 3.4.4); "UNK" is not one of them.
DIVISIONS = frozenset({"PRI", "ROD", "MAM", "VRT", "INV", "PLN", "BCT", "VRL", "PHG", "SYN",
                       "UNA", "EST", "PAT", "STS", "GSS", "HTG", "HTC", "ENV", "CON", "TSA"})
NAME_QUALIFIERS = ("label", "name", "gene", "product", "locus_tag", "note")


def _first(qualifiers: dict[str, Any], key: str) -> str | None:
    value = qualifiers.get(key)
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        return str(value[0]) if value else None
    return str(value)


def feature_name(feature: SeqFeature) -> str:
    for key in NAME_QUALIFIERS:
        value = _first(feature.qualifiers, key)
        if value:
            return value
    return feature.type


def _parts_to_ranges(feature: SeqFeature, seq_len: int, circular: bool) -> list[tuple[int, int]]:
    """Return inclusive (start, end) ranges in transcription order, merging origin wraps."""
    parts = list(feature.location.parts)
    if feature.location.strand == -1:
        # Biopython lists parts of minus-strand joins in transcription order (high to low);
        # OVE wants locations sorted by position.
        parts = sorted(parts, key=lambda p: int(p.start))
    ranges = [(int(p.start), int(p.end) - 1) for p in parts]
    if circular and len(ranges) >= 2:
        merged: list[tuple[int, int]] = []
        for rng in ranges:
            if merged and merged[-1][1] == seq_len - 1 and rng[0] == 0:
                merged[-1] = (merged[-1][0], rng[1])
            else:
                merged.append(rng)
        # A wrap can also be expressed as the last part ending at the origin.
        if len(merged) >= 2 and merged[0][0] == 0 and merged[-1][1] == seq_len - 1:
            merged = [(merged[-1][0], merged[0][1])] + merged[1:-1]
        ranges = merged
    return ranges


def _annotation(feature: SeqFeature, seq_len: int, circular: bool) -> dict[str, Any]:
    ranges = _parts_to_ranges(feature, seq_len, circular)
    strand = feature.location.strand or 1
    ann: dict[str, Any] = {
        "id": uuid.uuid4().hex[:12],
        "name": feature_name(feature),
        "type": feature.type,
        "strand": strand,
        "forward": strand != -1,
        "start": ranges[0][0],
        "end": ranges[-1][1],
        "notes": {k: [str(x) for x in v] if isinstance(v, (list, tuple)) else [str(v)]
                  for k, v in feature.qualifiers.items()},
    }
    if len(ranges) > 1:
        ann["locations"] = [{"start": s, "end": e} for s, e in ranges]
    if bases := _first(feature.qualifiers, "primer_sequence"):
        ann["bases"] = bases
    color = _first(feature.qualifiers, "ApEinfo_fwdcolor") or _first(feature.qualifiers, "color")
    if color:
        ann["color"] = color
    return ann


def record_to_ove(record: SeqRecord, fallback_name: str = "Untitled") -> dict[str, Any]:
    """Convert a SeqRecord to OVE ``sequenceData``."""
    sequence = str(record.seq)
    circular = str(record.annotations.get("topology", "linear")).lower() == "circular"
    name = record.name if record.name and record.name not in {"<unknown name>", "."} else ""
    features, primers = [], []
    for feature in record.features:
        if feature.type == "source":
            continue
        ann = _annotation(feature, len(sequence), circular)
        (primers if feature.type in PRIMER_TYPES else features).append(ann)
    doc = {
        "name": name or fallback_name,
        "description": record.description if record.description != "<unknown description>" else "",
        "sequence": sequence,
        "circular": circular,
        "features": features,
        "primers": primers,
    }
    if (division := record.annotations.get("data_file_division")) in DIVISIONS:
        doc["division"] = division       # written back when the document is saved as GenBank
    return doc


def locus_name(name: str) -> str:
    """A GenBank LOCUS name: no white space, at most 16 characters, shortened at a word or
    hyphen boundary when the name is longer (the full name goes in the DEFINITION line)."""
    compact = "_".join(name.split()) or "Untitled"
    if len(compact) <= LOCUS_LENGTH:
        return compact
    cut = max(compact.rfind(mark, 0, LOCUS_LENGTH + 1) for mark in "_-.")
    return compact[:cut] if cut >= LOCUS_LENGTH // 2 else compact[:LOCUS_LENGTH]


def _is_position(value: Any, length: int) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and 0 <= value < length


def annotation_ranges(ann: Any, length: int, circular: bool) -> list[tuple[int, int]] | None:
    """The inclusive (start, end) ranges of an OVE feature or primer, or None when its position
    cannot be placed on the sequence (missing, not a number, outside the sequence, or across
    the end of a linear sequence)."""
    if not isinstance(ann, dict):
        return None
    locations = ann.get("locations")
    if isinstance(locations, list) and locations:
        raw = [(loc.get("start"), loc.get("end")) if isinstance(loc, dict) else (None, None)
               for loc in locations]
    else:
        raw = [(ann.get("start"), ann.get("end"))]
    for start, end in raw:
        if not (_is_position(start, length) and _is_position(end, length)):
            return None
        if start > end and not circular:
            return None
    return raw


def annotation_label(ann: Any, default: str) -> str:
    name = ann.get("name") if isinstance(ann, dict) else None
    return str(name or (ann.get("type") if isinstance(ann, dict) else None) or default)


def unplaceable(data: dict[str, Any]) -> dict[str, list[str]]:
    """Names of the features and primers whose position cannot be placed on the sequence,
    by kind ("features", "primers"); conversion leaves them out."""
    sequence = data.get("sequence")
    length = len(sequence) if isinstance(sequence, str) else 0
    circular = bool(data.get("circular"))
    out: dict[str, list[str]] = {}
    for kind, default in (("features", "feature"), ("primers", "primer")):
        bad = [annotation_label(ann, default) for ann in _as_list(data.get(kind))
               if annotation_ranges(ann, length, circular) is None]
        if bad:
            out[kind] = bad
    return out


def _location(ann: dict[str, Any], ranges: list[tuple[int, int]],
              seq_len: int) -> SimpleLocation | CompoundLocation:
    strand = -1 if (ann.get("strand") == -1 or ann.get("forward") is False) else 1
    parts: list[SimpleLocation] = []
    for start, end in ranges:
        if start > end:  # crosses the origin
            parts.append(SimpleLocation(start, seq_len, strand))
            parts.append(SimpleLocation(0, end + 1, strand))
        else:
            parts.append(SimpleLocation(start, end + 1, strand))
    if strand == -1:
        parts = parts[::-1]
    return parts[0] if len(parts) == 1 else CompoundLocation(parts)


def _qualifiers(notes: Any) -> dict[str, list[str]]:
    if not isinstance(notes, dict):
        return {}
    return {str(k): [str(x) for x in v] if isinstance(v, (list, tuple)) else [str(v)]
            for k, v in notes.items()}


def ove_to_record(data: dict[str, Any]) -> SeqRecord:
    """Convert OVE ``sequenceData`` back to a SeqRecord (GenBank-ready). Features and primers
    whose position cannot be placed on the sequence are left out (see ``unplaceable``)."""
    sequence = data.get("sequence")
    sequence = sequence if isinstance(sequence, str) else ""
    circular = bool(data.get("circular"))
    name = str(data.get("name") or "Untitled").strip() or "Untitled"
    locus = locus_name(name)
    description = str(data.get("description") or "").strip() or name
    if locus != "_".join(name.split()) and name not in description:
        description = f"{name}: {description}"     # the LOCUS name was shortened
    record = SeqRecord(
        Seq(sequence),
        id=locus,
        name=locus,
        description=description,
        annotations={
            "molecule_type": "DNA",
            "topology": "circular" if data.get("circular") else "linear",
        },
    )
    for kind, default_type in (("features", "misc_feature"), ("primers", "primer_bind")):
        for ann in _as_list(data.get(kind)):
            ranges = annotation_ranges(ann, len(sequence), circular)
            if ranges is None:
                continue
            qualifiers = _qualifiers(ann.get("notes"))
            qualifiers["label"] = [str(ann.get("name") or default_type)]
            if ann.get("color"):
                qualifiers["ApEinfo_fwdcolor"] = [str(ann["color"])]
            if kind == "primers" and ann.get("bases"):
                qualifiers["primer_sequence"] = [str(ann["bases"])]
            record.features.append(
                SeqFeature(_location(ann, ranges, len(sequence)),
                           type=str(ann.get("type") or default_type), qualifiers=qualifiers)
            )
    return record


def _as_list(value: Any) -> list[Any]:
    """OVE stores annotations either as a list or as a dict keyed by id; anything else holds
    none."""
    if isinstance(value, dict):
        return list(value.values())
    if isinstance(value, list):
        return value
    return []
