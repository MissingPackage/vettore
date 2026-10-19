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
    return {
        "name": name or fallback_name,
        "description": record.description if record.description != "<unknown description>" else "",
        "sequence": sequence,
        "circular": circular,
        "features": features,
        "primers": primers,
    }


def _location(ann: dict[str, Any], seq_len: int) -> SimpleLocation | CompoundLocation:
    strand = -1 if (ann.get("strand") == -1 or ann.get("forward") is False) else 1
    ranges = [(loc["start"], loc["end"]) for loc in ann.get("locations") or []] or [
        (ann["start"], ann["end"])
    ]
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


def ove_to_record(data: dict[str, Any]) -> SeqRecord:
    """Convert OVE ``sequenceData`` back to a SeqRecord (GenBank-ready)."""
    sequence = data.get("sequence", "")
    name = (data.get("name") or "Untitled").strip()
    record = SeqRecord(
        Seq(sequence),
        id=name.replace(" ", "_")[:16] or "Untitled",
        name=name.replace(" ", "_")[:16] or "Untitled",
        description=data.get("description") or name,
        annotations={
            "molecule_type": "DNA",
            "topology": "circular" if data.get("circular") else "linear",
        },
    )
    for kind, default_type in (("features", "misc_feature"), ("primers", "primer_bind")):
        for ann in _as_list(data.get(kind)):
            qualifiers = {k: list(v) for k, v in (ann.get("notes") or {}).items()}
            qualifiers["label"] = [ann.get("name") or default_type]
            if ann.get("color"):
                qualifiers["ApEinfo_fwdcolor"] = [ann["color"]]
            if kind == "primers" and ann.get("bases"):
                qualifiers["primer_sequence"] = [ann["bases"]]
            record.features.append(
                SeqFeature(_location(ann, len(sequence)), type=ann.get("type") or default_type,
                           qualifiers=qualifiers)
            )
    return record


def _as_list(value: Any) -> list[dict[str, Any]]:
    """OVE stores annotations either as a list or as a dict keyed by id."""
    if value is None:
        return []
    if isinstance(value, dict):
        return list(value.values())
    return list(value)
