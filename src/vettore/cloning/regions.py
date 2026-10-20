"""Regions on linear and circular molecules, and feature transfer between molecules.

A region is a pair of 0-based inclusive top-strand coordinates. On a circular molecule
``start > end`` means the region spans the origin.
"""

from __future__ import annotations

from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from pydna.dseqrecord import Dseqrecord

from .base import CloningError


def region_sequence(sequence: str, start: int, end: int) -> str:
    if start <= end:
        return sequence[start : end + 1]
    return sequence[start:] + sequence[: end + 1]


def region_length(length: int, start: int, end: int) -> int:
    return end - start + 1 if start <= end else length - start + end + 1


def check_region(record: Dseqrecord, start: int, end: int, what: str) -> None:
    n = len(record)
    if not (0 <= start < n and 0 <= end < n):
        raise CloningError(f"The {what} region 1..{n} does not contain "
                           f"{start + 1}..{end + 1}")
    if start > end and not record.circular:
        raise CloningError(f"The {what} “{record.name}” is linear, so the region cannot span "
                           "the origin")


def _mapper(length: int, start: int, end: int, offset: int):
    """Map a template position into the extracted region, or None if outside."""
    def mapped(pos: int) -> int | None:
        if start <= end:
            return offset + pos - start if start <= pos <= end else None
        if pos >= start:
            return offset + pos - start
        if pos <= end:
            return offset + length - start + pos
        return None

    return mapped


def features_in_region(record: Dseqrecord, start: int, end: int, offset: int) -> list[SeqFeature]:
    """Features of ``record`` entirely inside the region, shifted to start at ``offset``."""
    mapped = _mapper(len(record), start, end, offset)
    out: list[SeqFeature] = []
    for feature in record.features:
        if feature.type == "source":
            continue
        parts = []
        for part in feature.location.parts:
            first = mapped(int(part.start))
            last = mapped(int(part.end) - 1)
            if first is None or last is None or last < first:
                parts = []
                break
            parts.append(SimpleLocation(first, last + 1, part.strand))
        if not parts:
            continue
        if feature.location.strand == -1:
            parts = parts[::-1]
        location = parts[0] if len(parts) == 1 else CompoundLocation(parts)
        out.append(SeqFeature(location, type=feature.type, qualifiers=dict(feature.qualifiers)))
    return out


def transfer_features(source: Dseqrecord, target: Dseqrecord, start: int, end: int,
                      offset: int) -> None:
    target.features.extend(features_in_region(source, start, end, offset))


def annotate(record: Dseqrecord, start: int, end: int, name: str, kind: str = "misc_feature",
             strand: int = 1) -> None:
    """Add a feature over an inclusive range (which may wrap the origin)."""
    n = len(record)
    if start <= end:
        location = SimpleLocation(start, end + 1, strand)
    else:
        parts = [SimpleLocation(start, n, strand), SimpleLocation(0, end + 1, strand)]
        location = CompoundLocation(parts[::-1] if strand == -1 else parts)
    record.features.append(SeqFeature(location, type=kind, qualifiers={"label": [name]}))
