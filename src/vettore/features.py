"""Detect common features in a sequence.

The database in ``data/common_features.json`` holds real subsequences of public GenBank
records (see ``scripts/build_feature_database.py``), so every annotation can be traced back to
an accession and a coordinate range.

Matching is deliberately not exact. Promoters and ITRs differ between vectors by a few bases,
and an ITR that is a few nucleotides shorter is still an ITR, so long features are found by
seeding on shared k-mers and then aligning locally, keeping matches above an identity
threshold. Short features (sites and tags) are found by exact search, on both strands, because
a single wrong base in a loxP site changes its behaviour.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources
from typing import Any

from Bio import Align
from Bio.Seq import Seq

SHORT_FEATURE = 50          # at or below this length, require an exact match
SEED = 18                   # k-mer length used to find candidate positions
DEFAULT_MIN_IDENTITY = 95.0
DEFAULT_MIN_COVERAGE = 90.0  # how much of the catalogue feature must be found

#: ITRs are routinely shortened in vectors (the common "130 bp" ITR is a 14 nt truncation of
#: the 145 nt wild type, and shorter variants exist), so they are matched permissively.
COVERAGE_BY_TYPE = {"repeat_region": 70.0}

#: A local alignment stops at the first mismatched base, which can leave a start codon out of
#: the match. A few bases are added back; a longer trim is a genuinely shorter feature (a
#: truncated ITR, say) and is left alone.
MAX_EXTENSION = 6


@dataclass(frozen=True)
class CatalogueFeature:
    name: str
    type: str
    sequence: str
    note: str
    source: str

    @property
    def length(self) -> int:
        return len(self.sequence)


@lru_cache(maxsize=1)
def catalogue() -> tuple[CatalogueFeature, ...]:
    text = resources.files("vettore.data").joinpath("common_features.json").read_text("utf-8")
    payload = json.loads(text)
    return tuple(CatalogueFeature(name=f["name"], type=f["type"], sequence=f["sequence"].upper(),
                                  note=f.get("note", ""), source=f["source"])
                 for f in payload["features"])


@lru_cache(maxsize=1)
def _aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner(scoring="blastn")
    aligner.mode = "local"
    return aligner


def _reverse_complement(sequence: str) -> str:
    return str(Seq(sequence).reverse_complement())


def _seeds(sequence: str, count: int = 6) -> list[tuple[int, str]]:
    """A few evenly spaced k-mers, used to locate candidate positions cheaply."""
    if len(sequence) <= SEED:
        return [(0, sequence)]
    positions = {round(i * (len(sequence) - SEED) / max(1, count - 1)) for i in range(count)}
    return [(p, sequence[p : p + SEED]) for p in sorted(positions)]


def _candidate_offsets(target: str, feature: str) -> set[int]:
    """Where the feature could start in the target, from exact seed hits."""
    offsets: set[int] = set()
    for position, seed in _seeds(feature):
        start = target.find(seed)
        while start != -1:
            offsets.add(max(0, start - position))
            start = target.find(seed, start + 1)
    return offsets


@dataclass
class Match:
    name: str
    type: str
    start: int           # 0-based inclusive on the original sequence
    end: int             # inclusive; end < start means the match crosses the origin
    strand: int
    identity: float
    note: str
    source: str

    def span(self, length: int) -> set[int]:
        if self.start <= self.end:
            return set(range(self.start, self.end + 1))
        return set(range(self.start, length)) | set(range(self.end + 1))

    def as_dict(self, length: int) -> dict[str, Any]:
        return {"name": self.name, "type": self.type, "start": self.start, "end": self.end,
                "strand": self.strand, "forward": self.strand == 1,
                "identity": round(self.identity, 1), "note": self.note, "source": self.source,
                "length": len(self.span(length))}


def _match_long(target: str, feature: CatalogueFeature, query: str, strand: int,
                min_identity: float, min_coverage: float) -> list[tuple[int, int, float]]:
    results: list[tuple[int, int, float]] = []
    aligner = _aligner()
    for offset in _candidate_offsets(target, query):
        margin = max(30, len(query) // 10)
        window_start = max(0, offset - margin)
        window = target[window_start : offset + len(query) + margin]
        if len(window) < SEED:
            continue
        alignments = aligner.align(window, query)
        if not len(alignments):
            continue
        best = alignments[0]
        aligned_ref, aligned_query = best[0], best[1]
        matches = sum(1 for a, b in zip(aligned_ref, aligned_query) if a == b and a != "-")
        coverage = 100.0 * sum(1 for b in aligned_query if b != "-") / len(query)
        identity = 100.0 * matches / max(1, len(aligned_ref))
        if identity >= min_identity and coverage >= min_coverage:
            start = window_start + int(best.coordinates[0][0])
            end = window_start + int(best.coordinates[0][-1]) - 1
            head = int(best.coordinates[1][0])
            tail = len(query) - int(best.coordinates[1][-1])
            if head <= MAX_EXTENSION and tail <= MAX_EXTENSION:
                extended_start = max(0, start - head)
                extended_end = min(len(target) - 1, end + tail)
                if extended_end - extended_start + 1 == len(query):
                    start, end = extended_start, extended_end
                    # Report identity over everything that is annotated, flanks included.
                    span = target[start : end + 1]
                    identity = 100.0 * sum(1 for a, b in zip(span, query) if a == b) / len(query)
            results.append((start, end, identity))
    del strand, feature
    return results


def _match_short(target: str, query: str) -> list[tuple[int, int, float]]:
    results = []
    start = target.find(query)
    while start != -1:
        results.append((start, start + len(query) - 1, 100.0))
        start = target.find(query, start + 1)
    return results


def detect(sequence: str, circular: bool = False, min_identity: float = DEFAULT_MIN_IDENTITY,
           min_coverage: float = DEFAULT_MIN_COVERAGE,
           features: tuple[CatalogueFeature, ...] | None = None) -> list[Match]:
    """Find catalogue features in ``sequence`` on both strands."""
    sequence = sequence.upper()
    length = len(sequence)
    if length < SEED:
        return []
    # Searching a doubled sequence finds features that cross the origin of a plasmid.
    target = sequence + sequence[: min(length, 20_000)] if circular else sequence
    found: list[Match] = []

    for feature in features if features is not None else catalogue():
        if feature.length > length:
            continue
        for strand in (1, -1):
            query = feature.sequence if strand == 1 else _reverse_complement(feature.sequence)
            coverage = COVERAGE_BY_TYPE.get(feature.type, min_coverage)
            hits = (_match_short(target, query) if feature.length <= SHORT_FEATURE
                    else _match_long(target, feature, query, strand, min_identity, coverage))
            for start, end, identity in hits:
                if start >= length:
                    continue  # the same hit also appears in the first copy
                found.append(Match(name=feature.name, type=feature.type, start=start % length,
                                   end=end % length if circular else min(end, length - 1),
                                   strand=strand, identity=identity, note=feature.note,
                                   source=feature.source))
    return _dedupe(found, length)


def _dedupe(matches: list[Match], length: int) -> list[Match]:
    """Drop matches that describe the same region: the same feature found twice by different
    seeds, or two catalogue entries (for example the 5′ and 3′ ITR) that are the same sequence.
    A short feature inside a longer one of a different kind is kept, because a tag inside a CDS
    is worth showing."""
    ordered = sorted(matches, key=lambda m: (-len(m.span(length)), -m.identity, m.name))
    kept: list[Match] = []
    for candidate in ordered:
        span = candidate.span(length)
        duplicate = False
        for existing in kept:
            other = existing.span(length)
            overlap = len(span & other)
            same_region = overlap >= 0.9 * len(span) and overlap >= 0.9 * len(other)
            # Variants of one feature (a full-length and a truncated ITR, say) are one match.
            variant = (existing.type == candidate.type
                       and overlap >= 0.9 * min(len(span), len(other)))
            if same_region or variant:
                duplicate = True
                break
        if not duplicate:
            kept.append(candidate)
    return sorted(kept, key=lambda m: (m.start, -len(m.span(length))))


def detect_payload(sequence: str, circular: bool, min_identity: float = DEFAULT_MIN_IDENTITY,
                   existing: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    """Matches, flagged against the features a document already carries."""
    length = len(sequence)
    matches = detect(sequence, circular, min_identity)
    already = []
    for feature in existing or []:
        start, end = feature.get("start"), feature.get("end")
        if isinstance(start, int) and isinstance(end, int):
            already.append((feature.get("name", "").lower(), start, end))

    rows = []
    for match in matches:
        row = match.as_dict(length)
        row["isNew"] = not any(
            name == match.name.lower() and abs(start - match.start) <= 5
            for name, start, _ in already)
        rows.append(row)
    return {"features": rows, "catalogue_size": len(catalogue()),
            "new_count": sum(1 for r in rows if r["isNew"])}
