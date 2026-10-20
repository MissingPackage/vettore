"""Sanger traces: reading .ab1, quality trimming and alignment to a reference.

Verification of a clone is the daily use: import the reads, align them to the plasmid you
expect, and look at the discrepancies. The trace is read with Biopython, trimmed with Mott's
algorithm (re-implemented here so the stringency is adjustable, which Biopython's ``abi-trim``
does not allow) and aligned with ``Bio.Align.PairwiseAligner`` in local mode, which is the
right mode because a read covers only part of the plasmid.

The chromatogram is converted to one trace bundle per called base, the shape Open Vector
Editor needs to draw it (``baseTraces``); doing it here rather than in the browser means a
reverse-strand read can be reverse-complemented together with its peaks.
"""

from __future__ import annotations

import io
from dataclasses import dataclass, field
from itertools import pairwise
from typing import Any

from Bio import Align, SeqIO
from Bio.Seq import Seq

#: Mott's algorithm: a base scores ``cutoff - 10**(-Q/10)``; the best-scoring run is kept.
DEFAULT_CUTOFF = 0.05          # ≈ Q13, Biopython's own default
MIN_SEGMENT = 20
COMPLEMENT = {"aTrace": "tTrace", "tTrace": "aTrace", "gTrace": "cTrace", "cTrace": "gTrace"}


class TraceError(ValueError):
    """A .ab1 file that cannot be used; the message is shown to the user."""


@dataclass
class Trace:
    name: str
    sequence: str
    quality: list[int]
    base_traces: list[dict[str, list[int]]]
    strand: int = 1
    annotations: dict[str, Any] = field(default_factory=dict)

    def slice(self, start: int, end: int) -> Trace:
        """Keep bases [start, end) together with their peaks."""
        return Trace(self.name, self.sequence[start:end], self.quality[start:end],
                     self.base_traces[start:end], self.strand, self.annotations)

    def reverse_complement(self) -> Trace:
        traces = [{COMPLEMENT[channel]: list(reversed(values)) for channel, values in base.items()}
                  for base in reversed(self.base_traces)]
        return Trace(self.name, str(Seq(self.sequence).reverse_complement()),
                     list(reversed(self.quality)), traces, -self.strand, self.annotations)

    def chromatogram(self) -> dict[str, Any]:
        """``chromatogramData`` in Open Vector Editor's shape."""
        position = 0
        base_pos = []
        for base in self.base_traces:
            width = len(base["aTrace"])
            base_pos.append(position + width // 2)
            position += width
        return {
            "baseTraces": self.base_traces,
            "basePos": base_pos,
            "baseCalls": list(self.sequence),
            "qualNums": self.quality,
            "aTrace": [v for base in self.base_traces for v in base["aTrace"]],
            "tTrace": [v for base in self.base_traces for v in base["tTrace"]],
            "gTrace": [v for base in self.base_traces for v in base["gTrace"]],
            "cTrace": [v for base in self.base_traces for v in base["cTrace"]],
        }


def _per_base_traces(channels: dict[str, list[int]], base_pos: list[int]) -> list[dict]:
    """Split each channel at the midpoints between peaks, one bundle per called base
    (the same binning Open Vector Editor's own converter uses)."""
    edges = [0]
    for left, right in pairwise(base_pos):
        edges.append(-(-(left + right) // 2))  # ceil of the midpoint
    edges.append(len(channels["A"]))
    if len(edges) > 3:  # make the first and last bins as wide as their neighbours
        if edges[1] - edges[0] > edges[2] - edges[1]:
            edges[0] = edges[1] - (edges[2] - edges[1])
        if edges[-1] - edges[-2] > edges[-2] - edges[-3]:
            edges[-1] = edges[-2] + (edges[-2] - edges[-3]) + 1
    return [{"aTrace": channels["A"][start:end], "tTrace": channels["T"][start:end],
             "gTrace": channels["G"][start:end], "cTrace": channels["C"][start:end]}
            for start, end in pairwise(edges)]


def read_trace(filename: str, content: bytes) -> Trace:
    try:
        record = SeqIO.read(io.BytesIO(content), "abi")
    except Exception as exc:
        raise TraceError(f"{filename} is not a readable .ab1 trace ({exc})") from exc
    raw = record.annotations.get("abif_raw", {})
    try:
        order = raw["FWO_1"].decode()                       # channel order, usually GATC
        channels = {base: list(raw[f"DATA{9 + index}"]) for index, base in enumerate(order)}
        base_pos = list(raw["PLOC2"])
    except KeyError as exc:
        raise TraceError(f"{filename} has no usable chromatogram (missing {exc})") from exc
    sequence = str(record.seq).upper()
    quality = list(record.letter_annotations.get("phred_quality") or [0] * len(sequence))
    if not sequence:
        raise TraceError(f"{filename} contains no base calls")
    if len(base_pos) != len(sequence):
        raise TraceError(f"{filename} has {len(base_pos)} peaks for {len(sequence)} bases")
    friendly = {key: str(record.annotations.get(key, "")) for key in
                ("machine_model", "dye", "polymer", "run_start", "sample_well")}
    return Trace(name=record.id or filename, sequence=sequence, quality=quality,
                 base_traces=_per_base_traces(channels, base_pos), annotations=friendly)


def mott_trim(quality: list[int], cutoff: float = DEFAULT_CUTOFF,
              min_segment: int = MIN_SEGMENT) -> tuple[int, int]:
    """Richard Mott's trimming algorithm: the highest-scoring run of good-quality bases.
    Returns [start, end); an empty range means nothing passed the threshold."""
    if len(quality) <= min_segment:
        return 0, len(quality)
    score = 0.0
    best = 0.0
    best_range = (0, 0)
    start = 0
    for index, value in enumerate(quality):
        score += cutoff - 10 ** (value / -10.0)
        if score < 0:
            score = 0.0
            start = index + 1
        elif score > best:
            best = score
            best_range = (start, index + 1)
    return best_range if best_range[1] - best_range[0] >= min_segment else (0, 0)


def _aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner(scoring="blastn")
    aligner.mode = "local"
    return aligner


def _gapped(target: str, query: str, coordinates) -> tuple[str, str, int, int, int, int]:
    """Gapped reference and read over the aligned span, plus the span's coordinates."""
    columns_ref: list[str] = []
    columns_read: list[str] = []
    for index in range(coordinates.shape[1] - 1):
        t0, t1 = int(coordinates[0][index]), int(coordinates[0][index + 1])
        q0, q1 = int(coordinates[1][index]), int(coordinates[1][index + 1])
        if t1 > t0 and q1 > q0:
            columns_ref.append(target[t0:t1])
            columns_read.append(query[q0:q1])
        elif t1 > t0:                                  # missing in the read
            columns_ref.append(target[t0:t1])
            columns_read.append("-" * (t1 - t0))
        else:                                          # extra bases in the read
            columns_ref.append("-" * (q1 - q0))
            columns_read.append(query[q0:q1])
    return ("".join(columns_ref), "".join(columns_read), int(coordinates[0][0]),
            int(coordinates[0][-1]), int(coordinates[1][0]), int(coordinates[1][-1]))


def _discrepancies(reference_gapped: str, read_gapped: str, offset: int,
                   reference_length: int) -> list[dict[str, Any]]:
    """Differences, with 1-based coordinates on the reference."""
    out: list[dict[str, Any]] = []
    position = offset
    for ref_base, read_base in zip(reference_gapped, read_gapped):
        if ref_base == "-":
            out.append({"kind": "insertion", "position": (position % reference_length) + 1,
                        "reference": "", "read": read_base})
            continue
        if read_base != "-" and read_base != ref_base:
            out.append({"kind": "mismatch", "position": (position % reference_length) + 1,
                        "reference": ref_base, "read": read_base})
        elif read_base == "-":
            out.append({"kind": "deletion", "position": (position % reference_length) + 1,
                        "reference": ref_base, "read": ""})
        position += 1
    return _merge_runs(out)


def _merge_runs(items: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Consecutive insertions or deletions are one event, as a biologist reads them."""
    merged: list[dict[str, Any]] = []
    for item in items:
        previous = merged[-1] if merged else None
        if (previous and previous["kind"] == item["kind"] != "mismatch"
                and item["position"] - previous["position"] <= 1):
            previous["reference"] += item["reference"]
            previous["read"] += item["read"]
            if item["kind"] == "deletion":
                previous["position"] = min(previous["position"], item["position"])
        else:
            merged.append(dict(item))
    return merged


@dataclass
class AlignedRead:
    trace: Trace
    reference_gapped: str
    read_gapped: str
    start: int                 # 0-based inclusive on the reference
    end: int                   # 0-based inclusive
    identity: float
    discrepancies: list[dict[str, Any]]
    trimmed: tuple[int, int]
    original_length: int
    spans_origin: bool
    offset: int = 0            # first reference position of the alignment, before wrapping


def align_trace(reference: str, circular: bool, trace: Trace, cutoff: float = DEFAULT_CUTOFF,
                min_quality_bases: int = MIN_SEGMENT) -> AlignedRead:
    start, end = mott_trim(trace.quality, cutoff, min_quality_bases)
    if end - start < min_quality_bases:
        raise TraceError(f"“{trace.name}” has no stretch of at least {min_quality_bases} bases "
                         "above the quality threshold; lower the trim stringency to see it")
    trimmed = trace.slice(start, end)
    target = reference + reference if circular else reference
    aligner = _aligner()
    candidates = []
    for candidate in (trimmed, trimmed.reverse_complement()):
        try:
            alignments = aligner.align(target, candidate.sequence)
        except ValueError as exc:
            raise TraceError(f"“{trace.name}” could not be aligned ({exc})") from exc
        if len(alignments):
            candidates.append((alignments[0].score, candidate, alignments[0]))
    if not candidates:
        raise TraceError(f"“{trace.name}” does not align to the reference")
    _, candidate, alignment = max(candidates, key=lambda item: item[0])

    ref_mid, read_mid, t_start, t_end, q_start, q_end = _gapped(
        target, candidate.sequence, alignment.coordinates)
    length = len(reference)
    spans_origin = circular and t_start < length < t_end
    if circular and t_start >= length:
        t_start -= length
        t_end -= length

    aligned = candidate.slice(q_start, q_end)
    matches = sum(1 for a, b in zip(ref_mid, read_mid) if a == b and a != "-")
    identity = 100.0 * matches / max(1, len(ref_mid))
    discrepancies = _discrepancies(ref_mid, read_mid, t_start, length)
    return AlignedRead(trace=aligned, reference_gapped=ref_mid, read_gapped=read_mid,
                       start=t_start % length, end=(t_end - 1) % length, identity=identity,
                       discrepancies=discrepancies, trimmed=(start, end),
                       original_length=len(trace.sequence), spans_origin=spans_origin,
                       offset=t_start)


def _columns(reference_gapped: str, read_gapped: str, offset: int,
             length: int) -> tuple[dict[int, str], dict[int, str]]:
    """Per reference position: the read base aligned to it, and any bases inserted before it."""
    bases: dict[int, str] = {}
    inserts: dict[int, str] = {}
    position = offset
    for ref_base, read_base in zip(reference_gapped, read_gapped):
        if ref_base == "-":
            key = position % length
            inserts[key] = inserts.get(key, "") + read_base
        else:
            bases[position % length] = read_base
            position += 1
    return bases, inserts


def merge_tracks(reference: str, reads: list[AlignedRead]) -> tuple[str, list[str]]:
    """Put the reference and every read on one coordinate system: an insertion found in any
    read opens a column in all tracks, so the reads can be stacked under the reference."""
    return merge_gapped(reference, [(r.reference_gapped, r.read_gapped, r.offset) for r in reads])


def merge_gapped(reference: str,
                 alignments: list[tuple[str, str, int]]) -> tuple[str, list[str]]:
    """``merge_tracks`` for any pairwise alignment to the reference, each given as
    (gapped reference span, gapped other sequence, first reference position of the span).
    Spans may run past the end of a circular reference; positions wrap around the origin."""
    length = len(reference)
    columns = [_columns(ref, other, offset, length) for ref, other, offset in alignments]
    reference_gapped: list[str] = []
    tracks: list[list[str]] = [[] for _ in alignments]
    for position in range(length):
        width = max((len(inserts.get(position, "")) for _, inserts in columns), default=0)
        if width:
            reference_gapped.append("-" * width)
            for index, (_, inserts) in enumerate(columns):
                tracks[index].append(inserts.get(position, "").ljust(width, "-"))
        reference_gapped.append(reference[position])
        for index, (bases, _) in enumerate(columns):
            tracks[index].append(bases.get(position, "-"))
    return "".join(reference_gapped), ["".join(track) for track in tracks]


def alignment_payload(reference_name: str, reference: str, circular: bool,
                      reads: list[tuple[str, bytes]], cutoff: float = DEFAULT_CUTOFF) -> dict:
    """Everything the UI needs: one stacked alignment plus a per-read summary."""
    aligned_reads: list[AlignedRead] = []
    summaries: list[dict[str, Any]] = []
    warnings: list[str] = []
    length = len(reference)

    for filename, content in reads:
        try:
            trace = read_trace(filename, content)
            aligned = align_trace(reference, circular, trace, cutoff)
        except TraceError as exc:
            warnings.append(str(exc))
            continue
        aligned_reads.append(aligned)
        summaries.append({
            "name": aligned.trace.name,
            "file": filename,
            "strand": aligned.trace.strand,
            "start": aligned.start + 1,
            "end": aligned.end + 1,
            "coverage": len(aligned.read_gapped.replace("-", "")),
            "identity": round(aligned.identity, 2),
            "discrepancies": aligned.discrepancies,
            "trimmed_bases": aligned.original_length - (aligned.trimmed[1] - aligned.trimmed[0]),
            "original_length": aligned.original_length,
            "spans_origin": aligned.spans_origin,
            "machine": aligned.trace.annotations.get("machine_model", ""),
        })

    reference_gapped, read_tracks = merge_tracks(reference, aligned_reads)
    tracks: list[dict[str, Any]] = []
    if aligned_reads:
        tracks.append({
            "sequenceData": {"id": "reference", "name": reference_name, "sequence": reference,
                             "circular": False},
            "alignmentData": {"id": "reference-alignment", "sequence": reference_gapped},
        })
        for index, (aligned, gapped) in enumerate(zip(aligned_reads, read_tracks)):
            chromatogram = aligned.trace.chromatogram()
            tracks.append({
                # Open Vector Editor reads the chromatogram from the track in the alignment
                # view and from sequenceData in the row view, so it is given in both places.
                "sequenceData": {"id": f"read-{index}", "name": aligned.trace.name,
                                 "sequence": aligned.trace.sequence,
                                 "chromatogramData": chromatogram},
                "chromatogramData": chromatogram,
                "alignmentData": {"id": f"read-alignment-{index}", "sequence": gapped},
            })

    covered: set[int] = set()
    for summary in summaries:
        start, end = summary["start"] - 1, summary["end"] - 1
        covered.update(range(start, end + 1) if start <= end
                       else list(range(start, length)) + list(range(end + 1)))
    first_covered = min((s["start"] for s in summaries), default=1)
    return {
        "alignmentTracks": tracks,
        "reads": summaries,
        "warnings": warnings,
        "reference": {"name": reference_name, "length": length, "circular": circular},
        "coverage_percent": round(100.0 * len(covered) / max(1, length), 1),
        "first_covered": first_covered,
    }
