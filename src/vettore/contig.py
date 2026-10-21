"""Sanger contig: several reads of one clone, one consensus over the expected sequence.

Reference-guided only: every read is first aligned to the reference on its own
(``sanger.align_trace``), so the reads already share the reference's coordinates and a read
that comes from another plasmid never gets here. At each reference position (and in each slot
between two positions, for insertions) every read that covers it casts one call with the Phred
quality behind it; the consensus takes the call supported by the highest-quality read, as a
biologist reading the chromatograms side by side would.

A difference from the reference is **confirmed** when two reads show it, or one read at
Q ≥ 30 (one wrong call in 1,000); otherwise it is to be checked in the trace. An ``N`` is no
call: it neither votes nor counts against the others.
"""

from __future__ import annotations

import uuid
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from .convert import annotation_ranges

if TYPE_CHECKING:
    from .sanger import AlignedRead

#: Phred Q30 (decided with the user on 2026-10-21): one read this good is enough to trust a
#: difference on its own.
CONFIRMING_QUALITY = 30
CONFIRMING_READS = 2


@dataclass(frozen=True)
class Call:
    read: int
    value: str          # a base, "-" for a deletion, or the bases inserted in a slot ("" = none)
    quality: int


def _no_call(value: str) -> bool:
    return bool(value) and set(value) == {"N"}


def _read_calls(read: AlignedRead, length: int) -> tuple[dict[int, Call], dict[int, Call]]:
    """One read's call at each reference position it covers, and its insertions by slot (the
    slot before reference position p is keyed p)."""
    from .sanger import _column_qualities

    qualities = _column_qualities(read.read_gapped, read.trace.quality)
    bases: dict[int, Call] = {}
    inserts: dict[int, tuple[str, int]] = {}
    position = read.offset
    for ref_base, read_base, quality in zip(read.reference_gapped, read.read_gapped, qualities):
        if ref_base == "-":
            key = position % length
            text, low = inserts.get(key, ("", quality))
            inserts[key] = (text + read_base, min(low, quality))
        else:
            bases[position % length] = Call(-1, read_base, quality)
            position += 1
    return bases, {key: Call(-1, text, quality) for key, (text, quality) in inserts.items()}


def _vote(calls: list[Call], reference_value: str) -> str:
    """The call backed by the best read; more reads, then the reference, break a tie."""
    groups: dict[str, list[Call]] = {}
    for call in calls:
        groups.setdefault(call.value, []).append(call)
    return max(groups, key=lambda value: (max(c.quality for c in groups[value]),
                                          len(groups[value]), value == reference_value))


def _runs(values: list[int], circular: bool) -> list[tuple[int, int, int]]:
    """Runs of equal values as (start, end, value), 0-based inclusive; on a circular sequence a
    run across the origin is one run (start > end)."""
    runs: list[list[int]] = []
    for index, value in enumerate(values):
        if runs and runs[-1][2] == value:
            runs[-1][1] = index
        else:
            runs.append([index, index, value])
    if circular and len(runs) > 1 and runs[0][2] == runs[-1][2]:
        runs[0][0] = runs.pop()[0]
    return [(start, end, value) for start, end, value in runs]


def _in_run(position: int, start: int, end: int) -> bool:
    return start <= position <= end if start <= end else position >= start or position <= end


def _merge(items: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Adjacent deleted bases are one deletion, confirmed only if each of its bases is."""
    merged: list[dict[str, Any]] = []
    for item in items:
        previous = merged[-1] if merged else None
        if (previous and previous["kind"] == item["kind"] == "deletion"
                and item["position"] == previous["position"] + len(previous["reference"])):
            previous["reference"] += item["reference"]
            previous["reads"] = min(previous["reads"], item["reads"])
            previous["quality"] = min(previous["quality"], item["quality"])
            previous["confirmed"] = previous["confirmed"] and item["confirmed"]
            previous["_end"] = item["_end"]
        else:
            merged.append(dict(item))
    return merged


def contig_payload(reference: str, circular: bool, reads: list[AlignedRead],
                   document: dict[str, Any] | None = None) -> dict[str, Any]:
    """Consensus, per-strand coverage and the verdict on the consensus. With ``document`` (the
    reference's ``sequenceData``) it also carries the consensus as a new document."""
    length = len(reference)
    per_read = [_read_calls(read, length) for read in reads]
    names = [read.trace.name for read in reads]
    strands = [read.trace.strand for read in reads]

    coverage = [0] * length
    consensus_base: dict[int, str] = {}
    consensus_insert: dict[int, str] = {}
    differences: list[dict[str, Any]] = []
    disagreements: list[dict[str, Any]] = []

    def listed(calls: list[Call]) -> list[dict[str, Any]]:
        return [{"read": names[c.read], "strand": strands[c.read], "call": c.value,
                 "quality": c.quality} for c in calls]

    def difference(kind: str, position: int, ref: str, value: str,
                   calls: list[Call]) -> dict[str, Any]:
        """``value`` is the winning call; a deletion ("-") is shown with no bases."""
        support = [c for c in calls if c.value == value]
        best = max(c.quality for c in support)
        return {"kind": kind, "position": position + 1, "reference": ref,
                "consensus": "" if value == "-" else value,
                "reads": len(support), "quality": best,
                "confirmed": len(support) >= CONFIRMING_READS or best >= CONFIRMING_QUALITY,
                "_end": position}

    for position in range(length):
        covering = [i for i, (bases, _) in enumerate(per_read) if position in bases]
        coverage[position] = len({strands[i] for i in covering})
        if not covering:
            continue
        # The slot before this base: a read votes on it when it covers both of its sides.
        previous = (position - 1) % length
        if position > 0 or circular:
            slot = [Call(i, per_read[i][1][position].value if position in per_read[i][1] else "",
                         per_read[i][1][position].quality if position in per_read[i][1]
                         else min(per_read[i][0][previous].quality,
                                  per_read[i][0][position].quality))
                    for i in covering if previous in per_read[i][0]]
            votes = [c for c in slot if not _no_call(c.value)]
            if votes:
                inserted = _vote(votes, "")
                if len({c.value for c in votes}) > 1:
                    disagreements.append({"kind": "insertion", "position": position + 1,
                                          "reference": "", "consensus": inserted,
                                          "calls": listed(slot)})
                if inserted:
                    consensus_insert[position] = inserted
                    differences.append(difference("insertion", position, "", inserted, votes))
        calls = [Call(i, per_read[i][0][position].value, per_read[i][0][position].quality)
                 for i in covering]
        votes = [c for c in calls if not _no_call(c.value)]
        ref_base = reference[position]
        if not votes:
            disagreements.append({"kind": "no call", "position": position + 1,
                                  "reference": ref_base, "consensus": ref_base,
                                  "calls": listed(calls)})
            continue
        base = _vote(votes, ref_base)
        consensus_base[position] = base
        if len({c.value for c in votes}) > 1:
            disagreements.append({"kind": "base", "position": position + 1,
                                  "reference": ref_base, "consensus": base,
                                  "calls": listed(calls)})
        if base != ref_base:
            differences.append(difference("deletion" if base == "-" else "mismatch", position,
                                          ref_base, base, votes))
    differences = _merge(differences)

    # The consensus: the reference where no read reaches, the reads' verdict elsewhere.
    pieces: list[str] = []
    start_of = [0] * length          # where each reference base lands in the consensus
    kept = [True] * length
    cursor = 0
    for position in range(length):
        inserted = consensus_insert.get(position, "")
        pieces.append(inserted)
        cursor += len(inserted)
        start_of[position] = cursor
        base = consensus_base.get(position, reference[position])
        if base == "-":
            kept[position] = False
        else:
            pieces.append(base)
            cursor += 1
    consensus = "".join(pieces)

    spans = [set(bases) for bases, _ in per_read]
    covered = [(start, end, sum(1 for span in spans
                                if any(_in_run(p, start, end) for p in span)))
               for start, end, value in _runs([min(1, c) for c in coverage], circular) if value]
    out: dict[str, Any] = {
        "reads": len(reads),
        "consensus_length": len(consensus),
        "coverage": [{"start": start + 1, "end": end + 1, "strands": value}
                     for start, end, value in _runs(coverage, circular)],
        "covered": [{"start": start + 1, "end": end + 1, "reads": n} for start, end, n in covered],
        "differences": [{k: v for k, v in d.items() if k != "_end"} for d in differences],
        "disagreements": disagreements,
    }
    if document is not None:
        out["document"] = _consensus_document(document, reference, circular, consensus,
                                              start_of, kept, covered, differences)
    return out


def _new_range(start: int, end: int, start_of: list[int], kept: list[bool],
               new_length: int, circular: bool) -> tuple[int, int] | None:
    """Where the reference bases [start, end] (inclusive) lie in the consensus; None when the
    reads show all of them deleted."""
    new_start = start_of[start]
    new_end = start_of[end] + (1 if kept[end] else 0) - 1
    if start <= end:
        if new_end < new_start:
            return None
        return new_start % new_length, new_end
    return new_start % new_length, new_end % new_length if circular else new_end


def _shifted(annotations: Any, length: int, circular: bool, start_of: list[int],
             kept: list[bool], new_length: int) -> list[dict[str, Any]]:
    values = (list(annotations.values()) if isinstance(annotations, dict)
              else annotations if isinstance(annotations, list) else [])
    out: list[dict[str, Any]] = []
    for ann in values:
        ranges = annotation_ranges(ann, length, circular)
        if ranges is None:
            continue
        moved = [_new_range(s, e, start_of, kept, new_length, circular) for s, e in ranges]
        moved = [r for r in moved if r is not None]
        if not moved:
            continue
        copy = {**ann, "start": moved[0][0], "end": moved[-1][1]}
        if "locations" in ann:
            copy["locations"] = [{"start": s, "end": e} for s, e in moved]
        out.append(copy)
    return out


def _note(text: str) -> dict[str, list[str]]:
    return {"note": [text]}


def _plural(count: int, noun: str) -> str:
    return f"{count} {noun}{'' if count == 1 else 's'}"


def _consensus_document(document: dict[str, Any], reference: str, circular: bool,
                        consensus: str, start_of: list[int], kept: list[bool],
                        covered: list[tuple[int, int, int]],
                        differences: list[dict[str, Any]]) -> dict[str, Any]:
    """The consensus as a new document: the reference's features moved over any insertion or
    deletion, the stretches the reads cover annotated as verified, the confirmed differences
    annotated where they are."""
    length = len(reference)
    new_length = max(1, len(consensus))
    features = _shifted(document.get("features"), length, circular, start_of, kept, new_length)
    primers = _shifted(document.get("primers"), length, circular, start_of, kept, new_length)

    def add(name: str, kind: str, span: tuple[int, int] | None, note: str) -> None:
        if span is None:
            return
        features.append({"id": uuid.uuid4().hex[:12], "name": name, "type": kind, "strand": 1,
                         "forward": True, "start": span[0], "end": span[1], "notes": _note(note)})

    for start, end, n in covered:
        add(f"Sanger verified ({_plural(n, 'read')})", "misc_feature",
            _new_range(start, end, start_of, kept, new_length, circular),
            f"Consensus of {_plural(n, 'Sanger read')} over reference bases "
            f"{start + 1:,}–{end + 1:,}")
    for d in differences:
        if not d["confirmed"]:
            continue
        evidence = (f"{d['reads']} reads" if d["reads"] > 1
                    else f"1 read at Q{d['quality']}")
        position = d["position"] - 1
        if d["kind"] == "mismatch":
            add(f"{d['reference']}→{d['consensus']} (Sanger)", "variation",
                (start_of[position], start_of[position]),
                f"Reference {d['reference']} at {d['position']:,}, read {d['consensus']} "
                f"({evidence})")
        elif d["kind"] == "insertion":
            end = start_of[position] - 1
            add(f"+{d['consensus']} (Sanger)", "variation",
                (end - len(d["consensus"]) + 1, end),
                f"{len(d['consensus'])} bp not in the reference, before its base "
                f"{d['position']:,} ({evidence})")
        else:
            # A deleted stretch has no bases left: the two bases around the gap carry it.
            after = start_of[position]
            span = (((after - 1) % new_length, after % new_length) if circular
                    else (max(0, after - 1), min(after, new_length - 1)))
            add(f"Δ{d['reference']} (Sanger)", "variation", span,
                f"Reference bases {d['position']:,}–{d['_end'] + 1:,} ({d['reference']}) are "
                f"missing between these two bases ({evidence})")

    name = str(document.get("name") or "Reference")
    out = {key: value for key, value in document.items()
           if key not in {"features", "primers", "parts", "translations", "sequence", "name",
                          "description", "id"}}
    out.update({
        "name": f"{name} consensus",
        "description": f"Sanger consensus over {name}",
        "sequence": consensus,
        "circular": circular,
        "features": features,
        "primers": primers,
    })
    return out
