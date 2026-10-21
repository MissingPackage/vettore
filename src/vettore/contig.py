"""Sanger contig: several reads of one clone, one consensus over the expected sequence.

Reference-guided only: every read is first aligned to the reference on its own
(``sanger.align_trace``), so the reads already share the reference's coordinates and a read
that comes from another plasmid never gets here.

The reads vote **event by event**. A stretch where any read has a deletion is one unit, and
each read covering all of it votes for what it shows across the whole stretch, so a read never
wins half of its own indel (a 2 bp deletion cannot become a 1 bp frameshift no read shows).
Every other base, and every slot between two bases for insertions, is a unit of its own. Each
call weighs its Phred quality, and the weights of reads showing the same thing add up, as in
phrap: three Q35 reads outvote one Q40 read.

A difference from the reference is **confirmed** when two reads show it at Q20 or better, or one
read at Q30 or better (one wrong call in 1,000, decided with the user on 2026-10-21), and no
read at Q20 or better shows something else there; otherwise it is to be checked in the trace.
A call that is not A, C, G or T (N, or a mixed IUPAC call) is no call: it neither votes nor
counts against the others.
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
#: Phred Q20: two reads at this quality confirm a difference together, and one read at this
#: quality that shows something else is enough to doubt it.
SUPPORTING_QUALITY = 20
CONFIRMING_READS = 2
BASES = set("ACGT")


@dataclass(frozen=True)
class Call:
    read: int
    value: str          # what the read shows over a unit (see ``_unit_call``)
    quality: int


@dataclass(frozen=True)
class Column:
    base: str           # the read's base at a reference position, "-" for a deletion
    quality: int


def _clear(text: str) -> bool:
    """Bases a read really called: A, C, G or T only (an empty insertion is clear too)."""
    return set(text) <= BASES


def _read_calls(read: AlignedRead,
                length: int) -> tuple[dict[int, Column], dict[int, tuple[str, int]]]:
    """One read's call at each reference position it covers, and its insertions by slot (the
    slot before reference position p is keyed p), each with its lowest quality."""
    from .sanger import _column_qualities

    qualities = _column_qualities(read.read_gapped, read.trace.quality)
    bases: dict[int, Column] = {}
    inserts: dict[int, tuple[str, int]] = {}
    position = read.offset
    for ref_base, read_base, quality in zip(read.reference_gapped, read.read_gapped, qualities):
        if ref_base == "-":
            key = position % length
            text, low = inserts.get(key, ("", quality))
            inserts[key] = (text + read_base, min(low, quality))
        else:
            bases[position % length] = Column(read_base, quality)
            position += 1
    return bases, inserts


def _vote(calls: list[Call], reference_value: str) -> tuple[str, bool]:
    """The call with the largest summed quality, and whether another call has as much (a tie:
    more reads, then the reference, pick the call shown, which is then not trusted)."""
    groups: dict[str, int] = {}
    for call in calls:
        groups[call.value] = groups.get(call.value, 0) + call.quality
    counts = {value: sum(1 for c in calls if c.value == value) for value in groups}
    winner = max(groups, key=lambda value: (groups[value], counts[value],
                                            value == reference_value))
    tied = any(value != winner and score == groups[winner] for value, score in groups.items())
    return winner, tied


def _confirmed(support: list[Call], against: list[Call]) -> bool:
    if any(c.quality >= SUPPORTING_QUALITY for c in against):
        return False
    return (sum(1 for c in support if c.quality >= SUPPORTING_QUALITY) >= CONFIRMING_READS
            or any(c.quality >= CONFIRMING_QUALITY for c in support))


def _runs(values: list[Any], circular: bool) -> list[tuple[int, int, Any]]:
    """Runs of equal values as (start, end, value), 0-based inclusive; on a circular sequence a
    run across the origin is one run (start > end)."""
    runs: list[list[Any]] = []
    for index, value in enumerate(values):
        if runs and runs[-1][2] == value:
            runs[-1][1] = index
        else:
            runs.append([index, index, value])
    if circular and len(runs) > 1 and runs[0][2] == runs[-1][2]:
        runs[0][0] = runs.pop()[0]
    return [(start, end, value) for start, end, value in runs]


def _positions(start: int, end: int, length: int) -> list[int]:
    return list(range(start, end + 1)) if start <= end else [*range(start, length),
                                                              *range(end + 1)]


def _units(per_read: list[tuple[dict[int, Column], dict]], length: int,
           circular: bool) -> list[list[int]]:
    """Reference positions grouped into voting units: each stretch where any read has a
    deletion is one unit, every other position a unit of its own."""
    deleted = [any(position in bases and bases[position].base == "-"
                   for bases, _ in per_read) for position in range(length)]
    units: list[list[int]] = []
    for start, end, is_deleted in _runs(deleted, circular):
        positions = _positions(start, end, length)
        units.extend([positions] if is_deleted else [[p] for p in positions])
    return units


def _merge(items: list[dict[str, Any]], length: int) -> list[dict[str, Any]]:
    """Adjacent deleted bases are one deletion, confirmed only if each of its bases is; on a
    circular reference this holds across the origin too."""
    merged: list[dict[str, Any]] = []
    for item in items:
        previous = merged[-1] if merged else None
        if (previous and previous["kind"] == item["kind"] == "deletion"
                and item["_start"] == (previous["_end"] + 1) % length):
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

    coverage = [len({strands[i] for i, (bases, _) in enumerate(per_read) if p in bases})
                for p in range(length)]
    depth = [sum(1 for bases, _ in per_read if p in bases) for p in range(length)]
    consensus_base: dict[int, str] = {}
    consensus_insert: dict[int, str] = {}
    differences: list[dict[str, Any]] = []
    disagreements: list[dict[str, Any]] = []
    unclear: set[int] = set()        # covered positions without a clear call (reference kept)
    doubtful: set[int] = set()       # positions where good reads conflict, or the vote is tied
    support = [0] * length           # reads whose clear call agrees with the consensus
    strong = [False] * length        # ... one of them at Q20 or better

    def listed(calls: list[Call]) -> list[dict[str, Any]]:
        return [{"read": names[c.read], "strand": strands[c.read], "call": c.value,
                 "quality": c.quality} for c in calls]

    def evidence(value: str, votes: list[Call]) -> dict[str, Any]:
        support = [c for c in votes if c.value == value]
        against = [c for c in votes if c.value != value]
        return {"reads": len(support), "quality": max(c.quality for c in support),
                "confirmed": _confirmed(support, against)}

    def slot_vote(position: int, touching: list[int]) -> None:
        """The slot before ``position``: a read votes on it when it covers both of its sides.
        No insertion weighs the weaker of the two bases around the slot."""
        previous = (position - 1) % length
        slot = []
        for i in touching:
            bases, inserts = per_read[i]
            if previous not in bases or position not in bases:
                continue
            text, quality = inserts.get(position, ("", min(bases[previous].quality,
                                                           bases[position].quality)))
            slot.append(Call(i, text, quality))
        votes = [c for c in slot if _clear(c.value)]
        if not votes:
            return
        inserted, tied = _vote(votes, "")
        if len({c.value for c in votes}) > 1:
            disagreements.append({"kind": "insertion", "position": position + 1,
                                  "reference": "", "consensus": inserted,
                                  "calls": listed(slot)})
        conflict = tied or any(c.value != inserted and c.quality >= SUPPORTING_QUALITY
                               for c in votes)
        if conflict:
            doubtful.add(position)
        if inserted:
            consensus_insert[position] = inserted
            found = evidence(inserted, votes)
            found["confirmed"] = found["confirmed"] and not conflict
            differences.append({"kind": "insertion", "position": position + 1,
                                "reference": "", "consensus": inserted, **found,
                                "_start": position, "_end": position})

    def column(i: int, unit: list[int], index: int) -> tuple[str, int, bool, bool]:
        """Read i at one position of a unit: its call (with what it inserts before it, when it
        also covers the previous position of the unit), the quality, whether the call is
        clear, and whether the inserted bases are part of the call."""
        bases, inserts = per_read[i]
        position = unit[index]
        here = bases[position]
        joined = index > 0 and unit[index - 1] in bases
        inserted, low = inserts.get(position, ("", 100)) if joined else ("", 100)
        clear = _clear(inserted) and (here.base == "-" or here.base in BASES)
        return inserted + here.base, min(here.quality, low), clear, joined

    for unit in _units(per_read, length, circular):
        first = unit[0]
        touching = [i for i, (bases, _) in enumerate(per_read) if any(p in bases for p in unit)]
        if not touching:
            continue
        if first > 0 or circular:
            slot_vote(first, touching)
        # The reads covering the whole unit vote for what they show across it, one string per
        # position ("" or the inserted bases, then the base or "-").
        covering = [i for i in touching if all(p in per_read[i][0] for p in unit)]
        votes: list[Call] = []
        for i in covering:
            columns = [column(i, unit, index) for index in range(len(unit))]
            if all(clear for _, _, clear, _ in columns):
                votes.append(Call(i, "|".join(call for call, _, _, _ in columns),
                                  min(quality for _, quality, _, _ in columns)))
        shown = [Call(i, "".join(column(i, unit, index)[0] for index, p in enumerate(unit)
                                 if p in per_read[i][0]),
                      min(column(i, unit, index)[1] for index, p in enumerate(unit)
                          if p in per_read[i][0]))
                 for i in touching]
        ref_value = "|".join(reference[p] for p in unit)
        if not votes:
            unclear.update(unit)
            disagreements.append({"kind": "no call", "position": first + 1,
                                  "reference": ref_value.replace("|", ""),
                                  "consensus": ref_value.replace("|", ""),
                                  "calls": listed(shown)})
            continue
        winner, tied = _vote(votes, ref_value)
        parts = winner.split("|")
        # Every clear call of every read touching the unit, position by position, counts: a
        # read that covers only part of it, or makes no clear call somewhere in it, still
        # speaks for or against the consensus where it does call.
        conflict = tied
        disputed = len({c.value for c in votes}) > 1
        for i in touching:
            for index, p in enumerate(unit):
                if p not in per_read[i][0]:
                    continue
                call, quality, clear, joined = column(i, unit, index)
                if not clear:
                    continue
                if call == (parts[index] if joined else parts[index][-1]):
                    support[p] += 1
                    strong[p] = strong[p] or quality >= SUPPORTING_QUALITY
                else:
                    disputed = True
                    conflict = conflict or quality >= SUPPORTING_QUALITY
        if disputed:
            disagreements.append({"kind": "base", "position": first + 1,
                                  "reference": ref_value.replace("|", ""),
                                  "consensus": winner.replace("|", ""),
                                  "calls": listed(shown)})
        if conflict:
            doubtful.update(unit)
        found = evidence(winner, votes)
        found["confirmed"] = found["confirmed"] and not conflict
        for p, part in zip(unit, parts):
            inserted, base = part[:-1], part[-1]
            consensus_base[p] = base
            if inserted:
                consensus_insert[p] = inserted
                differences.append({"kind": "insertion", "position": p + 1, "reference": "",
                                    "consensus": inserted, **found, "_start": p, "_end": p})
            if base != reference[p]:
                differences.append({"kind": "deletion" if base == "-" else "mismatch",
                                    "position": p + 1, "reference": reference[p],
                                    "consensus": "" if base == "-" else base, **found,
                                    "_start": p, "_end": p})
    differences = _merge(differences, length)

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

    # Verified: covered, a clear call that no good read contradicts, backed by a read at Q20 or
    # better, and no difference left to check. The count is the fewest reads agreeing with the
    # consensus over any base of the stretch.
    for p in range(length):
        if depth[p] and p not in unclear and not strong[p]:
            doubtful.add(p)          # only weak calls (under Q20) back this base
    for d in differences:
        if not d["confirmed"]:
            doubtful.update(_positions(d["_start"], d["_end"], length))
    doubtful -= unclear
    state = [0 if not depth[p] else 1 if p in doubtful or p in unclear else 2
             for p in range(length)]
    verified = [(start, end, min(support[p] for p in _positions(start, end, length)))
                for start, end, value in _runs(state, circular) if value == 2]
    covered = [(start, end) for start, end, value in _runs([min(1, d) for d in depth], circular)
               if value]
    out: dict[str, Any] = {
        "reads": len(reads),
        "consensus_length": len(consensus),
        "coverage": [{"start": start + 1, "end": end + 1, "strands": value}
                     for start, end, value in _runs(coverage, circular)],
        "covered": [{"start": start + 1, "end": end + 1} for start, end in covered],
        "verified": [{"start": start + 1, "end": end + 1, "reads": n}
                     for start, end, n in verified],
        "unclear": len(unclear),
        "doubtful": len(doubtful),
        "differences": [{k: v for k, v in d.items() if not k.startswith("_")}
                        for d in differences],
        "disagreements": disagreements,
    }
    if document is not None:
        def runs_of(positions: set[int]) -> list[tuple[int, int]]:
            return [(start, end) for start, end, value
                    in _runs([p in positions for p in range(length)], circular) if value]

        out["document"] = _consensus_document(document, reference, circular, consensus,
                                              start_of, kept, verified, runs_of(unclear),
                                              runs_of(doubtful), differences)
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
                        verified: list[tuple[int, int, int]],
                        unclear: list[tuple[int, int]],
                        doubtful: list[tuple[int, int]],
                        differences: list[dict[str, Any]]) -> dict[str, Any]:
    """The consensus as a new document: the reference's features moved over any insertion or
    deletion, the verified stretches annotated, every difference annotated where it is (those
    still to check say so), and the positions no read called clearly marked."""
    length = len(reference)
    new_length = max(1, len(consensus))
    features = _shifted(document.get("features"), length, circular, start_of, kept, new_length)
    primers = _shifted(document.get("primers"), length, circular, start_of, kept, new_length)

    def add(name: str, kind: str, span: tuple[int, int] | None, note: str) -> None:
        if span is None:
            return
        features.append({"id": uuid.uuid4().hex[:12], "name": name, "type": kind, "strand": 1,
                         "forward": True, "start": span[0], "end": span[1], "notes": _note(note)})

    for start, end, n in verified:
        add(f"Sanger verified ({_plural(n, 'read')})", "misc_feature",
            _new_range(start, end, start_of, kept, new_length, circular),
            f"Consensus of at least {_plural(n, 'Sanger read')} at every base, over reference "
            f"bases {start + 1:,}–{end + 1:,}")
    for start, end in unclear:
        add("Sanger: no clear call", "misc_feature",
            _new_range(start, end, start_of, kept, new_length, circular),
            f"No read makes a clear call over reference bases {start + 1:,}–{end + 1:,}: the "
            "reference is kept; check the traces")
    for start, end in doubtful:
        add("Sanger: check the trace", "misc_feature",
            _new_range(start, end, start_of, kept, new_length, circular),
            f"Reference bases {start + 1:,}–{end + 1:,}: the reads conflict here, or only "
            "calls under Q20 back the consensus; check the traces")
    for d in differences:
        evidence = (f"{d['reads']} reads" if d["reads"] > 1
                    else f"1 read at Q{d['quality']}")
        if not d["confirmed"]:
            evidence += "; not confirmed: check the trace"
        tag = "Sanger" if d["confirmed"] else "Sanger, check the trace"
        position = d["position"] - 1
        if d["kind"] == "mismatch":
            add(f"{d['reference']}→{d['consensus']} ({tag})", "variation",
                (start_of[position], start_of[position]),
                f"Reference {d['reference']} at {d['position']:,}, read {d['consensus']} "
                f"({evidence})")
        elif d["kind"] == "insertion":
            end = start_of[position] - 1
            add(f"+{d['consensus']} ({tag})", "variation",
                (end - len(d["consensus"]) + 1, end),
                f"{len(d['consensus'])} bp not in the reference, before its base "
                f"{d['position']:,} ({evidence})")
        else:
            # A deleted stretch has no bases left: the two bases around the gap carry it.
            after = start_of[position]
            span = (((after - 1) % new_length, after % new_length) if circular
                    else (max(0, after - 1), min(after, new_length - 1)))
            add(f"Δ{d['reference']} ({tag})", "variation", span,
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
