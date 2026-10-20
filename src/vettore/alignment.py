"""Multiple alignment of open documents against a reference. Routes under /api/alignment.

The question asked of several related sequences (capsid serotypes, variants of a construct, a
clone against its design) is "where do they differ from this one?". This is therefore a *star*
alignment: every sequence is aligned pairwise to the chosen reference with
``Bio.Align.PairwiseAligner``, and the pairwise alignments are stacked on the reference's
coordinates, an insertion in any sequence opening a column in every track (the same merge the
Sanger view uses, ``sanger.merge_gapped``). Differences are always reported in reference
coordinates. The other sequences are not aligned to each other, as a progressive multiple
alignment would do.

For each sequence:

1. **Strand and origin.** Every 11-base word shared with the reference votes for a diagonal
   (the k-tuple method of FASTA: Pearson & Lipman 1988, PNAS 85:2444). Both strands are
   tried; the one whose best diagonal collects more votes wins. On a circular sequence that
   diagonal says where its origin falls relative to the reference, so the sequence is rotated
   to start where the reference starts. When the reference is circular and the sequence is
   linear (a fragment), the reference is rotated internally instead, so a fragment that
   crosses its origin aligns in one piece; positions are reported on the original numbering.
2. **Mode, per pair.** The shared words are chained (longest collinear chain); runs of
   nearby chained words mark homologous stretches, and how much of the shorter sequence
   they span decides:

   - *semi-global* (global, with end gaps free on both sequences) when the runs span at
     least half of the shorter sequence. The two are related end to end: serotypes, variants,
     or a fragment contained in a plasmid. Every base of the overlap is aligned and internal
     indels are scored as in a global alignment, while a shorter sequence can sit anywhere in
     a longer one and ragged ends cost nothing;
   - *global* instead of semi-global when both sequences are circular: circles have no ends,
     so both are cut inside the same shared word and aligned end to end with every gap
     scored. With free end gaps a difference next to the cut would pass for an overhang;
   - *local* otherwise: only part of the sequence is related (a plasmid that shares one
     cassette with the reference). A global alignment would pair unrelated bases with each
     other and report them as differences; the local alignment reports the shared region
     only, and the coverage says how much of the reference it is.
3. Columns where one sequence runs past the other's end are not part of the comparison: they
   are neither differences nor counted in the identity, and are left out of the view.

Scores: NCBI BLAST+ ``blastn`` task defaults (reward 2, penalty -3, gap open 5, extend 2;
BLAST+ user manual, NBK279684), Biopython's ``blastn`` preset. Identity is computed as BLAST
does: identical columns over all columns of the aligned region, gaps included.

Limits: the dynamic programming matrix grows with the product of the two lengths (about
2 bytes per cell: 15 kb x 15 kb is ~460 MB and ~4 s on a laptop core), hence
``MAX_LENGTH``; the run time grows with the number of sequences, hence ``MAX_SEQUENCES``.
"""

from __future__ import annotations

from bisect import bisect_left
from collections import Counter
from dataclasses import dataclass
from typing import Any

from Bio import Align
from Bio.Seq import Seq
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from .sanger import merge_gapped

router = APIRouter(prefix="/api/alignment")

MAX_SEQUENCES = 10           # reference included
MAX_LENGTH = 15_000          # bp, each sequence
MIN_LENGTH = 15              # bp; shorter than this nothing meaningful can be anchored
WORD = 11                    # word length for the strand / origin vote
MAX_WORD_COPIES = 4          # words repeated more often in the reference (poly-A, repeats) do not vote
DIAGONAL_SLACK = 3           # votes from diagonals this close count together (small indels)
MIN_VOTES = 2                # fewer shared words than this on both strands: unrelated
SEMI_GLOBAL_MIN_CHAIN = 0.5  # fraction of the shorter sequence the chain must span
RUN_GAP = 500                # chained words further apart than this start a new run
RUN_DRIFT = 50               # ... and so does a diagonal shift larger than this (an indel)

IUPAC_DNA = set("ACGTRYSWKMBDHVN")


class AlignmentError(ValueError):
    """A user-facing error: the message is shown in the dialog as is."""


@dataclass
class Placement:
    strand: int              # 1, or -1 when the sequence is reverse-complemented
    rotation: int            # circular sequence: its base (oriented, 0-based) facing reference base 1
    reference_shift: int     # internal rotation of a circular reference (positions map back)
    mode: str                # "global", "semi-global" or "local"
    votes: int


@dataclass
class PairResult:
    name: str
    length: int
    placement: Placement
    reference_gapped: str    # over the aligned region only
    other_gapped: str
    offset: int              # first reference position (0-based) of the aligned region
    reference_span: int      # reference bases in the aligned region
    other_aligned: str       # the sequence's bases in the aligned region (oriented)
    identity: float
    differences: list[dict[str, Any]]


# ---- input -----------------------------------------------------------------------------

def clean_sequence(sequence: str, name: str) -> str:
    seq = (sequence or "").upper().replace("U", "T")
    bad = sorted({c for c in seq if c not in IUPAC_DNA})
    if bad:
        shown = ", ".join(repr(c) for c in bad[:5])
        raise AlignmentError(f"“{name}” contains characters that are not DNA ({shown}). "
                             "Remove them in the Sequence view and align again.")
    return seq


def check_limits(reference_name: str, reference: str, others: list[tuple[str, str]]) -> None:
    if len(reference) < MIN_LENGTH:
        raise AlignmentError(f"The reference “{reference_name}” is too short to align to "
                             f"(at least {MIN_LENGTH} bp are needed).")
    if not others:
        raise AlignmentError("Choose at least one sequence to align to the reference.")
    if len(others) + 1 > MAX_SEQUENCES:
        raise AlignmentError(f"{len(others) + 1} sequences were chosen, but at most "
                             f"{MAX_SEQUENCES} can be aligned at once (reference included). "
                             "Untick some and align again.")
    for name, seq in [(reference_name, reference), *others]:
        if len(seq) > MAX_LENGTH:
            raise AlignmentError(f"“{name}” is {len(seq):,} bp; sequences longer than "
                                 f"{MAX_LENGTH:,} bp cannot be aligned here. Align the region "
                                 "you need instead (Actions → New File from Selection).")


def rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def rotate(seq: str, start: int) -> str:
    start %= max(1, len(seq))
    return seq[start:] + seq[:start]


# ---- step 1: strand and origin from shared words ---------------------------------------

def _words(seq: str, circular: bool) -> list[tuple[int, str]]:
    if len(seq) < WORD:
        return []
    text = seq + seq[:WORD - 1] if circular else seq
    count = len(seq) if circular else len(seq) - WORD + 1
    return [(i, text[i:i + WORD]) for i in range(count)]


def _reference_index(reference: str, circular: bool) -> dict[str, list[int]]:
    index: dict[str, list[int]] = {}
    for position, word in _words(reference, circular):
        if set(word) <= set("ACGT"):
            index.setdefault(word, []).append(position)
    return {w: p for w, p in index.items() if len(p) <= MAX_WORD_COPIES}


def _hits(index: dict[str, list[int]], query: str, circular: bool) -> list[tuple[int, int]]:
    """(query position, reference position) of every shared word."""
    return [(q, r) for q, word in _words(query, circular) for r in index.get(word, ())]


def _diagonal(q: int, r: int, ref_len: int, query_len: int, ref_circular: bool,
              query_circular: bool) -> int:
    if query_circular:               # rotation of the query that puts q opposite r
        return (q - r) % query_len
    if ref_circular:                 # reference position opposite the query's first base
        return (r - q) % ref_len
    return r - q


def _best_diagonal(diagonals: list[int], modulus: int | None) -> tuple[int, int]:
    """Diagonal with the most votes, counting votes within DIAGONAL_SLACK of it."""
    if not diagonals:
        return 0, 0
    counts = Counter(diagonals)
    best, best_votes = 0, -1
    for key in counts:
        votes = 0
        for delta in range(-DIAGONAL_SLACK, DIAGONAL_SLACK + 1):
            other = key + delta
            votes += counts.get(other % modulus if modulus else other, 0)
        if votes > best_votes or (votes == best_votes and key < best):
            best, best_votes = key, votes
    return best, best_votes


def _near(key: int, diagonal: int, modulus: int | None) -> bool:
    distance = abs(key - diagonal)
    if modulus:
        distance = min(distance, modulus - distance)
    return distance <= DIAGONAL_SLACK


def _chain(points: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Longest chain of shared words increasing on both sequences (patience sorting)."""
    if not points:
        return []
    points = sorted(points, key=lambda p: (p[0], -p[1]))
    tails: list[int] = []              # smallest reference end of a chain of each length
    tail_index: list[int] = []
    parent = [-1] * len(points)
    for i, (_, r) in enumerate(points):
        k = bisect_left(tails, r)
        if k == len(tails):
            tails.append(r)
            tail_index.append(i)
        else:
            tails[k] = r
            tail_index[k] = i
        parent[i] = tail_index[k - 1] if k else -1
    chain = []
    i = tail_index[-1]
    while i != -1:
        chain.append(points[i])
        i = parent[i]
    return chain[::-1]


def _run_coverage(chain: list[tuple[int, int]], on_query: bool) -> int:
    """Bases of one sequence spanned by runs of the chain: words that follow each other within
    RUN_GAP bases on a nearly constant diagonal. A word shared by chance is isolated and adds
    nothing, however far apart it sits; homologous regions give runs, and a large indel only
    splits a run in two."""
    covered = 0
    run_start = 0
    for i in range(1, len(chain) + 1):
        if i < len(chain):
            (q0, r0), (q1, r1) = chain[i - 1], chain[i]
            if q1 - q0 <= RUN_GAP and r1 - r0 <= RUN_GAP and abs((r1 - q1) - (r0 - q0)) <= RUN_DRIFT:
                continue
        if i - run_start >= 2:
            first, last = chain[run_start], chain[i - 1]
            axis = 0 if on_query else 1
            covered += last[axis] - first[axis] + WORD
        run_start = i
    return covered


def place(reference: str, ref_circular: bool, query: str, query_circular: bool,
          index: dict[str, list[int]] | None = None) -> tuple[Placement, str, str] | None:
    """Strand, rotation and mode for one sequence. Returns the placement together with the
    reference and the sequence as they go into the aligner, or None when the two share no
    word on either strand."""
    index = index if index is not None else _reference_index(reference, ref_circular)
    ref_len, query_len = len(reference), len(query)
    modulus = query_len if query_circular else ref_len if ref_circular else None
    best: tuple[int, int, int, list[tuple[int, int]]] | None = None   # votes, strand, key, hits
    for strand, oriented in ((1, query), (-1, rc(query))):
        hits = _hits(index, oriented, query_circular)
        diagonals = [_diagonal(q, r, ref_len, query_len, ref_circular, query_circular)
                     for q, r in hits]
        key, votes = _best_diagonal(diagonals, modulus)
        if best is None or votes > best[0]:
            best = (votes, strand, key, hits)
    assert best is not None
    votes, strand, key, hits = best
    if votes < MIN_VOTES:
        return None
    oriented = query if strand == 1 else rc(query)

    rotation = shift = 0
    if query_circular:
        rotation = key
        oriented = rotate(oriented, rotation)
        points = [((q - rotation) % query_len, r) for q, r in hits]
    elif ref_circular:
        on_diagonal = [q for q, r in hits
                       if _near(_diagonal(q, r, ref_len, query_len, True, False), key, modulus)]
        if query_len < ref_len:      # centre the fragment in the rotated reference
            shift = (key - (ref_len - query_len) // 2) % ref_len
        else:                        # the reference lies inside the fragment
            shift = (key + min(on_diagonal, default=0)) % ref_len
        points = [(q, (r - shift) % ref_len) for q, r in hits]
    else:
        points = hits

    chain = _chain(points)
    covered = _run_coverage(chain, on_query=query_len <= ref_len)
    mode = "semi-global" if covered >= SEMI_GLOBAL_MIN_CHAIN * min(query_len, ref_len) else "local"
    if mode == "semi-global" and query_circular and ref_circular:
        # Two circles have no ends: cut both inside the same shared word, so that both start
        # and end in homology, and score the ends like any other gap. With free end gaps a
        # difference next to the cut would pass for an overhang and go unreported.
        cut_query, cut_reference = chain[len(chain) // 2]
        oriented = rotate(oriented, cut_query)
        shift = cut_reference
        mode = "global"
    target = rotate(reference, shift) if shift else reference
    return Placement(strand, rotation, shift, mode, votes), target, oriented


# ---- step 2: the pairwise alignment ----------------------------------------------------

def _aligner(mode: str) -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner(scoring="blastn")
    if mode == "local":
        aligner.mode = "local"
    else:
        aligner.mode = "global"
        if mode == "semi-global":
            aligner.end_gap_score = 0.0  # free end gaps on both sequences
    return aligner


def _columns(target: str, query: str, coordinates) -> list[tuple[int | None, int | None]]:
    """The alignment as (target index or None, query index or None) columns."""
    out: list[tuple[int | None, int | None]] = []
    for i in range(coordinates.shape[1] - 1):
        t0, t1 = int(coordinates[0][i]), int(coordinates[0][i + 1])
        q0, q1 = int(coordinates[1][i]), int(coordinates[1][i + 1])
        if t1 > t0 and q1 > q0:
            out.extend(zip(range(t0, t1), range(q0, q1)))
        elif t1 > t0:
            out.extend((t, None) for t in range(t0, t1))
        else:
            out.extend((None, q) for q in range(q0, q1))
    # Overhangs (one sequence running past the other's end) are not compared.
    paired = [i for i, (t, q) in enumerate(out) if t is not None and q is not None]
    return out[paired[0]:paired[-1] + 1] if paired else []


def differences(columns: list[tuple[str, str, int | None]], ref_len: int) -> list[dict[str, Any]]:
    """Mismatches, deletions and insertions in 1-based reference coordinates.

    ``columns`` holds (reference base or "-", other base or "-", reference position or None).
    A run of deleted or inserted bases is one event. An insertion is placed between the two
    reference bases that flank it (``start`` and ``end``), the way HGVS writes ``123_124ins``.
    """
    out: list[dict[str, Any]] = []
    previous_position: int | None = None
    for ref_base, other_base, position in columns:
        last = out[-1] if out else None
        if ref_base == "-":
            assert previous_position is not None   # overhangs were trimmed
            if last and last["kind"] == "insertion" and last["start"] == previous_position + 1:
                last["sequence"] += other_base
            else:
                out.append({"kind": "insertion", "start": previous_position + 1,
                            "end": (previous_position + 1) % ref_len + 1,
                            "reference": "", "sequence": other_base})
            continue
        assert position is not None
        if other_base == "-":
            if (last and last["kind"] == "deletion"
                    and last["end"] % ref_len == position % ref_len):
                last["end"] = position + 1
                last["reference"] += ref_base
            else:
                out.append({"kind": "deletion", "start": position + 1, "end": position + 1,
                            "reference": ref_base, "sequence": ""})
        elif other_base != ref_base:
            out.append({"kind": "mismatch", "start": position + 1, "end": position + 1,
                        "reference": ref_base, "sequence": other_base})
        previous_position = position
    return out


def align_pair(reference: str, ref_circular: bool, name: str, sequence: str, circular: bool,
               index: dict[str, list[int]] | None = None) -> PairResult | None:
    """Align one sequence to the reference; None when the two have nothing in common."""
    placed = place(reference, ref_circular, sequence, circular, index)
    if placed is None:
        return None
    placement, target, query = placed
    alignment = _aligner(placement.mode).align(target, query)[0]
    pairs = _columns(target, query, alignment.coordinates)
    if not pairs:
        return None
    ref_len = len(reference)
    shift = placement.reference_shift
    columns: list[tuple[str, str, int | None]] = []
    for t, q in pairs:
        position = None if t is None else (t + shift) % ref_len
        columns.append((target[t] if t is not None else "-",
                        query[q] if q is not None else "-", position))
    first_t = next(t for t, _ in pairs if t is not None)
    last_t = next(t for t, _ in reversed(pairs) if t is not None)
    first_q = next(q for _, q in pairs if q is not None)
    last_q = next(q for _, q in reversed(pairs) if q is not None)
    reference_gapped = "".join(c[0] for c in columns)
    other_gapped = "".join(c[1] for c in columns)
    matches = sum(1 for a, b, _ in columns if a == b != "-")
    return PairResult(
        name=name, length=len(sequence), placement=placement,
        reference_gapped=reference_gapped, other_gapped=other_gapped,
        offset=(first_t + shift) % ref_len, reference_span=last_t - first_t + 1,
        other_aligned=query[first_q:last_q + 1],
        identity=100.0 * matches / len(columns),
        # A circular alignment starts at its internal cut; list in reference order.
        differences=sorted(differences(columns, ref_len), key=lambda d: d["start"]),
    )


# ---- step 3: everything the dialog needs -----------------------------------------------

def _counts(diffs: list[dict[str, Any]]) -> dict[str, int]:
    return {
        "mismatches": sum(1 for d in diffs if d["kind"] == "mismatch"),
        "insertions": sum(1 for d in diffs if d["kind"] == "insertion"),
        "deletions": sum(1 for d in diffs if d["kind"] == "deletion"),
        "inserted_bases": sum(len(d["sequence"]) for d in diffs if d["kind"] == "insertion"),
        "deleted_bases": sum(len(d["reference"]) for d in diffs if d["kind"] == "deletion"),
    }


def alignment_payload(reference_doc: dict[str, Any], docs: list[dict[str, Any]]) -> dict:
    reference_name = reference_doc.get("name") or "Reference"
    reference = clean_sequence(reference_doc.get("sequence") or "", reference_name)
    ref_circular = bool(reference_doc.get("circular"))
    others = [(d.get("name") or f"Sequence {i + 1}",
               clean_sequence(d.get("sequence") or "", d.get("name") or f"Sequence {i + 1}"))
              for i, d in enumerate(docs)]
    check_limits(reference_name, reference, others)

    index = _reference_index(reference, ref_circular)
    results: list[tuple[int, PairResult]] = []
    warnings: list[str] = []
    for i, ((name, sequence), doc) in enumerate(zip(others, docs)):
        if len(sequence) < MIN_LENGTH:
            warnings.append(f"“{name}” was left out: it is shorter than {MIN_LENGTH} bp.")
            continue
        result = align_pair(reference, ref_circular, name, sequence, bool(doc.get("circular")),
                            index)
        if result is None:
            warnings.append(f"“{name}” was left out: it shares no stretch of {WORD} identical "
                            f"bases with “{reference_name}” on either strand, so the two "
                            "cannot be aligned meaningfully.")
            continue
        results.append((i, result))

    ref_len = len(reference)
    reference_gapped, gapped = merge_gapped(
        reference, [(r.reference_gapped, r.other_gapped, r.offset) for _, r in results])
    tracks: list[dict[str, Any]] = []
    if results:
        tracks.append({
            "sequenceData": {"id": "reference", "name": reference_name, "sequence": reference,
                             "circular": False,
                             "features": _linear_features(reference_doc.get("features"), ref_len)},
            "alignmentData": {"id": "reference-alignment", "sequence": reference_gapped},
        })
        for (i, result), track in zip(results, gapped):
            suffix = " (reverse complement)" if result.placement.strand == -1 else ""
            tracks.append({
                "sequenceData": {"id": f"sequence-{i}", "name": result.name + suffix,
                                 "sequence": result.other_aligned},
                "alignmentData": {"id": f"sequence-alignment-{i}", "sequence": track},
            })

    sequences = []
    for i, result in results:
        whole = result.reference_span == ref_len
        start = 0 if whole else result.offset
        end = ref_len - 1 if whole else (result.offset + result.reference_span - 1) % ref_len
        placement = result.placement
        sequences.append({
            "index": i,
            "name": result.name,
            "length": result.length,
            "strand": placement.strand,
            # Circular sequence: which of its own bases now faces the reference's base 1.
            "origin_moved_to": (None if not placement.rotation else
                                placement.rotation + 1 if placement.strand == 1
                                else result.length - placement.rotation),
            "mode": placement.mode,
            "start": start + 1,
            "end": end + 1,
            "spans_origin": not whole and result.offset + result.reference_span > ref_len,
            "coverage": round(100.0 * result.reference_span / ref_len, 1),
            "aligned_length": len(result.other_aligned),
            "identity": round(result.identity, 2),
            "differences": result.differences,
            "counts": _counts(result.differences),
        })
    return {
        "alignmentTracks": tracks,
        "reference": {"name": reference_name, "length": ref_len, "circular": ref_circular},
        "sequences": sequences,
        "warnings": warnings,
        "first_covered": min((s["start"] for s in sequences), default=1),
    }


def _linear_features(features: Any, length: int) -> list[dict[str, Any]]:
    """Reference features for the alignment view, which draws the reference as linear: a
    feature across the origin is left out rather than drawn backwards."""
    items = list(features.values()) if isinstance(features, dict) else list(features or [])
    out = []
    for f in items:
        start, end = f.get("start"), f.get("end")
        if isinstance(start, int) and isinstance(end, int) and 0 <= start <= end < length:
            out.append({k: f[k] for k in ("id", "name", "type", "start", "end", "strand",
                                          "forward", "color") if k in f})
    return out


# ---- HTTP ------------------------------------------------------------------------------

class AlignRequest(BaseModel):
    reference: dict[str, Any]
    sequences: list[dict[str, Any]]


@router.get("/limits")
def limits() -> dict[str, int]:
    return {"max_sequences": MAX_SEQUENCES, "max_length": MAX_LENGTH, "min_length": MIN_LENGTH}


@router.post("/align")
def align(body: AlignRequest) -> dict[str, Any]:
    try:
        return alignment_payload(body.reference, body.sequences)
    except AlignmentError as exc:
        raise HTTPException(422, str(exc)) from exc
