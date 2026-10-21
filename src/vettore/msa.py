"""Progressive multiple alignment: every sequence aligned to every other. Routes under
/api/alignment/multiple.

The star alignment (``vettore.alignment``) answers "where do these differ from the
reference?". This module answers "how do these line up with each other?": variants of a
construct or capsid genes of several serotypes, side by side, with no sequence privileged.
It is a progressive alignment in the manner of ClustalW (Thompson, Higgins & Gibson 1994,
NAR 22:4673) and MUSCLE's first stage (Edgar 2004, NAR 32:1792):

1. **Strand and origin.** Each sequence is turned to the strand and, when circular, to the
   origin of the first sequence, with the shared-word vote of the star alignment
   (``alignment.place``). A sequence unrelated to the first is tried against the others
   already placed; one related to none of them is left out with the reason.
2. **Guide tree.** Distances are shared-word distances: the fractional common k-mer count of
   MUSCLE (Edgar 2004, BMC Bioinformatics 5:113), d = 1 - F with F the shared words of
   length ``DISTANCE_WORD`` over the words of the shorter sequence. UPGMA joins the closest
   pair first, so the most similar sequences are aligned together before distant ones.
3. **Progressive merges.** Following the tree, two groups (profiles) are aligned to each
   other column by column with ``Bio.Align.PairwiseAligner`` and the scores of EMBOSS
   needle (Rice, Longden & Bleasby 2000, Trends Genet 16:276: EDNAFULL reward 5, penalty -4,
   gap open 10, extend 0.5). Two columns score the
   average of the scores of every pair of bases across them, a base facing a gap scoring 0
   (the profile score of ClustalW), so a column most sequences skip is cheap to pass. Gaps
   already in a group stay ("once a gap, always a gap"). A stretch longer than
   ``PROFILE_MAX`` columns between two anchors, which would need a score table too large for
   a laptop, is aligned on summaries instead: each column as the IUPAC code of its bases
   (those under ``CONSENSUS_SHARE`` of the column ignored), codes that share a base scoring
   as a match.
4. **Anchors.** Words of ``ANCHOR_WORD`` bases found once in each summary and lined up on a
   diagonal for at least ``MIN_ANCHOR`` bases are kept as they are (the longest collinear
   chain of them), and only the stretches between anchors go through dynamic programming.
   Related sequences share long identical stretches, so a merge of two 15 kb groups costs
   a fraction of a second instead of the ~4 s and ~460 MB of a full matrix. End gaps are
   free on the outer stretches, so a fragment can sit anywhere along a full-length sequence.

The view gets the stacked alignment in tree order (similar sequences next to each other), a
consensus row, the conservation of every column, and the differences between any two
sequences in the first one's numbering (``compare``).

Limits: ``MAX_SEQUENCES`` sequences of up to ``MAX_LENGTH`` bp and ``MAX_TOTAL`` bp in all
(10 sequences of 15 kb, or 50 of 3 kb). The alignment runs as a background job whose
progress the dialog shows.
"""

from __future__ import annotations

import logging
import threading
import uuid
from collections import Counter
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from .alignment import (
    MIN_LENGTH,
    AlignmentError,
    _chain,
    _counts,
    _name_and_sequence,
    _reference_index,
    differences,
    place,
    rc,
)

log = logging.getLogger(__name__)
router = APIRouter(prefix="/api/alignment/multiple")

MAX_SEQUENCES = 50
MAX_LENGTH = 15_000          # bp, each sequence
MAX_TOTAL = 150_000          # bp, all sequences together
DISTANCE_WORD = 8            # word length of the guide-tree distances
ANCHOR_WORD = 12             # word length of the anchors between two groups
MIN_ANCHOR = 24              # bp of identical sequence on one diagonal to anchor a merge
CONSENSUS_SHARE = 0.25       # bases rarer than this in a column do not shape its summary
MIN_OVERLAP = 0.5            # of the shorter group's columns, else redone with end gaps scored
# Scores of EMBOSS needle's defaults (EDNAFULL matrix, gap open 10, extend 0.5). On the
# simulated families of tests/test_msa.py they recover more of the true alignment than the
# blastn scores of the star alignment (0.92 against 0.88 sum-of-pairs): a progressive
# alignment must open the gaps of every member, and blastn's costly gaps resist that.
MATCH, MISMATCH, PARTIAL = 5.0, -4.0, 1.0
GAP_OPEN, GAP_EXTEND = -10.0, -0.5
PROFILE_MAX = 3_000          # columns; longer stretches between anchors use the summaries

IUPAC = "ACGTRYSWKMBDHVN"
CODE_BASES = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
BASES_CODE = {frozenset(bases): code for code, bases in CODE_BASES.items()}


def _substitution_matrix() -> substitution_matrices.Array:
    """MATCH and MISMATCH extended to IUPAC codes: MATCH for the same base, PARTIAL for codes
    that share a base, MISMATCH for codes that share none; N scores 0 (it says nothing)."""
    matrix = substitution_matrices.Array(IUPAC, dims=2)
    for a in IUPAC:
        for b in IUPAC:
            if "N" in (a, b):
                score = 0.0
            elif set(CODE_BASES[a]) & set(CODE_BASES[b]):
                score = MATCH if a == b and len(CODE_BASES[a]) == 1 else PARTIAL
            else:
                score = MISMATCH
            matrix[a, b] = score
    return matrix


MATRIX = _substitution_matrix()
#: Scores between the four bases, for profile scores.
BASE_SCORES = np.full((4, 4), MISMATCH) + np.eye(4) * (MATCH - MISMATCH)
#: Each character as its share of A, C, G and T (an ambiguity code spreads one base over the
#: bases it stands for; a gap is none).
WEIGHTS = np.zeros((256, 4))
for _code, _bases in CODE_BASES.items():
    for _base in _bases:
        WEIGHTS[ord(_code), "ACGT".index(_base)] = 1.0 / len(_bases)


class Cancelled(Exception):
    """The user closed the dialog or started another alignment."""


@dataclass
class Member:
    """A sequence as it enters the alignment."""
    index: int               # position in the request
    name: str
    sequence: str            # as given
    circular: bool
    strand: int = 1
    rotation: int = 0        # oriented sequence = rotate(sequence or its reverse complement, rotation)
    oriented: str = ""


@dataclass
class Job:
    progress: float = 0.0
    message: str = "Starting…"
    result: dict[str, Any] | None = None
    error: Exception | None = None
    cancelled: bool = False
    done: threading.Event = field(default_factory=threading.Event)


# ---- input -----------------------------------------------------------------------------

def check_limits(members: list[tuple[str, str]]) -> None:
    if len(members) < 2:
        raise AlignmentError("Choose at least two sequences to align to each other.")
    if len(members) > MAX_SEQUENCES:
        raise AlignmentError(f"{len(members)} sequences were chosen, but at most "
                             f"{MAX_SEQUENCES} can be aligned to each other at once. Untick "
                             "some and align again.")
    for name, seq in members:
        if len(seq) > MAX_LENGTH:
            raise AlignmentError(f"“{name}” is {len(seq):,} bp; sequences longer than "
                                 f"{MAX_LENGTH:,} bp cannot be aligned here. Align the region "
                                 "you need instead (Actions → New File from Selection).")
    total = sum(len(seq) for _, seq in members)
    if total > MAX_TOTAL:
        raise AlignmentError(f"The chosen sequences add up to {total:,} bp; at most "
                             f"{MAX_TOTAL:,} bp can be aligned to each other at once. Untick "
                             "some, or align the regions you need.")


# ---- step 1: strand and origin ---------------------------------------------------------

def orient(members: list[Member], tick: Callable[[str], None]) -> list[str]:
    """Turn every sequence to the strand and origin of the first one (or of another already
    turned, when it shares nothing with the first). Returns warnings for those left out."""
    first = members[0]
    first.oriented = first.sequence
    placed = [first]
    indexes = {first.index: _reference_index(first.oriented, first.circular)}
    pending = members[1:]
    # A sequence unrelated to the first may be related to one placed later: go round again
    # until a round places nothing new.
    progress = True
    while progress:
        progress = False
        for member in pending:
            if member.oriented:
                continue
            tick(f"Finding the strand of sequence {member.index + 1} of {len(members)}")
            for anchor in placed:
                index = indexes.setdefault(anchor.index,
                                           _reference_index(anchor.oriented, anchor.circular))
                found = place(anchor.oriented, anchor.circular, member.sequence,
                              member.circular, index, local=True)
                if found is not None:
                    placement, _, oriented = found
                    member.strand, member.rotation, member.oriented = (
                        placement.strand, placement.rotation, oriented)
                    placed.append(member)
                    progress = True
                    break
    warnings = []
    for member in pending:
        if not member.oriented:
            warnings.append(f"“{member.name}” was left out: it has nothing recognisably in "
                            "common with the other sequences on either strand, so it cannot "
                            "be aligned to them meaningfully. Check that it is the sequence "
                            "you meant.")
    members[:] = [m for m in members if m.oriented]
    return warnings


# ---- step 2: guide tree ----------------------------------------------------------------

def _word_counts(seq: str) -> Counter[str]:
    return Counter(seq[i:i + DISTANCE_WORD] for i in range(len(seq) - DISTANCE_WORD + 1))


def distance(a: Counter[str], b: Counter[str], length_a: int, length_b: int) -> float:
    """1 - the fraction of shared words (Edgar 2004): 0 for identical sequences, 1 for
    sequences that share no word."""
    words = max(1, min(length_a, length_b) - DISTANCE_WORD + 1)
    small, large = (a, b) if len(a) <= len(b) else (b, a)
    shared = sum(min(count, large[w]) for w, count in small.items() if w in large)
    return 1.0 - min(1.0, shared / words)


def guide_tree(seqs: list[str]) -> list[tuple[int, int]]:
    """UPGMA joins as (cluster, cluster) pairs. Leaves are 0..n-1; the k-th join creates
    cluster n + k. Ties are broken by the lower indexes, so the tree is deterministic."""
    n = len(seqs)
    counts = [_word_counts(s) for s in seqs]
    dist: dict[tuple[int, int], float] = {}
    for i in range(n):
        for j in range(i + 1, n):
            dist[i, j] = distance(counts[i], counts[j], len(seqs[i]), len(seqs[j]))
    size = {i: 1 for i in range(n)}
    joins = []
    for k in range(n - 1):
        (i, j), _ = min(dist.items(), key=lambda item: (item[1], item[0]))
        new = n + k
        joins.append((i, j))
        for other in list(size):
            if other in (i, j):
                continue
            d_i = dist.pop((min(i, other), max(i, other)))
            d_j = dist.pop((min(j, other), max(j, other)))
            dist[other, new] = (d_i * size[i] + d_j * size[j]) / (size[i] + size[j])
        del dist[i, j]
        size[new] = size.pop(i) + size.pop(j)
    return joins


# ---- step 3: aligning two groups -------------------------------------------------------

def column_codes(rows: list[str]) -> str:
    """One IUPAC code per column: the bases that make up at least CONSENSUS_SHARE of the
    column's bases. A column of ambiguous codes only keeps the first; N when nothing else."""
    if len(rows) == 1:
        return rows[0].replace("-", "N")
    out = []
    for column in zip(*rows):
        counts = [column.count(b) for b in "ACGT"]
        total = sum(counts)
        if not total:
            out.append(next((c for c in column if c != "-"), "N"))
            continue
        floor = max(1, CONSENSUS_SHARE * total)
        present = frozenset(b for b, c in zip("ACGT", counts) if c >= floor)
        out.append(BASES_CODE[present or frozenset("ACGT")])
    return "".join(out)


def _aligner(left_free: bool, right_free: bool) -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = MATRIX
    aligner.open_gap_score, aligner.extend_gap_score = GAP_OPEN, GAP_EXTEND
    if left_free:
        aligner.left_gap_score = 0.0
    if right_free:
        aligner.right_gap_score = 0.0
    return aligner


def _columns_of(coordinates) -> list[tuple[int | None, int | None]]:
    out: list[tuple[int | None, int | None]] = []
    for k in range(coordinates.shape[1] - 1):
        x0, x1 = int(coordinates[0][k]), int(coordinates[0][k + 1])
        y0, y1 = int(coordinates[1][k]), int(coordinates[1][k + 1])
        if x1 > x0 and y1 > y0:
            out.extend(zip(range(x0, x1), range(y0, y1)))
        elif x1 > x0:
            out.extend((i, None) for i in range(x0, x1))
        else:
            out.extend((None, j) for j in range(y0, y1))
    return out


def _segment(x: str, y: str, fx: np.ndarray, fy: np.ndarray, left_free: bool,
             right_free: bool) -> list[tuple[int | None, int | None]]:
    """Global alignment of two stretches as (x index or None, y index or None) columns.
    ``fx`` and ``fy`` hold each column's share of A, C, G and T (see ``frequencies``)."""
    if not x or not y:
        return [(i, None) for i in range(len(x))] + [(None, j) for j in range(len(y))]
    aligner = _aligner(left_free, right_free)
    side = max(len(x), len(y))
    if side > PROFILE_MAX:
        return _columns_of(aligner.align(x, y)[0].coordinates)
    # Profile scores: the aligner looks up score[i, j] for column i of x facing column j of y.
    scores = np.zeros((side, side))
    scores[:len(x), :len(y)] = fx @ BASE_SCORES @ fy.T
    aligner.substitution_matrix = scores
    return _columns_of(aligner.align(np.arange(len(x), dtype=np.int32),
                                     np.arange(len(y), dtype=np.int32))[0].coordinates)


def anchors(x: str, y: str) -> list[tuple[int, int, int]]:
    """Identical stretches shared by x and y as (x start, y start, length): words found once
    in each, joined along their diagonal, at least MIN_ANCHOR long, in one collinear chain
    that does not overlap itself."""
    def unique_words(s: str) -> dict[str, int]:
        seen: dict[str, int] = {}
        repeated = set()
        for i in range(len(s) - ANCHOR_WORD + 1):
            word = s[i:i + ANCHOR_WORD]
            if word in seen:
                repeated.add(word)
            else:
                seen[word] = i
        return {w: i for w, i in seen.items() if w not in repeated and set(w) <= set("ACGT")}

    in_y = unique_words(y)
    hits = sorted((j - i, i) for w, i in unique_words(x).items() if (j := in_y.get(w)) is not None)
    blocks = []
    for diagonal, i in hits:
        if blocks and blocks[-1][1] - blocks[-1][0] == diagonal and blocks[-1][3] == i - 1:
            blocks[-1][3] = i
        else:
            blocks.append([i, i + diagonal, i, i])   # x start, y start, first word, last word
    blocks = [(x0, y0, last - first + ANCHOR_WORD) for x0, y0, first, last in blocks
              if last - first + ANCHOR_WORD >= MIN_ANCHOR]
    length = {(x0, y0): n for x0, y0, n in blocks}
    chained = []
    end_x = end_y = 0
    for x0, y0 in _chain([(x0, y0) for x0, y0, _ in blocks]):
        n = length[x0, y0]
        overlap = max(end_x - x0, end_y - y0, 0)   # trim the start so blocks do not overlap
        x0, y0, n = x0 + overlap, y0 + overlap, n - overlap
        if n >= ANCHOR_WORD:
            chained.append((x0, y0, n))
            end_x, end_y = x0 + n, y0 + n
    return chained


def frequencies(rows: list[str]) -> np.ndarray:
    """Each column's share of A, C, G and T over all rows (gaps count as nothing)."""
    table = np.frombuffer("".join(rows).encode("ascii"), dtype=np.uint8).reshape(len(rows), -1)
    return WEIGHTS[table].sum(axis=0) / len(rows)


def align_columns(x: str, y: str, fx: np.ndarray, fy: np.ndarray, free_ends: bool = True
                  ) -> list[tuple[int | None, int | None]]:
    """Columns of x paired with columns of y: anchors kept, stretches between them aligned,
    end gaps free on the outer stretches unless ``free_ends`` is false."""
    out: list[tuple[int | None, int | None]] = []
    i = j = 0

    def stretch(x1: int, y1: int, left_free: bool, right_free: bool) -> None:
        out.extend((a + i if a is not None else None, b + j if b is not None else None)
                   for a, b in _segment(x[i:x1], y[j:y1], fx[i:x1], fy[j:y1],
                                        left_free, right_free))

    for k, (x0, y0, n) in enumerate(anchors(x, y)):
        stretch(x0, y0, free_ends and k == 0, False)
        out.extend(zip(range(x0, x0 + n), range(y0, y0 + n)))
        i, j = x0 + n, y0 + n
    stretch(len(x), len(y), free_ends and i == 0 and j == 0, free_ends)
    return out


def merge(a: list[str], b: list[str]) -> list[str]:
    """Two groups of gapped rows aligned to each other, as one group."""
    codes_a, codes_b, fa, fb = column_codes(a), column_codes(b), frequencies(a), frequencies(b)
    pairs = align_columns(codes_a, codes_b, fa, fb)
    paired = sum(1 for i, j in pairs if i is not None and j is not None)
    if paired < MIN_OVERLAP * min(len(codes_a), len(codes_b)):
        # Every sequence here is related to another (step 1), so two groups that slid apart
        # along their free end gaps did so because distant groups score little per column.
        pairs = align_columns(codes_a, codes_b, fa, fb, free_ends=False)
    rows = ["".join(row[i] if i is not None else "-" for i, _ in pairs) for row in a]
    rows += ["".join(row[j] if j is not None else "-" for _, j in pairs) for row in b]
    return rows


# ---- step 4: what the view shows -------------------------------------------------------

def _span(row: str) -> tuple[int, int]:
    """First and last column holding a base of this row."""
    first = next(i for i, c in enumerate(row) if c != "-")
    last = len(row) - 1 - next(i for i, c in enumerate(reversed(row)) if c != "-")
    return first, last


def consensus_and_conservation(rows: list[str]) -> tuple[str, list[int]]:
    """The consensus row (the most common base of each column, an IUPAC code for a tie, a gap
    where most sequences have one) and the conservation of each column: the percentage of
    the sequences covering it that hold its most common base. A sequence does not cover the
    columns before its first base or after its last (a fragment's ends)."""
    spans = [_span(r) for r in rows]
    consensus, conservation = [], []
    for k, column in enumerate(zip(*rows)):
        covering = [c for c, (first, last) in zip(column, spans) if first <= k <= last]
        counts = Counter(c for c in covering if c != "-")
        gaps = len(covering) - sum(counts.values())
        if not counts:
            consensus.append("-")
            conservation.append(0)
            continue
        top = max(counts.values())
        if gaps > top:
            consensus.append("-")
        else:
            tied = [c for c, n in counts.items() if n == top]
            if len(tied) == 1:
                consensus.append(tied[0])
            else:
                bases = frozenset("".join(CODE_BASES.get(c, "ACGT") for c in tied))
                consensus.append(BASES_CODE.get(bases, "N"))
        conservation.append(round(100 * max(top, gaps) / len(covering)))
    return "".join(consensus), conservation


def _identity(row: str, consensus: str) -> float:
    first, last = _span(row)
    columns = [(a, b) for a, b in zip(row[first:last + 1], consensus[first:last + 1])
               if not (a == "-" and b == "-")]
    same = sum(1 for a, b in columns if a == b)
    return round(100.0 * same / max(1, len(columns)), 2)


def leaf_order(joins: list[tuple[int, int]], n: int) -> list[int]:
    children = {n + k: pair for k, pair in enumerate(joins)}
    order: list[int] = []
    stack = [n + len(joins) - 1] if joins else [0]
    while stack:
        node = stack.pop()
        if node < n:
            order.append(node)
        else:
            left, right = children[node]
            stack.extend((right, left))
    return order


def multiple_alignment(docs: list[Any], report: Callable[[float, str], None] | None = None,
                       cancelled: Callable[[], bool] = lambda: False) -> dict[str, Any]:
    """Everything the dialog needs. ``report`` hears (fraction done, what is happening);
    ``cancelled`` is asked between steps."""
    named = [_name_and_sequence(d, f"Sequence {i + 1}") for i, d in enumerate(docs)]
    docs = [d if isinstance(d, dict) else {} for d in docs]
    check_limits(named)
    warnings = []
    members = []
    for i, ((name, sequence), doc) in enumerate(zip(named, docs)):
        if len(sequence) < MIN_LENGTH:
            warnings.append(f"“{name}” was left out: it is shorter than {MIN_LENGTH} bp.")
        else:
            members.append(Member(i, name, sequence, bool(doc.get("circular"))))
    if len(members) < 2:
        raise AlignmentError(f"At least two sequences of {MIN_LENGTH} bp or more are needed.")

    # Steps, for the progress bar: orienting is cheap, merges dominate (weighted by length).
    steps = {"done": 0.0}
    n = len(members)
    total = (n - 1) * 0.2 + 1 + (n - 1) * 3

    def tick(message: str, weight: float = 0.2) -> None:
        if cancelled():
            raise Cancelled()
        if report:
            report(min(0.99, steps["done"] / total), message)
        steps["done"] += weight

    warnings += orient(members, tick)
    if len(members) < 2:
        raise AlignmentError("None of the chosen sequences has anything recognisably in "
                             "common with the others, so there is nothing to align. Check "
                             "that they are the sequences you meant.")
    n = len(members)
    tick(f"Comparing the {n} sequences with each other", 1)
    joins = guide_tree([m.oriented for m in members])
    groups: dict[int, list[int]] = {i: [i] for i in range(n)}
    rows: dict[int, list[str]] = {i: [m.oriented] for i, m in enumerate(members)}
    for k, (a, b) in enumerate(joins):
        tick(f"Aligning: step {k + 1} of {len(joins)}", 3)
        rows[n + k] = merge(rows.pop(a), rows.pop(b))
        groups[n + k] = groups.pop(a) + groups.pop(b)
    (root,) = groups
    gapped = dict(zip(groups[root], rows[root]))
    if report:
        report(0.99, "Preparing the view")

    order = leaf_order(joins, n)
    consensus, conservation = consensus_and_conservation([gapped[i] for i in order])
    tracks = [{
        "sequenceData": {"id": "consensus", "name": "Consensus",
                         "sequence": consensus.replace("-", "")},
        "alignmentData": {"id": "consensus-alignment", "sequence": consensus},
    }]
    sequences = []
    for i in order:
        m = members[i]
        suffix = " (reverse complement)" if m.strand == -1 else ""
        tracks.append({
            "sequenceData": {"id": f"sequence-{m.index}", "name": m.name + suffix,
                             "sequence": m.oriented},
            "alignmentData": {"id": f"sequence-alignment-{m.index}", "sequence": gapped[i]},
        })
        first, last = _span(gapped[i])
        sequences.append({
            "index": m.index,
            "name": m.name,
            "length": len(m.sequence),
            "circular": m.circular,
            "strand": m.strand,
            "rotation": m.rotation,
            # Circular sequence: which of its own bases now comes first.
            "origin_moved_to": (None if not m.rotation else
                                m.rotation + 1 if m.strand == 1
                                else len(m.sequence) - m.rotation),
            "row": gapped[i],
            "first_column": first + 1,
            "last_column": last + 1,
            "identity": _identity(gapped[i], consensus),
        })
    return {
        "alignmentTracks": tracks,
        "sequences": sequences,
        "consensus": consensus,
        "conservation": conservation,
        "columns": len(consensus),
        "warnings": warnings,
    }


# ---- differences between two rows ------------------------------------------------------

def compare(a: dict[str, Any], b: dict[str, Any]) -> dict[str, Any]:
    """Differences of row b from row a, in a's own numbering (1-based, on a's strand as
    stored in its file), over the columns both cover."""
    row_a, row_b = str(a["row"]), str(b["row"])
    if len(row_a) != len(row_b):
        raise AlignmentError("The two rows come from different alignments. Align again.")
    length = int(a["length"])
    rotation = int(a.get("rotation") or 0)
    if int(a.get("strand") or 1) == -1:
        # Turn both rows over so that a reads on its own strand; see rc(rotate(s, r)).
        row_a, row_b = rc(row_a), rc(row_b)
        rotation = (length - rotation) % length
    pairs = [(x, y) for x, y in zip(row_a, row_b) if not (x == "-" and y == "-")]
    both = [k for k, (x, y) in enumerate(pairs) if x != "-" and y != "-"]
    if not both:
        return {"differences": [], "counts": _counts([]), "identity": 0.0, "start": None,
                "end": None, "overlap": 0}
    before = sum(1 for x, _ in pairs[:both[0]] if x != "-")
    columns: list[tuple[str, str, int | None]] = []
    position = before
    for x, y in pairs[both[0]:both[-1] + 1]:
        if x == "-":
            columns.append((x, y, None))
        else:
            columns.append((x, y, (position + rotation) % length))
            position += 1
    diffs = sorted(differences(columns, length), key=lambda d: d["start"])
    same = sum(1 for x, y, _ in columns if x == y)
    first = next(p for _, _, p in columns if p is not None)
    last = next(p for _, _, p in reversed(columns) if p is not None)
    return {
        "differences": diffs,
        "counts": _counts(diffs),
        "identity": round(100.0 * same / len(columns), 2),
        "start": first + 1,
        "end": last + 1,
        "overlap": position - before,
    }


# ---- background jobs and HTTP ----------------------------------------------------------

JOBS: dict[str, Job] = {}
JOBS_LOCK = threading.Lock()
MAX_JOBS = 8                 # finished jobs nobody collected are dropped beyond this


def start(docs: list[Any]) -> str:
    job = Job()
    job_id = uuid.uuid4().hex
    with JOBS_LOCK:
        finished = [k for k, j in JOBS.items() if j.done.is_set()]
        for k in finished[:max(0, len(JOBS) - MAX_JOBS + 1)]:
            del JOBS[k]
        JOBS[job_id] = job

    def report(fraction: float, message: str) -> None:
        job.progress, job.message = fraction, message

    def run() -> None:
        try:
            job.result = multiple_alignment(docs, report, lambda: job.cancelled)
            job.progress, job.message = 1.0, "Done"
        except Exception as exc:  # noqa: BLE001  (handed to the request that collects the job)
            job.error = exc
        finally:
            job.done.set()

    threading.Thread(target=run, name=f"msa-{job_id[:8]}", daemon=True).start()
    return job_id


class MultipleRequest(BaseModel):
    sequences: list[dict[str, Any]]


class CompareRequest(BaseModel):
    a: dict[str, Any]
    b: dict[str, Any]


@router.get("/limits")
def limits() -> dict[str, int]:
    return {"max_sequences": MAX_SEQUENCES, "max_length": MAX_LENGTH, "max_total": MAX_TOTAL,
            "min_length": MIN_LENGTH}


@router.post("/start")
def start_job(body: MultipleRequest) -> dict[str, str]:
    # Refuse what is over the limits at once, rather than from the job a moment later.
    try:
        check_limits([_name_and_sequence(d, f"Sequence {i + 1}")
                      for i, d in enumerate(body.sequences)])
    except AlignmentError as exc:
        raise HTTPException(422, str(exc)) from exc
    return {"job": start(body.sequences)}


@router.get("/job/{job_id}")
def job_state(job_id: str) -> dict[str, Any]:
    with JOBS_LOCK:
        job = JOBS.get(job_id)
    if job is None:
        raise HTTPException(404, "This alignment is no longer running. Align again.")
    if not job.done.is_set():
        return {"state": "running", "progress": round(job.progress, 3), "message": job.message}
    with JOBS_LOCK:
        JOBS.pop(job_id, None)
    if isinstance(job.error, Cancelled):
        return {"state": "cancelled"}
    if isinstance(job.error, AlignmentError):
        raise HTTPException(422, str(job.error))
    if job.error is not None:
        raise job.error                   # the server's handler logs it and explains
    return {"state": "done", "result": job.result}


@router.post("/job/{job_id}/cancel")
def cancel(job_id: str) -> dict[str, bool]:
    with JOBS_LOCK:
        job = JOBS.pop(job_id, None)
    if job is not None:
        job.cancelled = True
    return {"cancelled": job is not None}


@router.post("/compare")
def compare_rows(body: CompareRequest) -> dict[str, Any]:
    try:
        return compare(body.a, body.b)
    except AlignmentError as exc:
        raise HTTPException(422, str(exc)) from exc
