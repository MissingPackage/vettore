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
    DIAGONAL_SLACK,
    MIN_LENGTH,
    MIN_VOTES,
    WORD,
    AlignmentError,
    _counts,
    _hits,
    _name_and_sequence,
    _reference_index,
    differences,
    place,
    rc,
    rotate,
)

log = logging.getLogger(__name__)
router = APIRouter(prefix="/api/alignment/multiple")

MAX_SEQUENCES = 50
MAX_LENGTH = 15_000          # bp, each sequence
MAX_TOTAL = 150_000          # bp, all sequences together
DISTANCE_WORD = 8            # word length of the guide-tree distances
ANCHOR_WORD = 12             # word length of the anchors between two groups
MIN_ANCHOR = 24              # bp of identical sequence on one diagonal to anchor a merge
REFINE_WINDOW = 500          # bp at the start of a circular anchor that decide where a circle starts
JOIN_WINDOW = 300            # columns each side of two circles' cut aligned again by compare
WRAP_IDENTITY = 0.8          # a fragment's overhang this identical to the circle's other end crosses it
SPAN_MARGIN = 50             # bp kept clear on each side of a fragment when cutting circles
CONSENSUS_SHARE = 0.25       # bases rarer than this in a column do not shape its summary
MIN_OVERLAP = 0.5            # of the shorter group's columns, else redone with end gaps scored
# Scores of EMBOSS needle's defaults (EDNAFULL matrix, gap open 10, extend 0.5). On the
# simulated families of tests/test_msa.py they recover more of the true alignment than the
# blastn scores of the star alignment (0.92 against 0.88 sum-of-pairs): a progressive
# alignment must open the gaps of every member, and blastn's costly gaps resist that.
MATCH, MISMATCH, PARTIAL = 5.0, -4.0, 1.0
GAP_OPEN, GAP_EXTEND = -10.0, -0.5
#: A segment of a pairwise comparison must score this much (blastn scores, see
#: ``_column_scores``) to count as homology: about 30 identical bases.
SEGMENT_MIN = 60.0
#: A homologous stretch ends where its score falls this far below its best: unrelated inserts
#: of 100 bp cost 104-158 in tests, 200 bp 199-300; the most divergent stretch of AAV2 against
#: AAV9 VP1 costs 90. Shorter unrelated stretches are reported as differences.
X_DROP = 150.0
#: A stretch with bases on one side only, give or take this many (or this share), is an
#: insertion or deletion, reported as such however long.
INDEL_SLACK, INDEL_SHARE = 5, 0.1
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
    if _code == "N":          # an unknown base scores 0 against anything, as in MATRIX
        continue
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

def _refine_circular(oriented: str, index: dict[str, list[int]]) -> int:
    """Extra turn for a circle already turned to a circular anchor by the majority vote.

    The majority diagonal is that of the longer arc on either side of the largest insertion
    or deletion between the two (a different insert in the same backbone), so the shorter
    arc, up to half the plasmid, would be out of frame and hang over the ends. The circle
    must start at the base that faces the anchor's base 1: the diagonal most voted for by
    the words shared at the start of the anchor (the first REFINE_WINDOW bases that hold
    at least MIN_VOTES words on one diagonal), however far it is from the majority one."""
    length = len(oriented)
    hits = sorted((r, (q - r + length // 2) % length - length // 2)
                  for q, r in _hits(index, oriented, True))
    start = 0
    while start < (hits[-1][0] + 1 if hits else 0):
        window = Counter(shift for r, shift in hits if start <= r < start + REFINE_WINDOW)
        best, votes = 0, 0
        for shift in sorted(window):
            near = sum(window.get(shift + d, 0) for d in range(-DIAGONAL_SLACK, DIAGONAL_SLACK + 1))
            if near > votes:
                best, votes = shift, near
        if votes >= MIN_VOTES:
            return best
        start += REFINE_WINDOW // 2
    return 0


def _orient_from(root: Member, members: list[Member], tick: Callable[[str], None]) -> None:
    """Place every member on the frame of ``root`` (already oriented). Each round takes the
    pairing of a pending member with a placed one that shares the most words, so a
    sequence is turned by its closest relative, not by whichever was placed first: two
    constructs sharing only an inverted element (a FLEx cassette, a loxP, a promoter) do not
    flip each other while a sequence matching end to end is there."""
    for member in members:
        if member is not root:
            member.oriented = ""
    placed = [root]
    pending = [m for m in members if m is not root]
    indexes: dict[int, dict[str, list[int]]] = {}
    tried: dict[tuple[int, int], Any] = {}
    while pending:
        best = None
        for member in pending:
            for anchor in placed:
                key = (member.index, anchor.index)
                if key not in tried:
                    tick(f"Finding the strand of “{member.name}”", 0)
                    if anchor.index not in indexes:
                        indexes[anchor.index] = _reference_index(anchor.oriented,
                                                                 anchor.circular)
                    index = indexes[anchor.index]
                    tried[key] = place(anchor.oriented, anchor.circular, member.sequence,
                                       member.circular, index, local=True)
                found = tried[key]
                if found is not None and (best is None or found[0].votes > best[0]):
                    best = (found[0].votes, member, anchor, found)
        if best is None:
            break
        _, member, anchor, (placement, _, oriented) = best
        member.strand, member.rotation, member.oriented = (placement.strand,
                                                           placement.rotation, oriented)
        if member.circular and anchor.circular:
            turn = _refine_circular(oriented, indexes[anchor.index])
            if turn:
                member.rotation = (member.rotation + turn) % len(oriented)
                member.oriented = rotate(oriented, turn)
        placed.append(member)
        pending.remove(member)
        tick(f"Placed “{member.name}”")


def _uncovered_middle(spans: list[tuple[int, int]], circle: int) -> int | None:
    """The middle of the longest stretch of the circle that no span covers, or None when
    the spans cover all of it."""
    covered = [False] * circle
    for start, length in spans:
        for k in range(start - SPAN_MARGIN, start + length + SPAN_MARGIN):
            covered[k % circle] = True
    if all(covered):
        return None
    first = covered.index(False)
    best = (0, 0)       # length, start
    run_start, run = first, 0
    for step in range(circle):
        k = (first + step) % circle
        if covered[k]:
            if run > best[0]:
                best = (run, run_start)
            run = 0
        else:
            if run == 0:
                run_start = k
            run += 1
    if run > best[0]:
        best = (run, run_start)
    return (best[1] + best[0] // 2) % circle


def orient(members: list[Member], tick: Callable[[str], None], root_turn: int = 0
           ) -> tuple[list[str], Member]:
    """Turn every sequence to the strand and origin of the first circular one ticked (the
    first one, when none is circular), that root itself turned by ``root_turn``. Returns
    warnings for the sequences left out, and the root."""
    root = next((m for m in members if m.circular), members[0])
    root.strand, root.rotation = 1, root_turn
    root.oriented = rotate(root.sequence, root_turn)
    _orient_from(root, members, tick)
    warnings = []
    for member in members:
        if not member.oriented:
            warnings.append(f"“{member.name}” was left out: it has nothing recognisably in "
                            "common with the other sequences on either strand (not even two "
                            f"stretches of {WORD} identical bases in line), so it cannot be "
                            "aligned to them meaningfully. Check that it is the sequence you "
                            "meant; sequences less than about 70 % identical cannot be lined "
                            "up here.")
    members[:] = [m for m in members if m.oriented]
    return warnings, root


def _wrap_share(overhang: str, root_end: str) -> float:
    """Identity of a fragment's overhang with the root's other end (ungapped, same length)."""
    return sum(a == b for a, b in zip(overhang, root_end)) / max(1, len(overhang))


def crossing_cut(root: Member, members: list[Member], gapped: dict[int, str]
                 ) -> tuple[int | None, list[str]]:
    """A circle has no ends, but the rows of an alignment do. When a linear sequence (an
    amplicon, a read) crosses the root circle's base 1, part of it hangs past one end of the
    circle's row while it continues the circle's other end. Returns where to open the root
    instead (a turn of its current row: the middle of the longest stretch no linear sequence
    covers) and the names of the crossing sequences; the turn is None when nothing crosses,
    or when the linear sequences cover the whole circle between them."""
    ri = members.index(root)
    row = gapped[ri]
    first, last = _span(row)
    length = len(root.oriented)
    position, k = [], 0
    for c in row:
        position.append(k if c != "-" else None)
        k += c != "-"
    spans, crossing = [], []
    for i, member in enumerate(members):
        if member.circular:
            continue
        other = gapped[i]
        paired = [c for c in range(len(row)) if row[c] != "-" and other[c] != "-"]
        if not paired:
            continue
        after = other[last + 1:].replace("-", "")
        before = other[:first].replace("-", "")
        start = position[paired[0]]
        if after and _wrap_share(after, root.oriented[:len(after)]) >= WRAP_IDENTITY:
            crossing.append(member.name)
        elif before and _wrap_share(before, root.oriented[-len(before):]) >= WRAP_IDENTITY:
            crossing.append(member.name)
            start = (start - len(before)) % length
        spans.append((start, len(member.sequence)))
    if not crossing:
        return None, []
    return _uncovered_middle(spans, length), crossing


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
    return _heaviest_chain(blocks)


def _heaviest_chain(blocks: list[tuple[int, int, int]]) -> list[tuple[int, int, int]]:
    """The collinear chain of blocks that anchors the most bases. A block may overlap the
    one before it on either sequence (a repeat around an indel); the overlap is trimmed from
    its start and not counted. Weighing bases, not blocks, keeps a 3 kb identical stretch
    rather than four short blocks of a rearranged cassette."""
    if not blocks:
        return []
    blocks = sorted(blocks)
    xs = np.array([b[0] for b in blocks])
    ys = np.array([b[1] for b in blocks])
    ns = np.array([b[2] for b in blocks])
    score = np.zeros(len(blocks))
    parent = np.full(len(blocks), -1)
    for i in range(len(blocks)):
        score[i] = ns[i]
        if i:
            overlap = np.maximum(np.maximum(xs[:i] + ns[:i] - xs[i], ys[:i] + ns[:i] - ys[i]), 0)
            usable = (xs[:i] < xs[i]) & (ys[:i] < ys[i]) & (overlap <= ns[i] - ANCHOR_WORD)
            if usable.any():
                gain = np.where(usable, score[:i] + ns[i] - overlap, -1.0)
                j = int(np.argmax(gain))
                if gain[j] > score[i]:
                    score[i], parent[i] = gain[j], j
    chain = []
    i = int(np.argmax(score))
    while i != -1:
        chain.append(i)
        i = int(parent[i])
    out = []
    end_x = end_y = 0
    for i in reversed(chain):
        x0, y0, n = int(xs[i]), int(ys[i]), int(ns[i])
        overlap = max(end_x - x0, end_y - y0, 0)
        x0, y0, n = x0 + overlap, y0 + overlap, n - overlap
        out.append((x0, y0, n))
        end_x, end_y = x0 + n, y0 + n
    return out


def frequencies(rows: list[str]) -> np.ndarray:
    """Each column's share of A, C, G and T over all rows (gaps count as nothing)."""
    table = np.frombuffer("".join(rows).encode("ascii"), dtype=np.uint8).reshape(len(rows), -1)
    return WEIGHTS[table].sum(axis=0) / len(rows)


def align_columns(x: str, y: str, fx: np.ndarray, fy: np.ndarray, free_ends: bool = True,
                  check: Callable[[], None] = lambda: None
                  ) -> list[tuple[int | None, int | None]]:
    """Columns of x paired with columns of y: anchors kept, stretches between them aligned,
    end gaps free on the outer stretches unless ``free_ends`` is false."""
    out: list[tuple[int | None, int | None]] = []
    i = j = 0

    def stretch(x1: int, y1: int, left_free: bool, right_free: bool) -> None:
        check()
        out.extend((a + i if a is not None else None, b + j if b is not None else None)
                   for a, b in _segment(x[i:x1], y[j:y1], fx[i:x1], fy[j:y1],
                                        left_free, right_free))

    for k, (x0, y0, n) in enumerate(anchors(x, y)):
        stretch(x0, y0, free_ends and k == 0, False)
        out.extend(zip(range(x0, x0 + n), range(y0, y0 + n)))
        i, j = x0 + n, y0 + n
    stretch(len(x), len(y), free_ends and i == 0 and j == 0, free_ends)
    return out


def merge(a: list[str], b: list[str], check: Callable[[], None] = lambda: None) -> list[str]:
    """Two groups of gapped rows aligned to each other, as one group. ``check`` is called
    between stretches (it raises to stop)."""
    codes_a, codes_b, fa, fb = column_codes(a), column_codes(b), frequencies(a), frequencies(b)
    pairs = align_columns(codes_a, codes_b, fa, fb, check=check)
    paired = sum(1 for i, j in pairs if i is not None and j is not None)
    if paired < MIN_OVERLAP * min(len(codes_a), len(codes_b)):
        # Every sequence here is related to another (step 1), so two groups that slid apart
        # along their free end gaps did so because distant groups score little per column.
        pairs = align_columns(codes_a, codes_b, fa, fb, free_ends=False, check=check)
    rows = ["".join(row[i] if i is not None else "-" for i, _ in pairs) for row in a]
    rows += ["".join(row[j] if j is not None else "-" for _, j in pairs) for row in b]
    return rows


# ---- step 4: what the view shows -------------------------------------------------------

def _span(row: str) -> tuple[int, int]:
    """First and last column holding a base of this row."""
    first = next(i for i, c in enumerate(row) if c != "-")
    last = len(row) - 1 - next(i for i, c in enumerate(reversed(row)) if c != "-")
    return first, last


def consensus_and_conservation(rows: list[str]) -> tuple[str, list[int | None]]:
    """The consensus row (the most common base of each column, an IUPAC code for a tie, a gap
    where most sequences have one) and the conservation of each column: the percentage of
    the sequences covering it that hold its most common base. A gap or an N is never
    agreement, so an insertion carried by one sequence in ten is 10 % conserved, not 90 %.
    A sequence does not cover the columns before its first base or after its last (a
    fragment's ends); a column covered by fewer than two sequences compares nothing and
    has no conservation (None)."""
    spans = [_span(r) for r in rows]
    consensus, conservation = [], []
    for k, column in enumerate(zip(*rows)):
        covering = [c for c, (first, last) in zip(column, spans) if first <= k <= last]
        counts = Counter(c for c in covering if c != "-")
        gaps = len(covering) - sum(counts.values())
        if not counts:
            consensus.append("-")
            conservation.append(0 if len(covering) > 1 else None)
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
        known = {c: n for c, n in counts.items() if c != "N"}
        # Nothing but N: unknown bases say nothing about conservation.
        conservation.append(round(100 * max(known.values()) / len(covering))
                            if len(covering) > 1 and known else
                            0 if len(covering) > 1 and not counts else None)
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
    steps = {"done": 0.0, "total": (len(members) - 1) * 3.2 + 1}

    def tick(message: str, weight: float = 0.2) -> None:
        if cancelled():
            raise Cancelled()
        if report:
            report(min(0.99, steps["done"] / steps["total"]), message)
        steps["done"] += weight

    def stop_if_cancelled() -> None:
        if cancelled():
            raise Cancelled()

    def progressive() -> tuple[list[tuple[int, int]], dict[int, str]]:
        n = len(members)
        tick(f"Comparing the {n} sequences with each other", 1)
        joins = guide_tree([m.oriented for m in members])
        groups: dict[int, list[int]] = {i: [i] for i in range(n)}
        rows: dict[int, list[str]] = {i: [m.oriented] for i, m in enumerate(members)}
        for k, (a, b) in enumerate(joins):
            tick(f"Aligning: step {k + 1} of {len(joins)}", 3)
            rows[n + k] = merge(rows.pop(a), rows.pop(b), stop_if_cancelled)
            groups[n + k] = groups.pop(a) + groups.pop(b)
        (top,) = groups
        return joins, dict(zip(groups[top], rows[top]))

    left_out, root = orient(members, tick)
    warnings += left_out
    if len(members) < 2:
        raise AlignmentError("None of the chosen sequences has anything recognisably in "
                             "common with the others, so there is nothing to align. Check "
                             "that they are the sequences you meant.")
    joins, gapped = progressive()
    if root.circular:
        cut, crossing = crossing_cut(root, members, gapped)
        if crossing and cut is None:
            warnings.append(f"{', '.join(crossing)} cross{'es' if len(crossing) == 1 else ''} "
                            f"base 1 of “{root.name}”, and the linear sequences cover all of "
                            "it between them, so the circle cannot be opened anywhere else: "
                            "the part past base 1 is shown in columns of its own at the end "
                            "of the circle.")
        elif crossing:
            # Open every circle elsewhere and align again, so the fragment is in one piece.
            steps["total"] += 1 + (len(members) - 1) * 3.2
            orient(members, tick, (root.rotation + cut) % len(root.sequence))
            joins, gapped = progressive()
    n = len(members)
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
            # A reverse-complemented circle starts at its own last base even unturned.
            "origin_moved_to": (None if not m.circular or (m.strand == 1 and not m.rotation)
                                else m.rotation + 1 if m.strand == 1
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

def _column_scores(columns: list[tuple[str, str, int | None]]) -> list[float]:
    """The blastn score of each column, as in the star alignment (reward 2, penalty -3, gap
    open 5 + extend 2 on the first gap of a run, 2 on the others; N scores 0, it says
    nothing). Not the alignment's own scores: their cheap gaps, right for lining up many
    sequences, let two unrelated stretches score above zero once aligned."""
    out = []
    gap = None
    for x, y, _ in columns:
        if x == "-" or y == "-":
            side = "a" if x == "-" else "b"
            out.append(-2.0 if gap == side else -7.0)
            gap = side
            continue
        gap = None
        out.append(0.0 if "N" in (x, y) else 2.0 if x == y else -3.0)
    return out


def maximal_segments(scores: list[float]) -> list[tuple[int, int, float]]:
    """All maximal scoring segments (Ruzzo & Tompa 1999, ISMB 7:234) as (start, end
    exclusive, score): the stretches a local alignment would report one after another."""
    segments: list[list[float]] = []      # [start, end, left cumulative, right cumulative]
    total = 0.0
    for k, score in enumerate(scores):
        before = total
        total += score
        if score <= 0:
            continue
        current = [k, k + 1, before, total]
        while True:
            j = next((i for i in range(len(segments) - 1, -1, -1)
                      if segments[i][2] < current[2]), None)
            if j is None or segments[j][3] >= current[3]:
                segments.append(current)
                break
            current = [segments[j][0], current[1], segments[j][2], current[3]]
            del segments[j:]
    return [(int(a), int(b), r - lft) for a, b, lft, r in segments]


def homologous_segments(scores: list[float], offset: int = 0) -> list[tuple[int, int]]:
    """The stretches of a pairwise comparison that are homologous, as (start, end exclusive):
    maximal scoring segments of at least SEGMENT_MIN, each split wherever its running score
    drops by more than X_DROP below the best it reached before (BLAST's X-drop rule)."""
    out = []
    for i, j, _ in maximal_segments(scores):
        part = scores[i:j]
        if sum(part) < SEGMENT_MIN:
            continue
        total = peak = 0.0
        peak_at, drop, cut = 0, 0.0, None
        for k, score in enumerate(part):
            total += score
            if total > peak:
                peak, peak_at = total, k + 1
            elif peak - total > drop:
                drop, cut = peak - total, (peak_at, k + 1)
        if drop <= X_DROP:
            out.append((offset + i, offset + j))
            continue
        p, t = cut
        out += homologous_segments(part[:p], offset + i)
        out += homologous_segments(part[t:], offset + i + t)
    return sorted(out)


def _realign_join(columns: list[tuple[str, str, int | None]], length: int
                  ) -> list[tuple[str, str, int | None]]:
    """Two circles' columns with the stretch around the rows' cut aligned again in one
    piece. The rows of a multiple alignment are linear: an insertion or deletion next to
    where the circles were opened is split over the two ends, half at each, which read
    round would be two events instead of one."""
    window = min(JOIN_WINDOW, len(columns) // 4)
    if window < ANCHOR_WORD:
        return columns
    join = columns[-window:] + columns[:window]
    a = [(x, p) for x, _, p in join if x != "-"]
    b = "".join(y for _, y, _ in join if y != "-")
    if not a or not b:
        return columns
    aligner = Align.PairwiseAligner(scoring="blastn")
    pairs = _columns_of(aligner.align("".join(x for x, _ in a), b)[0].coordinates)
    redone = [(a[i][0] if i is not None else "-", b[j] if j is not None else "-",
               a[i][1] if i is not None else None) for i, j in pairs]
    return redone + columns[window:-window]


def _longest_identical_run(columns: list[tuple[str, str, int | None]]) -> int:
    """The middle column of the longest run of identical bases (0 when there is none)."""
    best, best_at, run = 0, 0, 0
    for k, (x, y, _) in enumerate(columns):
        run = run + 1 if x == y != "-" and x != "N" else 0
        if run > best:
            best, best_at = run, k - run // 2
    return best_at


def _stretch(columns: list[tuple[str, str, int | None]]) -> dict[str, Any]:
    """A stretch left out of a comparison, in the first sequence's numbering."""
    on_a = [p for _, _, p in columns if p is not None]
    return {"start": on_a[0] + 1 if on_a else None, "end": on_a[-1] + 1 if on_a else None,
            "a_bases": len(on_a), "b_bases": sum(1 for _, y, _ in columns if y != "-")}


def compare(a: dict[str, Any], b: dict[str, Any]) -> dict[str, Any]:
    """Differences of row b from row a, in a's own numbering (1-based, on a's strand as
    stored in its file), over the columns both cover.

    A multiple alignment lines up every column, related or not: two plasmids that share a
    backbone but carry different inserts face each other insert against insert. Such a
    stretch is not a list of mutations. Only homologous segments are compared
    (``homologous_segments``); a stretch between or around them is reported as unrelated
    when it costs more than X_DROP and both sequences carry bases there. A divergent but
    homologous loop costs less (the variable regions of AAV2 and AAV9 VP1 cost 90 at most)
    and stays compared; so does an insertion or deletion of any length (bases on one side
    only). Unrelated sequence on both sides is recognised from about 150 bp. N facing a
    base is not a difference; such positions are counted apart."""
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
    empty = {"differences": [], "counts": _counts([]), "identity": 0.0, "start": None,
             "end": None, "overlap": 0, "unrelated": [], "unknown": 0}
    columns: list[tuple[str, str, int | None]] = []
    position = 0
    for x, y in pairs:
        if x == "-":
            columns.append((x, y, None))
        else:
            columns.append((x, y, (position + rotation) % length))
            position += 1
    both = [k for k, (x, y, _) in enumerate(columns) if x != "-" and y != "-"]
    if not both:
        return empty
    if a.get("circular") and b.get("circular"):
        # Two circles: the alignment's last column is followed by its first on both, so read
        # it round, starting inside the longest identical run. A difference where the rows
        # happen to be cut is then reported whole, never split over the two ends.
        columns = rotate(_realign_join(columns, length), _longest_identical_run(columns))
    else:
        # Columns where only one of the two has begun or still goes on (a fragment's ends,
        # an overhang) are not compared, as in the star alignment.
        columns = columns[both[0]:both[-1] + 1]
    scores = _column_scores(columns)
    segments = [list(seg) for seg in homologous_segments(scores)]
    if not segments:
        return {**empty, "unrelated": [_stretch(columns)]}
    # Ends: a stretch before the first homologous segment or after the last one is still
    # compared unless it costs more than X_DROP (a difference next to the cut of a circle,
    # or at the first base, must not vanish).
    if sum(scores[:segments[0][0]]) > -X_DROP:
        segments[0][0] = 0
    if sum(scores[segments[-1][1]:]) > -X_DROP:
        segments[-1][1] = len(columns)
    # Between two segments: an insertion or deletion (bases on one side only, give or take a
    # few) is a difference however long; so is a stretch that costs less than X_DROP.
    merged = [segments[0]]
    for seg in segments[1:]:
        between = columns[merged[-1][1]:seg[0]]
        a_bases = sum(1 for x, _, _ in between if x != "-")
        b_bases = sum(1 for _, y, _ in between if y != "-")
        if (min(a_bases, b_bases) <= max(INDEL_SLACK, INDEL_SHARE * max(a_bases, b_bases))
                or sum(scores[merged[-1][1]:seg[0]]) > -X_DROP):
            merged[-1][1] = seg[1]
        else:
            merged.append(seg)
    diffs: list[dict[str, Any]] = []
    compared: list[tuple[str, str, int | None]] = []
    unrelated = []
    edges = [0] + [x for seg in merged for x in seg] + [len(columns)]
    for k in range(0, len(edges), 2):
        if edges[k] < edges[k + 1]:
            unrelated.append(_stretch(columns[edges[k]:edges[k + 1]]))
    for i, j in merged:
        part = columns[i:j]
        compared += part
        diffs += differences(part, length)
    unknown = [d for d in diffs if d["kind"] == "mismatch" and "N" in (d["reference"],
                                                                      d["sequence"])]
    diffs = sorted((d for d in diffs if d not in unknown), key=lambda d: d["start"])
    same = sum(1 for x, y, _ in compared if x == y)
    counted = len(compared) - len(unknown)
    first = next(p for _, _, p in compared if p is not None)
    last = next(p for _, _, p in reversed(compared) if p is not None)
    if a.get("circular") and len({p for _, _, p in compared if p is not None}) == length:
        first, last = 0, length - 1          # all the way round
    return {
        "differences": diffs,
        "counts": _counts(diffs),
        "identity": round(100.0 * same / max(1, counted), 2),
        "start": first + 1,
        "end": last + 1,
        "overlap": sum(1 for _, _, p in compared if p is not None),
        "unrelated": unrelated,
        "unknown": len(unknown),
    }


# ---- background jobs and HTTP ----------------------------------------------------------

JOBS: dict[str, Job] = {}
JOBS_LOCK = threading.Lock()
MAX_JOBS = 8                 # finished jobs nobody collected are dropped beyond this


def start(docs: list[Any]) -> str:
    job = Job()
    job_id = uuid.uuid4().hex
    with JOBS_LOCK:
        # One window, one alignment at a time: a new one stops any still running (a
        # 10 x 15 kb alignment holds close to 1 GB while it runs).
        for other in JOBS.values():
            if not other.done.is_set():
                other.cancelled = True
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
