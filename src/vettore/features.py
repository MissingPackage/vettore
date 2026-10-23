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
import math
import re
from collections import deque
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
    """Where the feature could start in the target, from exact seed hits.

    `str.find` is load-bearing for the speed criterion, and not only for tidiness. It is a C
    scan over the whole document for every seed of every catalogue entry, which is most of what
    detection does on a plasmid whose elements are not diverged — the common case, since an
    element present verbatim is then taken whole below without being aligned at all. Pulling
    that loop into Python costs, measured as medians on the fixtures of `tests/test_features.py`
    and on the records in `tests/data`:

      * 40x on the verbatim AAV vector, 0.040 s to 1.60 s. This is the case the user meets most
        often, and the one the figures used to miss, because the diverged fixtures dilute the
        scan with alignment;
      * on the real records, 19x on `cloning-gateway-destination.gb` (0.065 s to 1.21 s), 9.0x
        on `PZ629150.1.gb`, and 2.8x on `cloning-multisite-entry-l4r1.gb`, which is the dearest
        record there is and so the one with the most alignment to dilute it. A review measured
        others in that directory above 30x; the factor is widest wherever least of the document
        has to be aligned, which is the case the user meets most;
      * 4.5x on the AAV fixture at three deleted sites, 0.47 s to 2.12 s;
      * 1.4x on the packed stress plasmid, where alignment dominates everything.

    Which tests see it depends on what the rewrite changes, and the distinction matters because
    the wrong reading of a failure here is to raise a bound:

      * a rewrite that keeps the scan **exact** — the same seeds, found in the same places — is
        seen by nothing else. Annotations are identical, every cell `ALIGNMENT_WORK` pins is
        identical, and every other test passes. The test named for this case,
        `test_the_document_is_not_walked_in_python_once_per_element`, charges the document for
        what is asked of it rather than timing the scan: 267 things for the engine against
        3,962,927 of its bases for the rewrite, and no clock;
      * a scan that tolerates a **mismatch** finds *more* starting points, so more guesses are
        aligned, and `ALIGNMENT_WORK` fails: on the AAV fixture 15.49 M cells become 42.20 M,
        on the verbatim one 0.34 M become 27.04 M, on the packed stress plasmid 111.75 M become
        152.15 M, with the annotations unchanged throughout. Those failures are the alarm
        working. The bounds are not stale and must not be raised to accommodate it;
      * a scan that reads **IUPAC ambiguity codes** in the document changes nothing today only
        because no catalogue entry carries an ambiguity code — 0 of the 36 — so no seed does
        either. Add one and it behaves like the mismatch case.

    So if tolerance is ever wanted here, the tolerant scan belongs in C or in a compiled
    matcher: a regex-based reader of ambiguity codes takes the document into C and is charged
    nothing at all, finds the same offsets, and passes both guards. And whatever route is taken,
    the cost has to be measured by hand on a *verbatim* plasmid before it is believed, because
    that is where the factor is largest and where no timed test can see it.
    """
    offsets: set[int] = set()
    for position, seed in _seeds(feature):
        start = target.find(seed)
        while start != -1:
            offsets.add(max(0, start - position))
            start = target.find(seed, start + 1)
    return offsets


#: (start, end, identity, first, last, frameshift) of a hit: start and end on the target,
#: inclusive; first and last the bases of the query it covers, as a 0-based half-open range;
#: frameshift whether the indels inside a coding sequence break its reading frame.
Hit = tuple[int, int, float, int, int, bool]

Column = tuple[int | None, int | None, bool]

STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})


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
    entry_length: int    # length of the catalogue entry
    part: tuple[int, int]  # the bases of the entry the match covers, 1-based inclusive
    frameshift: bool = False  # a coding sequence whose indels break its reading frame

    def span(self, length: int) -> set[int]:
        if self.start <= self.end:
            return set(range(self.start, self.end + 1))
        return set(range(self.start, length)) | set(range(self.end + 1))

    def as_dict(self, length: int) -> dict[str, Any]:
        return {"name": self.name, "type": self.type, "start": self.start, "end": self.end,
                "strand": self.strand, "forward": self.strand == 1,
                "identity": shown_identity(self.identity), "note": self.note,
                "source": self.source, "length": len(self.span(length)),
                "catalogueLength": self.entry_length, "catalogueFrom": self.part[0],
                "catalogueTo": self.part[1], "frameshift": self.frameshift}


def shown_identity(identity: float) -> float:
    """The identity to one decimal, rounded down: one mismatch in SaCas9 is 99.97 %, and it
    must not read as 100.0 %."""
    return math.floor(identity * 10 + 1e-9) / 10


def _columns(alignment: Any, offset: int = 0) -> list[Column]:
    """(target index, query index, identical) for each column of a pairwise alignment, with
    None for the sequence that has a gap there; target indices are shifted by ``offset``."""
    columns: list[Column] = []
    t, q = offset + int(alignment.coordinates[0][0]), int(alignment.coordinates[1][0])
    for a, b in zip(alignment[0], alignment[1]):
        columns.append((t if a != "-" else None, q if b != "-" else None, a == b and a != "-"))
        t += a != "-"
        q += b != "-"
    return columns


def _whole_codons(columns: list[Column]) -> list[Column]:
    """The columns from the first base of a codon of the catalogue entry to the last base of
    one, both present in the target, so a coding sequence is annotated in its own frame."""
    first = next((i for i, (t, q, _) in enumerate(columns)
                  if t is not None and q is not None and q % 3 == 0), None)
    last = next((i for i in range(len(columns) - 1, -1, -1)
                 if columns[i][0] is not None and columns[i][1] is not None
                 and columns[i][1] % 3 == 2), None)
    if first is None or last is None or first > last:
        return []
    return columns[first : last + 1]


def _codon_ends(columns: list[Column], target: str, query: str, strand: int) -> list[Column]:
    """A coding sequence cut to whole codons, with its own start and stop codons put back when
    the vector has them in frame next to the match: a TAG where the catalogue has TAA is still
    the stop, and an ATG one codon before the match is the start, after a codon put in to make
    an NcoI site (ATG GGG ATT) too. The codon put in is counted as an insertion, and is only
    accepted in that context (it starts with G, as CCATGG requires, and is not a stop), so a
    Met in the linker of a fusion is not taken for the start. On the minus strand the query is
    the entry's reverse complement, so the stop codon is at its head."""
    starts, stops = {"ATG"}, STOP_CODONS
    if strand == -1:
        starts, stops = {_reverse_complement(c) for c in stops}, {_reverse_complement("ATG")}
    head_reach, tail_reach = (2, 1) if strand == 1 else (1, 2)
    length = len(query)

    first_t = next(t for t, _, _ in columns if t is not None)
    first_q = next(q for _, q, _ in columns if q is not None)
    if first_q == 3:
        for k in range(1, head_reach + 1):
            begin = first_t - 3 * k
            put_in = target[begin + 3 : first_t]
            if begin >= 0 and target[begin : begin + 3] in starts and _ncoi_codon(put_in, strand):
                columns = ([(begin + i, i, target[begin + i] == query[i]) for i in range(3)]
                           + [(begin + 3 + i, None, False) for i in range(3 * (k - 1))]
                           + columns)
                break

    last_t = next(t for t, _, _ in reversed(columns) if t is not None)
    last_q = next(q for _, q, _ in reversed(columns) if q is not None)
    if last_q == length - 4:
        for k in range(1, tail_reach + 1):
            begin = last_t + 1 + 3 * (k - 1)
            put_in = target[last_t + 1 : begin]
            if target[begin : begin + 3] in stops and _ncoi_codon(put_in, strand):
                columns = (columns
                           + [(last_t + 1 + i, None, False) for i in range(3 * (k - 1))]
                           + [(begin + i, length - 3 + i,
                               target[begin + i] == query[length - 3 + i]) for i in range(3)])
                break
    return columns


def _ncoi_codon(codon: str, strand: int) -> bool:
    """Whether a codon put in after the ATG (none at all is fine too) is the one an NcoI site
    makes, read on the coding strand: it starts with G and is not a stop."""
    if not codon:
        return True
    if strand == -1:
        codon = _reverse_complement(codon)
    return codon[0] == "G" and codon not in STOP_CODONS


def _frameshift(columns: list[Column]) -> bool:
    """Whether the gaps inside an alignment break its reading frame: overall, or over a codon
    or more between an insertion and a deletion that make up for each other."""
    shift = 0
    out_of_frame = 0
    for t, q, _ in columns:
        shift += (q is None) - (t is None)
        if shift % 3 and q is not None:
            out_of_frame += 1
            if out_of_frame >= 3:
                return True
        elif shift % 3 == 0:
            out_of_frame = 0
    return shift % 3 != 0


#: Two candidate positions this close together are the same copy of the element read in two
#: registers, not two copies: an indel inside an element shifts every seed after it by the
#: length of the indel, so one deleted base puts the seeds before it and the seeds after it a
#: base apart, and a lost codon three. Aligning such a run once, over a window that holds all
#: of it, costs one pass where it used to cost one per register, and the window grows by the
#: distance between the first and the last of them rather than by its own length again.
MERGE_GAP = 4
#: The furthest a run may reach, however many registers line up: beyond a seed's length the
#: registers no longer come from one indel inside one copy. `_merge_limit` narrows both figures
#: for an entry that repeats itself, which is what makes the merge safe; see there.
MERGE_SPAN = SEED


@lru_cache(maxsize=4096)
def _self_similarity(query: str, shift: int) -> float:
    """How much of the entry matches itself shifted ``shift`` bases along, without gaps."""
    if shift >= len(query):
        return 0.0
    same = sum(a == b for a, b in zip(query[shift:], query[:-shift]))
    return 100.0 * same / (len(query) - shift)


def _merge_limit(query: str, min_identity: float, min_coverage: float = 100.0,
                 span: int = MERGE_SPAN) -> int:
    """How far apart two guesses may be before they could be two copies rather than one.

    This is what makes the merge safe, and it is not a constant. A reported copy of an entry
    `L` bases long agrees with the entry at all but ``slack`` of its bases: it may leave
    ``100 - min_coverage`` per cent of them out of the alignment altogether, and differ at
    ``100 - min_identity`` per cent of the rest. Two copies ``shift`` bases apart overlap over
    ``L - shift`` bases, and every base of that overlap where both copies agree with the entry
    is a base where the entry agrees with *itself* shifted ``shift`` along. So the entry has to
    be at least ``(L - shift - 2 * slack) / (L - shift)`` the same as itself at that shift for
    two copies to sit there at all, and below that figure two guesses that close together
    cannot be two copies — aligning them once cannot lose one.

    Both terms matter, and the coverage one is the one that bites. A `repeat_region` entry only
    has to be 70 % covered, so a copy may ignore a third of the entry and still be reported at
    100 % identity: a (CAG)16 polyQ tract with unique flanks is reported on its tract alone,
    nine times over in a longer tract, and an identity-only bound merges six of them away.

    The bound is ungapped, which is why the catalogue invariant is a test of its own rather
    than something to infer from here.
    """
    length = len(query)
    slack = length * (1.0 - min_identity * min_coverage / 10_000.0)
    for shift in range(1, span + 1):
        if shift >= length:
            return span
        allowed = 100.0 * (1.0 - 2.0 * slack / (length - shift))
        if _self_similarity(query, shift) >= allowed:
            return shift - 1
    return span


def _offset_runs(offsets: list[int], gap: int = MERGE_GAP,
                 span: int = MERGE_SPAN) -> list[tuple[int, int]]:
    """``offsets``, in order, as (first, last) runs of neighbours: one alignment each."""
    runs: list[tuple[int, int]] = []
    for offset in sorted(offsets):
        if runs and offset - runs[-1][1] <= gap and offset - runs[-1][0] <= span:
            runs[-1] = (runs[-1][0], offset)
        else:
            runs.append((offset, offset))
    return runs


def _match_long(target: str, query: str, length: int, min_identity: float,
                min_coverage: float, strand: int = 1, coding: bool = False) -> list[Hit]:
    """Matches of ``query`` (the catalogue entry, or its reverse complement on the minus
    ``strand``) in ``target``. With ``coding``, each is annotated in whole codons of its own
    frame, from its start to its stop codon when the vector has them (the catalogue's CDS
    entries run from ATG to stop, so the frame is the same on both strands)."""
    results: list[Hit] = []
    aligner = _aligner()
    margin = max(30, len(query) // 10)
    inexact: list[int] = []
    for offset in sorted(_candidate_offsets(target, query)):
        # A guess at or past the end of the sequence itself can only be the copy a circular
        # search appends, and the same seed hit also appears in the sequence a length earlier,
        # so that copy is found there instead. Aligning this one would be wasted work.
        if offset >= length:
            continue
        # An element present verbatim is the best local alignment of itself under this
        # scoring, in full and at no offset, so it needs no alignment: most catalogue elements
        # sit in a vector unchanged, and this is what keeps detection quick on a plasmid full
        # of them. Every guess is tested, because a repeated element (a 6×HA tag is two 3×HA
        # in a row) gives one guess per copy and per register.
        if 100.0 >= min_identity and target[offset : offset + len(query)] == query:
            results.append((offset, offset + len(query) - 1, 100.0, 0, len(query), False))
            continue
        inexact.append(offset)
    # Four is the narrow part of the rule and is meant to be: it is as far apart as the
    # registers of one copy get when the indel that split them is a lost codon or two. A cut
    # landing inside a seed window pushes the surviving registers further apart, and widening
    # the gap to the span cap would merge those too. Measured, it helps only where the cuts
    # are artificially even and costs almost everywhere else. The AAV fixture's four-site
    # trough halves, 29.6 M cells of dynamic programming to 15.3 M, and the all-in-one
    # fixture's goes 87.3 M to 62.1 M; of the points the tests pin, aav 1/3/5 and all-in-one
    # 3/5 move by under a per cent, the verbatim AAV vector does not move at all, and the
    # remaining five get *worse* — all-in-one 0 by 5.6 %, all-in-one 1 by 2.4 %, the packed
    # stress plasmid by 18.4 % with no deleted site, 4.4 % with one and 4.6 % with three. A
    # wider run is a wider window, a wider window drifts off its own guesses more often, and
    # each one it drifts off is paid for separately below. The
    # troughs are a property of evenly spaced cuts coinciding with the seed positions rather
    # than of divergence, so there is nothing here worth the extra base of distance at which
    # one alignment would have to stand in for two guesses.
    #
    # Which is the thing this loop checks rather than assumes. One window yields one
    # alignment, and the best alignment over a run's window is not always the best alignment
    # over the window of each guess in it: a wider window can hold a better-scoring placement
    # elsewhere and drift to it, and then the copy a narrower window would have found is gone.
    # A 6×HA tag with four bases inserted does exactly that — the run's window reaches back
    # over the first copy and reports it twice instead of reporting the second. So after the
    # alignment, the guesses it actually passes through (its diagonals, one per register it
    # reads) are the ones it stands in for; any other guess in the run is put back on the
    # worklist on its own, which is what it would have cost before the merge. The common case,
    # the registers of one copy, is wholly accounted for by its own alignment and still costs
    # one pass.
    limit = _merge_limit(query, min_identity, min_coverage) if len(inexact) > 1 else 0
    pending = deque((first, last, tuple(o for o in inexact if first <= o <= last))
                    for first, last in _offset_runs(inexact, MERGE_GAP, limit))
    while pending:
        first, last, group = pending.popleft()
        window_start = max(0, first - margin)
        window = target[window_start : last + len(query) + margin]
        hit: Hit | None = None
        diagonals: set[int] = set()
        if len(window) >= SEED:
            alignments = aligner.align(window, query)
            if len(alignments):
                # The annotation is the aligned region and nothing else, so the identity tested
                # here is the identity the user reads. Ends of the catalogue entry that the
                # local alignment left out are not added back: under this scoring a dropped end
                # is at least 40 % mismatches, and in real vectors it is the junction the
                # element was cloned into (a BamHI site before an SV40 polyA, the linker of a
                # GFP fusion), not the element. The match says which bases of the entry it
                # covers, so a shorter match is never silent.
                columns = _columns(alignments[0], window_start)
                diagonals = {t - q for t, q, _ in columns if t is not None and q is not None}
                choices = [columns]
                if coding:
                    cut = _whole_codons(columns)
                    # The start and stop codons count towards the identity; if they would take
                    # it under the threshold, the match is shown without them rather than not
                    # at all.
                    choices = [_codon_ends(cut, target, query, strand), cut] if cut else []
                for chosen in choices:
                    in_target = [t for t, _, _ in chosen if t is not None]
                    in_query = [q for _, q, _ in chosen if q is not None]
                    coverage = 100.0 * len(in_query) / len(query)
                    identity = 100.0 * sum(same for _, _, same in chosen) / len(chosen)
                    if identity >= min_identity and coverage >= min_coverage:
                        hit = (in_target[0], in_target[-1], identity, in_query[0],
                               in_query[-1] + 1, coding and _frameshift(chosen))
                        break
        # A guess is stood in for only by an answer that was accepted and that the alignment
        # passes through. When the alignment over the wider window scores highest on a
        # placement that then fails the thresholds, no guess in the run has been answered,
        # however many of them lie on it, so they all go back: identity and coverage are not
        # monotone in alignment score, so a window that can only score higher can still come
        # back with nothing where a narrower one would have reported a match. Belt and braces:
        # a review looked for a construct where this half of the condition is the only thing
        # holding — the run rejected, its diagonals over a guess, that guess accepted on its
        # own — and found none in 5262 runs of more than one guess, so no test isolates it.
        # `len(group) > 1` is what makes this terminate as well as what makes it cheap: a
        # re-queued guess is a group of one, so it can never put anything back, and the
        # worklist is at most the runs plus the guesses in them.
        missed = ([o for o in group if hit is None or o not in diagonals]
                  if len(group) > 1 else [])
        # The alignment is this run's answer when it stands in for at least one of the guesses
        # it was made from; the guesses it never passes through are aligned on their own,
        # next, so that the results stay in the order the guesses are in.
        if hit is not None and len(missed) < len(group):
            results.append(hit)
        # When it stands in for none of them it has drifted to a placement some other guess
        # will find on its own, so its answer is dropped rather than reported from here: that
        # leaves the results exactly the ones the guesses give one at a time.
        # At the front and in order, so that a guess put back is aligned where it would have
        # been had the run never been made: `_dedupe` keeps whichever of two overlapping
        # copies it meets first, and an earlier version of this loop reordered the results and
        # silently changed which of nine overlapping copies of a repeat survived. Nothing in
        # the suite distinguishes this from appending, because a run either answers every
        # guess in it or none, so the ones put back are the whole run and come back adjacent
        # either way; it is written this way because the partial case is the one that would
        # not, and it costs nothing to be right about it in advance.
        pending.extendleft((o, o, (o,)) for o in reversed(missed))
    return results


def _match_short(target: str, query: str) -> list[Hit]:
    results = []
    start = target.find(query)
    while start != -1:
        results.append((start, start + len(query) - 1, 100.0, 0, len(query), False))
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
                    else _match_long(target, query, length, min_identity, coverage,
                                     strand=strand,
                                     coding=feature.type == "CDS" and feature.length % 3 == 0))
            for start, end, identity, first, last, frameshift in hits:
                if start >= length:
                    continue  # the same hit also appears in the first copy
                # The query is the entry's reverse complement on the minus strand.
                part = ((first + 1, last) if strand == 1
                        else (feature.length - last + 1, feature.length - first))
                found.append(Match(name=feature.name, type=feature.type, start=start % length,
                                   end=end % length if circular else min(end, length - 1),
                                   strand=strand, identity=identity, note=feature.note,
                                   source=feature.source, entry_length=feature.length,
                                   part=part, frameshift=frameshift))
    return _dedupe(found, length)


#: Types whose elements are folded into a longer element the catalogue says holds them (see
#: ``HOLDS``): the CMV enhancer at the head of a full-length CAG promoter is the CAG's own
#: first third, not a second element. A tag inside a CDS (3×HA in SaCas9) stays a match of its
#: own, which is why this is a short list rather than the whole of ``HOLDS``.
PART_OF_TYPES = frozenset({"promoter", "enhancer"})


def _dedupe(matches: list[Match], length: int) -> list[Match]:
    """Drop matches that describe the same region: the same feature found twice by different
    seeds, or two catalogue entries (for example the 5′ and 3′ ITR) that are the same sequence,
    or a regulatory element the catalogue lists inside a longer one that was found here too.
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
            # Part of the longer element, by the catalogue's own account.
            part_of = (candidate.name in HOLDS.get(existing.name, ())
                       and existing.type in PART_OF_TYPES
                       and candidate.type in PART_OF_TYPES
                       and overlap >= 0.9 * len(span))
            if same_region or variant or part_of:
                duplicate = True
                break
        if not duplicate:
            kept.append(candidate)
    return sorted(kept, key=lambda m: (m.start, -len(m.span(length))))


# A match is already annotated when an existing feature covers this much of it, whatever its
# strand (an origin is often annotated without a direction), and either
# - the match is at least ANNOTATED_SHARE_OF_EXISTING of that feature, whatever its name (so a
#   feature spanning half the plasmid does not hide the elements inside it), or
# - the feature's name names the element: vendors annotate one block, "CMV promoter", over the
#   CMV enhancer and the CMV promoter.
ANNOTATED_SHARE_OF_MATCH = 0.9
ANNOTATED_SHARE_OF_EXISTING = 0.5

#: A feature named after one element whose sequence is mostly another (nearly exact matches,
#: at least LABEL_CHECK_MIN_IDENTITY, covering at least LABEL_CHECK_SHARE_OF_EXISTING of it)
#: probably carries a wrong label: a CMV enhancer and promoter labelled "PGK promoter". A
#: short element at the edge of a longer feature (a polyA signal inside an "EGFP" block) is
#: not enough to doubt the label.
LABEL_CHECK_MIN_IDENTITY = 98.0
LABEL_CHECK_SHARE_OF_EXISTING = 0.5

#: Words that do not tell elements apart: "SV40 polyA" and "bGH polyA" share "polyA" but are
#: different elements. Words of one character and plain numbers ("5'", "3'") never count.
GENERIC_WORDS = frozenset({
    "promoter", "enhancer", "polya", "pa", "signal", "ori", "origin", "region", "sequence",
    "gene", "late", "early", "of", "the", "and", "with", "for", "from", "in", "to",
    "replication", "polyadenylation", "terminator", "site", "element", "truncated", "partial",
    "full", "length", "minimal", "core", "fragment", "intron", "exon", "cds", "orf", "protein",
    "tag", "repeat", "inverted", "terminal", "human", "mouse", "rat", "chicken", "codon",
    "optimized", "optimised", "resistance", "marker", "cassette", "insert", "left", "right",
    "forward", "reverse", "fwd", "rev", "start", "stop", "binding", "misc", "feature", "variant",
    "short", "long", "copy", "upstream", "downstream", "fusion", "wild", "type", "wt", "new",
    "rabbit",
})

#: Other spellings of a word, as found in vendors' and labs' labels.
WORD_ALIASES = {"hcmv": "cmv", "hef1a": "ef1a", "camkii": "camkiia", "hsyn1": "hsyn",
                "hpgk": "pgk", "actb": "bactin"}

#: Catalogue elements that hold another catalogue element, so a block named after the first is
#: not a wrong label for the second (checked against the catalogue sequences in tests):
#: - SaCas9: the catalogue entry carries a C-terminal 3×HA tag (PZ629150.1:3842..3922, inside
#:   581..3925).
#: - CAG promoter: the CMV enhancer, the chicken β-actin promoter and the chimeric chicken
#:   β-actin/rabbit β-globin intron (Niwa H, Yamamura K, Miyazaki J, Gene 1991;108:193-199).
#:   The full-length entry (ON012604.1:4799..6471) opens with the whole CMV enhancer; the
#:   chicken β-actin promoter that follows it is not a catalogue element of its own.
#: - CAG promoter (truncated): the same elements, but it starts 131 bp inside the enhancer, so
#:   only part of it is there and a search does not find the enhancer in it.
#: - EF-1α promoter (human): the EFS promoter is its intron-less core — EFS bases 7..239 are
#:   EF-1α bases 6..238, 233 of them — so a search finds EFS inside it and a block named
#:   "EF-1α" is not a wrong label for the core.
#: - WPRE: WPRE3 is the WPRE with its middle taken out, in two pieces (WPRE3 3..160 is WPRE
#:   1..158, WPRE3 166..246 is WPRE 302..382) and seven bases in neither, so it is not a
#:   contiguous stretch of the WPRE and a search does not find it inside one. Listed here by
#:   hand, like the truncated CAG below: a cassette carrying the short element and labelled
#:   "WPRE" by its author is labelled loosely, not wrongly.
HOLDS = {"SaCas9": frozenset({"3×HA"}),
         "CAG promoter": frozenset({"CMV enhancer"}),
         "CAG promoter (truncated)": frozenset({"CMV enhancer"}),
         "EF-1α promoter (human)": frozenset({"EFS promoter (human)"}),
         "WPRE": frozenset({"WPRE3"})}

_LETTERS = str.maketrans({"α": "a", "β": "b", "γ": "g", "δ": "d", "Δ": "d", "κ": "k",
                          "×": "x", "′": "'", "’": "'"})


def name_words(name: str, pieces: bool = True) -> frozenset[str]:
    """The words of a feature name that can identify an element, lower case: "EF-1α
    promoter" gives "ef1a", "CMV-promoter" gives "cmv". With ``pieces`` the parts of a
    hyphenated word count too ("cmv" from "CMV-promoter"); catalogue names are read without,
    so that "EF-1α" is not matched by any "1a"."""
    text = name.translate(_LETTERS).lower().replace("alpha", "a").replace("beta", "b")
    words: set[str] = set()
    for token in re.split(r"[\s/(),;:_\[\]{}|+=]+", text):
        parts = [p for p in re.split(r"[-.']+", token) if p]
        for word in {"".join(parts), *(parts if pieces else [])}:
            word = WORD_ALIASES.get(word, word)
            words.add(word)
            if tagged := re.fullmatch(r"\d+x(\D\w*)", word):
                words.add(tagged.group(1))       # "3×HA" names the HA tag too
    words = {w for w in words if len(w) > 1 and not w.isdigit() and w not in GENERIC_WORDS}
    return frozenset(words)


@lru_cache(maxsize=1)
def _element_words() -> tuple[tuple[str, frozenset[str]], ...]:
    """(element name, its identifying words) for every element the check knows. They are the
    catalogue's own names, so a label that names one is always read against a real sequence:
    "PGK promoter" over a CMV promoter is doubted because the catalogue carries the two PGK
    promoters and they are not this sequence."""
    names = [f.name for f in catalogue()]
    return tuple((name, words) for name in dict.fromkeys(names)
                 if (words := name_words(name, pieces=False)))


def _label_words(label: str) -> frozenset[str]:
    """A label's identifying words, with a species prefix also read off ("mU6" names U6).
    What is left after the prefix is read like any other word, so a spelling ``WORD_ALIASES``
    knows is understood there too ("hACTB" names the β-actin promoter)."""
    words = name_words(label)
    stripped = {w[1:] for w in words if w[:1] in "hm" and len(w) > 2}
    return words | {WORD_ALIASES.get(w, w) for w in stripped}


def named_elements(label: str) -> set[str]:
    """The known elements a label names: "CMV promoter" names the catalogue's CMV enhancer and
    CMV promoter; "my insert" names none."""
    words = _label_words(label)
    return {name for name, element in _element_words() if words & element}


def _names_element(label: str, element: str) -> bool:
    """Whether a feature called ``label`` annotates the catalogue ``element``: the label
    shares an identifying word with it ("CMV" and "CMV enhancer")."""
    return bool(_label_words(label) & name_words(element, pieces=False))


def _names_another_element(label: str, element: str) -> bool:
    """Whether ``label`` names a known element other than ``element`` and not one that holds
    it (a "CAG" block holds the CMV enhancer, so that label is right).

    An element the label names but contradicts does not count, because the label is then
    not claiming to be that element at all: it only shares a word with it. "chicken β-actin
    promoter" shares "β-actin" with the human β-actin promoter and names the other species;
    "rabbit β-globin intron" shares "β-globin" with the rabbit β-globin polyA signal and
    names another kind of element. Both are right labels for the CAG promoter they sit on,
    and before this both raised a wrong-label warning on every CAG a user annotated by hand.
    """
    others = {name for name in named_elements(label)
              if not _species_conflict(label, name) and not _kind_conflict(label, name)}
    return bool(others) and not any(element in HOLDS.get(name, ()) for name in others)


#: Kinds a label can name outright, where naming one rules the others out: a feature called
#: "…intron" is not a polyA signal, whatever name it shares with one. Promoters and enhancers
#: are deliberately absent, because a CMV enhancer-and-promoter block is labelled either way
#: and neither label is wrong.
#: A word in a label and the catalogue types it rules out, so that a label sharing a name with
#: an entry of a kind it cannot be is not read as claiming to be that entry: "rabbit β-globin
#: intron", the piece a CAG promoter is made of, shares "β-globin" with the rabbit β-globin
#: polyA signal and is not it.
#:
#: Stated as what each word excludes, not as what it allows. Reading it the other way around
#: made a word the catalogue has no type for exclude everything at once: nothing here is typed
#: ``intron``, so "EF-1α promoter with intron A" and "hPGK (intron-less)" contradicted all
#: sixty entries and switched the wrong-label warning off altogether.
#:
#: Promoters and enhancers do not exclude each other on purpose, because a CMV
#: enhancer-and-promoter block is labelled either way and neither label is wrong.
_KIND_WORDS = {"intron": frozenset({"polyA_signal"}),
               "polya": frozenset({"promoter", "enhancer"}),
               "polyadenylation": frozenset({"promoter", "enhancer"})}


def _kind_conflict(label: str, element: str) -> bool:
    """Whether ``label`` says outright that it is a kind of element the catalogue ``element``
    is not."""
    types = _element_types()
    if element not in types:
        return False
    tokens = {t for t in re.split(r"[^a-z0-9]+", label.translate(_LETTERS).lower()) if t}
    kinds = tokens & _KIND_WORDS.keys()
    if not kinds:
        return False
    # A label naming more than one kind ("β-globin intron and polyA") rules out only what every
    # one of them rules out, because it describes a composite and claims to be none of them.
    ruled_out = frozenset.intersection(*(_KIND_WORDS[kind] for kind in kinds))
    return types[element] in ruled_out


@lru_cache(maxsize=1)
def _element_types() -> dict[str, str]:
    """Each catalogue entry's name and type."""
    return {feature.name: feature.type for feature in catalogue()}


#: How labels say which species an element comes from: a word ("human", "murine"), or a
#: one-letter prefix on the element's own word ("hPGK", "mU6"). The species that are not
#: human or mouse are here because a label naming one contradicts a "(human)" or "(mouse)"
#: entry: "chicken β-actin promoter" is not the human β-actin promoter.
_SPECIES_WORDS = {"human": "human", "homo": "human", "mouse": "mouse", "murine": "mouse",
                  "mus": "mouse", "chicken": "chicken", "gallus": "chicken",
                  "rabbit": "rabbit", "oryctolagus": "rabbit", "rat": "rat",
                  "rattus": "rat", "bovine": "bovine", "woodchuck": "woodchuck"}
_PREFIXES = {"h": "human", "m": "mouse"}


def entry_species(element: str) -> str | None:
    """The species a catalogue entry names, from its "(human)" or "(mouse)" suffix."""
    found = re.search(r"\((human|mouse)\)", element)
    return found.group(1) if found else None


def label_species(label: str, element: str) -> str | None:
    """The species a feature label gives for ``element``, or None when it gives none."""
    tokens = [t for t in re.split(r"[^a-z0-9]+", label.translate(_LETTERS).lower()) if t]
    for token in tokens:
        if token in _SPECIES_WORDS:
            return _SPECIES_WORDS[token]
    words = name_words(element, pieces=False)
    for token in tokens:
        rest = token[1:]
        if token[:1] in _PREFIXES and len(rest) > 1 and WORD_ALIASES.get(rest, rest) in words:
            return _PREFIXES[token[:1]]
    return None


def _species_conflict(label: str, element: str) -> bool:
    """Whether a label names the element of the other species ("PGK promoter (human)" over
    the mouse PGK promoter, or "mU6" over the human U6)."""
    ours, theirs = entry_species(element), label_species(label, element)
    return bool(ours and theirs and ours != theirs)


def _identity_text(identity: float) -> str:
    shown = shown_identity(identity)
    return "100 %" if shown >= 100 else f"{shown:.1f} %"


def _int(value: Any) -> int | None:
    return value if isinstance(value, int) and not isinstance(value, bool) else None


def _existing_spans(existing: Any, length: int) -> list[tuple[str, set[int]]]:
    """(name, positions) of the document's features; anything malformed is skipped."""
    items = list(existing.values()) if isinstance(existing, dict) else existing
    out = []
    for feature in items if isinstance(items, list) else []:
        if not isinstance(feature, dict):
            continue
        locations = feature.get("locations")
        ranges = ([(_int(loc.get("start")), _int(loc.get("end"))) for loc in locations
                   if isinstance(loc, dict)] if isinstance(locations, list) and locations
                  else [(_int(feature.get("start")), _int(feature.get("end")))])
        span: set[int] = set()
        for start, end in ranges:
            if start is None or end is None or not (0 <= start < length and 0 <= end < length):
                continue
            span |= (set(range(start, end + 1)) if start <= end
                     else set(range(start, length)) | set(range(end + 1)))
        if span:
            out.append((str(feature.get("name") or feature.get("type") or "feature"), span))
    return out


def _label_warning(existing: str, match: Match) -> str:
    return (f"Annotated as “{existing}”, but the sequence is the {match.name} "
            f"({_identity_text(match.identity)} identity). Check the label.")


def detect_payload(sequence: str, circular: bool, min_identity: float = DEFAULT_MIN_IDENTITY,
                   existing: Any = None) -> dict[str, Any]:
    """Matches, flagged against the features a document already carries.

    ``isNew`` is false and ``annotatedAs`` names the existing feature when one already covers
    the match. ``labelWarning`` says when that feature is named after a different element
    while its sequence is mostly this one: the label is probably wrong (a CMV promoter
    labelled "PGK promoter"). ``label_warnings`` counts the existing features concerned.
    """
    length = len(sequence)
    matches = detect(sequence, circular, min_identity)
    already = _existing_spans(existing or [], length)

    found = []
    unexplained: dict[int, set[int]] = {}     # existing feature -> bases of other elements
    for match in matches:
        span = match.span(length)
        covering = [(index, name, len(span & other), len(other))
                    for index, (name, other) in enumerate(already)
                    if len(span & other) >= ANNOTATED_SHARE_OF_MATCH * len(span)]
        named = [item for item in covering if _names_element(item[1], match.name)]
        # A label naming this element but the other species is a wrong label, not an
        # annotation of it.
        other_species = [item for item in named if _species_conflict(item[1], match.name)]
        named = [item for item in named if item not in other_species]
        sized = [item for item in covering if item[2] >= ANNOTATED_SHARE_OF_EXISTING * item[3]]
        suspect = ([] if named or match.identity < LABEL_CHECK_MIN_IDENTITY
                   else [item for item in covering
                         if item in other_species
                         or _names_another_element(item[1], match.name)])
        for item in suspect:
            unexplained.setdefault(item[0], set()).update(span & already[item[0]][1])
        found.append((match, named, sized, suspect))
    doubtful = {index for index, bases in unexplained.items()
                if len(bases) >= LABEL_CHECK_SHARE_OF_EXISTING * len(already[index][1])}

    rows = []
    flagged: set[int] = set()
    for match, named, sized, suspect in found:
        row = match.as_dict(length)
        mislabelled = [item for item in suspect if item[0] in doubtful]
        # The feature that covers most of the match and, of those, fits it most closely.
        best = max(named or mislabelled or sized, key=lambda item: (item[2], -item[3]),
                   default=None)
        row["annotatedAs"] = best[1] if best else None
        row["isNew"] = best is None
        row["labelWarning"] = _label_warning(best[1], match) if best and mislabelled else None
        if best and mislabelled:
            flagged.add(best[0])
        rows.append(row)
    return {"features": rows, "catalogue_size": len(catalogue()),
            "new_count": sum(1 for r in rows if r["isNew"]),
            "label_warnings": len(flagged)}
