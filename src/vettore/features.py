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
import re
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


def _loci(offsets: set[int], margin: int) -> list[tuple[int, int]]:
    """``offsets`` grouped into the (first, last) of each locus to align.

    Several seeds of one feature usually land on one copy of it, each giving a slightly
    different guess at where it starts. Aligning every guess separately costs a full dynamic
    programming pass each and ends in the same match, which ``_dedupe`` then drops. Guesses
    less than a margin apart cannot be two copies of the feature — two copies are at least a
    feature length apart — so they are one locus, aligned once in a window wide enough to hold
    whichever guess is right.
    """
    loci: list[tuple[int, int]] = []
    for offset in sorted(offsets):
        if loci and offset - loci[-1][1] < margin:
            loci[-1] = (loci[-1][0], offset)
        else:
            loci.append((offset, offset))
    return loci


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


def _match_long(target: str, query: str, length: int,
                min_identity: float, min_coverage: float) -> list[tuple[int, int, float]]:
    results: list[tuple[int, int, float]] = []
    aligner = _aligner()
    margin = max(30, len(query) // 10)
    for first, last in _loci(_candidate_offsets(target, query), margin):
        window_start = max(0, first - margin)
        # A window that opens at or after the end of the sequence itself only ever matches in
        # the copy a circular search appends, where every hit is a duplicate of one in the
        # sequence and is dropped below. Aligning it would be wasted work.
        if window_start >= length:
            continue
        window = target[window_start : last + len(query) + margin]
        if len(window) < SEED:
            continue
        # An element that is present verbatim is the best local alignment of itself, so it
        # needs no alignment: most catalogue elements sit in a vector unchanged, and this is
        # what keeps detection quick on a plasmid full of them.
        exact = {offset for offset in (first, last)
                 if target[offset : offset + len(query)] == query}
        if exact:
            results.extend((offset, offset + len(query) - 1, 100.0) for offset in exact)
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
                    else _match_long(target, query, length, min_identity, coverage))
            for start, end, identity in hits:
                if start >= length:
                    continue  # the same hit also appears in the first copy
                found.append(Match(name=feature.name, type=feature.type, start=start % length,
                                   end=end % length if circular else min(end, length - 1),
                                   strand=strand, identity=identity, note=feature.note,
                                   source=feature.source))
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
})

#: Other spellings of a word, as found in vendors' and labs' labels.
WORD_ALIASES = {"hcmv": "cmv", "hef1a": "ef1a", "camkii": "camkiia", "hsyn1": "hsyn",
                "hpgk": "pgk"}

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
HOLDS = {"SaCas9": frozenset({"3×HA"}),
         "CAG promoter": frozenset({"CMV enhancer"}),
         "CAG promoter (truncated)": frozenset({"CMV enhancer"})}

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
    """A label's identifying words, with a species prefix also read off ("mU6" names U6)."""
    words = name_words(label)
    return words | {w[1:] for w in words if w[:1] in "hm" and len(w) > 2}


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
    it (a "CAG" block holds the CMV enhancer, so that label is right)."""
    others = named_elements(label)
    return bool(others) and not any(element in HOLDS.get(name, ()) for name in others)


#: How labels say which species an element comes from: a word ("human", "murine"), or a
#: one-letter prefix on the element's own word ("hPGK", "mU6").
_SPECIES_WORDS = {"human": "human", "homo": "human", "mouse": "mouse", "murine": "mouse",
                  "mus": "mouse"}
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
    return "100 %" if round(identity, 1) >= 100 else f"{identity:.1f} %"


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
