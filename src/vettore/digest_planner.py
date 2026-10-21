"""Diagnostic digest planner: the digest whose gel tells the expected clone apart from the
usual wrong ones (insert in the other orientation, empty vector, uncut parent vector...).

Candidates are single enzymes and pairs of enzymes from the set the user chose. Each candidate
is digested *in silico* on every molecule and run through the gel simulator's migration model
(``gel.py``). A candidate is **diagnostic** when, for every wrong clone, at least one visible
band of one lane runs apart from every band of the other lane by at least the gel's resolution
(``gel.RESOLUTION``, the distance at which the simulator draws two bands as one). Bands under
``gel.MIN_VISIBLE_BP`` or above the resolving range of the agarose percentage do not count.

Cut positions come from Biopython's REBASE data (the same enzymes the cloning engine and the
gel simulator use), so the sizes listed here are the sizes the gel simulator draws. Methylation
is decided by the interface, which holds REBASE's methylation table: it sends the enzymes that
Dam or Dcm methylation would block on these molecules as ``excluded``.

Ranking, deterministic: diagnostic only; then clear results (every wrong clone at least twice
the resolution away, every lane digested) before results at the limit of the gel; then fewer
enzymes; then an expected lane of 2 to 5 visible bands; then larger separations; then names.
Buffer compatibility is not ranked: the REBASE tables bundled with Vettore carry no buffer
activity.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations, pairwise
from typing import Any

from Bio import Restriction
from Bio.Seq import Seq

from . import gel
from .cloning.base import CloningError, dna_sequence
from .enzymes import PREFERRED_NAMES, build_database

#: Enzymes that cut the expected clone more often than this give smears of bands, not a pattern.
MAX_CUTS = 8
#: How many results the dialog lists.
TOP = 10
#: A result is "clear" when every wrong clone differs by at least this many resolutions.
CLEAR_FACTOR = 2

_DATABASE: dict | None = None


def _database() -> dict:
    global _DATABASE
    if _DATABASE is None:
        _DATABASE = build_database()
    return _DATABASE


def commercial_enzymes() -> list[str]:
    """Every enzyme of the nonredundant commercial set (the editor's database)."""
    return sorted((entry["name"] for entry in _database().values()), key=str.lower)


def _suppliers(name: str) -> int:
    """How many suppliers sell the enzyme: a tie goes to the one most labs already have."""
    return len(_database().get(name.lower(), {}).get("suppliers", []))


@dataclass
class Molecule:
    name: str          # how the dialog calls it: "Expected clone", "empty vector", a file name
    sequence: str
    circular: bool
    cuts: dict[str, list[int]] = field(default_factory=dict)

    @property
    def length(self) -> int:
        return len(self.sequence)


def molecule(name: str, data: dict[str, Any]) -> Molecule:
    sequence = data.get("sequence") if isinstance(data, dict) else None
    label = str((data or {}).get("name") or name)
    if not isinstance(sequence, str) or not sequence:
        raise CloningError(f"“{label}” has no sequence. Remove it from the list or choose "
                           "another file.")
    return Molecule(name, dna_sequence(label, sequence).upper(), bool(data.get("circular")))


def _cut_positions(molecules: list[Molecule], names: list[str]) -> None:
    """Fill ``cuts``: per enzyme, the top-strand cut positions (Biopython's, 1-based)."""
    enzymes = [getattr(Restriction, n) for n in names if hasattr(Restriction, n)]
    batch = Restriction.RestrictionBatch(enzymes)
    for mol in molecules:
        found = batch.search(Seq(mol.sequence), linear=not mol.circular)
        mol.cuts = {str(e): sorted(set(sites)) for e, sites in found.items() if sites}


def fragment_sizes(mol: Molecule, enzymes: tuple[str, ...]) -> list[int]:
    """Fragment lengths between top-strand cuts, largest first; the whole molecule if none."""
    positions = sorted({p for e in enzymes for p in mol.cuts.get(e, [])})
    if not positions:
        return [mol.length]
    if mol.circular:
        pairs = zip(positions, positions[1:] + [positions[0] + mol.length])
        sizes = [b - a for a, b in pairs]
    else:
        bounds = [1, *positions, mol.length + 1]
        sizes = [b - a for a, b in pairwise(bounds)]
    return sorted((s for s in sizes if s > 0), reverse=True)


@dataclass
class Lane:
    name: str
    sizes: list[int]
    uncut: bool          # no enzyme of the candidate cuts it: a supercoiled plasmid or the
                         # whole linear molecule
    bands: list[gel.Band]
    visible: list[gel.Band]


def _lane(mol: Molecule, enzymes: tuple[str, ...], percent: float) -> Lane:
    sizes = fragment_sizes(mol, enzymes)
    uncut = not any(mol.cuts.get(e) for e in enzymes)
    bands = gel.lane_bands(sizes, percent, supercoiled=uncut and mol.circular)
    _, high = gel.resolving_range(percent)
    visible = [b for b in bands if b.size >= gel.MIN_VISIBLE_BP and b.apparent_size <= high]
    return Lane(mol.name, sizes, uncut, bands, visible)


def separation(a: Lane, b: Lane) -> float:
    """How far apart the two lanes are: the largest distance from a visible band of one lane
    to the nearest band of the other (0 when every visible band has a partner)."""
    best = 0.0
    for one, other in ((a, b), (b, a)):
        for band in one.visible:
            nearest = min((abs(band.rf - o.rf) for o in other.bands), default=1.0)
            best = max(best, nearest)
    return best


def _kb(size: int) -> str:
    # A no-break space keeps a size whole when the sentence wraps.
    return f"{size:,}\u00a0bp" if size < 1000 else f"{size / 1000:.1f}\u00a0kb"


def lane_text(lane: Lane) -> str:
    """'3.2 kb + 1.1 kb', '4.3 kb uncut', with the bands too small to see named apart."""
    if lane.uncut:
        return f"{_kb(lane.sizes[0])} uncut"
    shown = [s for s in lane.sizes if s >= gel.MIN_VISIBLE_BP]
    hidden = [s for s in lane.sizes if s < gel.MIN_VISIBLE_BP]
    text = " + ".join(_kb(s) for s in shown) or "no visible band"
    if hidden:
        text += f" (and {' + '.join(_kb(s) for s in hidden)}, too small to see)"
    return text


@dataclass
class Candidate:
    enzymes: tuple[str, ...]
    lanes: list[Lane]
    separations: list[float]      # per wrong clone

    @property
    def diagnostic(self) -> bool:
        return all(s >= gel.RESOLUTION for s in self.separations)

    @property
    def min_separation(self) -> float:
        return min(self.separations)

    @property
    def clear(self) -> bool:
        return (self.min_separation >= CLEAR_FACTOR * gel.RESOLUTION
                and not any(lane.uncut for lane in self.lanes))

    def rank(self) -> tuple:
        expected_bands = len(self.lanes[0].visible)
        return (not self.diagnostic, not self.clear, len(self.enzymes),
                not 2 <= expected_bands <= 5, -round(self.min_separation, 4),
                sum(e not in PREFERRED_NAMES for e in self.enzymes),
                -sum(_suppliers(e) for e in self.enzymes),
                tuple(e.lower() for e in self.enzymes))

    def pattern(self) -> tuple:
        """The fragments of every lane: candidates with the same pattern are one result."""
        return tuple(tuple(lane.sizes) + (lane.uncut,) for lane in self.lanes)

    def sentence(self) -> str:
        return "; ".join(f"{lane.name}: {lane_text(lane)}" for lane in self.lanes)

    def as_dict(self, same: list[Candidate] | None = None) -> dict[str, Any]:
        return {
            "enzymes": list(self.enzymes),
            "same_pattern": [list(c.enzymes) for c in same or []],
            "sentence": self.sentence(),
            "clear": self.clear,
            "min_separation": round(self.min_separation, 4),
            "lanes": [{"name": lane.name, "sizes": lane.sizes, "uncut": lane.uncut,
                       "bands": [b.as_dict() for b in lane.bands],
                       "text": lane_text(lane)} for lane in self.lanes],
        }


def _usable(mols: list[Molecule], name: str) -> bool:
    expected = mols[0]
    if len(expected.cuts.get(name, [])) > MAX_CUTS:
        return False
    return any(m.cuts.get(name) for m in mols)


def candidates(mols: list[Molecule], names: list[str], percent: float) -> list[Candidate]:
    usable = [n for n in names if _usable(mols, n)]
    combos: list[tuple[str, ...]] = [(n,) for n in usable]
    for a, b in combinations(usable, 2):
        cuts = len(mols[0].cuts.get(a, [])) + len(mols[0].cuts.get(b, []))
        if cuts <= MAX_CUTS + 2:
            combos.append((a, b))
    out = []
    for combo in combos:
        lanes = [_lane(m, combo, percent) for m in mols]
        seps = [separation(lanes[0], other) for other in lanes[1:]]
        out.append(Candidate(combo, lanes, seps))
    return out


def _why_not(mols: list[Molecule], found: list[Candidate], set_name: str,
             percent: float) -> str:
    """Why no digest in the set works, for the wrong clone that no candidate tells apart."""
    expected = mols[0]
    for i, other in enumerate(mols[1:]):
        if any(c.separations[i] >= gel.RESOLUTION for c in found):
            continue
        same = _same_molecule(expected, other)
        if same:
            return (f"“{other.name}” is the same molecule as the expected clone, so no digest "
                    "can tell them apart. Remove it from the list.")
        diff = abs(expected.length - other.length)
        changes = _point_changes(expected, other)
        if changes is not None:
            return (f"“{other.name}” differs from the expected clone by only {changes} base"
                    f"{'' if changes == 1 else 's'}: a digest cannot show that on a gel. Use "
                    "colony PCR or sequencing instead.")
        if diff == 0:
            return (f"“{other.name}” has the same length as the expected clone: only an enzyme "
                    f"that cuts inside the insert, off its centre, shows the difference, and "
                    f"none in {set_name} does. Try all commercial enzymes.")
        if diff < 60:
            return (f"“{other.name}” differs in length from the expected clone by only "
                    f"{diff} bp, too little for a {percent:.1f} % gel. Use colony PCR or "
                    "sequencing instead.")
        return (f"No enzyme in {set_name} gives “{other.name}” a band that runs apart from the "
                f"expected clone on a {percent:.1f} % gel. Try all commercial enzymes or "
                "another agarose percentage.")
    return f"No digest with {set_name} tells every wrong clone apart."


def _rotations(mol: Molecule) -> str:
    return mol.sequence + mol.sequence if mol.circular else mol.sequence


def _same_molecule(a: Molecule, b: Molecule) -> bool:
    if a.length != b.length or a.circular != b.circular:
        return False
    rc = str(Seq(b.sequence).reverse_complement())
    if a.circular:
        doubled = _rotations(a)
        return b.sequence in doubled or rc in doubled
    return a.sequence in (b.sequence, rc)


def _point_changes(a: Molecule, b: Molecule, limit: int = 30) -> int | None:
    """Number of substituted bases when two molecules of the same length line up base for
    base (as written, or one reverse-complemented); None when they are not that close."""
    if a.length != b.length:
        return None
    best = None
    for seq in (b.sequence, str(Seq(b.sequence).reverse_complement())):
        count = sum(1 for x, y in zip(a.sequence, seq) if x != y)
        if count <= limit and (best is None or count < best):
            best = count
    return best


def plan(expected: tuple[str, dict], alternatives: list[tuple[str, dict]],
         enzymes: list[str] | None, excluded: list[str], percent: float,
         ladder: str) -> dict[str, Any]:
    """Rank the diagnostic digests. ``enzymes`` None means the whole commercial set."""
    if ladder not in gel.LADDERS:
        raise ValueError(gel.unknown_ladder(ladder))
    if not alternatives:
        raise CloningError("Add at least one wrong clone to tell apart from the expected clone.")
    mols = [molecule(*expected)] + [molecule(name, data) for name, data in alternatives]
    known = {n.lower(): n for n in commercial_enzymes()}
    set_name = "your enzymes" if enzymes is not None else "the commercial set"
    chosen = known.values() if enzymes is None else [known[n.lower()] for n in enzymes
                                                     if n.lower() in known]
    left_out = {n.lower() for n in excluded}
    names = sorted({n for n in chosen if n.lower() not in left_out}, key=str.lower)
    if not names:
        raise CloningError("No enzyme left to try: every enzyme of the set is blocked by "
                           "methylation on these plasmids. Tick “DNA from a dam−/dcm− strain” "
                           "if yours is, or try all commercial enzymes.")
    _cut_positions(mols, names)
    found = candidates(mols, names, percent)
    ranked = sorted((c for c in found if c.diagnostic), key=Candidate.rank)
    # Digests that give the very same bands are one result, listed with the best enzymes.
    groups: dict[tuple, list[Candidate]] = {}
    for c in ranked:
        groups.setdefault(c.pattern(), []).append(c)
    best = list(groups.values())
    return {
        "results": [group[0].as_dict(group[1:]) for group in best[:TOP]],
        "diagnostic_count": len(best),
        "tried": len(found),
        "enzymes_tried": len(names),
        "reason": None if ranked else _why_not(mols, found, set_name, percent),
        "ladder": gel.ladder_bands(ladder, percent),
        "ladder_label": gel.LADDERS[ladder]["label"],
        "percent": percent,
    }
