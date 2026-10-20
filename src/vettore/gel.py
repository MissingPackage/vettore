"""Agarose gel simulation.

Migration uses the semi-logarithmic relation between mobility and fragment length
(Helling, Goodman & Boyer 1974, J Virol 14:1235), scaled to the resolving range of the chosen
agarose percentage. That is the model behind most quick gel tools: good enough to predict the
*pattern* of a diagnostic digest and to read sizes against a ladder, but not a physical
simulation — it ignores field strength, buffer, run time, temperature and base composition.
A fragment outside the resolving range is drawn compressed near the well or near the dye
front, as it would run on a real gel.

Ladder band sizes are taken from the vendors' own product documentation: NEB product pages for
N3232, N3200 and N3231, and the Thermo Fisher manuals MAN0013004 (SM0311), MAN0013047 (SM1331),
MAN0013010 (SM0321 family) and MAN0000898 (Invitrogen 10787018).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# bp, high to low. Bands printed in bold by the vendor are marked as reference bands.
LADDERS: dict[str, dict] = {
    "neb-1kb": {
        "label": "NEB 1 kb (N3232)",
        "bands": [10002, 8001, 6001, 5001, 4001, 3001, 2000, 1500, 1000, 517, 500],
        "reference": [3001],
    },
    "neb-1kb-plus": {
        "label": "NEB 1 kb Plus (N3200)",
        "bands": [10002, 8001, 6001, 5001, 4001, 3001, 2017, 1517, 1200, 1000, 900, 800, 700,
                  600, 517, 500, 400, 300, 200, 100],
        "reference": [3001, 1000, 500],
    },
    "neb-100bp": {
        "label": "NEB 100 bp (N3231)",
        "bands": [1517, 1200, 1000, 900, 800, 700, 600, 517, 500, 400, 300, 200, 100],
        "reference": [1000, 500],
    },
    "generuler-1kb": {
        "label": "Thermo GeneRuler 1 kb (SM0311)",
        "bands": [10000, 8000, 6000, 5000, 4000, 3500, 3000, 2500, 2000, 1500, 1000, 750, 500,
                  250],
        "reference": [6000, 3000, 1000],
    },
    "generuler-1kb-plus": {
        "label": "Thermo GeneRuler 1 kb Plus (SM1331)",
        "bands": [20000, 10000, 7000, 5000, 4000, 3000, 2000, 1500, 1000, 700, 500, 400, 300,
                  200, 75],
        "reference": [5000, 1500, 500],
    },
    "generuler-100bp-plus": {
        "label": "Thermo GeneRuler 100 bp Plus (SM0321)",
        "bands": [3000, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
        "reference": [1000, 500],
    },
    "invitrogen-1kb-plus": {
        "label": "Invitrogen 1 Kb Plus (10787018)",
        "bands": [15000, 10000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 850, 650,
                  500, 400, 300, 200, 100],
        "reference": [1500],
    },
}
DEFAULT_LADDER = "neb-1kb-plus"

# Resolving range (bp) per agarose percentage, from standard electrophoresis tables.
RESOLVING_RANGE: dict[float, tuple[int, int]] = {
    0.5: (1000, 30000),
    0.7: (800, 12000),
    0.8: (600, 10000),
    1.0: (400, 10000),
    1.2: (300, 7000),
    1.5: (200, 3000),
    2.0: (50, 2000),
    3.0: (20, 1000),
}
DEFAULT_PERCENT = 1.0

#: A supercoiled plasmid runs ahead of a linear molecule of the same length. The factor is an
#: approximation used only to place the "uncut" lane; it is labelled as apparent size.
SUPERCOILED_FACTOR = 0.65


def resolving_range(percent: float) -> tuple[int, int]:
    nearest = min(RESOLVING_RANGE, key=lambda p: abs(p - percent))
    return RESOLVING_RANGE[nearest]


def migration(size_bp: float, percent: float = DEFAULT_PERCENT) -> float:
    """Relative migration distance, 0 at the well and 1 at the dye front."""
    low, high = resolving_range(percent)
    if size_bp <= 0:
        return 1.0
    span = math.log10(high) - math.log10(low)
    rf = 0.08 + 0.84 * (math.log10(high) - math.log10(size_bp)) / span
    # Outside the resolving range bands pile up at the ends instead of leaving the gel.
    if rf < 0.08:
        rf = 0.08 - 0.06 * math.tanh((0.08 - rf) * 3)
    if rf > 0.92:
        rf = 0.92 + 0.06 * math.tanh((rf - 0.92) * 3)
    return min(max(rf, 0.01), 0.99)


@dataclass
class Band:
    size: int
    apparent_size: float
    rf: float
    intensity: float
    label: str

    def as_dict(self) -> dict:
        return {"size": self.size, "apparent_size": round(self.apparent_size), "rf": self.rf,
                "intensity": self.intensity, "label": self.label}


def _merge(bands: list[Band], resolution: float = 0.012) -> list[Band]:
    """Bands that would run together on a real gel are drawn as one, brighter band."""
    merged: list[Band] = []
    for band in sorted(bands, key=lambda b: b.rf):
        if merged and abs(band.rf - merged[-1].rf) < resolution:
            previous = merged[-1]
            total = previous.intensity + band.intensity
            previous.intensity = min(1.0, total)
            previous.label = (f"{previous.label} + {band.label}"
                              if band.label not in previous.label else previous.label)
        else:
            merged.append(band)
    return merged


def lane_bands(sizes: list[int], percent: float, supercoiled: bool = False) -> list[Band]:
    """Bands for a lane. Intensity follows mass, which for an equimolar digest is the
    fragment length (a 4 kb fragment stains eight times as strongly as a 500 bp one)."""
    if not sizes:
        return []
    largest = max(sizes)
    bands = []
    for size in sizes:
        apparent = size * SUPERCOILED_FACTOR if supercoiled else size
        label = f"{size:,} bp" + (" (supercoiled)" if supercoiled else "")
        bands.append(Band(size=size, apparent_size=apparent, rf=migration(apparent, percent),
                          intensity=round(max(0.18, min(1.0, size / largest)), 3), label=label))
    return _merge(bands)


def ladder_bands(ladder: str, percent: float) -> list[dict]:
    if ladder not in LADDERS:
        raise ValueError(unknown_ladder(ladder))
    entry = LADDERS[ladder]
    reference = set(entry["reference"])
    return [{"size": size, "rf": migration(size, percent),
             "intensity": 0.95 if size in reference else 0.55,
             "label": f"{size:,}", "reference": size in reference}
            for size in entry["bands"]]


def ladders() -> list[dict]:
    return [{"id": key, "label": value["label"], "bands": value["bands"]}
            for key, value in LADDERS.items()]


def unknown_ladder(ladder: str) -> str:
    return (f"The ladder “{ladder}” is not one Vettore knows. Choose a ladder from the list and "
            "simulate the gel again.")


def _lane(name: str, data: object, enzyme_names: list[str], percent: float) -> dict:
    """One lane: digested, uncut, or empty, with a note saying why when it is not what was
    asked. A non-cutting enzyme is a normal control on a diagnostic gel, so it never fails the
    whole gel."""
    from pydna.dseqrecord import Dseqrecord

    from .cloning.base import CloningError, dna_sequence
    from .cloning.restriction import cut_counts, enzyme_objects, fragments

    data = data if isinstance(data, dict) else {}
    label = str(data.get("name") or name or "Untitled")
    title = f"{name} / {' + '.join(enzyme_names)}" if enzyme_names else f"{name} / uncut"
    sequence = data.get("sequence")
    if not isinstance(sequence, str) or not sequence:
        return {"name": title, "bands": [], "sizes": [],
                "note": f"“{label}” has no sequence, so its lane is empty."}
    try:
        record = Dseqrecord(dna_sequence(label, sequence), circular=bool(data.get("circular")))
        enzymes = enzyme_objects(enzyme_names) if enzyme_names else []
    except CloningError as exc:
        return {"name": title, "bands": [], "sizes": [], "note": f"{exc} Its lane is empty."}
    record.name = label
    counts = cut_counts(record, enzymes)
    cutting = [e for e in enzymes if counts[str(e)]]
    idle = " and ".join(str(e) for e in enzymes if not counts[str(e)])
    verb = "does" if " and " not in idle else "do"
    note = None
    if cutting:
        sizes = [f.length for f in fragments(record, cutting)]
        supercoiled = False
        if idle:
            note = (f"{idle} {verb} not cut “{label}”; the lane shows the "
                    f"{' + '.join(map(str, cutting))} digest.")
    else:
        # Incubated with an enzyme that does not cut, a plasmid stays supercoiled.
        sizes = [len(record)]
        supercoiled = record.circular
        if idle:
            shape = "uncut plasmid" if record.circular else "uncut molecule"
            note = f"{idle} {verb} not cut “{label}”, so the lane shows the {shape}."
    bands = lane_bands(sizes, percent, supercoiled)
    return {"name": title, "bands": [b.as_dict() for b in bands],
            "sizes": sorted(sizes, reverse=True), "note": note}


def simulate(lanes: list[tuple[str, object, list[str]]], ladder: str, percent: float) -> dict:
    """The whole gel: (lane name, sequenceData, enzyme names) per lane, then the ladder."""
    if ladder not in LADDERS:
        raise ValueError(unknown_ladder(ladder))
    return {"lanes": [_lane(name, data, enzymes, percent) for name, data, enzymes in lanes],
            "ladder": ladder_bands(ladder, percent), "ladder_label": LADDERS[ladder]["label"],
            "percent": percent}
