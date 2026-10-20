"""Primer melting temperature.

Nearest-neighbour Tm from Biopython (``Bio.SeqUtils.MeltingTemp.Tm_NN``) with explicit,
named conditions. The same primer differs by several degrees between tools, so every value is
reported together with the preset that produced it.

Presets:

- ``standard``: SantaLucia & Hicks 2004 nearest-neighbour parameters (Annu Rev Biophys
  33:415), 50 mM Na+, 250 nM primer, SantaLucia 1998 entropy salt correction — the textbook
  conditions for a PCR primer.
- ``neb_q5``: tracks NEB's Tm calculator for Q5 / high-fidelity polymerases within ~1.5 °C on
  a 10-primer benchmark against NEB's own API (Owczarzy 2008 with 2 mM Mg2+, 500 nM primer).
- ``neb_taq``: tracks NEB's calculator for Taq within ~1.2 °C (Owczarzy 2004, 50 mM Na+).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from Bio.SeqUtils import MeltingTemp as mt

MIN_LENGTH = 2      # the nearest-neighbour model needs at least one pair of bases
MAX_LENGTH = 100


@dataclass(frozen=True)
class TmPreset:
    key: str
    label: str
    description: str
    params: dict[str, Any]


PRESETS: dict[str, TmPreset] = {
    p.key: p
    for p in (
        TmPreset(
            "standard",
            "Standard (SantaLucia & Hicks 2004)",
            "SantaLucia & Hicks 2004; Na⁺ 50 mM; primer 250 nM",
            {"nn_table": mt.DNA_NN4, "dnac1": 250, "dnac2": 0, "Na": 50, "saltcorr": 5},
        ),
        TmPreset(
            "neb_q5",
            "NEB Q5 / high-fidelity",
            "SantaLucia 1998; K⁺ 50 mM, Tris 10 mM, Mg²⁺ 2 mM, dNTP 0.2 mM (Owczarzy 2008); "
            "primer 500 nM",
            {"nn_table": mt.DNA_NN3, "dnac1": 500, "dnac2": 0, "Na": 0, "K": 50, "Tris": 10,
             "Mg": 2, "dNTPs": 0.2, "saltcorr": 7},
        ),
        TmPreset(
            "neb_taq",
            "NEB Taq",
            "SantaLucia 1998; Na⁺ 50 mM (Owczarzy 2004); primer 500 nM",
            {"nn_table": mt.DNA_NN3, "dnac1": 500, "dnac2": 0, "Na": 50, "saltcorr": 6},
        ),
    )
}
DEFAULT_PRESET = "standard"


def primer_tm(sequence: str, preset: str = DEFAULT_PRESET) -> float:
    seq = sequence.upper().replace("U", "T")
    if not seq or any(c not in "ACGT" for c in seq):
        raise ValueError("Tm needs an unambiguous DNA sequence (A, C, G, T)")
    if len(seq) < MIN_LENGTH:
        raise ValueError(f"A Tm needs at least {MIN_LENGTH} bases")
    if len(seq) > MAX_LENGTH:
        raise ValueError(f"Tm is computed for oligos up to {MAX_LENGTH} nt")
    if preset not in PRESETS:
        raise ValueError(f"Unknown Tm preset: {preset}")
    return float(mt.Tm_NN(seq, **PRESETS[preset].params))


def all_presets(sequence: str) -> list[dict[str, Any]]:
    return [
        {"key": p.key, "label": p.label, "description": p.description,
         "tm": primer_tm(sequence, p.key)}
        for p in PRESETS.values()
    ]
