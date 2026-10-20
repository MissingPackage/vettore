"""Primer design with primer3. Routes under /api/primer-design.

Picks PCR primer pairs that flank a target region, so that the product contains the whole
region, with primer3 (primer3-py ``design_primers``; Untergasser et al. 2012, Nucleic Acids
Res 40:e115). The best pairs are returned ranked by primer3's penalty, those without a
warning (below) first.

Conditions. primer3 is run in the conditions of Vettore's Standard Tm preset
(``tm.PRESETS["standard"]``): SantaLucia nearest-neighbour Tm and salt correction, 50 mM Na⁺,
no Mg²⁺ or dNTPs, 250 nM primer. primer3's Tm formula divides the oligo concentration by four
(two complementary strands at equal concentration, SantaLucia 1998 PNAS 95:1460), whereas the
preset, like Biopython's ``Tm_NN`` with ``dnac2=0``, uses the primer concentration itself
(primer in excess over template); passing 4 x 250 nM to primer3 makes the two agree within
about 0.1 °C, which ``tests/test_primer_design.py`` checks. The Tm shown is the preset's, so it
is the same number the status bar shows for that primer.

Defaults for a standard PCR: length 18–24 nt (optimum 20), GC 40–60 %, pair within 5 °C of
each other (Addgene, "How to Design a Primer", addgene.org/protocols/primer-design); Tm
57–63 °C, optimum 60 °C (primer3's own defaults, ``libprimer3.c``
``pr_set_default_global_args_1``).

Secondary structures. Self-dimer, hairpin and hetero-dimer ΔG come from primer3's
thermodynamic alignment (``calc_homodimer``, ``calc_hairpin``, ``calc_heterodimer``) at 37 °C
in the same conditions; a structure more stable than −9.0 kcal/mol is flagged (IDT, "Designing
PCR primers and probes": "the ΔG value of any self-dimers, hairpins, and heterodimers should be
weaker (more positive) than –9.0 kcal/mol"). primer3 already rejects primers whose structures
melt above 47 °C (its ``PRIMER_MAX_*_TH`` defaults), so most returned pairs carry no flag.

Circular templates are rotated so that the target sits in the middle; primers and product
may then cross the origin. Positions are 0-based and inclusive, like the editor's.

Packaging. primer3-py loads its thermodynamic parameter files from
``primer3/src/libprimer3/primer3_config/`` next to its compiled module, and swallows the
error at ``import primer3`` if they are missing. This module imports ``primer3.bindings``
explicitly, which raises, and turns that into a message the user can act on.
"""

from __future__ import annotations

import importlib
import logging
import os
import re
from pathlib import Path
from typing import Any

from Bio.Seq import Seq
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from . import tm

router = APIRouter(prefix="/api/primer-design")
log = logging.getLogger(__name__)

NUM_RETURN = 5
# primer3 is asked for more candidates than are shown, so that pairs free of the warnings
# below can be listed first (primer3's own structure limits are looser, see the docstring).
NUM_CANDIDATES = 25
TM_PRESET = tm.PRESETS[tm.DEFAULT_PRESET]
# Standard preset conditions (see the module docstring for the factor of four).
CONDITIONS = {"mv_conc": 50.0, "dv_conc": 0.0, "dntp_conc": 0.0, "dna_conc": 250.0}
PRIMER3_DNA_CONC = 4 * CONDITIONS["dna_conc"]
DG_TEMPERATURE = 37.0          # °C, the usual reference temperature for structure ΔG
DG_WARNING = -9.0              # kcal/mol, IDT guideline (module docstring)
THREE_PRIME_CHECK = 15         # 3' bases searched elsewhere on the template
MAX_PRIMER_LENGTH = 36         # primer3's built-in limit (MAX_PRIMER_LENGTH in libprimer3.c)
PARAMETER_FILES = ("stack.ds", "stack.dh", "loops.ds", "dangle.ds", "tstack2.ds")


class PrimerDesignError(ValueError):
    """A user-facing error: the message is shown in the dialog as is."""


class DesignRequest(BaseModel):
    template: dict[str, Any]
    start: int                                 # target, 0-based inclusive
    end: int
    product_min: int | None = None             # default: the smallest product that fits
    product_max: int | None = None             # default: product_min + 300
    size_min: int = Field(18, ge=8, le=MAX_PRIMER_LENGTH)
    size_opt: int = Field(20, ge=8, le=MAX_PRIMER_LENGTH)
    size_max: int = Field(24, ge=8, le=MAX_PRIMER_LENGTH)
    tm_min: float = 57.0
    tm_opt: float = 60.0
    tm_max: float = 63.0
    gc_min: float = Field(40.0, ge=0, le=100)
    gc_max: float = Field(60.0, ge=0, le=100)
    max_tm_difference: float = Field(5.0, ge=0)


# ---- primer3 and its parameter files ---------------------------------------------------

def parameters_dir() -> Path:
    """Where primer3-py reads its thermodynamic parameters (see ``thermoanalysis.pyx``)."""
    import primer3

    return Path(primer3.LOCAL_DIR) / "src" / "libprimer3" / "primer3_config"


def _unavailable_message(directory: Path) -> str:
    message = ("Primer design cannot run: primer3's thermodynamic data files are missing or "
               "unreadable in this copy of Vettore. Unzip a fresh copy of Vettore and try again.")
    if os.name == "nt" and not str(directory).isascii():
        # primer3 opens the files with the C library, which on Windows cannot open a path
        # with characters outside the system code page.
        message += (" Vettore is in a folder whose name has accented or non-Latin characters; "
                    "move the Vettore folder to a plain path such as C:\\Vettore.")
    return message


def primer3_bindings():
    """``primer3.bindings``, or a PrimerDesignError explaining why primer3 cannot run."""
    try:
        # Imported by name: when the files are missing ``import primer3`` has already
        # swallowed the error, and only importing the submodule again raises it.
        bindings = importlib.import_module("primer3.bindings")
    except (ImportError, OSError) as exc:
        try:
            directory = parameters_dir()
        except ImportError:
            directory = Path("primer3")
        log.error("primer3 unavailable (parameters expected in %s): %s", directory, exc)
        raise PrimerDesignError(_unavailable_message(directory)) from exc
    return bindings


def status() -> dict[str, Any]:
    """Whether primer3 and its parameter files can be used (also checked by packaging)."""
    try:
        directory = parameters_dir()
    except ImportError as exc:
        return {"ready": False, "message": f"primer3 is not installed: {exc}", "parameters": ""}
    missing = [name for name in PARAMETER_FILES if not (directory / name).is_file()]
    try:
        primer3_bindings()
    except PrimerDesignError as exc:
        return {"ready": False, "message": str(exc), "parameters": str(directory),
                "missing": missing}
    return {"ready": not missing, "message": "" if not missing else
            f"Missing primer3 parameter files: {', '.join(missing)}",
            "parameters": str(directory), "missing": missing}


# ---- helpers -----------------------------------------------------------------------------

def rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def _clean_template(sequence: str, name: str) -> str:
    seq = (sequence or "").upper().replace("U", "T")
    bad = sorted({c for c in seq if c not in "ACGTRYSWKMBDHVN"})
    if bad:
        raise PrimerDesignError(f"“{name}” contains characters that are not DNA "
                                f"({', '.join(repr(c) for c in bad[:5])}).")
    # primer3 avoids primers over N; other ambiguity codes are treated the same way.
    return re.sub(r"[^ACGT]", "N", seq)


def _gc(seq: str) -> float:
    return 100.0 * sum(1 for c in seq if c in "GC") / max(1, len(seq))


def _occurrences(template: str, circular: bool, word: str) -> int:
    """Exact matches of ``word`` on both strands of the template."""
    text = template + template[:len(word) - 1] if circular else template
    count = 0
    for probe in {word, rc(word)}:
        start = text.find(probe)
        while start != -1:
            count += 1
            start = text.find(probe, start + 1)
    return count


def _structures(bindings, forward: str, reverse: str) -> dict[str, float]:
    """ΔG (kcal/mol, 37 °C) of the most stable structure; 0 when none forms."""
    kwargs = {**CONDITIONS, "temp_c": DG_TEMPERATURE}

    def dg(result) -> float:
        return round(result.dg / 1000.0, 2) if result.structure_found else 0.0

    return {
        "forward_hairpin": dg(bindings.calc_hairpin(forward, **kwargs)),
        "forward_self_dimer": dg(bindings.calc_homodimer(forward, **kwargs)),
        "reverse_hairpin": dg(bindings.calc_hairpin(reverse, **kwargs)),
        "reverse_self_dimer": dg(bindings.calc_homodimer(reverse, **kwargs)),
        "hetero_dimer": dg(bindings.calc_heterodimer(forward, reverse, **kwargs)),
    }


STRUCTURE_LABELS = {
    "forward_hairpin": "The forward primer can fold into a hairpin",
    "forward_self_dimer": "The forward primer can pair with itself (self-dimer)",
    "reverse_hairpin": "The reverse primer can fold into a hairpin",
    "reverse_self_dimer": "The reverse primer can pair with itself (self-dimer)",
    "hetero_dimer": "The two primers can pair with each other (primer-dimer)",
}

def pair_warnings(structures: dict[str, float], forward: str, reverse: str, template: str,
                  circular: bool, name: str) -> list[str]:
    """Structures more stable than DG_WARNING, and 3' ends that also match elsewhere."""
    warnings = [f"{STRUCTURE_LABELS[key]} (ΔG {value:.1f} kcal/mol at 37 °C; IDT advises "
                f"keeping it above {DG_WARNING:.0f})"
                for key, value in structures.items() if value < DG_WARNING]
    for label, primer in (("forward", forward), ("reverse", reverse)):
        tail = primer[-THREE_PRIME_CHECK:]
        extra = _occurrences(template, circular, tail) - 1
        if extra > 0:
            warnings.append(f"The last {len(tail)} bases of the {label} primer also match "
                            f"{extra} other place{'s' if extra > 1 else ''} in “{name}”: "
                            "it may prime there too.")
    return warnings


EXPLAIN_WORDS = {
    "low tm": "Tm too low",
    "high tm": "Tm too high",
    "GC content failed": "GC outside the range",
    "GC clamp failed": "no G/C at the 3′ end",
    "high any compl": "pairs with itself",
    "high end compl": "3′ end pairs with itself",
    "high hairpin stability": "forms a hairpin",
    "long poly-x seq": "run of one base too long",
    "high end stability": "3′ end too stable",
    "contains N": "covers unknown bases (N)",
    "unacceptable product size": "product size outside the range",
    "tm diff too large": "Tm difference too large",
    "high any compl TH": "the primers pair with each other",
    "high end compl TH": "the 3′ ends pair with each other",
    "not in any ok left region": "left primer outside the allowed region",
    "not in any ok right region": "right primer outside the allowed region",
}


def _explain(text: str) -> str:
    """primer3's 'considered 1285, low tm 363, ..., ok 0' in plain words."""
    parts = []
    for item in (text or "").split(","):
        item = item.strip()
        match = re.match(r"(.+?) (\d+)$", item)
        if not match or match.group(1) in {"considered", "ok"}:
            continue
        label, count = match.groups()
        parts.append(f"{EXPLAIN_WORDS.get(label, label)} ({count})")
    return ", ".join(parts)


def _no_pairs_message(result: dict[str, Any]) -> str:
    reasons = []
    for key, side in (("PRIMER_LEFT_EXPLAIN", "forward candidates rejected"),
                      ("PRIMER_RIGHT_EXPLAIN", "reverse candidates rejected"),
                      ("PRIMER_PAIR_EXPLAIN", "pairs rejected")):
        if explained := _explain(result.get(key, "")):
            reasons.append(f"{side}: {explained}")
    detail = f" ({'; '.join(reasons)})" if reasons else ""
    return ("primer3 found no primer pair with these settings" + detail + ". Widen the Tm, "
            "length or GC range, or allow a larger product, then design again.")


# ---- design ------------------------------------------------------------------------------

def design(body: DesignRequest) -> dict[str, Any]:
    data = body.template
    name = data.get("name") or "Template"
    template = _clean_template(data.get("sequence") or "", name)
    circular = bool(data.get("circular"))
    length = len(template)
    if length < 2 * body.size_min + 1:
        raise PrimerDesignError(f"“{name}” is too short to design PCR primers on.")
    if not (0 <= body.start < length and 0 <= body.end < length):
        raise PrimerDesignError(f"The region must lie between 1 and {length:,}.")
    if body.start > body.end and not circular:
        raise PrimerDesignError("The region's start is after its end. On a linear sequence "
                                "the start must come first.")
    if not body.size_min <= body.size_opt <= body.size_max:
        raise PrimerDesignError("Primer length: the minimum must not exceed the optimum, nor "
                                "the optimum the maximum.")
    if not body.tm_min <= body.tm_opt <= body.tm_max:
        raise PrimerDesignError("Tm: the minimum must not exceed the optimum, nor the optimum "
                                "the maximum.")
    if body.gc_min > body.gc_max:
        raise PrimerDesignError("GC: the minimum must not exceed the maximum.")

    target_length = (body.end - body.start) % length + 1
    smallest = target_length + 2 * body.size_min
    product_min = body.product_min if body.product_min is not None else smallest
    product_max = body.product_max if body.product_max is not None else product_min + 300
    if product_max < smallest:
        raise PrimerDesignError(f"The product must be at least {smallest:,} bp to contain the "
                                f"{target_length:,} bp region and both primers. Raise the "
                                "maximum product size.")
    if product_min > product_max:
        raise PrimerDesignError("Product size: the minimum is larger than the maximum.")
    if product_min > length:
        raise PrimerDesignError(f"The smallest product allowed ({product_min:,} bp) is longer "
                                f"than “{name}” ({length:,} bp).")

    # A circular template is rotated so that the target sits in the middle: primers can then
    # be placed on either side of it, across the origin if need be.
    rotation = 0
    if circular:
        rotation = (body.start + target_length // 2 - length // 2) % length
    rotated = template[rotation:] + template[:rotation]
    target_start = (body.start - rotation) % length
    if target_start < body.size_min or target_start + target_length + body.size_min > length:
        where = ("the sequence is too short around it" if circular
                 else "it is too close to the start or the end of this linear sequence")
        raise PrimerDesignError(f"Primers need at least {body.size_min} bases on each side of "
                                f"the region, but {where}. Choose a smaller region.")

    bindings = primer3_bindings()
    global_args = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_NUM_RETURN": NUM_CANDIDATES,
        "PRIMER_MIN_SIZE": body.size_min,
        "PRIMER_OPT_SIZE": body.size_opt,
        "PRIMER_MAX_SIZE": body.size_max,
        "PRIMER_MIN_TM": body.tm_min,
        "PRIMER_OPT_TM": body.tm_opt,
        "PRIMER_MAX_TM": body.tm_max,
        "PRIMER_MIN_GC": body.gc_min,
        "PRIMER_MAX_GC": body.gc_max,
        "PRIMER_PAIR_MAX_DIFF_TM": body.max_tm_difference,
        "PRIMER_PRODUCT_SIZE_RANGE": [[product_min, min(product_max, length)]],
        "PRIMER_TM_FORMULA": 1,            # SantaLucia 1998
        "PRIMER_SALT_CORRECTIONS": 1,      # SantaLucia 1998
        "PRIMER_SALT_MONOVALENT": CONDITIONS["mv_conc"],
        "PRIMER_SALT_DIVALENT": CONDITIONS["dv_conc"],
        "PRIMER_DNTP_CONC": CONDITIONS["dntp_conc"],
        "PRIMER_DNA_CONC": PRIMER3_DNA_CONC,
        "PRIMER_EXPLAIN_FLAG": 1,
    }
    seq_args = {"SEQUENCE_ID": name, "SEQUENCE_TEMPLATE": rotated,
                "SEQUENCE_TARGET": [target_start, target_length]}
    try:
        result = bindings.design_primers(seq_args, global_args)
    except (OSError, ValueError) as exc:
        log.warning("primer3 rejected the request: %s", exc)
        raise PrimerDesignError(f"primer3 could not run with these settings: {exc}") from exc

    count = int(result.get("PRIMER_PAIR_NUM_RETURNED", 0))
    if not count:
        raise PrimerDesignError(_no_pairs_message(result))

    def place(start: int, size: int) -> tuple[int, int]:
        return (start + rotation) % length, (start + size - 1 + rotation) % length

    pairs = []
    for i in range(count):
        left_start, left_size = result[f"PRIMER_LEFT_{i}"]
        right_end, right_size = result[f"PRIMER_RIGHT_{i}"]     # primer3 gives the 5' end
        forward = result[f"PRIMER_LEFT_{i}_SEQUENCE"]
        reverse = result[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        f_start, f_end = place(left_start, left_size)
        r_start, r_end = place(right_end - right_size + 1, right_size)
        structures = _structures(bindings, forward, reverse)
        warnings = pair_warnings(structures, forward, reverse, template, circular, name)
        pairs.append({
            "primer3_rank": i + 1,
            "penalty": round(float(result[f"PRIMER_PAIR_{i}_PENALTY"]), 3),
            "product_size": int(result[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]),
            "product_start": f_start,
            "product_end": r_end,
            "forward": _primer(forward, f_start, f_end, result, f"PRIMER_LEFT_{i}"),
            "reverse": _primer(reverse, r_start, r_end, result, f"PRIMER_RIGHT_{i}"),
            "structures": structures,
            "warnings": warnings,
        })
    # Ranked by primer3's penalty, pairs without warnings first.
    pairs.sort(key=lambda p: (bool(p["warnings"]), p["penalty"]))
    pairs = pairs[:NUM_RETURN]
    for rank, pair in enumerate(pairs, 1):
        pair["rank"] = rank
    return {
        "pairs": pairs,
        "candidates": count,
        "template": {"name": name, "length": length, "circular": circular},
        "target": {"start": body.start, "end": body.end, "length": target_length},
        "product_range": [product_min, product_max],
        # e.g. "SantaLucia & Hicks 2004; Na⁺ 50 mM; primer 250 nM, Standard preset"
        "tm_method": f"{TM_PRESET.description}, {TM_PRESET.label.split(' (')[0]} preset",
        "dg_conditions": "ΔG at 37 °C, Na⁺ 50 mM (primer3 thermodynamic alignment)",
    }


def _primer(sequence: str, start: int, end: int, result: dict[str, Any], key: str) -> dict:
    return {
        "sequence": sequence,
        "start": start,
        "end": end,
        "length": len(sequence),
        "tm": round(tm.primer_tm(sequence), 1),
        "tm_primer3": round(float(result[f"{key}_TM"]), 1),
        "tm_presets": tm.all_presets(sequence),
        "gc": round(_gc(sequence), 1),
    }


# ---- HTTP --------------------------------------------------------------------------------

@router.get("/status")
def primer3_status() -> dict[str, Any]:
    return status()


@router.post("/design")
def design_route(body: DesignRequest) -> dict[str, Any]:
    try:
        return design(body)
    except PrimerDesignError as exc:
        raise HTTPException(422, str(exc)) from exc
