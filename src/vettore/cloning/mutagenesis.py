"""Site-directed mutagenesis with back-to-back primers (NEB Q5 SDM style).

The whole plasmid is amplified with two primers that point away from each other at the edit
site; the linear product is then circularised by KLD (kinase, ligase, DpnI). The change is
carried on the primers' 5' ends:

- substitutions and insertions up to :data:`SINGLE_PRIMER_LIMIT` bases sit on the 5' end of the
  forward primer;
- longer changes are split between both primers, the 3' half on the forward primer and the
  reverse complement of the 5' half on the reverse primer (NEBaseChanger's rule);
- deletions need no tail at all.
"""

from __future__ import annotations

from typing import Literal

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from fastapi import APIRouter
from pydantic import BaseModel
from pydna.dseqrecord import Dseqrecord

from ..tm import primer_tm
from .base import (
    CloningError,
    CloningResult,
    Oligo,
    Product,
    SequenceDoc,
    from_dseqrecord,
    to_dseqrecord,
)
from .pcr import clean_bases
from .regions import region_length, region_sequence

router = APIRouter()

SINGLE_PRIMER_LIMIT = 6
MIN_ANNEALING = 16
MAX_ANNEALING = 36
DEFAULT_TARGET_TM = 60.0

Kind = Literal["substitution", "insertion", "deletion"]


def annealing_part(sequence: str, target_tm: float) -> str:
    """Shortest stretch from the 5' end of ``sequence`` whose Tm reaches ``target_tm``."""
    best = sequence[:MIN_ANNEALING]
    for length in range(MIN_ANNEALING, min(MAX_ANNEALING, len(sequence)) + 1):
        best = sequence[:length]
        try:
            if primer_tm(best) >= target_tm:
                return best
        except ValueError:
            continue
    return best


def _remap_features(record: Dseqrecord, start: int, removed: int, added: int) -> list[SeqFeature]:
    """Shift features across an edit; features spanning the edit keep their boundaries and
    change length, which is what a bench scientist expects after a small edit."""
    end = start + removed  # exclusive
    delta = added - removed

    def move(pos: int, is_end: bool) -> int | None:
        if pos <= start:
            return pos
        if pos >= end:
            return pos + delta
        return None if removed and not is_end else start + added

    out: list[SeqFeature] = []
    for feature in record.features:
        if feature.type == "source":
            continue
        parts = []
        for part in feature.location.parts:
            first = move(int(part.start), False)
            last = move(int(part.end), True)
            if first is None or last is None or last <= first:
                parts = []
                break
            parts.append(SimpleLocation(first, last, part.strand))
        if not parts:
            continue
        if feature.location.strand == -1:
            parts = parts[::-1]
        location = parts[0] if len(parts) == 1 else CompoundLocation(parts)
        out.append(SeqFeature(location, type=feature.type, qualifiers=dict(feature.qualifiers)))
    return out


def mutagenesis(doc: SequenceDoc, kind: Kind, start: int, end: int, new_sequence: str = "",
                target_tm: float = DEFAULT_TARGET_TM) -> CloningResult:
    record = to_dseqrecord(doc)
    sequence = str(record.seq).upper()
    n = len(sequence)
    if not record.circular:
        raise CloningError(f"“{doc.name}” is linear; this tool mutates circular plasmids")
    if not (0 <= start < n and 0 <= end < n):
        raise CloningError(f"The position must be inside 1..{n}")
    new_bases = clean_bases(new_sequence, "The new bases") if new_sequence else ""

    if kind == "insertion":
        end = start
        removed = 0
        old = ""
    else:
        if start > end:
            raise CloningError("The region must not span the origin; rotate the plasmid first")
        removed = region_length(n, start, end)
        old = region_sequence(sequence, start, end)
    if kind == "deletion":
        new_bases = ""
    elif not new_bases:
        raise CloningError("Enter the new bases")

    # the edit replaces sequence[prefix_end:suffix_start] with the new bases
    prefix_end = start + 1 if kind == "insertion" else start      # exclusive
    suffix_start = start + 1 if kind == "insertion" else end + 1

    product_sequence = sequence[:prefix_end] + new_bases + sequence[suffix_start:]
    if len(product_sequence) < 200:
        raise CloningError("The product would be too short for whole-plasmid PCR")

    doubled = sequence + sequence
    forward_anneal = annealing_part(doubled[suffix_start : suffix_start + MAX_ANNEALING],
                                    target_tm)
    upstream = doubled[n + prefix_end - MAX_ANNEALING : n + prefix_end]
    reverse_anneal = annealing_part(str(Seq(upstream).reverse_complement()), target_tm)

    if len(new_bases) <= SINGLE_PRIMER_LIMIT:
        f_tail, r_tail = new_bases, ""
    else:
        half = len(new_bases) // 2
        f_tail, r_tail = new_bases[half:], str(Seq(new_bases[:half]).reverse_complement())

    oligos = [
        Oligo(name=f"{doc.name}_SDM_F", sequence=f_tail + forward_anneal,
              tm=primer_tm(forward_anneal), annealing_length=len(forward_anneal)),
        Oligo(name=f"{doc.name}_SDM_R", sequence=r_tail + reverse_anneal,
              tm=primer_tm(reverse_anneal), annealing_length=len(reverse_anneal)),
    ]

    product = Dseqrecord(product_sequence, circular=True)
    product.features.extend(_remap_features(record, prefix_end, removed, len(new_bases)))

    if kind == "substitution":
        details = (f"Substituted {removed} bp at {start + 1}..{end + 1} "
                   f"({old} → {new_bases})")
        suffix = "mut"
    elif kind == "insertion":
        details = f"Inserted {len(new_bases)} bp after position {start + 1} ({new_bases})"
        suffix = "ins"
    else:
        shown = f" ({old})" if removed <= 20 else ""
        details = f"Deleted {removed} bp at {start + 1}..{end + 1}{shown}"
        suffix = "del"

    warnings = []
    if len(new_bases) > SINGLE_PRIMER_LIMIT:
        warnings.append(f"The {len(new_bases)} bp change is split between both primers' 5' ends")
    if abs(oligos[0].tm - oligos[1].tm) > 5:
        warnings.append("The two annealing Tm values differ by more than 5 °C")
    warnings.append("After PCR, treat with KLD (kinase, ligase, DpnI) to circularise the product "
                    "and remove the template")

    name = f"{doc.name} {suffix}"
    return CloningResult(
        products=[Product(doc=from_dseqrecord(product, name), operation="Mutagenesis",
                          parents=[doc.name], details=details.strip(), oligos=oligos)],
        warnings=warnings, oligos=oligos)


class MutagenesisRequest(BaseModel):
    doc: SequenceDoc
    kind: Kind
    start: int
    end: int
    new_sequence: str = ""
    target_tm: float = DEFAULT_TARGET_TM


@router.post("/mutagenesis")
def mutagenesis_route(body: MutagenesisRequest) -> CloningResult:
    return mutagenesis(body.doc, body.kind, body.start, body.end, body.new_sequence,
                       body.target_tm)
