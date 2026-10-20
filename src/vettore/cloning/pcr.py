"""PCR and primer design.

Primer Tm is always the Tm of the **annealing part** (5' tails excluded), computed with
Vettore's default preset (:mod:`vettore.tm`), the convention most bench tools follow.
"""

from __future__ import annotations

from Bio.Seq import Seq
from fastapi import APIRouter
from pydantic import BaseModel
from pydna.amplify import Anneal
from pydna.design import primer_design
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer

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
from .regions import annotate, check_region, features_in_region, region_length, region_sequence

router = APIRouter()

TM_DIFFERENCE_WARNING = 5.0
MIN_PRIMER_LENGTH = 12
TARGET_TM_RANGE = (40.0, 80.0)   # °C: annealing temperatures a PCR is run at


class PrimerIn(BaseModel):
    name: str = "Primer"
    sequence: str


def tm_of(sequence: str) -> float | None:
    try:
        return primer_tm(str(sequence))
    except ValueError:
        return None


def clean_bases(sequence: str, what: str) -> str:
    seq = "".join(sequence.split()).upper().replace("U", "T")
    if not seq:
        raise CloningError(f"{what} is empty")
    if any(c not in "ACGTRYSWKMBDHVN" for c in seq):
        raise CloningError(f"{what} contains characters that are not DNA")
    return seq


def anneal(template: Dseqrecord, forward: PrimerIn, reverse: PrimerIn, limit: int = 13):
    f = Primer(clean_bases(forward.sequence, f"Forward primer “{forward.name}”"),
               name=forward.name)
    r = Primer(clean_bases(reverse.sequence, f"Reverse primer “{reverse.name}”"),
               name=reverse.name)
    for primer, label in ((f, "Forward"), (r, "Reverse")):
        if len(primer) < MIN_PRIMER_LENGTH:
            raise CloningError(f"{label} primer “{primer.name}” is only {len(primer)} nt; "
                               f"primers must be at least {MIN_PRIMER_LENGTH} nt")
    annealed = Anneal([f, r], template, limit=limit)
    if not annealed.forward_primers or not annealed.reverse_primers:
        raise CloningError(_binding_problem(annealed, f, r, template.name))
    return annealed


def _matches(primers, primer: Primer | PrimerIn) -> bool:
    """Is that primer (by sequence) among the primers pydna found on one strand?"""
    bases = "".join(str(primer.sequence if isinstance(primer, PrimerIn) else primer.seq).split())
    bases = bases.upper().replace("U", "T")
    return any(str(p.seq).upper() == bases for p in primers)


def _binding_problem(annealed, forward: Primer, reverse: Primer, template: str) -> str:
    """Why no product forms, with the fix: a primer written on the wrong strand (the commonest
    slip, a reverse primer copied from the top strand) is named with its reverse complement."""
    top, bottom = annealed.forward_primers, annealed.reverse_primers
    f_top, f_bottom = _matches(top, forward), _matches(bottom, forward)
    r_top, r_bottom = _matches(top, reverse), _matches(bottom, reverse)
    check = "Check that it is written 5′→3′ and that the right template is selected."
    if r_top and not r_bottom:
        return (f"Reverse primer “{reverse.name}” matches the top strand of “{template}”, so it "
                "points the same way as the forward primer. A reverse primer is written on the "
                "other strand: use its reverse complement, "
                f"5′-{reverse.seq.reverse_complement()}-3′.")
    if f_bottom and not f_top:
        return (f"Forward primer “{forward.name}” matches the bottom strand of “{template}”, so "
                "it points the same way as the reverse primer. A forward primer is written on "
                "the top strand: use its reverse complement, "
                f"5′-{forward.seq.reverse_complement()}-3′.")
    if not f_top:
        return f"Forward primer “{forward.name}” does not bind “{template}”. {check}"
    return f"Reverse primer “{reverse.name}” does not bind “{template}”. {check}"


def amplified_region(annealed) -> tuple[int, int, int, int]:
    """(template start, template end inclusive, forward tail length, reverse tail length)."""
    fwd = annealed.forward_primers[0]
    rev = annealed.reverse_primers[0]
    n = len(annealed.template)
    f_footprint = len(fwd.footprint)
    r_footprint = len(rev.footprint)
    start = (fwd.position - f_footprint) % n
    end = (rev.position + r_footprint - 1) % n
    return start, end, len(fwd) - f_footprint, len(rev) - r_footprint


def run_pcr(template_doc: SequenceDoc, forward: PrimerIn, reverse: PrimerIn,
            mismatches: int = 0) -> CloningResult:
    if mismatches:
        raise CloningError("Primers with mismatches are not supported yet; "
                           "use primers that match the template exactly")
    template = to_dseqrecord(template_doc)
    annealed = anneal(template, forward, reverse)
    swapped: list[str] = []
    if (_matches(annealed.forward_primers, reverse)
            and not _matches(annealed.forward_primers, forward)):
        # Both bind, each on the other's strand: the product is the same, and it is drawn with
        # each primer where it really binds.
        forward, reverse = reverse, forward
        swapped.append(f"“{forward.name}” binds as the forward primer and “{reverse.name}” as "
                       "the reverse primer, the other way round from how they were entered. "
                       "The product is the same; keep it in mind when you name and order them.")
    start, end, f_tail, r_tail = amplified_region(annealed)
    if not template.circular and start > end:
        raise CloningError(f"The primers face away from each other on the linear template "
                           f"“{template_doc.name}”; no product is formed")

    fwd_seq = clean_bases(forward.sequence, "Forward primer")
    rev_seq = clean_bases(reverse.sequence, "Reverse primer")
    core = region_sequence(str(template.seq).upper(), start, end)
    sequence = fwd_seq[:f_tail] + core + str(Seq(rev_seq[:r_tail]).reverse_complement())

    product = Dseqrecord(sequence, circular=False)
    product.features.extend(features_in_region(template, start, end, f_tail))
    annotate(product, 0, len(fwd_seq) - 1, forward.name, "primer_bind", 1)
    annotate(product, len(sequence) - len(rev_seq), len(sequence) - 1, reverse.name,
             "primer_bind", -1)

    f_annealing = str(annealed.forward_primers[0].footprint)
    r_annealing = str(annealed.reverse_primers[0].footprint)
    oligos = [
        Oligo(name=forward.name, sequence=fwd_seq, tm=tm_of(f_annealing),
              annealing_length=len(f_annealing)),
        Oligo(name=reverse.name, sequence=rev_seq, tm=tm_of(r_annealing),
              annealing_length=len(r_annealing)),
    ]
    warnings = swapped + binding_warnings(annealed, oligos)

    name = f"{template_doc.name} PCR product"
    details = (f"Amplified {region_length(len(template), start, end)} bp of "
               f"{template_doc.name} ({start + 1}..{end + 1}) with {forward.name} and "
               f"{reverse.name}")
    return CloningResult(
        products=[Product(doc=from_dseqrecord(product, name), operation="PCR",
                          parents=[template_doc.name], details=details, oligos=oligos)],
        warnings=warnings, oligos=oligos)


def binding_warnings(annealed, oligos: list[Oligo]) -> list[str]:
    warnings: list[str] = []
    for label, primers in (("Forward", annealed.forward_primers),
                           ("Reverse", annealed.reverse_primers)):
        if len(primers) > 1:
            positions = ", ".join(str(p.position) for p in primers)
            warnings.append(f"{label} primer “{primers[0].name}” binds at {len(primers)} "
                            f"positions ({positions}); the first one is used")
    tms = [o.tm for o in oligos if o.tm is not None]
    if len(tms) == 2 and abs(tms[0] - tms[1]) > TM_DIFFERENCE_WARNING:
        warnings.append(f"The annealing Tm values differ by {abs(tms[0] - tms[1]):.1f} °C "
                        f"({tms[0]:.1f} / {tms[1]:.1f} °C); consider redesigning one primer")
    return warnings


def design_primers(template_doc: SequenceDoc, start: int, end: int, target_tm: float = 60.0,
                   forward_tail: str = "", reverse_tail: str = "") -> CloningResult:
    low, high = TARGET_TM_RANGE
    if not low <= target_tm <= high:
        raise CloningError(f"The target Tm must be between {low:.0f} and {high:.0f} °C; "
                           f"{target_tm:g} °C is not a temperature a primer anneals at")
    template = to_dseqrecord(template_doc)
    check_region(template, start, end, "template")
    core = region_sequence(str(template.seq).upper(), start, end)
    if len(core) < 2 * MIN_PRIMER_LENGTH:
        raise CloningError(f"The region is only {len(core)} bp; choose at least "
                           f"{2 * MIN_PRIMER_LENGTH} bp")
    f_tail = clean_bases(forward_tail, "Forward tail") if forward_tail else ""
    r_tail = clean_bases(reverse_tail, "Reverse tail") if reverse_tail else ""
    amplicon = design_amplicon(core, target_tm, f_tail, r_tail,
                               name=f"{template_doc.name}_{start + 1}-{end + 1}")

    product = Dseqrecord(amplicon.sequence, circular=False)
    product.features.extend(features_in_region(template, start, end, len(f_tail)))
    annotate(product, 0, len(amplicon.forward.sequence) - 1, amplicon.forward.name,
             "primer_bind", 1)
    annotate(product, len(amplicon.sequence) - len(amplicon.reverse.sequence),
             len(amplicon.sequence) - 1, amplicon.reverse.name, "primer_bind", -1)

    oligos = [amplicon.forward, amplicon.reverse]
    name = f"{template_doc.name} {start + 1}..{end + 1} PCR product"
    details = (f"Designed primers to amplify {len(core)} bp of {template_doc.name} "
               f"({start + 1}..{end + 1}) at about {target_tm:.0f} °C")
    return CloningResult(
        products=[Product(doc=from_dseqrecord(product, name), operation="PCR",
                          parents=[template_doc.name], details=details, oligos=oligos)],
        warnings=binding_warnings_for_tm(oligos), oligos=oligos)


def binding_warnings_for_tm(oligos: list[Oligo]) -> list[str]:
    tms = [o.tm for o in oligos if o.tm is not None]
    if len(tms) == 2 and abs(tms[0] - tms[1]) > TM_DIFFERENCE_WARNING:
        return [f"The annealing Tm values differ by {abs(tms[0] - tms[1]):.1f} °C"]
    return []


class DesignedAmplicon:
    def __init__(self, sequence: str, forward: Oligo, reverse: Oligo):
        self.sequence = sequence
        self.forward = forward
        self.reverse = reverse


def design_amplicon(core: str, target_tm: float, forward_tail: str = "",
                    reverse_tail: str = "", name: str = "product") -> DesignedAmplicon:
    """Primers whose annealing parts reach ``target_tm`` (Vettore's Tm preset), plus 5' tails."""
    try:
        designed = primer_design(Dseqrecord(core, circular=False), target_tm=target_tm,
                                 tm_func=lambda seq: primer_tm(str(seq)))
    except ValueError as exc:      # a primer grew past the longest oligo a Tm is computed for
        raise CloningError(f"No primer of up to 100 nt reaches {target_tm:g} °C on this region "
                           "(it may be too AT-rich, or contain N); choose a lower target Tm or "
                           "another region") from exc
    f_anneal = str(designed.forward_primer.seq).upper()
    r_anneal = str(designed.reverse_primer.seq).upper()
    return DesignedAmplicon(
        sequence=forward_tail + core + str(Seq(reverse_tail).reverse_complement()),
        forward=Oligo(name=f"{name}_F", sequence=forward_tail + f_anneal,
                      tm=tm_of(f_anneal), annealing_length=len(f_anneal)),
        reverse=Oligo(name=f"{name}_R", sequence=reverse_tail + r_anneal,
                      tm=tm_of(r_anneal), annealing_length=len(r_anneal)),
    )


class PcrRequest(BaseModel):
    template: SequenceDoc
    forward: PrimerIn
    reverse: PrimerIn
    mismatches: int = 0


class DesignRequest(BaseModel):
    template: SequenceDoc
    start: int
    end: int
    target_tm: float = 60.0
    forward_tail: str = ""
    reverse_tail: str = ""


@router.post("/pcr")
def pcr_route(body: PcrRequest) -> CloningResult:
    return run_pcr(body.template, body.forward, body.reverse, body.mismatches)


@router.post("/design-pcr-primers")
def design_route(body: DesignRequest) -> CloningResult:
    return design_primers(body.template, body.start, body.end, body.target_tm,
                          body.forward_tail, body.reverse_tail)
