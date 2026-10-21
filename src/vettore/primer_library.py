"""Lab primer library (S5.3): the primers a lab already keeps, so they are reused instead of
re-ordered. Routes under /api/primer-library.

The library is **one .csv file** the lab keeps in a shared folder (OneDrive, Dropbox, a lab
drive) so colleagues share it; Vettore reads it from the path the user chose and writes it back
there. The file's columns are ``name, sequence, notes, location, added`` (a header row, so a
file edited by hand or in Excel stays readable). Nothing is locked: the window re-reads the
file when its modification time changes, and if it changed under an edit the dialog offers to
reload.

Binding — *the one rule used everywhere* (here, on the map, and in primer design). A primer
binds a template where its 3′ end anneals: the last ``ANCHOR`` (15) nucleotides match the
template exactly, and the rest of the primer, read 5′ from that anchor, keeps matching with at
most ``MAX_MISMATCHES`` (2) mismatches. A 5′ tail is allowed: the bases beyond the longest
annealing stretch (anchor + up to two mismatches) are the tail and are not scored. Both strands
are searched; a circular template is searched across its origin. The annealing footprint is
reported as a 0-based inclusive span on the template, like the editor's own primers.
"""

from __future__ import annotations

import csv
import io
import logging
from dataclasses import dataclass
from datetime import date
from typing import Any

from Bio.Seq import Seq
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from . import settings, tm

router = APIRouter(prefix="/api/primer-library")
log = logging.getLogger(__name__)

#: The 3′ bases that must match the template exactly for a primer to prime there.
ANCHOR = 15
#: Mismatches allowed in the annealing stretch 5′ of the anchor; more starts the tail.
MAX_MISMATCHES = 2
#: Columns of the shared .csv, in order. ``name`` and ``sequence`` are required.
COLUMNS = ("name", "sequence", "notes", "location", "added")
#: How far upstream a sequencing primer reads from before the region it covers [ASSUMED].
SEQ_PRIMER_MIN, SEQ_PRIMER_MAX = 50, 100
_IUPAC = set("ACGTRYSWKMBDHVN")


class PrimerLibraryError(ValueError):
    """A user-facing error: the message is shown in the dialog as is."""


@dataclass(frozen=True)
class Primer:
    """One primer of the library. ``sequence`` is upper-case DNA, 5′→3′."""

    name: str
    sequence: str
    notes: str = ""
    location: str = ""
    added: str = ""

    def as_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "sequence": self.sequence,
            "notes": self.notes,
            "location": self.location,
            "added": self.added,
            "length": len(self.sequence),
            "tm": round(tm.primer_tm(self.sequence), 1) if _is_plain_dna(self.sequence) else None,
        }


def _is_plain_dna(sequence: str) -> bool:
    return bool(sequence) and all(c in "ACGT" for c in sequence)


def clean_sequence(raw: str) -> str:
    """A primer sequence in upper case, U read as T. Raises on anything but IUPAC DNA codes."""
    seq = (raw or "").strip().upper().replace("U", "T")
    seq = "".join(seq.split())  # spaces are common when a sequence is copied in blocks
    bad = sorted({c for c in seq if c not in _IUPAC})
    if not seq:
        raise PrimerLibraryError("the sequence is empty.")
    if bad:
        raise PrimerLibraryError(
            "it contains characters that are not DNA (" + ", ".join(repr(c) for c in bad[:5]) + ").")
    return seq


# ---- the shared .csv file ----------------------------------------------------------------

def _normalise_header(fields: list[str]) -> list[str] | None:
    """Map a file's header to our column names, or None when the first row is not a header."""
    lowered = [f.strip().lower() for f in fields]
    if "name" in lowered and ("sequence" in lowered or "seq" in lowered):
        return ["sequence" if f in ("seq",) else f for f in lowered]
    return None


def parse_csv(text: str) -> tuple[list[Primer], list[dict[str, str]]]:
    """Primers read from .csv/.tsv text (or columns pasted from a spreadsheet), and the rows
    left out. Each rejected row carries its ``name``, ``sequence`` and the ``reason``.

    The delimiter is sniffed (comma or tab). A header naming the columns is used when present;
    otherwise the first two columns are name and sequence, the third notes, the fourth location.
    """
    text = text.lstrip("﻿")
    if not text.strip():
        return [], []
    delimiter = "\t" if text.count("\t") >= text.count(",") and "\t" in text else ","
    rows = list(csv.reader(io.StringIO(text), delimiter=delimiter))
    rows = [r for r in rows if any(cell.strip() for cell in r)]
    if not rows:
        return [], []
    header = _normalise_header(rows[0])
    index = {name: header.index(name) for name in COLUMNS if header and name in header} if header else {}
    body = rows[1:] if header else rows

    def cell(row: list[str], key: str, default: int) -> str:
        pos = index.get(key, default)
        return row[pos].strip() if 0 <= pos < len(row) else ""

    primers: list[Primer] = []
    rejected: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in body:
        name = cell(row, "name", 0)
        raw_seq = cell(row, "sequence", 1)
        if not raw_seq:
            continue  # a blank row, or a stray name with no sequence: nothing to import
        try:
            sequence = clean_sequence(raw_seq)
        except PrimerLibraryError as exc:
            rejected.append({"name": name, "sequence": raw_seq, "reason": str(exc)})
            continue
        if sequence in seen:  # the same oligo twice in the file: keep the first, note the rest
            rejected.append({"name": name, "sequence": raw_seq,
                             "reason": "it repeats a sequence already in the file."})
            continue
        seen.add(sequence)
        primers.append(Primer(name=name or "(unnamed)", sequence=sequence,
                              notes=cell(row, "notes", 2), location=cell(row, "location", 3),
                              added=cell(row, "added", 4)))
    return primers, rejected


def to_csv(primers: list[Primer]) -> str:
    """The library as .csv text, with a header row; written with \\r\\n so Excel opens it cleanly."""
    out = io.StringIO()
    writer = csv.writer(out, lineterminator="\r\n")
    writer.writerow(COLUMNS)
    for p in primers:
        writer.writerow([p.name, p.sequence, p.notes, p.location, p.added])
    return out.getvalue()


def read_library() -> dict[str, Any]:
    """The library at the chosen path: its primers, the rows that could not be read, the path
    and its modification time (so the window can tell when a colleague changed the file)."""
    path = settings.primer_library_path()
    if not path.is_file():
        return {"path": str(path), "exists": False, "mtime": None, "primers": [], "rejected": []}
    try:
        text = path.read_text(encoding="utf-8-sig")
        mtime = path.stat().st_mtime
    except OSError as exc:
        log.error("Could not read the primer library %s: %r", path, exc)
        raise PrimerLibraryError(
            f"The primer library file “{path.name}” could not be read: {_os_reason(exc)}") from exc
    primers, rejected = parse_csv(text)
    return {"path": str(path), "exists": True, "mtime": mtime,
            "primers": [p.as_dict() for p in primers], "rejected": rejected}


def write_library(primers: list[Primer]) -> dict[str, Any]:
    path = settings.primer_library_path()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(to_csv(primers), encoding="utf-8")
        mtime = path.stat().st_mtime
    except OSError as exc:
        log.error("Could not write the primer library %s: %r", path, exc)
        raise PrimerLibraryError(
            f"The primer library file “{path.name}” could not be saved: {_os_reason(exc)}. "
            "Choose another location for the library in the Primer Library dialog.") from exc
    return {"path": str(path), "mtime": mtime, "count": len(primers)}


def _os_reason(exc: OSError) -> str:
    import errno
    if isinstance(exc, PermissionError):
        return "Vettore is not allowed to use it there, or another program has it open"
    if isinstance(exc, FileNotFoundError):
        return "the folder does not exist any more"
    if exc.errno == errno.ENOSPC:
        return "the disk is full"
    if exc.errno == errno.EROFS:
        return "the disk or folder is read-only"
    return "the system reported a problem"


# ---- binding: the one rule used everywhere -----------------------------------------------

@dataclass(frozen=True)
class Binding:
    """Where a primer anneals on a template. ``start``/``end`` are 0-based inclusive on the top
    strand; ``strand`` is +1 (the primer reads along the top strand) or −1 (the bottom strand).
    ``footprint`` bases anneal; ``tail`` bases at the 5′ end do not; ``mismatches`` fall inside
    the footprint, outside the 3′ anchor."""

    start: int
    end: int
    strand: int
    footprint: int
    tail: int
    mismatches: int

    def as_dict(self) -> dict[str, Any]:
        return {"start": self.start, "end": self.end, "strand": self.strand,
                "footprint": self.footprint, "tail": self.tail, "mismatches": self.mismatches}


def _rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def _anneal(primer: str, template_window: str) -> tuple[int, int] | None:
    """Score a primer 3′-flush against ``template_window`` (same length, primer orientation).
    Returns (footprint, mismatches) or None when the 3′ anchor does not match exactly.

    The footprint is the longest 3′ stretch (the anchor and the bases 5′ of it) that holds at
    most MAX_MISMATCHES mismatches and whose own 5′ end is a match; everything 5′ of it is the
    tail. So scattered mismatches inside an otherwise-annealing primer count against it, but a
    5′ run that no longer anneals is dropped as a tail rather than scored."""
    anchor = min(ANCHOR, len(primer))
    if primer[-anchor:] != template_window[-anchor:]:
        return None
    body_end = len(primer) - anchor          # the non-anchor region is primer[0:body_end]
    for start in range(body_end + 1):
        if start < body_end and primer[start] != template_window[start]:
            continue                         # a footprint never begins on a mismatch (that's tail)
        mismatches = sum(primer[i] != template_window[i] for i in range(start, body_end))
        if mismatches <= MAX_MISMATCHES:
            return len(primer) - start, mismatches
    return anchor, 0                          # unreachable: start == body_end always matches


def _find_anchor(text: str, probe: str) -> list[int]:
    """Start indices of every exact occurrence of ``probe`` in ``text`` (overlapping)."""
    hits, start = [], text.find(probe)
    while start != -1:
        hits.append(start)
        start = text.find(probe, start + 1)
    return hits


def find_bindings(primer: str, sequence: str, circular: bool) -> list[Binding]:
    """Every place the primer anneals on either strand of the template (see the module
    docstring for the rule). Deterministic: sorted by position, then strand."""
    primer = clean_sequence(primer)
    top = sequence.upper()
    n = len(top)
    if n == 0 or len(primer) < 1:
        return []
    anchor = min(ANCHOR, len(primer))
    # A circular template is searched on a copy extended by one primer length, so a primer that
    # straddles the origin is found once; positions are taken modulo n.
    extend = min(len(primer) - 1, n) if circular else 0
    text = top + top[:extend]
    results: dict[tuple[int, int], Binding] = {}

    # Forward: the primer's 3′ anchor matches the top strand; it reads along the top strand (+1).
    for hit in _find_anchor(text, primer[-anchor:]):
        anchor_start = hit                      # 5′ base of the anchor on the top strand
        window_start = anchor_start - (len(primer) - anchor)
        if window_start < 0:
            continue
        window = text[window_start:window_start + len(primer)]
        scored = _anneal(primer, window)
        if scored is None:
            continue
        footprint, mism = scored
        end = (anchor_start + anchor - 1) % n   # 3′ base
        start = (end - footprint + 1) % n       # 5′ base of the annealing footprint
        results[(start, 1)] = Binding(start, end, 1, footprint, len(primer) - footprint, mism)

    # Reverse: the primer reads along the bottom strand (−1); its reverse complement matches the
    # top strand, with the primer's 3′ end at the low-coordinate end of the match.
    rc = _rc(primer)
    for hit in _find_anchor(text, rc[:anchor]):
        window = text[hit:hit + len(primer)]
        if len(window) < len(primer):
            continue
        scored = _anneal(primer, _rc(window))
        if scored is None:
            continue
        footprint, mism = scored
        end = hit % n                           # 3′ base of the primer (low coordinate)
        start = (hit + footprint - 1) % n       # 5′ base (high coordinate)
        results[(end, -1)] = Binding(start, end, -1, footprint, len(primer) - footprint, mism)

    return [results[k] for k in sorted(results)]


def _span(low: int, high: int, n: int) -> int:
    """Length of the inclusive span from ``low`` to ``high`` on a sequence of length n, wrapping
    across the origin when high is before low."""
    return (high - low) % n + 1


def library_on_sequence(primers: list[Primer], sequence: str, circular: bool) -> list[dict[str, Any]]:
    """Every library primer that binds the sequence, with its binding sites (the map overlay and
    the Primer Library dialog use this)."""
    seq = sequence.upper()
    found = []
    for p in primers:
        try:
            bindings = find_bindings(p.sequence, seq, circular)
        except PrimerLibraryError:
            continue  # a malformed stored sequence never breaks the map
        if bindings:
            found.append({**p.as_dict(), "bindings": [b.as_dict() for b in bindings]})
    return found


def amplified_products(primers: list[Primer], sequence: str, circular: bool,
                       max_product: int = 10000) -> list[dict[str, Any]]:
    """Pairs of library primers that make a PCR product on the template: a forward binding and a
    reverse binding downstream of it, within ``max_product`` bp. Each product carries both primer
    names, the amplicon size (tails included) and the span it covers. Smallest product first."""
    seq = sequence.upper()
    n = len(seq)
    forwards, reverses = [], []
    for p in primers:
        try:
            for b in find_bindings(p.sequence, seq, circular):
                (forwards if b.strand == 1 else reverses).append((p, b))
        except PrimerLibraryError:
            continue
    products = []
    for fp, fb in forwards:
        for rp, rb in reverses:
            # The amplicon runs from the forward primer's 5′ end to the reverse primer's 5′ end.
            length = _span(fb.start, rb.start, n) if circular else rb.start - fb.start + 1
            if length <= 0 or length > max_product:
                continue
            if fb.end == rb.end:  # the very same site read both ways: not a product
                continue
            products.append({
                "forward": fp.name, "reverse": rp.name,
                "forward_sequence": fp.sequence, "reverse_sequence": rp.sequence,
                "size": length + fb.tail + rb.tail,
                "product_start": fb.start, "product_end": rb.start,
            })
    products.sort(key=lambda x: (x["size"], x["forward"], x["reverse"]))
    return products


def product_over_region(primers: list[Primer], sequence: str, circular: bool,
                        start: int, end: int) -> dict[str, Any] | None:
    """The smallest library-primer product that already contains the region [start, end]
    (0-based inclusive), or None. Used to tell the user "you already have F and R"."""
    n = len(sequence)
    if n == 0:
        return None
    for product in amplified_products(primers, sequence, circular):
        p_start, p_end = product["product_start"], product["product_end"]
        if circular:
            if _span(p_start, start, n) <= _span(p_start, p_end, n) and \
               _span(p_start, end, n) <= _span(p_start, p_end, n):
                return product
        elif p_start <= start and end <= p_end:
            return product
    return None


def sequencing_primers(primers: list[Primer], sequence: str, circular: bool,
                      start: int, end: int) -> list[dict[str, Any]]:
    """Library primers that read across the region [start, end]: a forward primer whose 3′ end
    sits SEQ_PRIMER_MIN–SEQ_PRIMER_MAX bp upstream of the region reads into it; a reverse primer
    the same distance downstream reads back across it. Listed in reading order."""
    seq = sequence.upper()
    n = len(seq)
    if n == 0:
        return []
    out = []
    for p in primers:
        try:
            bindings = find_bindings(p.sequence, seq, circular)
        except PrimerLibraryError:
            continue
        for b in bindings:
            if b.strand == 1:
                gap = _span(b.end, start, n) - 1 if circular else start - b.end - 1
                reads_into = gap >= 0
            else:
                gap = _span(end, b.end, n) - 1 if circular else b.end - end - 1
                reads_into = gap >= 0
            if reads_into and SEQ_PRIMER_MIN <= gap + 1 <= SEQ_PRIMER_MAX:
                out.append({"name": p.name, "sequence": p.sequence, "strand": b.strand,
                            "start": b.start, "end": b.end, "reads_from": gap + 1})
    out.sort(key=lambda x: (x["strand"] != 1, x["end"]))
    return out


# ---- HTTP --------------------------------------------------------------------------------

class SetPath(BaseModel):
    path: str


class SaveLibrary(BaseModel):
    primers: list[dict[str, Any]]


class ImportText(BaseModel):
    text: str


class OnSequence(BaseModel):
    """A document, and optionally the primers to look for on it: the Primer Library dialog sends
    its working list (with changes not yet saved); without one, the saved shared file is used."""

    sequenceData: dict[str, Any]
    primers: list[dict[str, Any]] | None = None


class SequencingRequest(BaseModel):
    sequenceData: dict[str, Any]
    start: int
    end: int


def _primers_from_payload(rows: list[dict[str, Any]]) -> list[Primer]:
    primers, seen = [], set()
    for row in rows:
        try:
            sequence = clean_sequence(str(row.get("sequence", "")))
        except PrimerLibraryError as exc:
            raise HTTPException(422, f"“{row.get('name') or 'A primer'}”: {exc}") from exc
        if sequence in seen:
            continue
        seen.add(sequence)
        primers.append(Primer(
            name=str(row.get("name") or "(unnamed)").strip(), sequence=sequence,
            notes=str(row.get("notes") or "").strip(), location=str(row.get("location") or "").strip(),
            added=str(row.get("added") or "").strip()))
    return primers


def _sequence_of(data: dict[str, Any]) -> tuple[str, bool]:
    sequence = data.get("sequence")
    return (sequence.upper() if isinstance(sequence, str) else ""), bool(data.get("circular"))


@router.get("")
def get_library() -> dict[str, Any]:
    try:
        return read_library()
    except PrimerLibraryError as exc:
        raise HTTPException(422, str(exc)) from exc


@router.post("/path")
def set_path(body: SetPath) -> dict[str, Any]:
    raw = body.path.strip()
    if not raw:
        raise HTTPException(422, "Choose a .csv file for the primer library.")
    settings.remember_primer_library_path(raw)
    try:
        return read_library()
    except PrimerLibraryError as exc:
        raise HTTPException(422, str(exc)) from exc


@router.post("/save")
def save_library(body: SaveLibrary) -> dict[str, Any]:
    primers = _primers_from_payload(body.primers)
    try:
        return write_library(primers)
    except PrimerLibraryError as exc:
        raise HTTPException(422, str(exc)) from exc


@router.post("/parse")
def parse(body: ImportText) -> dict[str, Any]:
    """Parse pasted or file text without saving it, so the dialog can show what will be imported
    and what will be left out before the user keeps it."""
    primers, rejected = parse_csv(body.text)
    return {"primers": [p.as_dict() for p in primers], "rejected": rejected,
            "today": date.today().isoformat()}


@router.post("/bindings")
def bindings(body: OnSequence) -> dict[str, Any]:
    sequence, circular = _sequence_of(body.sequenceData)
    if body.primers is not None:
        rows = body.primers
    else:
        try:
            rows = read_library()["primers"]
        except PrimerLibraryError as exc:
            raise HTTPException(422, str(exc)) from exc
    primers = _primers_from_payload(rows)
    return {"primers": library_on_sequence(primers, sequence, circular)}


@router.post("/product")
def product(body: SequencingRequest) -> dict[str, Any]:
    """The smallest pair of library primers that already amplifies a product covering the region
    [start, end], or ``null`` — so primer design can say "you already have F and R"."""
    sequence, circular = _sequence_of(body.sequenceData)
    if not sequence:
        return {"product": None}
    try:
        library = read_library()
    except PrimerLibraryError as exc:
        raise HTTPException(422, str(exc)) from exc
    primers = _primers_from_payload(library["primers"])
    return {"product": product_over_region(primers, sequence, circular, body.start, body.end)}


@router.post("/sequencing")
def sequencing(body: SequencingRequest) -> dict[str, Any]:
    sequence, circular = _sequence_of(body.sequenceData)
    if not sequence:
        raise HTTPException(422, "Open a sequence first.")
    try:
        library = read_library()
    except PrimerLibraryError as exc:
        raise HTTPException(422, str(exc)) from exc
    primers = _primers_from_payload(library["primers"])
    return {"primers": sequencing_primers(primers, sequence, circular, body.start, body.end)}
