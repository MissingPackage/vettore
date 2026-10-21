"""Lab primer library (S5.3): the primers a lab already keeps, so they are reused instead of
re-ordered. Routes under /api/primer-library.

The library is **one .csv file** the lab keeps in a shared folder (OneDrive, Dropbox, a lab
drive) so colleagues share it; Vettore reads it from the path the user chose and writes it back
there. The columns are ``name, sequence, notes, location, added`` under a header row, so a file
edited by hand or in Excel stays readable. Files from Excel are read as they come: separated by
commas, semicolons (Excel set to Italian and other European languages) or tabs, in UTF-8 or the
Windows code page, with the usual header names ("Primer name", "Sequence (5'-3')", "Box"...).
Vettore writes the file back with the separator it found and a UTF-8 byte-order mark, so Excel
shows accents and µ correctly. A row Vettore cannot read (a base it does not know, the same oligo
twice) is kept in the file, unchanged, on every save: a shared file never loses a colleague's
line. Nothing is locked: the dialog checks the file's modification time before saving and
offers to reload when a colleague changed it.

Binding — *the one rule used everywhere* (the map, the Primer Library dialog, primer design and
sequencing primers). A primer primes where its 3′ end anneals:

- its last ``ANCHOR`` (15) bases match the template (no mismatch at the 3′ end);
- from there toward the 5′ end, the annealing footprint grows over matching bases and up to
  ``MAX_MISMATCHES`` (2) mismatches, and always ends on a match;
- whatever lies 5′ of the footprint is a tail (a restriction site, an overlap for assembly) and
  is not scored; it adds to the product.

IUPAC codes in a primer (a degenerate primer: N, R, Y...) match every base they stand for; an
ambiguous base in the template matches only a primer base that covers it. Primers shorter than
``MIN_LENGTH`` are not searched (they would prime almost anywhere), and a primer that binds more
than ``MAX_SITES`` places is reported with its count instead of being drawn. Overlapping hits
of one primer on one strand (a primer inside a repeat) count once, the best of them.

Positions are 0-based and inclusive, as in the editor: ``start`` is always the lower coordinate
and ``end`` the higher, on both strands, and ``start > end`` means the site crosses the origin
of a circular sequence. A forward primer's 3′ end is ``end``; a reverse primer's is ``start``.
"""

from __future__ import annotations

import csv
import io
import logging
import re
from dataclasses import dataclass
from datetime import date
from typing import Any

from Bio.Seq import Seq
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from . import settings, tm

router = APIRouter(prefix="/api/primer-library")
log = logging.getLogger(__name__)

#: The 3′ bases that must match the template for a primer to prime there.
ANCHOR = 15
#: Mismatches allowed in the annealing footprint 5′ of the anchor.
MAX_MISMATCHES = 2
#: Shorter primers are not searched: a 10-mer matches a plasmid by chance many times.
MIN_LENGTH = ANCHOR
#: A primer binding more places than this is listed with its count, not drawn.
MAX_SITES = 10
#: Two library primers make a product only when it is at least this long (shorter is a dimer).
MIN_PRODUCT = 40
#: Sequencing primers: Sanger leaves the first ~40 bases after the primer unreadable and reads
#: about 800 bases in all, so a primer reads a region whose first base is at least
#: SEQ_PRIMER_MIN from its 3′ end and whose last base is within SEQ_READ_LENGTH of it.
SEQ_PRIMER_MIN, SEQ_READ_LENGTH = 40, 800
#: Columns of the shared .csv, in order. ``name`` and ``sequence`` are required.
COLUMNS = ("name", "sequence", "notes", "location", "added")

#: IUPAC nucleotide codes and the bases each stands for.
_CODES = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

# Header cells as labs write them, mapped to our columns.
_HEADER_WORDS = {
    "name": ("name", "primer", "primer name", "primer_name", "oligo", "oligo name", "id", "nome"),
    "notes": ("notes", "note", "comment", "comments", "description", "note/commenti", "commenti"),
    "location": ("location", "box", "position", "where", "freezer", "posizione", "scatola"),
    "added": ("added", "date", "date added", "data"),
}


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
            "degenerate": not _is_plain_dna(self.sequence),
            "searched": len(self.sequence) >= MIN_LENGTH,
        }


def _is_plain_dna(sequence: str) -> bool:
    return bool(sequence) and all(c in "ACGT" for c in sequence)


def clean_sequence(raw: str) -> str:
    """A primer sequence in upper case, U read as T. Raises on anything but IUPAC DNA codes."""
    seq = "".join((raw or "").split()).upper().replace("U", "T")
    seq = seq.replace("-", "")       # 5'-ACGT-3' style dashes
    if seq.startswith("5'") and seq.endswith("3'"):
        seq = seq[2:-2]
    bad = sorted({c for c in seq if c not in _CODES})
    if not seq:
        raise PrimerLibraryError("it has no sequence.")
    if bad:
        raise PrimerLibraryError(
            "its sequence has characters that are not DNA ("
            + ", ".join(repr(c) for c in bad[:5]) + ").")
    return seq


# ---- the shared .csv file ----------------------------------------------------------------

def _column_of(cell: str) -> str | None:
    word = cell.strip().lower().strip("\"' ")
    if word.startswith("seq") or "sequence" in word or "sequenza" in word:
        return "sequence"
    for column, words in _HEADER_WORDS.items():
        if word in words:
            return column
    return None


def _header_index(cells: list[str]) -> dict[str, int] | None:
    """Columns of a header row, or None when the row is not a header (it has no name and
    sequence columns, or its "sequence" cell is itself DNA)."""
    index: dict[str, int] = {}
    for pos, cell in enumerate(cells):
        column = _column_of(cell)
        if column and column not in index:
            index[column] = pos
    if "name" not in index or "sequence" not in index:
        return None
    try:
        clean_sequence(cells[index["sequence"]])
        return None                      # "ACGT..." under the name column: data, not a header
    except PrimerLibraryError:
        return index


def sniff_delimiter(text: str) -> str:
    """Tab when the text has any (columns pasted from a spreadsheet), otherwise the commoner of
    ``;`` and ``,`` on the first line (Excel writes ``;`` when the system uses a decimal comma)."""
    if "\t" in text:
        return "\t"
    first = next((line for line in text.splitlines() if line.strip()), "")
    return ";" if first.count(";") > first.count(",") else ","


def parse_csv(text: str) -> tuple[list[Primer], list[dict[str, Any]]]:
    """Primers read from .csv/.tsv text (or columns pasted from a spreadsheet), and the rows
    left out. Each rejected row carries its ``name``, ``sequence``, the ``reason`` and the
    original ``row`` (its cells), so a save can write it back unchanged.

    A header naming the columns is used when present; otherwise the columns are taken in
    order: name, sequence, notes, location, added.
    """
    text = text.lstrip("﻿")
    if not text.strip():
        return [], []
    delimiter = sniff_delimiter(text)
    rows = [r for r in csv.reader(io.StringIO(text), delimiter=delimiter)
            if any(cell.strip() for cell in r)]
    if rows and len(rows[0]) == 1 and rows[0][0].strip().lower().startswith("sep="):
        rows = rows[1:]                   # Excel's "sep=;" first line
    if not rows:
        return [], []
    header = _header_index(rows[0])
    index = header or {name: pos for pos, name in enumerate(COLUMNS)}
    body = rows[1:] if header else rows

    def cell(row: list[str], key: str) -> str:
        pos = index.get(key, -1)
        return row[pos].strip() if 0 <= pos < len(row) else ""

    primers: list[Primer] = []
    rejected: list[dict[str, Any]] = []
    seen: dict[str, str] = {}
    for row in body:
        name = cell(row, "name")
        raw_seq = cell(row, "sequence")
        try:
            sequence = clean_sequence(raw_seq)
        except PrimerLibraryError as exc:
            rejected.append({"name": name, "sequence": raw_seq, "reason": str(exc), "row": row})
            continue
        if sequence in seen:
            rejected.append({"name": name, "sequence": raw_seq, "row": row,
                             "reason": f"it is the same oligo as “{seen[sequence]}”."})
            continue
        seen[sequence] = name or "(unnamed)"
        primers.append(Primer(name=name or "(unnamed)", sequence=sequence,
                              notes=cell(row, "notes"), location=cell(row, "location"),
                              added=cell(row, "added")))
    return primers, rejected


def to_csv(primers: list[Primer], delimiter: str = ",",
           kept_rows: list[list[str]] | None = None) -> str:
    """The library as .csv text: a header row, the primers, then ``kept_rows`` (rows of the
    file Vettore could not read, written back as they were). Lines end in \\r\\n, as Excel's."""
    out = io.StringIO()
    writer = csv.writer(out, delimiter=delimiter, lineterminator="\r\n")
    writer.writerow(COLUMNS)
    for p in primers:
        writer.writerow([p.name, p.sequence, p.notes, p.location, p.added])
    for row in kept_rows or []:
        writer.writerow(row)
    return out.getvalue()


def _decode(data: bytes) -> str:
    """The file's text: UTF-8 (with or without a byte-order mark), else the Windows code page
    Excel uses for .csv on Western European systems."""
    try:
        return data.decode("utf-8-sig")
    except UnicodeDecodeError:
        return data.decode("cp1252", errors="replace")


def _read_file() -> tuple[str | None, float | None]:
    path = settings.primer_library_path()
    if not path.is_file():
        return None, None
    try:
        return _decode(path.read_bytes()), path.stat().st_mtime
    except OSError as exc:
        log.error("Could not read the primer library %s: %r", path, exc)
        raise PrimerLibraryError(
            f"The primer library file “{path.name}” could not be read: {_os_reason(exc)}.") from exc


def read_library() -> dict[str, Any]:
    """The library at the chosen path: its primers, the rows that could not be read, the path
    and its modification time (so the window can tell when a colleague changed the file)."""
    path = settings.primer_library_path()
    text, mtime = _read_file()
    if text is None:
        return {"path": str(path), "exists": False, "mtime": None, "primers": [], "rejected": []}
    primers, rejected = parse_csv(text)
    return {"path": str(path), "exists": True, "mtime": mtime,
            "primers": [p.as_dict() for p in primers],
            "rejected": [{k: v for k, v in r.items() if k != "row"} for r in rejected]}


def write_library(primers: list[Primer]) -> dict[str, Any]:
    """Write the primers to the shared file. Rows of the file as it is now that Vettore could
    not read are kept at the end, unchanged, and the file keeps its separator."""
    path = settings.primer_library_path()
    text, _ = _read_file()
    delimiter, kept = ",", []
    if text is not None and text.strip():
        delimiter = sniff_delimiter(text)
        kept = [r["row"] for r in parse_csv(text)[1]]
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(f".{path.name}.vettore-tmp")
        tmp.write_text(to_csv(primers, delimiter, kept), encoding="utf-8-sig", newline="")
        tmp.replace(path)
        mtime = path.stat().st_mtime
    except OSError as exc:
        log.error("Could not write the primer library %s: %r", path, exc)
        raise PrimerLibraryError(
            f"The primer library file “{path.name}” could not be saved: {_os_reason(exc)}. "
            "Choose another location for the library in the Primer Library dialog.") from exc
    return {"path": str(path), "mtime": mtime, "count": len(primers), "kept": len(kept)}


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
    """Where a primer anneals. ``start``/``end``: 0-based inclusive, ``start`` the lower
    coordinate on both strands (``start > end`` crosses the origin). ``strand`` +1: the primer
    reads along the top strand, its 3′ end at ``end``; −1: along the bottom strand, its 3′ end
    at ``start``. ``footprint`` bases anneal, ``tail`` 5′ bases do not, and ``mismatches`` fall
    inside the footprint."""

    start: int
    end: int
    strand: int
    footprint: int
    tail: int
    mismatches: int

    @property
    def three_prime(self) -> int:
        return self.end if self.strand == 1 else self.start

    @property
    def five_prime(self) -> int:
        """The 5′ base of the annealing footprint (the tail, if any, lies beyond it)."""
        return self.start if self.strand == 1 else self.end

    def as_dict(self) -> dict[str, Any]:
        return {"start": self.start, "end": self.end, "strand": self.strand,
                "footprint": self.footprint, "tail": self.tail, "mismatches": self.mismatches}


def _rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def _pairs(primer_base: str, template_base: str | None) -> bool:
    """A primer base anneals to a template base: the same base, or a code covering it. An
    ambiguous template base pairs only with a primer code covering all it stands for."""
    if template_base is None:
        return False
    allowed = _CODES.get(primer_base, "")
    return all(b in allowed for b in _CODES.get(template_base, "?"))


def _anneal(primer: str, window: list[str | None]) -> tuple[int, int] | None:
    """Score a primer 3′-flush against ``window`` (the template bases it would pair with, in
    primer orientation; None past the end of a linear sequence). Returns (footprint,
    mismatches), or None when the 3′ anchor does not pair."""
    anchor = ANCHOR
    if not all(_pairs(p, t) for p, t in zip(primer[-anchor:], window[-anchor:], strict=True)):
        return None
    body_end = len(primer) - anchor
    for start in range(body_end + 1):
        if start < body_end and not _pairs(primer[start], window[start]):
            continue                         # a footprint never begins on a mismatch
        mismatches = sum(not _pairs(primer[i], window[i]) for i in range(start, body_end))
        if mismatches <= MAX_MISMATCHES:
            return len(primer) - start, mismatches
    return anchor, 0


def _anchor_pattern(probe: str) -> re.Pattern[str]:
    """Template stretches the anchor pairs with: each IUPAC code becomes the set of template
    letters it covers (an ambiguous template letter only where the code covers all of it)."""
    parts = []
    for code in probe:
        letters = "".join(t for t in _CODES if all(b in _CODES[code] for b in _CODES[t]))
        parts.append(f"[{letters}]" if len(letters) > 1 else letters)
    return re.compile("(?=(" + "".join(parts) + "))")


def _candidates(primer: str, top: str, n: int, circular: bool) -> list[Binding]:
    anchor = ANCHOR
    length = len(primer)
    text = top + top[:anchor - 1] if circular else top

    def base(i: int) -> str | None:
        if circular:
            return top[i % n]
        return top[i] if 0 <= i < n else None

    found: list[Binding] = []
    # Forward: the anchor pairs with the top strand, the primer's 3′ end at the anchor's end.
    for m in _anchor_pattern(primer[-anchor:]).finditer(text):
        three = m.start() + anchor - 1
        window = [base(three - length + 1 + k) for k in range(length)]
        if (scored := _anneal(primer, window)) is None:
            continue
        footprint, mism = scored
        found.append(Binding((three - footprint + 1) % n, three % n, 1, footprint,
                             length - footprint, mism))
    # Reverse: the reverse complement of the anchor lies on the top strand, the primer's 3′ end
    # at its low end; the primer pairs with the bases above it, read right to left.
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    for m in _anchor_pattern(_rc(primer[-anchor:])).finditer(text):
        three = m.start()
        window = []
        for k in range(length):          # primer base k pairs with top base three + length-1-k
            t = base(three + length - 1 - k)
            window.append(None if t is None else
                          complement.get(t) or _rc(t))
        if (scored := _anneal(primer, window)) is None:
            continue
        footprint, mism = scored
        found.append(Binding(three % n, (three + footprint - 1) % n, -1, footprint,
                             length - footprint, mism))
    return found


def _covers(b: Binding, pos: int, n: int) -> bool:
    return (pos - b.start) % n <= (b.end - b.start) % n


def find_bindings(primer: str, sequence: str, circular: bool) -> list[Binding]:
    """Every place the primer anneals on either strand (see the module docstring for the rule).
    Overlapping hits on one strand count once, the best. Sorted by position, then strand.
    Primers shorter than MIN_LENGTH, and templates shorter than the primer, give none."""
    primer = clean_sequence(primer)
    top = sequence.upper()
    n = len(top)
    if len(primer) < MIN_LENGTH or n < len(primer):
        return []
    hits = sorted(_candidates(primer, top, n, circular),
                  key=lambda b: (b.mismatches, -b.footprint, b.start, b.strand))
    kept: list[Binding] = []
    for b in hits:
        if any(k.strand == b.strand and _covers(k, b.three_prime, n) for k in kept):
            continue
        kept.append(b)
    return sorted(kept, key=lambda b: (b.start, b.strand))


def _span(low: int, high: int, n: int) -> int:
    """Length of the inclusive span from ``low`` to ``high`` on a circle of length n."""
    return (high - low) % n + 1


def library_on_sequence(primers: list[Primer], sequence: str, circular: bool) -> list[dict[str, Any]]:
    """Every library primer that binds the sequence, with its binding sites (the map overlay and
    the Primer Library dialog use this). A primer binding more than MAX_SITES places is listed
    with ``sites`` and no bindings, so the map is not flooded."""
    seq = sequence.upper()
    found = []
    for p in primers:
        try:
            bindings = find_bindings(p.sequence, seq, circular)
        except PrimerLibraryError:
            continue  # a malformed stored sequence never breaks the map
        if not bindings:
            continue
        entry = {**p.as_dict(), "sites": len(bindings)}
        entry["bindings"] = [] if len(bindings) > MAX_SITES else [b.as_dict() for b in bindings]
        found.append(entry)
    return found


def _bindings_by_strand(primers: list[Primer], seq: str, circular: bool):
    forwards, reverses = [], []
    for p in primers:
        try:
            bindings = find_bindings(p.sequence, seq, circular)
        except PrimerLibraryError:
            continue
        if len(bindings) > MAX_SITES:
            continue
        for b in bindings:
            (forwards if b.strand == 1 else reverses).append((p, b))
    return forwards, reverses


def amplified_products(primers: list[Primer], sequence: str, circular: bool,
                       max_product: int = 10000) -> list[dict[str, Any]]:
    """Pairs of library primers that make a PCR product: a forward binding and a reverse
    binding whose 3′ end lies beyond the forward's footprint, MIN_PRODUCT to ``max_product``
    bp apart. The amplicon runs from the forward footprint's 5′ base to the reverse's; its size
    includes both tails. Smallest product first."""
    seq = sequence.upper()
    n = len(seq)
    forwards, reverses = _bindings_by_strand(primers, seq, circular)
    products = []
    for fp, fb in forwards:
        for rp, rb in reverses:
            if circular:
                to_r3 = _span(fb.start, rb.start, n)
                length = _span(fb.start, rb.end, n)
                if to_r3 <= fb.footprint or length < to_r3:
                    continue                 # the reverse 3′ end is not past the forward one
            else:
                if rb.start <= fb.end:
                    continue
                length = rb.end - fb.start + 1
            if not MIN_PRODUCT <= length <= max_product:
                continue
            products.append({
                "forward": fp.name, "reverse": rp.name,
                "forward_sequence": fp.sequence, "reverse_sequence": rp.sequence,
                "size": length + fb.tail + rb.tail,
                "product_start": fb.start, "product_end": rb.end,
            })
    products.sort(key=lambda x: (x["size"], x["forward"], x["reverse"]))
    return products


def product_over_region(primers: list[Primer], sequence: str, circular: bool,
                        start: int, end: int) -> dict[str, Any] | None:
    """The smallest library-primer product that contains the whole region [start, end]
    (0-based inclusive; on a circular sequence ``start > end`` crosses the origin), or None.
    Used to tell the user "you already have F and R"."""
    n = len(sequence)
    if n == 0:
        return None
    for product in amplified_products(primers, sequence, circular):
        p_start, p_end = product["product_start"], product["product_end"]
        if circular:
            to_start, to_end = _span(p_start, start, n), _span(p_start, end, n)
            if to_start <= to_end <= _span(p_start, p_end, n):
                return product
        elif p_start <= start <= end <= p_end:
            return product
    return None


def sequencing_primers(primers: list[Primer], sequence: str, circular: bool,
                      start: int, end: int) -> list[dict[str, Any]]:
    """Library primers that read across the region [start, end]: a forward primer whose 3′ end
    is at least SEQ_PRIMER_MIN bases before the region, a reverse primer the same after it, and
    the whole region within SEQ_READ_LENGTH of the 3′ end. ``reads_from`` is the distance from
    the 3′ end to the region's nearer edge. Forward primers first, nearest first."""
    seq = sequence.upper()
    n = len(seq)
    if n == 0:
        return []
    region = _span(start, end, n) if circular else end - start + 1
    forwards, reverses = _bindings_by_strand(primers, seq, circular)
    out = []
    for p, b in forwards + reverses:
        if b.strand == 1:
            distance = _span(b.end, start, n) - 1 if circular else start - b.end
        else:
            distance = _span(end, b.start, n) - 1 if circular else b.start - end
        if distance >= SEQ_PRIMER_MIN and distance + region <= SEQ_READ_LENGTH:
            out.append({"name": p.name, "sequence": p.sequence, "strand": b.strand,
                        "start": b.start, "end": b.end, "reads_from": distance})
    out.sort(key=lambda x: (x["strand"] != 1, x["reads_from"]))
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
