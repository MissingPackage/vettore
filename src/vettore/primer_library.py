"""Lab primer library (S5.3): the primers a lab already keeps, so they are reused instead of
re-ordered. Routes under /api/primer-library.

The library is **one .csv file** the lab keeps in a shared folder (OneDrive, Dropbox, a lab
drive) so colleagues share it; Vettore reads it from the path the user chose and writes it back
there. The columns are ``name, sequence, notes, location, added`` under a header row, so a file
edited by hand or in Excel stays readable. Files from Excel are read as they come: separated by
commas, semicolons (Excel set to Italian and other European languages) or tabs, in UTF-8 or the
Windows code page, with the usual header names ("Primer name", "Sequence (5'-3')", "Box"...).
A save rewrites only the rows that changed: the header, extra columns (Tm, who ordered it),
a ``sep=;`` line, the separator and the way each sequence is written all stay, and the file is
written with a UTF-8 byte-order mark so Excel shows accents and µ correctly. A row Vettore cannot
read (a base it does not know, the same oligo twice) is kept, unchanged, on every save. Nothing
is locked: a save reads the file again and merges into it (see ``merged_text``), so a primer a
colleague added meanwhile is kept. A workbook or any file that is not a primer list is refused
and never overwritten.

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
from dataclasses import dataclass, field, replace
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

# Header cells as labs write them, mapped to our columns, best match first.
_HEADER_WORDS = {
    "name": (("name", "primer name", "primer_name", "oligo name", "nome", "nome primer"),
             ("primer", "oligo"), ("id",)),
    "notes": (("notes", "note", "comment", "comments", "description", "note/commenti", "commenti"),),
    "location": (("location", "box", "position", "where", "freezer", "posizione", "scatola"),),
    "added": (("added", "date", "date added", "data"),),
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
    seq = re.sub("[’′]", "'", seq).replace("-", "")   # 5'-ACGT-3' (or 5′-…-3′)
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

#: Extensions a library file may have: a spreadsheet saved as text.
LIBRARY_SUFFIXES = (".csv", ".tsv", ".txt")
_BOM = b"\xef\xbb\xbf"
#: Rows looked at for the header: a title line or two may sit above it.
_HEADER_SEARCH = 10
#: Words that make a "Sequence …" header cell something else than the sequence itself.
_NOT_THE_SEQUENCE = ("id", "length", "len", "lunghezza", "no", "no.", "number", "numero", "#", "nr")


def _column_of(cell: str) -> tuple[str, int] | None:
    """The column a header cell names and how sure that is (0 best): "Name" beats "Primer",
    which beats "ID", so a file with both an ID and a Name column takes names from Name.
    "Sequence name" is a name, and "Sequence ID" or "Sequence length" are not the sequence."""
    word = cell.strip().lower().strip("\"' ")
    if word in ("sequence", "seq", "sequenza"):
        return "sequence", 0
    if word.startswith("seq") or "sequence" in word or "sequenza" in word:
        parts = set(re.split(r"[\s_.()/-]+", word))
        if parts & {"name", "nome"}:
            return "name", 1
        if parts & set(_NOT_THE_SEQUENCE) or word.endswith("#"):
            return None
        return "sequence", 1
    for column, words in _HEADER_WORDS.items():
        for rank, group in enumerate(words):
            if word in group:
                return column, rank
    return None


def _header_index(cells: list[str]) -> tuple[dict[str, int], list[int]] | None:
    """Columns of a header row and every column that could hold the sequence (best first), or
    None when the row is not a header (it has no name and sequence columns, or its "sequence"
    cell is itself DNA)."""
    best: dict[str, tuple[int, int]] = {}
    sequences: list[tuple[int, int]] = []
    for pos, cell in enumerate(cells):
        found = _column_of(cell)
        if not found:
            continue
        if found[0] == "sequence":
            sequences.append((found[1], pos))
        if found[0] not in best or found[1] < best[found[0]][0]:
            best[found[0]] = (found[1], pos)
    index = {column: pos for column, (_, pos) in best.items()}
    if "name" not in index or "sequence" not in index:
        return None
    if _is_dna(cells[index["sequence"]]):
        return None                      # "ACGT..." under the name column: data, not a header
    return index, [pos for _, pos in sorted(sequences)]


def _split(line: str, delimiter: str) -> list[str]:
    try:
        return next(csv.reader([line], delimiter=delimiter), [])
    except csv.Error:
        return []


def _is_dna(cell: str) -> bool:
    try:
        clean_sequence(cell)
        return True
    except PrimerLibraryError:
        return False


def sniff_delimiter(lines: list[str]) -> str:
    """The separator of a library or a paste: the one that turns one of the first lines into a
    header; without a header, the one that puts DNA in the first or second column of most rows
    (pasted columns are name, sequence...); failing both, the one that splits the first line
    into most cells. Excel writes ``;`` when the system uses a decimal comma, and pasted
    columns use tabs."""
    rows = [line for line in lines if line.strip()][:20]
    if not rows:
        return ","
    candidates = (",", ";", "\t")
    for row in rows[:_HEADER_SEARCH]:
        headed = [d for d in candidates if _header_index(_split(row, d))]
        if headed:
            return max(headed, key=lambda d: len(_split(row, d)))

    def dna_cells(d: str) -> int:
        count = 0
        for r in rows:
            cells = _split(r, d)
            count += len(cells) > 1 and (_is_dna(cells[1]) or _is_dna(cells[0]))
        return count

    dna = {d: dna_cells(d) for d in candidates}
    if max(dna.values()) > 0:
        return max(candidates, key=lambda d: (dna[d], -candidates.index(d)))
    return max(candidates, key=lambda d: (len(_split(rows[0], d)), -candidates.index(d)))


@dataclass
class _Row:
    """One row of the file: its cells, its text exactly as in the file (written back as is
    while nobody changes it), and the cleaned sequence when the row is a primer."""

    cells: list[str]
    raw: str
    key: str | None = None


@dataclass
class LibraryText:
    """A library file (or a paste) as read: the separator, any ``sep=`` or title lines, the
    header and every row in order, plus the primers and the rows that could not be read."""

    delimiter: str
    newline: str = "\r\n"
    preamble: list[str] = field(default_factory=list)
    header: _Row | None = None
    index: dict[str, int] = field(default_factory=dict)
    rows: list[_Row] = field(default_factory=list)
    primers: list[Primer] = field(default_factory=list)
    rejected: list[dict[str, Any]] = field(default_factory=list)


def _lines(text: str) -> list[str]:
    """The text's lines, each with its own line ending (\\r\\n, \\n or \\r); the last one keeps
    none when the text does not end with a newline."""
    return re.findall(r"[^\r\n]*(?:\r\n|\r|\n)|[^\r\n]+$", text)


def _positional(body: list[_Row]) -> dict[str, int]:
    """Columns of a file or paste without a header: name, sequence, notes, location, added;
    or sequence first (then name) when it is the first column that holds DNA."""
    rows = [r.cells for r in body if any(c.strip() for c in r.cells)][:50]
    first = sum(bool(c) and _is_dna(c[0]) for c in rows)
    second = sum(len(c) > 1 and _is_dna(c[1]) for c in rows)
    if rows and first > second:
        return {"sequence": 0, "name": 1, "notes": 2, "location": 3, "added": 4}
    return {name: pos for pos, name in enumerate(COLUMNS)}


def parse_text(text: str) -> LibraryText:
    """Read .csv/.tsv text (or columns pasted from a spreadsheet). A header naming the columns
    is used when present (title lines above it are kept as they are); otherwise the columns are
    taken in order: name, sequence, notes, location, added. Raises PrimerLibraryError when the
    text is not readable as columns (a quote that is never closed, for one)."""
    text = text.lstrip("﻿")
    newline = next((nl for nl in ("\r\n", "\n", "\r") if nl in text), "\r\n")
    lines = _lines(text)
    preamble: list[str] = []
    forced = None
    if lines and (m := re.fullmatch(r"\s*sep=(.)\s*", lines[0], re.IGNORECASE)):
        forced = m.group(1)              # Excel's "sep=;" first line names the separator
        preamble.append(lines.pop(0))
    delimiter = forced or sniff_delimiter(lines)
    lib = LibraryText(delimiter=delimiter, newline=newline, preamble=preamble)

    rows: list[_Row] = []
    consumed: list[str] = []

    def feed():
        for line in lines:
            consumed.append(line)
            yield line

    try:
        for cells in csv.reader(feed(), delimiter=delimiter, strict=True):
            rows.append(_Row(cells, "".join(consumed)))
            consumed.clear()
    except csv.Error as exc:
        line = sum(len(_lines(r.raw)) for r in rows) + 1
        raise PrimerLibraryError(
            f"the row starting on line {line + len(preamble)} could not be read as spreadsheet "
            f"columns ({exc}); a quote (\") that is opened and never closed does this.") from exc

    found = None
    for pos, row in enumerate(rows):
        if sum(any(c.strip() for c in r.cells) for r in rows[:pos]) >= _HEADER_SEARCH:
            break
        if any(c.strip() for c in row.cells) and (found := _header_index(row.cells)):
            lib.header = row
            lib.preamble += [r.raw for r in rows[:pos]]
            lib.rows = rows[pos + 1:]
            break
    if found:
        lib.index, sequence_columns = found
        if len(sequence_columns) > 1:     # "Sequence (5'-3')" and "Sequence (RNA)": the DNA one
            lib.index["sequence"] = max(sequence_columns, key=lambda pos: sum(
                pos < len(r.cells) and _is_dna(r.cells[pos]) for r in lib.rows))
    else:
        lib.rows = rows
        lib.index = _positional(rows)

    seen: dict[str, str] = {}
    for row in lib.rows:
        if not any(c.strip() for c in row.cells):
            continue
        name = _cell(row.cells, lib.index, "name")
        raw_seq = _cell(row.cells, lib.index, "sequence")
        try:
            sequence = clean_sequence(raw_seq)
        except PrimerLibraryError as exc:
            lib.rejected.append({"name": name, "sequence": raw_seq, "reason": str(exc)})
            continue
        if sequence in seen:
            lib.rejected.append({"name": name, "sequence": raw_seq,
                                 "reason": f"it is the same oligo as “{seen[sequence]}”."})
            continue
        seen[sequence] = name or "(unnamed)"
        row.key = sequence
        lib.primers.append(Primer(
            name=name or "(unnamed)", sequence=sequence, notes=_cell(row.cells, lib.index, "notes"),
            location=_cell(row.cells, lib.index, "location"), added=_cell(row.cells, lib.index, "added")))
    return lib


def _cell(cells: list[str], index: dict[str, int], key: str) -> str:
    pos = index.get(key, -1)
    return cells[pos].strip() if 0 <= pos < len(cells) else ""


def parse_csv(text: str) -> tuple[list[Primer], list[dict[str, Any]]]:
    """Primers read from .csv/.tsv text or pasted columns, and the rows left out (each with
    its ``name``, ``sequence`` and the ``reason``)."""
    if not text.strip():
        return [], []
    lib = parse_text(text)
    return lib.primers, lib.rejected


def _write_row(cells: list[str], delimiter: str, newline: str) -> str:
    out = io.StringIO()
    csv.writer(out, delimiter=delimiter, lineterminator=newline).writerow(cells)
    return out.getvalue()


#: The fields of a primer the dialog can change.
EDITABLE = ("name", "notes", "location", "added")


def merged_text(current: LibraryText | None, added: list[Primer],
                edits: dict[str, dict[str, str]] | None = None,
                removed: set[str] = frozenset()) -> tuple[str, dict[str, int]]:
    """The file's new text after a save from the dialog. The file is read again at save time
    and only the dialog's own changes are applied to it, so two people editing it never undo
    each other's work:

    - ``edits`` (sequence → the fields changed in the dialog) change only those cells of that
      primer's row; its other cells, any extra columns and the sequence as written stay;
    - a row whose sequence is in ``removed`` (deleted in the dialog) goes;
    - ``added`` primers not in the file are added at the end (one already there is left as is);
    - every other row stays exactly as it is, byte for byte: primers a colleague added, edited
      or kept, and rows Vettore could not read. A primer a colleague deleted stays deleted.

    Returns the text and counts: ``count`` primers written (added or changed), ``kept``
    unreadable rows kept, and ``gone`` edits not saved because someone had deleted that primer
    meanwhile.
    """
    edits = edits or {}
    lib = current or LibraryText(delimiter=";")
    nl = lib.newline
    index = dict(lib.index)
    header = list(lib.header.cells) if lib.header else None
    if current is None or (header is None and not any(r.cells for r in lib.rows)):
        header = list(COLUMNS)
        index = {name: pos for pos, name in enumerate(COLUMNS)}
    in_file = {row.key for row in lib.rows if row.key is not None}
    fresh = [p for p in added if p.sequence not in in_file]
    if header is not None:                # a column the file lacks, but the dialog now fills
        width = max([len(header)] + [len(r.cells) for r in lib.rows])
        for column in COLUMNS:
            needed = any(getattr(p, column) for p in fresh) or any(
                e.get(column) for e in edits.values())
            if column not in index and needed:
                index[column] = width     # past every cell of every row: nothing is overwritten
                header += [""] * (width - len(header)) + [column]
                width += 1

    def apply(cells: list[str], values: dict[str, str]) -> list[str]:
        cells = list(cells)
        for column, value in values.items():
            pos = index.get(column)
            if pos is None:
                continue
            now = cells[pos].strip() if pos < len(cells) else ""
            if now == value or (column == "name" and value == "(unnamed)" and not now):
                continue                  # unchanged: the cell stays as written
            cells += [""] * (pos + 1 - len(cells))
            cells[pos] = value
        return cells

    out = list(lib.preamble)
    if header is not None:
        same = lib.header is not None and header == lib.header.cells
        out.append(lib.header.raw if same else _write_row(header, lib.delimiter, nl))
    counts = {"count": 0, "kept": len(lib.rejected), "gone": 0}
    seen_edits = set()
    for row in lib.rows:
        if row.key is not None and row.key in removed:
            continue
        if row.key is not None and row.key in edits:
            seen_edits.add(row.key)
            cells = apply(row.cells, {k: v for k, v in edits[row.key].items() if k in EDITABLE})
            if cells != row.cells:
                out.append(_write_row(cells, lib.delimiter, nl))
                counts["count"] += 1
                continue
        out.append(row.raw)
    counts["gone"] = len(set(edits) - seen_edits - removed)
    for p in fresh:
        values = {column: getattr(p, column) for column in COLUMNS}
        out.append(_write_row(apply([], values), lib.delimiter, nl))
        counts["count"] += 1
    for i, line in enumerate(out[:-1]):    # a last line without a newline, now followed by more
        if not line.endswith(("\n", "\r")):
            out[i] = line + nl
    return "".join(out), counts


def _decode(data: bytes) -> tuple[str, str]:
    """The file's text and how to write it back: UTF-16 (Excel's "Unicode text"), UTF-8 (with
    or without a byte-order mark), else the Windows code page Excel uses for .csv on Western
    European systems, line by line (a UTF-8 file with a line added in Notepad keeps both).
    Workbooks and other binary files are refused, never rewritten."""
    if data.startswith((b"\xff\xfe", b"\xfe\xff")):
        try:
            return data.decode("utf-16"), "utf-16"
        except UnicodeDecodeError as exc:
            raise PrimerLibraryError("it is not a text file.") from exc
    if data.startswith((b"PK\x03\x04", b"\xd0\xcf\x11\xe0")):
        raise PrimerLibraryError(
            "it is an Excel workbook, not a .csv. In Excel choose File → Save As → CSV, "
            "then choose that file.")
    if b"\x00" in data:
        raise PrimerLibraryError("it is not a text file.")
    body = data.removeprefix(_BOM)
    try:
        return body.decode("utf-8"), "utf-8-sig"
    except UnicodeDecodeError:
        pass

    def line(raw: bytes) -> str:
        for encoding in ("utf-8", "cp1252"):
            try:
                return raw.decode(encoding)
            except UnicodeDecodeError:
                continue
        return raw.decode("latin-1")      # decodes every byte

    return "".join(line(raw) for raw in body.splitlines(keepends=True)), "utf-8-sig"


def _read_file() -> tuple[LibraryText | None, float | None, str]:
    """The library file as read, its modification time and the encoding to write it back in;
    (None, None, "utf-8-sig") when there is no file yet. Raises PrimerLibraryError when the
    file cannot be read or is not a primer list (so it is never overwritten)."""
    path = settings.primer_library_path()
    if not path.is_file():
        return None, None, "utf-8-sig"
    try:
        data, mtime = path.read_bytes(), path.stat().st_mtime
    except OSError as exc:
        log.error("Could not read the primer library %s: %r", path, exc)
        raise PrimerLibraryError(
            f"The primer library file “{path.name}” could not be read: {_os_reason(exc)}.") from exc
    try:
        text, encoding = _decode(data)
        lib = parse_text(text)
    except PrimerLibraryError as exc:
        raise PrimerLibraryError(
            f"“{path.name}” cannot be used as the primer library: {exc} It is left as it is."
        ) from exc
    if lib.header is None and not lib.primers and any(r.cells for r in lib.rows):
        raise PrimerLibraryError(
            f"“{path.name}” is not a primer library: no row has a name and a DNA sequence, and "
            "there is no header naming those columns. Choose another file; this one is left "
            "as it is.")
    return lib, mtime, encoding


def check_library_path(raw: str) -> None:
    """Refuse a library path that is not a spreadsheet saved as text."""
    if not raw.lower().endswith(LIBRARY_SUFFIXES):
        raise PrimerLibraryError(
            "The primer library is a .csv file (a spreadsheet saved as CSV). Choose or name a "
            ".csv file; .tsv and .txt work as well.")


def read_library() -> dict[str, Any]:
    """The library at the chosen path: its primers, the rows that could not be read, the path,
    its separator and its modification time (so the window can tell when it changed)."""
    path = settings.primer_library_path()
    lib, mtime, _ = _read_file()
    if lib is None:
        return {"path": str(path), "exists": False, "mtime": None, "delimiter": ";",
                "primers": [], "rejected": []}
    return {"path": str(path), "exists": True, "mtime": mtime, "delimiter": lib.delimiter,
            "primers": [p.as_dict() for p in lib.primers], "rejected": lib.rejected}


def write_library(added: list[Primer], edits: dict[str, dict[str, str]] | None = None,
                  removed: list[str] | tuple[str, ...] = ()) -> dict[str, Any]:
    """Apply the dialog's changes to the shared file as it is now (see ``merged_text``): a
    colleague's work and the rows Vettore could not read stay, and the file keeps its
    separator, header and extra columns. Returns the library as read back, with the counts."""
    path = settings.primer_library_path()
    check_library_path(str(path))
    current, _, encoding = _read_file()
    gone = set()
    for raw in removed:
        try:
            gone.add(clean_sequence(raw))
        except PrimerLibraryError:
            continue
    text, counts = merged_text(current, added, edits, gone)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(f".{path.name}.vettore-tmp")
        tmp.write_text(text, encoding=encoding, newline="")
        tmp.replace(path)
    except OSError as exc:
        log.error("Could not write the primer library %s: %r", path, exc)
        raise PrimerLibraryError(
            f"The primer library file “{path.name}” could not be saved: {_os_reason(exc)}. "
            "Choose another location for the library in the Primer Library dialog.") from exc
    return {**read_library(), **counts}


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
    inside the footprint. ``positions`` > 1: the site is a repeat where the primer anneals,
    overlapping itself, with its 3′ end at that many positions (it would prime a ladder)."""

    start: int
    end: int
    strand: int
    footprint: int
    tail: int
    mismatches: int
    positions: int = 1

    @property
    def three_prime(self) -> int:
        return self.end if self.strand == 1 else self.start

    @property
    def five_prime(self) -> int:
        """The 5′ base of the annealing footprint (the tail, if any, lies beyond it)."""
        return self.start if self.strand == 1 else self.end

    def as_dict(self) -> dict[str, Any]:
        return {"start": self.start, "end": self.end, "strand": self.strand,
                "footprint": self.footprint, "tail": self.tail, "mismatches": self.mismatches,
                "positions": self.positions}


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
    mismatches), or None when the 3′ anchor does not pair. One pass from the anchor toward the
    5′ end: the footprint reaches the farthest matching base with at most MAX_MISMATCHES
    mismatches between it and the anchor."""
    body_end = len(primer) - ANCHOR
    if not all(_pairs(p, t) for p, t in zip(primer[body_end:], window[body_end:], strict=True)):
        return None
    start, best_mismatches, mismatches = body_end, 0, 0
    for i in range(body_end - 1, -1, -1):
        if _pairs(primer[i], window[i]):
            start, best_mismatches = i, mismatches   # a footprint always begins on a match
        else:
            mismatches += 1
            if mismatches > MAX_MISMATCHES:
                break
    return len(primer) - start, best_mismatches


def _anchor_pattern(probe: str) -> re.Pattern[str]:
    """Template stretches the anchor pairs with: each IUPAC code becomes the set of template
    letters it covers (an ambiguous template letter only where the code covers all of it)."""
    parts = []
    for code in probe:
        letters = "".join(t for t in _CODES if all(b in _CODES[code] for b in _CODES[t]))
        parts.append(f"[{letters}]" if len(letters) > 1 else letters)
    return re.compile("(?=(" + "".join(parts) + "))")


def _candidates(primer: str, top: str, n: int, circular: bool) -> list[tuple[int, Binding]]:
    """Every anchor hit on both strands, scored, with the footprint's low coordinate before it
    is wrapped around the origin (so overlapping hits can be found in one sweep)."""
    anchor = ANCHOR
    length = len(primer)
    text = top + top[:anchor - 1] if circular else top

    def base(i: int) -> str | None:
        if circular:
            return top[i % n]
        return top[i] if 0 <= i < n else None

    found: list[tuple[int, Binding]] = []
    # Forward: the anchor pairs with the top strand, the primer's 3′ end at the anchor's end.
    for m in _anchor_pattern(primer[-anchor:]).finditer(text):
        three = m.start() + anchor - 1
        window = [base(three - length + 1 + k) for k in range(length)]
        if (scored := _anneal(primer, window)) is None:
            continue
        footprint, mism = scored
        found.append((three - footprint + 1, Binding((three - footprint + 1) % n, three % n, 1,
                                                     footprint, length - footprint, mism)))
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
        found.append((three, Binding(three % n, (three + footprint - 1) % n, -1, footprint,
                                     length - footprint, mism)))
    return found


def _better(a: Binding, b: Binding) -> Binding:
    return min(a, b, key=lambda x: (x.mismatches, -x.footprint, x.start))


def find_bindings(primer: str, sequence: str, circular: bool) -> list[Binding]:
    """Every place the primer anneals on either strand (see the module docstring for the rule).
    Hits of one primer on one strand whose footprints overlap (a primer inside a repeat) are one
    site, the best of them, with ``positions`` saying how many 3′ ends it has there. Sorted by position, then strand. Primers shorter than MIN_LENGTH,
    and templates shorter than the primer, give none. Linear in the template's length."""
    primer = clean_sequence(primer)
    top = sequence.upper()
    n = len(top)
    if len(primer) < MIN_LENGTH or n < len(primer):
        return []
    candidates = _candidates(primer, top, n, circular)
    kept: list[Binding] = []
    for strand in (1, -1):
        hits = sorted((h for h in candidates if h[1].strand == strand), key=lambda h: h[0])
        clusters: list[list] = []          # [low, high, best, hits] in unwrapped coordinates
        for low, b in hits:
            high = low + b.footprint - 1
            if clusters and low <= clusters[-1][1]:
                clusters[-1][1] = max(clusters[-1][1], high)
                clusters[-1][2] = _better(clusters[-1][2], b)
                clusters[-1][3] += 1
            else:
                clusters.append([low, high, b, 1])
        if circular and len(clusters) > 1 and clusters[-1][1] >= clusters[0][0] + n:
            last = clusters.pop()          # a repeat running across the origin
            clusters[0][2] = _better(clusters[0][2], last[2])
            clusters[0][3] += last[3]
        kept += [replace(best, positions=count) for _, _, best, count in clusters]
    return sorted(kept, key=lambda b: (b.start, b.strand))


def _span(low: int, high: int, n: int) -> int:
    """Length of the inclusive span from ``low`` to ``high`` on a circle of length n."""
    return (high - low) % n + 1


def library_on_sequence(primers: list[Primer], sequence: str, circular: bool) -> list[dict[str, Any]]:
    """Every library primer that binds the sequence, with its binding sites (the map overlay and
    the Primer Library dialog use this). ``sites`` counts every 3′ position, those inside a
    repeat included; a primer with more than MAX_SITES is listed with no bindings, so the map
    is not flooded."""
    seq = sequence.upper()
    found = []
    for p in primers:
        try:
            bindings = find_bindings(p.sequence, seq, circular)
        except PrimerLibraryError:
            continue  # a malformed stored sequence never breaks the map
        if not bindings:
            continue
        sites = sum(b.positions for b in bindings)
        entry = {**p.as_dict(), "sites": sites}
        entry["bindings"] = [] if sites > MAX_SITES else [b.as_dict() for b in bindings]
        found.append(entry)
    return found


def _bindings_by_strand(primers: list[Primer], seq: str, circular: bool):
    """Forward and reverse sites of the primers that bind specifically: a primer with more than
    MAX_SITES sites, or one in a repeat (it would prime a ladder of products or mixed reads),
    is not offered for products or sequencing."""
    forwards, reverses = [], []
    for p in primers:
        try:
            bindings = find_bindings(p.sequence, seq, circular)
        except PrimerLibraryError:
            continue
        if len(bindings) > MAX_SITES or any(b.positions > 1 for b in bindings):
            continue
        for b in bindings:
            (forwards if b.strand == 1 else reverses).append((p, b))
    return forwards, reverses


def amplified_products(primers: list[Primer], sequence: str, circular: bool,
                       max_product: int = 10000) -> list[dict[str, Any]]:
    """Pairs of library primers that make a PCR product: a forward and a reverse binding, each
    extending over the other's site (the reverse footprint ends past the forward's and starts
    past the forward's 5′ base; overlapping footprints still make a product), MIN_PRODUCT to
    ``max_product`` bp apart. The amplicon runs from the forward footprint's 5′ base to the
    reverse's; its size includes both tails. Smallest product first."""
    seq = sequence.upper()
    n = len(seq)
    forwards, reverses = _bindings_by_strand(primers, seq, circular)
    products = []
    for fp, fb in forwards:
        for rp, rb in reverses:
            if circular:
                to_r3, to_r5 = (rb.start - fb.start) % n, (rb.end - fb.start) % n
                if not 0 < to_r3 <= to_r5 or to_r5 <= (fb.end - fb.start) % n:
                    continue
                length = to_r5 + 1
            else:
                if not (rb.start > fb.start and rb.end > fb.end):
                    continue
                length = rb.end - fb.start + 1
            if not MIN_PRODUCT <= length <= max_product:
                continue
            products.append({
                "forward": fp.name, "reverse": rp.name,
                "forward_sequence": fp.sequence, "reverse_sequence": rp.sequence,
                "size": length + fb.tail + rb.tail,
                "product_start": fb.start, "product_end": rb.end,
                "forward_three_prime": fb.end, "reverse_three_prime": rb.start,
            })
    products.sort(key=lambda x: (x["size"], x["forward"], x["reverse"]))
    return products


def product_over_region(primers: list[Primer], sequence: str, circular: bool,
                        start: int, end: int, min_size: int = MIN_PRODUCT,
                        max_size: int = 10000) -> dict[str, Any] | None:
    """The smallest library-primer product of ``min_size``–``max_size`` bp (tails included)
    whose template between the two primers holds the whole region [start, end] (0-based
    inclusive; on a circular sequence ``start > end`` crosses the origin), or None. Bases under
    a primer come from the primer, not the template, so the region must lie between the 3′
    ends. Used to tell the user "you already have F and R"."""
    n = len(sequence)
    if n == 0:
        return None
    for product in amplified_products(primers, sequence, circular, max_product=max_size):
        if not min_size <= product["size"] <= max_size:
            continue
        f3, r3 = product["forward_three_prime"], product["reverse_three_prime"]
        if circular:
            to_start, to_end, to_r3 = (start - f3) % n, (end - f3) % n, (r3 - f3) % n
            if 0 < to_start <= to_end < to_r3:
                return product
        elif f3 < start <= end < r3:
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
        if distance >= SEQ_PRIMER_MIN and distance + region - 1 <= SEQ_READ_LENGTH:
            out.append({"name": p.name, "sequence": p.sequence, "strand": b.strand,
                        "start": b.start, "end": b.end, "reads_from": distance})
    out.sort(key=lambda x: (x["strand"] != 1, x["reads_from"]))
    return out



# ---- HTTP --------------------------------------------------------------------------------

class SetPath(BaseModel):
    path: str


class SaveLibrary(BaseModel):
    """What the dialog changed since it read the file: primers it added, the fields it changed
    on primers it had (``edited``: each with its ``sequence`` and only the changed fields), and
    the sequences it removed. Only those are applied to the file as it is at save time."""

    added: list[dict[str, Any]] = []
    edited: list[dict[str, Any]] = []
    removed: list[str] = []


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


class ProductRequest(SequencingRequest):
    """A region, and the product sizes primer design was asked for (tails included)."""

    product_min: int | None = None
    product_max: int | None = None


def _check_region(sequence: str, circular: bool, start: int, end: int) -> None:
    n = len(sequence)
    if not (0 <= start < n and 0 <= end < n):
        raise HTTPException(422, f"The region must lie within the sequence (1–{n:,}).")
    if start > end and not circular:
        raise HTTPException(422, "The region's end comes before its start; on a linear "
                                 "sequence a region cannot cross the origin.")


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
    try:
        check_library_path(raw)
    except PrimerLibraryError as exc:
        raise HTTPException(422, str(exc)) from exc
    previous = settings.primer_library_path()
    settings.remember_primer_library_path(raw)
    try:
        return read_library()
    except PrimerLibraryError as exc:
        settings.remember_primer_library_path(str(previous))   # keep the library that worked
        raise HTTPException(422, str(exc)) from exc


@router.post("/save")
def save_library(body: SaveLibrary) -> dict[str, Any]:
    added = _primers_from_payload(body.added)
    edits: dict[str, dict[str, str]] = {}
    for row in body.edited:
        try:
            sequence = clean_sequence(str(row.get("sequence", "")))
        except PrimerLibraryError as exc:
            raise HTTPException(422, f"“{row.get('name') or 'A primer'}”: {exc}") from exc
        edits[sequence] = {k: str(row[k] or "").strip() for k in EDITABLE if k in row}
    try:
        return write_library(added, edits, body.removed)
    except PrimerLibraryError as exc:
        raise HTTPException(422, str(exc)) from exc


@router.post("/parse")
def parse(body: ImportText) -> dict[str, Any]:
    """Parse pasted or file text without saving it, so the dialog can show what will be imported
    and what will be left out before the user keeps it."""
    try:
        primers, rejected = parse_csv(body.text)
    except PrimerLibraryError as exc:
        raise HTTPException(422, f"The text could not be read: {exc}") from exc
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
def product(body: ProductRequest) -> dict[str, Any]:
    """The smallest pair of library primers that already amplifies a product covering the region
    [start, end], or ``null`` — so primer design can say "you already have F and R"."""
    sequence, circular = _sequence_of(body.sequenceData)
    if not sequence:
        return {"product": None}
    _check_region(sequence, circular, body.start, body.end)
    try:
        library = read_library()
    except PrimerLibraryError as exc:
        raise HTTPException(422, str(exc)) from exc
    primers = _primers_from_payload(library["primers"])
    hit = product_over_region(primers, sequence, circular, body.start, body.end,
                              body.product_min or MIN_PRODUCT, body.product_max or 10000)
    return {"product": hit}


@router.post("/sequencing")
def sequencing(body: SequencingRequest) -> dict[str, Any]:
    sequence, circular = _sequence_of(body.sequenceData)
    if not sequence:
        raise HTTPException(422, "Open a sequence first.")
    _check_region(sequence, circular, body.start, body.end)
    try:
        library = read_library()
    except PrimerLibraryError as exc:
        raise HTTPException(422, str(exc)) from exc
    primers = _primers_from_payload(library["primers"])
    return {"primers": sequencing_primers(primers, sequence, circular, body.start, body.end)}
