"""Reading and writing sequence files.

Opening is forgiving, because files exchanged between labs often carry the wrong extension or
none at all: the format is recognised from the content first (GenBank, EMBL, FASTA, plain
sequence text, and the binary .dna, .xdna, .gck and .ab1 formats by their headers) and from the
extension only when the content does not say. Every error raised here is a sentence a biologist
can act on; parser details go to the log file.

Saving writes GenBank (the default exchange format), FASTA or .dna, chosen by the extension.
"""

from __future__ import annotations

import html
import io
import logging
import re
import struct
import uuid
import warnings
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq, UndefinedSequenceError, reverse_complement
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from .convert import _as_list, ove_to_record, record_to_ove

log = logging.getLogger(__name__)

# Biopython's SeqIO name for the binary .dna format is fixed by that library.
DNA_FILE_FORMAT = "snapgene"
# Vettore's own reader for plain sequence text (raw or numbered sequence, FASTA without header).
RAW = "raw"

# Extension -> format. ``.txt``/``.seq`` say only "text": the content decides.
FORMATS: dict[str, str] = {
    ".dna": DNA_FILE_FORMAT,
    ".gb": "genbank",
    ".gbk": "genbank",
    ".gbff": "genbank",
    ".genbank": "genbank",
    ".ape": "genbank",
    ".embl": "embl",
    ".emb": "embl",
    ".fa": "fasta",
    ".fasta": "fasta",
    ".fna": "fasta",
    ".fas": "fasta",
    ".fsa": "fasta",
    ".xdna": "xdna",
    ".gck": "gck",
    ".ab1": "abi",
    ".abi": "abi",
    ".txt": RAW,
    ".seq": RAW,
}
BINARY_FORMATS = {DNA_FILE_FORMAT, "abi", "xdna", "gck"}

LABELS = {
    "genbank": "GenBank",
    "embl": "EMBL",
    "fasta": "FASTA",
    DNA_FILE_FORMAT: ".dna",
    "xdna": ".xdna",
    "gck": ".gck",
    "abi": "Sanger trace (.ab1)",
    RAW: "plain sequence text",
}

# Formats Vettore writes, by extension; the first extension of each format is its default.
SAVE_FORMATS: dict[str, str] = {
    ".gb": "genbank",
    ".gbk": "genbank",
    ".gbff": "genbank",
    ".genbank": "genbank",
    ".ape": "genbank",
    ".fasta": "fasta",
    ".fa": "fasta",
    ".fna": "fasta",
    ".fas": "fasta",
    ".fsa": "fasta",
    ".dna": DNA_FILE_FORMAT,
}
DEFAULT_EXTENSION = {"genbank": ".gb", "fasta": ".fasta", DNA_FILE_FORMAT: ".dna"}
# Formats whose files Vettore can overwrite without losing anything it did not read.
IN_PLACE_FORMATS = {"genbank", "fasta"}

SUPPORTED_SUMMARY = (
    "Vettore opens GenBank (.gb, .gbk, .ape), EMBL (.embl), FASTA (.fa, .fasta), .dna, .xdna, "
    ".gck, Sanger traces (.ab1) and plain sequence text (.txt, .seq)."
)

# IUPAC nucleotide codes (Cornish-Bowden 1985, NC-IUB recommendations), plus U for RNA.
DNA_LETTERS = frozenset("ACGTURYSWKMBDHVN")
# Letters that occur in protein sequences but are not nucleotide codes.
PROTEIN_ONLY = frozenset("EFILPQXZJO*")


class UnsupportedFile(ValueError):
    """The file cannot be opened; the message is shown to the user as is."""


class SaveError(ValueError):
    """The document cannot be saved as asked; the message is shown to the user as is."""


@dataclass
class ReadResult:
    documents: list[dict[str, Any]]
    format: str
    warnings: list[str] = field(default_factory=list)
    # True when Save may overwrite the file: one record, in a format Vettore writes, whose
    # extension matches its content.
    writable: bool = False

    @property
    def label(self) -> str:
        return LABELS.get(self.format, self.format)

    def as_dict(self) -> dict[str, Any]:
        return {"documents": self.documents, "format": self.format, "format_label": self.label,
                "warnings": self.warnings, "writable": self.writable}


# ---- format detection -------------------------------------------------------------------

def _decode(content: bytes) -> str:
    """Text of a file in any encoding labs produce: UTF-8 (with or without BOM), UTF-16 from
    Windows editors, or a legacy 8-bit code page. Line endings become ``\\n``."""
    if content.startswith((b"\xff\xfe", b"\xfe\xff")):
        text = content.decode("utf-16", errors="replace")
    else:
        try:
            text = content.decode("utf-8-sig")
        except UnicodeDecodeError:
            # Windows tools write cp1252 (e.g. a µ or ° in a note); latin-1 maps every byte.
            try:
                text = content.decode("cp1252")
            except UnicodeDecodeError:
                text = content.decode("latin-1")
    return text.replace("\r\n", "\n").replace("\r", "\n")


def _is_text(content: bytes) -> bool:
    head = content[:4096]
    return head.startswith((b"\xff\xfe", b"\xfe\xff")) or b"\x00" not in head


def _looks_like_xdna(content: bytes) -> bool:
    # 112-byte header: version 0, sequence type 0-4, topology 0-1, sequence length at 28,
    # comment length at 96, then the sequence (Bio.SeqIO.XdnaIO).
    if len(content) < 112 or content[0] != 0 or content[1] > 4 or content[2] > 1:
        return False
    (seq_len,) = struct.unpack(">I", content[28:32])
    (comment_len,) = struct.unpack(">I", content[96:100])
    if not seq_len or 112 + seq_len + comment_len > len(content):
        return False
    return content[112:112 + min(seq_len, 64)].isalpha()


def _looks_like_gck(content: bytes) -> bool:
    # 24-byte header whose 4th and 8th bytes are 12, then a length-prefixed packet holding the
    # sequence length and the sequence (Bio.SeqIO.GckIO).
    if len(content) < 32 or content[:8] != b"\x00\x00\x00\x0c\x00\x00\x00\x0c":
        return False
    packet_len, seq_len = struct.unpack(">II", content[24:32])
    if not seq_len or seq_len > packet_len - 4 or 28 + packet_len > len(content):
        return False
    return content[32:32 + min(seq_len, 64)].isalpha()


def sniff(content: bytes) -> str | None:
    """The format the content itself declares, or None for unknown binary data."""
    if content.startswith(b"\x09\x00\x00\x00\x0e"):  # cookie packet: type 9, 14 bytes
        return DNA_FILE_FORMAT
    if content.startswith(b"ABIF"):
        return "abi"
    if _looks_like_xdna(content):
        return "xdna"
    if _looks_like_gck(content):
        return "gck"
    if not _is_text(content):
        return None
    lines = [line for line in _decode(content[:8192]).split("\n") if line.strip()]
    if not lines:
        return RAW
    if any(line.startswith("LOCUS ") for line in lines[:50]):
        return "genbank"
    if re.match(r"ID {3}\S", lines[0]):
        return "embl"
    if lines[0].lstrip().startswith(">"):
        return "fasta"
    return RAW


def detect_format(filename: str, content: bytes = b"") -> str:
    """Format of a file: the content decides, the extension is the fallback."""
    by_extension = FORMATS.get(Path(filename).suffix.lower())
    found = sniff(content) if content else None
    if found == RAW and by_extension in BINARY_FORMATS:
        # Text inside a .dna/.ab1/... file: report it as a damaged file of that kind.
        found = None
    if found:
        if by_extension and found != by_extension:
            log.info("%s: content is %s, extension says %s", filename, found, by_extension)
        return found
    if by_extension:
        return by_extension
    raise UnsupportedFile(f"“{Path(filename).name}” is not a sequence file Vettore can read. "
                          f"{SUPPORTED_SUMMARY}")


# ---- reading ----------------------------------------------------------------------------

def _parse_raw(text: str, display: str) -> list[SeqRecord]:
    """A plain sequence: letters, with spaces, line breaks and position numbers ignored. Also
    accepts a GenBank ORIGIN block pasted on its own."""
    bases: list[str] = []
    lines = text.split("\n")
    for lineno, line in enumerate(lines, 1):
        if line.strip() in {"//", "ORIGIN"}:
            continue
        for col, ch in enumerate(line, 1):
            if ch.isspace() or ch.isdigit():
                continue
            if ch.upper() not in DNA_LETTERS:
                hint = ""
                if lineno == 1 and len(lines) > 1:
                    hint = (" If the first line is the sequence name, start it with “>” "
                            "to make it a FASTA record.")
                raise UnsupportedFile(
                    f"{display} could not be read as a DNA sequence: line {lineno}, column "
                    f"{col} has “{ch}”, which is not a DNA letter. A plain sequence may contain "
                    "only A, C, G, T (or U) and the ambiguity codes N, R, Y, S, W, K, M, B, D, "
                    f"H, V; spaces and numbers are ignored.{hint}")
            bases.append(ch)
    if not bases:
        raise UnsupportedFile(f"{display} contains no DNA sequence.")

    record = SeqRecord(Seq("".join(bases)), id="", name="", description="")
    record.annotations["molecule_type"] = "DNA"
    return [record]


def _prepare_text(text: str, fmt: str) -> str:
    if fmt == "fasta":
        # Biopython's FASTA reader rejects anything before the first record.
        start = text.find(">")
        return text[start:] if start > 0 else text
    if fmt in {"genbank", "embl"}:
        stripped = text.rstrip()
        if stripped and not stripped.endswith("//"):
            text = stripped + "\n//\n"  # the record terminator is often lost when copying
    return text


def _parse(content: bytes, fmt: str, display: str) -> list[SeqRecord]:
    if fmt == RAW:
        return _parse_raw(_decode(content), display)
    if fmt in BINARY_FORMATS:
        handle: io.IOBase = io.BytesIO(content)
    else:
        handle = io.StringIO(_prepare_text(_decode(content), fmt))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        records = list(SeqIO.parse(handle, fmt))
    for w in caught:
        log.info("%s: parser warning: %s", display, w.message)
    return records


def _sequence_text(record: SeqRecord) -> str | None:
    try:
        return str(record.seq)
    except UndefinedSequenceError:  # header-only record, e.g. a contig
        return None


def _record_label(record: SeqRecord) -> str:
    name = _record_name(record)
    return f"“{name}”" if name else "The sequence"


_NOT_DNA = re.compile(f"[^{''.join(sorted(DNA_LETTERS))}{''.join(sorted(DNA_LETTERS)).lower()}]")


def _check_letters(record: SeqRecord, sequence: str) -> str | None:
    """A sentence naming the first character that is not a DNA letter, or None."""
    bad = _NOT_DNA.search(sequence)
    if not bad:
        return None
    ch = bad.group()
    hint = ""
    if ch in "-.":
        hint = " Remove the alignment gaps first."
    elif sum(c.upper() in PROTEIN_ONLY for c in sequence) > len(sequence) / 20:
        hint = " It looks like a protein sequence; Vettore opens DNA sequences."
    return (f"{_record_label(record)} has “{ch}” at position {bad.start() + 1}, which is not a "
            f"DNA letter (A, C, G, T, U or an IUPAC ambiguity code).{hint}")


def _short(name: str, limit: int = 40) -> str:
    return name if len(name) <= limit else name[: limit - 1] + "…"


def _clean_features(record: SeqRecord) -> tuple[list[str], list[str]]:
    """Drop or trim features the editor cannot draw; return the (left out, shortened) ones."""
    from .convert import feature_name

    length = len(record.seq)
    dropped: list[str] = []
    trimmed: list[str] = []
    kept: list[SeqFeature] = []
    for feature in record.features:
        label = f"“{_short(feature_name(feature))}”"
        report = feature.type != "source"
        location = feature.location
        problem = None
        if location is None:
            problem = "location not understood"
        elif any(part.ref for part in location.parts):
            problem = "it refers to another sequence"
        elif any(int(part.end) <= int(part.start) for part in location.parts):
            problem = "between two bases"
        elif all(int(part.start) >= length for part in location.parts):
            problem = "beyond the end of the sequence"
        if problem:
            if report:
                dropped.append(f"{label} ({problem})")
            continue
        parts = list(location.parts)
        if any(int(p.end) > length for p in parts):
            inside = [SimpleLocation(int(p.start), min(int(p.end), length), p.strand)
                      for p in parts if int(p.start) < length]
            feature.location = inside[0] if len(inside) == 1 else sum(inside[1:], inside[0])
            if report:
                trimmed.append(label)
        kept.append(feature)
    record.features = kept
    return dropped, trimmed


def _feature_warnings(display: str, dropped: list[str], trimmed: list[str]) -> list[str]:
    out = []
    if dropped:
        noun = "feature" if len(dropped) == 1 else "features"
        out.append(f"{display}: {len(dropped)} {noun} could not be placed on the sequence and "
                   f"{'was' if len(dropped) == 1 else 'were'} left out: {', '.join(dropped)}.")
    if trimmed:
        noun = "feature ran" if len(trimmed) == 1 else "features ran"
        out.append(f"{display}: {len(trimmed)} {noun} past the end of the sequence and "
                   f"{'was' if len(trimmed) == 1 else 'were'} shortened to fit: "
                   f"{', '.join(trimmed)}.")
    return out


_PLACEHOLDER_NAMES = {"", ".", "<unknown name>", "<unknown id>"}


def _record_name(record: SeqRecord) -> str:
    for candidate in (record.name, record.id):
        if candidate and candidate not in _PLACEHOLDER_NAMES:
            return candidate
    return ""


def read_file(filename: str, content: bytes, *, pasted: bool = False) -> ReadResult:
    """Parse a file (or pasted text, when ``pasted``) into OVE ``sequenceData`` documents.

    Names: a single record takes the file name (what the user double-clicked); several records
    take their own names, made distinct. Pasted text takes the record's own name when it has
    one, else ``filename``.
    """
    display = "The pasted text" if pasted else f"“{Path(filename).name}”"
    stem = filename if pasted else Path(filename).stem
    if not content.strip():
        raise UnsupportedFile(f"{display} is empty.")
    fmt = detect_format("" if pasted else filename, content)
    try:
        records = _parse(content, fmt, display)
    except UnsupportedFile:
        raise
    except Exception as exc:
        log.exception("Could not parse %s as %s", filename, fmt)
        advice = ("Ask for the trace file again; it may have been damaged in transfer."
                  if fmt == "abi" else
                  "If a colleague sent it, ask them to export it again as GenBank.")
        raise UnsupportedFile(
            f"{display} looks like a {LABELS.get(fmt, fmt)} file but could not be read: it may "
            f"be damaged or cut short. {advice}") from exc

    problems: list[str] = []
    usable: list[tuple[SeqRecord, str]] = []
    for record in records:
        if record.annotations.get("molecule_type") == "protein":
            problems.append(f"{_record_label(record)} is a protein sequence; Vettore opens DNA "
                            "sequences.")
            continue
        sequence = _sequence_text(record)
        if not sequence:
            problems.append(f"{_record_label(record)} has no sequence (it may be a header-only "
                            "record).")
            continue
        if bad := _check_letters(record, sequence):
            problems.append(bad)
            continue
        usable.append((record, sequence))
    if not usable:
        detail = problems[0] if problems else "No DNA sequence was found in it."
        raise UnsupportedFile(f"{display} could not be opened. {detail}")

    result = ReadResult(documents=[], format=fmt)
    if problems:
        result.warnings.append(f"{display}: {len(problems)} of {len(records)} sequences could "
                               f"not be opened. {' '.join(problems)}")
    names: dict[str, int] = {}
    for i, (record, _sequence) in enumerate(usable):
        dropped, trimmed = _clean_features(record)
        result.warnings += _feature_warnings(display, dropped, trimmed)
        own = _record_name(record)
        if len(records) == 1:
            name = (own or stem or "Pasted sequence") if pasted else (stem or own)
        else:
            name = own or f"{stem} ({i + 1})"
        if name in names:  # several records sharing a name get distinct tabs
            names[name] += 1
            name = f"{name} ({names[name]})"
        else:
            names[name] = 1
        doc = record_to_ove(record, fallback_name=name)
        doc["name"] = name
        if fmt == DNA_FILE_FORMAT:
            _enrich_dna_file(doc, content)
        result.documents.append(doc)
    # Overwriting is safe only when the document holds everything the file held.
    by_extension = FORMATS.get(Path(filename).suffix.lower())
    result.writable = (not pasted and len(records) == 1 and not result.warnings
                       and fmt in IN_PLACE_FORMATS and by_extension == fmt)
    return result


def read_documents(filename: str, content: bytes) -> list[dict]:
    """Parse a file into OVE ``sequenceData`` documents (one per record)."""
    return read_file(filename, content).documents


def read_records(filename: str, content: bytes) -> list[SeqRecord]:
    """Parse a file into Biopython records (no conversion, no clean-up)."""
    fmt = detect_format(filename, content)
    records = _parse(content, fmt, f"“{Path(filename).name}”")
    if not records:
        raise UnsupportedFile(f"No sequence found in {Path(filename).name}")
    return records


_NEAR_WHITE = re.compile(r"^#?(f[a-f0-9]){3}$", re.IGNORECASE)


def _enrich_dna_file(doc: dict, content: bytes) -> None:
    """Add what Biopython does not read from .dna files: feature colours, full primer
    sequences (with 5' tails) and the description. Locations stay Biopython's, which already
    applies the file's primer hybridization parameters."""
    try:
        from sgffp import SgffReader

        sg = SgffReader.from_bytes(content)
    except Exception:  # enrichment is best effort; the Biopython parse already succeeded
        log.exception("sgffp could not read the .dna file; colours and primer tails skipped")
        return

    colors = {}
    for feature in sg.features:
        color = feature.color or (feature.segments[0].color if feature.segments else None)
        if feature.segments and color and not _NEAR_WHITE.match(color):
            colors[(feature.name, feature.segments[0].start)] = color
    for ann in doc["features"]:
        start = ann["locations"][0]["start"] if ann.get("locations") else ann["start"]
        if color := colors.get((ann["name"], start)):
            ann["color"] = color

    sequences = {p.name: p.sequence for p in sg.primers if p.sequence}
    for ann in doc["primers"]:
        if bases := sequences.get(ann["name"]):
            ann["bases"] = bases.upper()
            ann.setdefault("notes", {})["primer_sequence"] = [bases.upper()]

    description = sg.notes.get("Description") if hasattr(sg.notes, "get") else None
    if description and not doc.get("description"):
        doc["description"] = html.unescape(re.sub(r"<[^>]+>", " ", str(description))).strip()


# ---- writing ----------------------------------------------------------------------------

def write_genbank(data: dict) -> str:
    out = io.StringIO()
    SeqIO.write(ove_to_record(data), out, "genbank")
    return out.getvalue()


def write_fasta(data: dict) -> str:
    record = ove_to_record(data)
    # FASTA has no 16-character limit on the identifier, unlike the GenBank LOCUS name.
    full = (data.get("name") or "Untitled").strip().replace(" ", "_") or "Untitled"
    record.id = record.name = full
    out = io.StringIO()
    SeqIO.write(record, out, "fasta")
    return out.getvalue()


_CONTROL = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f]")
_COLOUR_QUALIFIERS = {"label", "ApEinfo_fwdcolor", "ApEinfo_revcolor", "color"}


def _xml_text(value: Any) -> str:
    return _CONTROL.sub(" ", str(value))


def _ranges(ann: dict[str, Any]) -> list[tuple[int, int]]:
    """OVE inclusive ranges in position order; a range with start > end crosses the origin."""
    return [(int(loc["start"]), int(loc["end"])) for loc in ann.get("locations") or []] or [
        (int(ann["start"]), int(ann["end"]))]


def _forward(ann: dict[str, Any]) -> bool:
    return not (ann.get("strand") == -1 or ann.get("forward") is False)


def write_dna(data: dict) -> bytes:
    """A .dna file with the sequence, topology, features (with colours) and primers.

    The block layout follows the format as read by Biopython (Bio.SeqIO.SnapGeneIO) and sgffp;
    metadata defaults are those found in the Biopython test-suite files in tests/data. History,
    alignments and traces are not written.
    """
    from sgffp import SgffBindingSite, SgffFeature, SgffObject, SgffPrimer, SgffSegment, SgffWriter

    sequence = data.get("sequence") or ""
    if not sequence:
        raise SaveError("There is no sequence to save: the document is empty.")
    circular = bool(data.get("circular"))
    sg = SgffObject.new(sequence, topology="circular" if circular else "linear",
                        strandedness="double")

    for index, ann in enumerate(_as_list(data.get("features"))):
        forward = _forward(ann)
        color = ann.get("color")
        segments = [SgffSegment(start=start, end=end + 1, color=color or None)
                    for start, end in _ranges(ann)]
        qualifiers = {k: [_xml_text(v) for v in (vals if isinstance(vals, list) else [vals])]
                      for k, vals in (ann.get("notes") or {}).items()
                      if k not in _COLOUR_QUALIFIERS}
        extras = {"recentID": str(index), "allowSegmentOverlaps": "0",
                  "consecutiveTranslationNumbering": "1"}
        if not forward:
            extras["swappedSegmentNumbering"] = "1"
        sg.features.add(SgffFeature(
            name=_xml_text(ann.get("name") or ann.get("type") or "feature"),
            type=ann.get("type") or "misc_feature", strand="+" if forward else "-",
            segments=segments, qualifiers=qualifiers, color=color or None, extras=extras))
    if sg.has_features:
        sg.features._wrapper_extras = {"nextValidID": str(len(sg.features))}
        sg.features._sync()

    for ann in _as_list(data.get("primers")):
        start, end = _ranges(ann)[0]
        site = sequence[start:end + 1] if start <= end else sequence[start:] + sequence[:end + 1]
        forward = _forward(ann)
        annealed = site if forward else reverse_complement(site)
        sg.primers.add(SgffPrimer(
            name=_xml_text(ann.get("name") or "primer"),
            sequence=(ann.get("bases") or annealed).upper(),
            binding_sites=[SgffBindingSite(start=start, end=end + 1,
                                           bound_strand="+" if forward else "-",
                                           annealed_bases=annealed.upper())]))

    now = datetime.now(UTC)
    stamp = {"UTC": f"{now.hour}:{now.minute}:{now.second}",
             "_text": f"{now.year}.{now.month}.{now.day}"}
    notes = {"UUID": str(uuid.uuid4()), "Type": "Synthetic", "ConfirmedExperimentally": "0",
             "Created": stamp, "LastModified": dict(stamp), "SequenceClass": "UNA",
             "TransformedInto": "unspecified"}
    if description := (data.get("description") or "").strip():
        # Stored as HTML, like the descriptions in the reference files.
        body = html.escape(_xml_text(description), quote=False)
        notes["Description"] = f"<html><body>{body}</body></html>"
    sg.notes._data = notes
    sg.notes._sync()
    sg.bset(8, {"AdditionalSequenceProperties": {
        "UpstreamStickiness": "0", "DownstreamStickiness": "0",
        "UpstreamModification": "Unmodified", "DownstreamModification": "Unmodified"}})
    return SgffWriter.to_bytes(sg)


def encode(data: dict, fmt: str) -> bytes:
    """The bytes of a saved document in ``fmt`` (genbank, fasta or the .dna format)."""
    if fmt == "genbank":
        return write_genbank(data).encode("utf-8")
    if fmt == "fasta":
        return write_fasta(data).encode("utf-8")
    if fmt == DNA_FILE_FORMAT:
        return write_dna(data)
    raise SaveError(f"Vettore cannot save in the {fmt} format; choose GenBank, FASTA or .dna.")


def save_target(path: str | Path, fmt: str | None = None) -> tuple[Path, str]:
    """Where to save and in which format: the extension decides; without a known extension the
    requested format (GenBank by default) is used and its extension appended."""
    target = Path(path)
    by_extension = SAVE_FORMATS.get(target.suffix.lower())
    if by_extension:
        return target, by_extension
    chosen = fmt if fmt in DEFAULT_EXTENSION else "genbank"
    return target.with_name(target.name + DEFAULT_EXTENSION[chosen]), chosen
