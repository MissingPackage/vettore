"""Conversion between Biopython SeqRecords and the Open Vector Editor (OVE) JSON model.

Coordinate conventions:
- Biopython: 0-based, end-exclusive (``SimpleLocation(start, end)``).
- OVE: 0-based, end-inclusive (``start``/``end``). A feature that crosses the origin of a
  circular sequence has ``start > end``.
"""

from __future__ import annotations

import re
import uuid
from typing import Any

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

PRIMER_TYPES = {"primer_bind", "primer"}
LOCUS_LENGTH = 16            # characters of a GenBank LOCUS name that every reader accepts
# GenBank divisions (NCBI GenBank release notes, section 3.4.4); "UNK" is not one of them.
DIVISIONS = frozenset({"PRI", "ROD", "MAM", "VRT", "INV", "PLN", "BCT", "VRL", "PHG", "SYN",
                       "UNA", "EST", "PAT", "STS", "GSS", "HTG", "HTC", "ENV", "CON", "TSA"})
NAME_QUALIFIERS = ("label", "name", "gene", "product", "locus_tag", "note")

# ---- feature types ------------------------------------------------------------------------

#: The feature types the editor colours and groups features by: the INSDC feature keys (and
#: the few older keys plasmid maps still use) in the editor component's own list
#: (``genbankFeatureTypes`` in @teselagen/sequence-utils). The editor matches a type to this
#: list without regard to case and shows anything else as misc_feature.
STANDARD_TYPES = (
    "-10_signal", "-35_signal", "3'clip", "3'UTR", "5'clip", "5'UTR", "D-loop", "assembly_gap",
    "centromere", "Het", "mobile_element", "ncRNA", "proprotein", "regulatory", "SecStr", "Site",
    "telomere", "tmRNA", "unsure", "V_segment", "allele", "attenuator", "C_region",
    "CAAT_signal", "CDS", "conserved", "D_segment", "enhancer", "exon", "gap", "GC_signal",
    "gene", "iDNA", "intron", "J_region", "LTR", "m_rna", "mat_peptide", "misc_binding",
    "misc_difference", "misc_feature", "misc_marker", "misc_part", "misc_recomb", "misc_RNA",
    "misc_signal", "misc_structure", "modified_base", "mRNA", "N_region", "old_sequence",
    "operator", "operon", "oriT", "plasmid", "polyA_signal", "polyA_site", "precursor_RNA",
    "prim_transcript", "primer_bind", "promoter", "protein_bind", "protein_domain", "protein",
    "RBS", "rep_origin", "repeat_region", "repeat_unit", "rRNA", "s_mutation", "S_region",
    "satellite", "scRNA", "sig_peptide", "snoRNA", "snRNA", "source", "start", "stem_loop",
    "stop", "STS", "tag", "TATA_signal", "terminator", "transit_peptide", "transposon", "tRNA",
    "V_region", "variation",
)
_STANDARD_BY_LOWER = {key.lower(): key for key in STANDARD_TYPES}

#: Feature keys that vector vendors and other tools write instead of the standard key that
#: means the same thing (compared without regard to case). A GenBank feature key is at most
#: 15 characters long, so "Transcription_terminator" arrives cut to "Transcription_t".
TYPE_SYNONYMS = {
    "transcription_t": "terminator",
    "transcription_terminator": "terminator",
    "miscellaneous": "misc_feature",
    "itr": "repeat_region",                      # inverted terminal repeat
    "inverted_terminal_repeat": "repeat_region",
    "long_terminal_repeat": "LTR",
    "polyadenylation_signal": "polyA_signal",
    "origin_of_replication": "rep_origin",
    "replication_origin": "rep_origin",
    "signal_peptide": "sig_peptide",
    "ribosome_binding_site": "RBS",
}

#: The note that keeps a feature's key from the file when the editor shows another type.
ORIGINAL_TYPE = "original_type"


def standard_type(key: str) -> str | None:
    """The editor's type for a feature key from a file, or None for a key it does not know."""
    lowered = key.strip().lower()
    return _STANDARD_BY_LOWER.get(lowered) or TYPE_SYNONYMS.get(lowered)


# ---- feature colours ----------------------------------------------------------------------

#: CSS colour names (CSS Color Module Level 4), as written in feature notes ("color: Blue").
_CSS_TABLE = """
aliceblue f0f8ff antiquewhite faebd7 aqua 00ffff aquamarine 7fffd4 azure f0ffff
beige f5f5dc bisque ffe4c4 black 000000 blanchedalmond ffebcd blue 0000ff
blueviolet 8a2be2 brown a52a2a burlywood deb887 cadetblue 5f9ea0 chartreuse 7fff00
chocolate d2691e coral ff7f50 cornflowerblue 6495ed cornsilk fff8dc crimson dc143c
cyan 00ffff darkblue 00008b darkcyan 008b8b darkgoldenrod b8860b darkgray a9a9a9
darkgreen 006400 darkgrey a9a9a9 darkkhaki bdb76b darkmagenta 8b008b
darkolivegreen 556b2f darkorange ff8c00 darkorchid 9932cc darkred 8b0000
darksalmon e9967a darkseagreen 8fbc8f darkslateblue 483d8b darkslategray 2f4f4f
darkslategrey 2f4f4f darkturquoise 00ced1 darkviolet 9400d3 deeppink ff1493
deepskyblue 00bfff dimgray 696969 dimgrey 696969 dodgerblue 1e90ff firebrick b22222
floralwhite fffaf0 forestgreen 228b22 fuchsia ff00ff gainsboro dcdcdc ghostwhite f8f8ff
gold ffd700 goldenrod daa520 gray 808080 green 008000 greenyellow adff2f grey 808080
honeydew f0fff0 hotpink ff69b4 indianred cd5c5c indigo 4b0082 ivory fffff0 khaki f0e68c
lavender e6e6fa lavenderblush fff0f5 lawngreen 7cfc00 lemonchiffon fffacd
lightblue add8e6 lightcoral f08080 lightcyan e0ffff lightgoldenrodyellow fafad2
lightgray d3d3d3 lightgreen 90ee90 lightgrey d3d3d3 lightpink ffb6c1 lightsalmon ffa07a
lightseagreen 20b2aa lightskyblue 87cefa lightslategray 778899 lightslategrey 778899
lightsteelblue b0c4de lightyellow ffffe0 lime 00ff00 limegreen 32cd32 linen faf0e6
magenta ff00ff maroon 800000 mediumaquamarine 66cdaa mediumblue 0000cd
mediumorchid ba55d3 mediumpurple 9370db mediumseagreen 3cb371 mediumslateblue 7b68ee
mediumspringgreen 00fa9a mediumturquoise 48d1cc mediumvioletred c71585
midnightblue 191970 mintcream f5fffa mistyrose ffe4e1 moccasin ffe4b5
navajowhite ffdead navy 000080 oldlace fdf5e6 olive 808000 olivedrab 6b8e23
orange ffa500 orangered ff4500 orchid da70d6 palegoldenrod eee8aa palegreen 98fb98
paleturquoise afeeee palevioletred db7093 papayawhip ffefd5 peachpuff ffdab9
peru cd853f pink ffc0cb plum dda0dd powderblue b0e0e6 purple 800080
rebeccapurple 663399 red ff0000 rosybrown bc8f8f royalblue 4169e1 saddlebrown 8b4513
salmon fa8072 sandybrown f4a460 seagreen 2e8b57 seashell fff5ee sienna a0522d
silver c0c0c0 skyblue 87ceeb slateblue 6a5acd slategray 708090 slategrey 708090
snow fffafa springgreen 00ff7f steelblue 4682b4 tan d2b48c teal 008080 thistle d8bfd8
tomato ff6347 turquoise 40e0d0 violet ee82ee wheat f5deb3 white ffffff
whitesmoke f5f5f5 yellow ffff00 yellowgreen 9acd32
"""
_CSS_WORDS = _CSS_TABLE.split()
_CSS_NAMES = dict(zip(_CSS_WORDS[::2], _CSS_WORDS[1::2], strict=True))
_HEX_COLOUR = re.compile(r"#(?:[0-9a-f]{3}|[0-9a-f]{6})", re.IGNORECASE)
#: White or nearly white: invisible on the map, so the type's colour is used instead.
NEAR_WHITE = re.compile(r"^#?(f[a-f0-9]){3}$", re.IGNORECASE)
#: One entry of the "key: value; key: value" notes some vendors write, e.g.
#: /note="color: #63b349; direction: RIGHT".
_COLOUR_ENTRY = re.compile(r"^\s*colou?r\s*[:=]\s*(?P<value>.*?)\s*$", re.IGNORECASE)
_NO_COLOUR = {"", "none", "null", "default", "auto"}


def css_colour(value: str) -> str | None:
    """A colour the map can draw: ``#rgb``/``#rrggbb`` as written, a CSS colour name as
    ``#rrggbb``; None for anything else."""
    text = value.strip()
    if _HEX_COLOUR.fullmatch(text):
        return text
    if hex_value := _CSS_NAMES.get(text.lower()):
        return f"#{hex_value}"
    return None


def split_colour_note(note: str) -> tuple[str | None, str | None]:
    """Take the colour out of a note in the vendors' "color: …; direction: …" form.

    Returns (colour, the rest of the note); the rest is None when the note has no colour
    entry, so the note is left as it is, and empty when the colour was all it said. An entry
    saying there is no colour ("color: None") is removed too; one whose value is not a colour
    stays in the note.
    """
    entries = note.split(";")
    colour = None
    kept: list[str] = []
    found = False
    for entry in entries:
        match = _COLOUR_ENTRY.match(entry)
        value = match.group("value") if match else None
        parsed = css_colour(value) if value is not None else None
        if parsed or (value is not None and value.lower() in _NO_COLOUR):
            found = True
            colour = colour or parsed
        elif entry.strip():
            kept.append(entry.strip())
    if not found:
        return None, None
    return colour, "; ".join(kept)


def _first(qualifiers: dict[str, Any], key: str) -> str | None:
    value = qualifiers.get(key)
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        return str(value[0]) if value else None
    return str(value)


def feature_name(feature: SeqFeature) -> str:
    for key in NAME_QUALIFIERS:
        value = _first(feature.qualifiers, key)
        if value and key == "note":
            _colour, rest = split_colour_note(value)
            value = value if rest is None else rest
        if value:
            return value
    return feature.type


def _parts_to_ranges(feature: SeqFeature, seq_len: int, circular: bool) -> list[tuple[int, int]]:
    """Return inclusive (start, end) ranges in transcription order, merging origin wraps."""
    parts = list(feature.location.parts)
    if feature.location.strand == -1:
        # Biopython lists parts of minus-strand joins in transcription order (high to low);
        # OVE wants locations sorted by position.
        parts = sorted(parts, key=lambda p: int(p.start))
    ranges = [(int(p.start), int(p.end) - 1) for p in parts]
    if circular and len(ranges) >= 2:
        merged: list[tuple[int, int]] = []
        for rng in ranges:
            if merged and merged[-1][1] == seq_len - 1 and rng[0] == 0:
                merged[-1] = (merged[-1][0], rng[1])
            else:
                merged.append(rng)
        # A wrap can also be expressed as the last part ending at the origin.
        if len(merged) >= 2 and merged[0][0] == 0 and merged[-1][1] == seq_len - 1:
            merged = [(merged[-1][0], merged[0][1])] + merged[1:-1]
        ranges = merged
    return ranges


def _type_and_notes(key: str, notes: dict[str, list[str]]) -> str:
    """The editor's type for a feature key. A key the editor would show under another type
    (a vendor's "ITR", or a key it does not know) is kept in the notes, so it is not lost."""
    standard = standard_type(key)
    if standard and standard.lower() == key.strip().lower():
        return standard                  # the same key, perhaps in other letter case
    if key not in PRIMER_TYPES:
        notes.setdefault(ORIGINAL_TYPE, [key])
    return standard or key


def _note_colour(notes: dict[str, list[str]]) -> str | None:
    """The colour written in the notes' "color: …" entries; the entries are removed and the
    rest of each note kept."""
    colour = None
    kept = []
    for note in notes.get("note", []):
        found, rest = split_colour_note(note)
        colour = colour or found
        if rest is None:
            kept.append(note)
        elif rest:
            kept.append(rest)
    if "note" in notes:
        if kept:
            notes["note"] = kept
        else:
            del notes["note"]
    return colour


def _annotation(feature: SeqFeature, seq_len: int, circular: bool) -> dict[str, Any]:
    ranges = _parts_to_ranges(feature, seq_len, circular)
    strand = feature.location.strand or 1
    notes = {k: [str(x) for x in v] if isinstance(v, (list, tuple)) else [str(v)]
             for k, v in feature.qualifiers.items()}
    kind = _type_and_notes(feature.type, notes)
    note_colour = _note_colour(notes)
    ann: dict[str, Any] = {
        "id": uuid.uuid4().hex[:12],
        "name": feature_name(feature),
        "type": kind,
        "strand": strand,
        "forward": strand != -1,
        "start": ranges[0][0],
        "end": ranges[-1][1],
        "notes": notes,
    }
    if len(ranges) > 1:
        ann["locations"] = [{"start": s, "end": e} for s, e in ranges]
    if bases := _first(feature.qualifiers, "primer_sequence"):
        ann["bases"] = bases
    # A colour qualifier is the more specific statement; the note convention comes after it.
    written = _first(feature.qualifiers, "ApEinfo_fwdcolor") or _first(feature.qualifiers, "color")
    if written:
        ann["color"] = css_colour(written) or written
    elif note_colour and not NEAR_WHITE.match(note_colour):
        ann["color"] = note_colour       # white would be invisible: the type's colour is used
    return ann


def record_to_ove(record: SeqRecord, fallback_name: str = "Untitled") -> dict[str, Any]:
    """Convert a SeqRecord to OVE ``sequenceData``."""
    sequence = str(record.seq)
    circular = str(record.annotations.get("topology", "linear")).lower() == "circular"
    name = record.name if record.name and record.name not in {"<unknown name>", "."} else ""
    features, primers = [], []
    for feature in record.features:
        if feature.type == "source":
            continue
        ann = _annotation(feature, len(sequence), circular)
        (primers if ann["type"] in PRIMER_TYPES else features).append(ann)
    doc = {
        "name": name or fallback_name,
        "description": record.description if record.description != "<unknown description>" else "",
        "sequence": sequence,
        "circular": circular,
        "features": features,
        "primers": primers,
    }
    if (division := record.annotations.get("data_file_division")) in DIVISIONS:
        doc["division"] = division       # written back when the document is saved as GenBank
    return doc


def locus_name(name: str) -> str:
    """A GenBank LOCUS name: no white space, at most 16 characters, shortened at a word or
    hyphen boundary when the name is longer (the full name goes in the DEFINITION line)."""
    compact = "_".join(name.split()) or "Untitled"
    if len(compact) <= LOCUS_LENGTH:
        return compact
    cut = max(compact.rfind(mark, 0, LOCUS_LENGTH + 1) for mark in "_-.")
    return compact[:cut] if cut >= LOCUS_LENGTH // 2 else compact[:LOCUS_LENGTH]


def _is_position(value: Any, length: int) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and 0 <= value < length


def annotation_ranges(ann: Any, length: int, circular: bool) -> list[tuple[int, int]] | None:
    """The inclusive (start, end) ranges of an OVE feature or primer, or None when its position
    cannot be placed on the sequence (missing, not a number, outside the sequence, or across
    the end of a linear sequence)."""
    if not isinstance(ann, dict):
        return None
    locations = ann.get("locations")
    if isinstance(locations, list) and locations:
        raw = [(loc.get("start"), loc.get("end")) if isinstance(loc, dict) else (None, None)
               for loc in locations]
    else:
        raw = [(ann.get("start"), ann.get("end"))]
    for start, end in raw:
        if not (_is_position(start, length) and _is_position(end, length)):
            return None
        if start > end and not circular:
            return None
    return raw


def annotation_label(ann: Any, default: str) -> str:
    name = ann.get("name") if isinstance(ann, dict) else None
    return str(name or (ann.get("type") if isinstance(ann, dict) else None) or default)


def unplaceable(data: dict[str, Any]) -> dict[str, list[str]]:
    """Names of the features and primers whose position cannot be placed on the sequence,
    by kind ("features", "primers"); conversion leaves them out."""
    sequence = data.get("sequence")
    length = len(sequence) if isinstance(sequence, str) else 0
    circular = bool(data.get("circular"))
    out: dict[str, list[str]] = {}
    for kind, default in (("features", "feature"), ("primers", "primer")):
        bad = [annotation_label(ann, default) for ann in _as_list(data.get(kind))
               if annotation_ranges(ann, length, circular) is None]
        if bad:
            out[kind] = bad
    return out


def _location(ann: dict[str, Any], ranges: list[tuple[int, int]],
              seq_len: int) -> SimpleLocation | CompoundLocation:
    strand = -1 if (ann.get("strand") == -1 or ann.get("forward") is False) else 1
    parts: list[SimpleLocation] = []
    for start, end in ranges:
        if start > end:  # crosses the origin
            parts.append(SimpleLocation(start, seq_len, strand))
            parts.append(SimpleLocation(0, end + 1, strand))
        else:
            parts.append(SimpleLocation(start, end + 1, strand))
    if strand == -1:
        parts = parts[::-1]
    return parts[0] if len(parts) == 1 else CompoundLocation(parts)


def _qualifiers(notes: Any) -> dict[str, list[str]]:
    if not isinstance(notes, dict):
        return {}
    return {str(k): [str(x) for x in v] if isinstance(v, (list, tuple)) else [str(v)]
            for k, v in notes.items()}


def annotation_qualifiers(ann: dict[str, Any], kind: str) -> dict[str, list[str]]:
    """The qualifiers an annotation's notes become in a saved file. Its colour is written once,
    by the caller, so a "color: …" entry still in a note is left out; a kept original key
    that is the same as the type ``kind`` is dropped."""
    qualifiers = _qualifiers(ann.get("notes"))
    if ann.get("color"):
        _note_colour(qualifiers)
    if qualifiers.get(ORIGINAL_TYPE) == [kind]:
        del qualifiers[ORIGINAL_TYPE]
    return qualifiers


def ove_to_record(data: dict[str, Any]) -> SeqRecord:
    """Convert OVE ``sequenceData`` back to a SeqRecord (GenBank-ready). Features and primers
    whose position cannot be placed on the sequence are left out (see ``unplaceable``)."""
    sequence = data.get("sequence")
    sequence = sequence if isinstance(sequence, str) else ""
    circular = bool(data.get("circular"))
    name = str(data.get("name") or "Untitled").strip() or "Untitled"
    locus = locus_name(name)
    description = str(data.get("description") or "").strip() or name
    if locus != "_".join(name.split()) and name not in description:
        description = f"{name}: {description}"     # the LOCUS name was shortened
    record = SeqRecord(
        Seq(sequence),
        id=locus,
        name=locus,
        description=description,
        annotations={
            "molecule_type": "DNA",
            "topology": "circular" if data.get("circular") else "linear",
        },
    )
    for kind, default_type in (("features", "misc_feature"), ("primers", "primer_bind")):
        for ann in _as_list(data.get(kind)):
            ranges = annotation_ranges(ann, len(sequence), circular)
            if ranges is None:
                continue
            feature_type = str(ann.get("type") or default_type)
            qualifiers = annotation_qualifiers(ann, feature_type)
            qualifiers["label"] = [str(ann.get("name") or default_type)]
            if ann.get("color"):
                qualifiers["ApEinfo_fwdcolor"] = [str(ann["color"])]
            if kind == "primers" and ann.get("bases"):
                qualifiers["primer_sequence"] = [str(ann["bases"])]
            record.features.append(
                SeqFeature(_location(ann, ranges, len(sequence)),
                           type=feature_type, qualifiers=qualifiers)
            )
    return record


def _as_list(value: Any) -> list[Any]:
    """OVE stores annotations either as a list or as a dict keyed by id; anything else holds
    none."""
    if isinstance(value, dict):
        return list(value.values())
    if isinstance(value, list):
        return value
    return []
