"""Files from vector vendors and other labs open as their authors meant them.

``fidelity-notes.gb`` is the sequence of NCBI PZ629150.1 with the record's own feature
coordinates, rewritten in the notation vector vendors use: the colour in a note
(``/note="color: #63b349; direction: RIGHT"``, also as a CSS colour name or ``None``) and
non-standard feature keys (``ITR``, ``Promoter``, ``PROMOTER``, ``Miscellaneous``,
``Transcription_t``, ``Rep_origin``) next to standard ones. See tests/data/README.md.
"""

import io
import re
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from vettore import fileio
from vettore.convert import (
    ORIGINAL_TYPE,
    STANDARD_TYPES,
    css_colour,
    ove_to_record,
    record_to_ove,
    split_colour_note,
    standard_type,
)

DATA = Path(__file__).parent / "data"
REPO = Path(__file__).resolve().parents[1]


def vendor_doc() -> dict:
    path = DATA / "fidelity-notes.gb"
    [doc] = fileio.read_documents(path.name, path.read_bytes())
    return doc


def by_name(doc: dict) -> dict[str, dict]:
    return {f["name"]: f for f in doc["features"]}


def test_the_fixture_is_the_published_sequence():
    published = SeqIO.read(DATA / "PZ629150.1.gb", "genbank")
    assert vendor_doc()["sequence"].upper() == str(published.seq).upper()


def test_colour_notes_give_features_the_authors_colours():
    features = by_name(vendor_doc())
    assert features["5' ITR"]["color"] == "#63b349"
    assert features["SaCas9"]["color"] == "#ce0084"
    assert features["terminator"]["color"] == "#FF378A"
    assert features["cTNT promoter"]["color"] == "#0000ff"        # "color: Blue"
    # "color: None" and a feature without a colour note: the map uses the type's colour.
    assert "color" not in features["3' ITR"]
    assert "color" not in features["AmpR promoter"]


def test_the_rest_of_a_colour_note_is_kept():
    features = by_name(vendor_doc())
    assert features["5' ITR"]["notes"]["note"] == ["direction: RIGHT"]
    assert features["3' ITR"]["notes"]["note"] == ["direction: RIGHT"]
    every_note = [n for f in vendor_doc()["features"] for n in f["notes"].get("note", [])]
    assert not any("color" in note.lower() for note in every_note)


def test_vendor_feature_keys_become_the_types_the_map_groups_by():
    features = by_name(vendor_doc())
    expected = {
        "5' ITR": ("repeat_region", "ITR"),
        "3' ITR": ("repeat_region", "ITR"),
        "Kozak": ("misc_feature", "Miscellaneous"),
        "terminator": ("terminator", "Transcription_t"),
        "cTNT promoter": ("promoter", None),          # "Promoter": only the letter case differs
        "U6 promoter": ("promoter", None),            # "PROMOTER"
        "AmpR promoter": ("promoter", None),
        "f1 ori": ("rep_origin", None),               # "Rep_origin"
        "SaCas9": ("CDS", None),
    }
    for name, (kind, original) in expected.items():
        feature = features[name]
        assert feature["type"] == kind, name
        assert feature["notes"].get(ORIGINAL_TYPE) == ([original] if original else None), name


def test_an_unknown_key_is_left_as_it_is_and_also_kept_in_the_notes():
    """The editor shows a key it does not know as misc_feature; the note keeps the key."""
    sgrna = by_name(vendor_doc())["Cd36 sgRNA"]
    assert sgrna["type"] == "sgRNA"
    assert sgrna["notes"][ORIGINAL_TYPE] == ["sgRNA"]


@pytest.mark.parametrize(("key", "expected"), [
    ("Promoter", "promoter"), ("PROMOTER", "promoter"), ("PolyA_signal", "polyA_signal"),
    ("polya_signal", "polyA_signal"), ("Rep_origin", "rep_origin"), ("misc_FEATURE", "misc_feature"),
    ("Transcription_t", "terminator"), ("transcription_terminator", "terminator"),
    ("Miscellaneous", "misc_feature"), ("ITR", "repeat_region"), ("itr", "repeat_region"),
    ("LTR", "LTR"), ("ltr", "LTR"), ("cds", "CDS"), ("Repeat_region", "repeat_region"),
    ("shRNA", None), ("Kozak", None), ("", None),
])
def test_standard_types_ignore_letter_case(key, expected):
    assert standard_type(key) == expected


def test_the_standard_types_are_the_editors_own_list():
    """The list in convert.py mirrors the editor component's; check it when the frontend's
    packages are installed."""
    source = (REPO / "frontend" / "node_modules" / "@teselagen" / "sequence-utils"
              / "index.js")
    if not source.exists():
        pytest.skip("frontend packages are not installed")
    text = source.read_text(encoding="utf-8")
    block = text[text.index("const genbankFeatureTypes = ["):]
    block = block[: block.index("];")]
    editor = [name for name in re.findall(r'name: "([^"]+)"', block) if name != "default"]
    assert sorted(editor) == sorted(STANDARD_TYPES)


@pytest.mark.parametrize(("note", "expected"), [
    ("color: #63b349; direction: RIGHT", ("#63b349", "direction: RIGHT")),
    ("Color: Blue", ("#0000ff", "")),
    ("colour=#abc", ("#abc", "")),
    ("color: None; direction: LEFT", (None, "direction: LEFT")),
    ("direction: RIGHT; color: #FF378A", ("#FF378A", "direction: RIGHT")),
    ("an ordinary note", (None, None)),
    ("the colour of the colony; blue/white screen", (None, None)),
    ("color: reddish", (None, None)),                   # not a colour: the note stays as it is
])
def test_colour_entries_are_taken_out_of_notes(note, expected):
    assert split_colour_note(note) == expected


def test_css_colour_names_are_drawn_as_hex():
    assert css_colour("Blue") == "#0000ff"
    assert css_colour("rebeccapurple") == "#663399"
    assert css_colour("grey") == css_colour("gray") == "#808080"
    assert css_colour("#AbC") == "#AbC"
    assert css_colour("#12345") is None and css_colour("None") is None


def record_with(feature_type: str, qualifiers: dict) -> SeqRecord:
    record = SeqRecord(Seq("ACGT" * 25), annotations={"molecule_type": "DNA",
                                                      "topology": "linear"})
    record.features.append(SeqFeature(SimpleLocation(10, 40, 1), type=feature_type,
                                      qualifiers=qualifiers))
    return record


def test_a_colour_qualifier_comes_before_the_note():
    record = record_with("CDS", {"label": ["x"], "ApEinfo_fwdcolor": ["#123456"],
                                 "note": ["color: #abcdef; direction: RIGHT"]})
    [feature] = record_to_ove(record)["features"]
    assert feature["color"] == "#123456"
    assert feature["notes"]["note"] == ["direction: RIGHT"]


def test_a_white_note_colour_is_not_used():
    [feature] = record_to_ove(record_with("CDS", {"label": ["x"],
                                                  "note": ["color: #ffffff"]}))["features"]
    assert "color" not in feature and "note" not in feature["notes"]


def test_a_feature_named_only_by_its_note_is_not_named_after_the_colour():
    [feature] = record_to_ove(record_with(
        "misc_feature", {"note": ["color: #63b349; direction: RIGHT"]}))["features"]
    assert feature["name"] == "direction: RIGHT"
    [bare] = record_to_ove(record_with("misc_feature", {"note": ["color: #63b349"]}))["features"]
    assert bare["name"] == "misc_feature"


def qualifiers_of(text: str) -> list[tuple[str, dict]]:
    record = SeqIO.read(io.StringIO(text), "genbank")
    return [(f.type, f.qualifiers) for f in record.features]


def test_saved_genbank_writes_each_colour_once_and_the_standard_keys():
    doc = vendor_doc()
    saved = qualifiers_of(fileio.write_genbank(doc))
    features = {q["label"][0]: (kind, q) for kind, q in saved}
    kind, itr = features["5' ITR"]
    assert kind == "repeat_region"
    assert itr["ApEinfo_fwdcolor"] == ["#63b349"]
    assert itr["note"] == ["direction: RIGHT"]
    assert itr[ORIGINAL_TYPE] == ["ITR"]
    assert features["terminator"][0] == "terminator"
    for _kind, qualifiers in saved:
        assert not any("color" in n.lower() for n in qualifiers.get("note", []))
    # A key the file kept as it was needs no second copy in the notes.
    kind, sgrna = features["Cd36 sgRNA"]
    assert kind == "sgRNA" and ORIGINAL_TYPE not in sgrna


def test_a_key_the_editor_showed_as_misc_feature_is_saved_with_its_original_key():
    doc = vendor_doc()
    for feature in doc["features"]:
        if feature["type"] == "sgRNA":
            feature["type"] = "misc_feature"          # what the editor does with unknown keys
    saved = {q["label"][0]: (kind, q) for kind, q in qualifiers_of(fileio.write_genbank(doc))}
    kind, sgrna = saved["Cd36 sgRNA"]
    assert kind == "misc_feature" and sgrna[ORIGINAL_TYPE] == ["sgRNA"]


def test_a_colour_note_left_in_a_document_is_not_written_twice():
    doc = {"name": "t", "sequence": "ACGT" * 25, "circular": False, "features": [
        {"name": "f", "type": "CDS", "start": 0, "end": 29, "strand": 1, "color": "#00ff00",
         "notes": {"note": ["color: #ff0000; direction: RIGHT"]}}]}
    [feature] = ove_to_record(doc).features
    assert feature.qualifiers["ApEinfo_fwdcolor"] == ["#00ff00"]
    assert feature.qualifiers["note"] == ["direction: RIGHT"]


@pytest.mark.parametrize("extension", [".gb", ".dna"])
def test_colours_types_and_notes_survive_saving_and_reopening(extension):
    doc = vendor_doc()
    content = fileio.encode(doc, "genbank" if extension == ".gb" else fileio.DNA_FILE_FORMAT)
    [again] = fileio.read_documents(f"again{extension}", content)

    def key(f):
        return (f["name"], f["type"], f["start"], f["end"], f.get("color"),
                tuple(f["notes"].get("note", [])), tuple(f["notes"].get(ORIGINAL_TYPE, [])))

    assert sorted(map(key, again["features"])) == sorted(map(key, doc["features"]))
