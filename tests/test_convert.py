from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from vettore import fileio
from vettore.convert import ove_to_record, record_to_ove

DATA = Path(__file__).parent / "data"


def read_doc(name: str) -> dict:
    path = DATA / name
    [doc] = fileio.read_documents(path.name, path.read_bytes())
    return doc


def test_dna_file_is_read_with_features_and_topology():
    doc = read_doc("pFA-KanMX4.dna")
    assert doc["name"] == "pFA-KanMX4"
    assert len(doc["sequence"]) == 3941
    assert doc["circular"] is True
    assert len(doc["features"]) == 9
    assert {f["type"] for f in doc["features"]} == {
        "CDS", "gene", "promoter", "rep_origin", "terminator"}


def test_dna_file_primers_become_ove_primers():
    doc = read_doc("sample-hybridization-params.dna")
    names = sorted(p["name"] for p in doc["primers"])
    assert names == ["XhoI-hht2(US)-Fwd", "XmaI-hht2(Cter)-Rev"]
    rev = next(p for p in doc["primers"] if p["name"].endswith("Rev"))
    # Biopython [1396:1411](-) -> OVE inclusive 1396..1410 on the minus strand
    assert (rev["start"], rev["end"], rev["strand"], rev["forward"]) == (1396, 1410, -1, False)


def _locations(record: SeqRecord) -> list[tuple[str, str]]:
    return [(f.type, str(f.location)) for f in record.features]


@pytest.mark.parametrize("name", ["pFA-KanMX4.dna", "sample-d.dna", "sample-e.dna",
                                  "sample-hybridization-params.dna"])
def test_genbank_round_trip_preserves_sequence_and_annotations(name):
    doc = read_doc(name)
    text = fileio.write_genbank(doc)
    [again] = fileio.read_documents("x.gb", text.encode())
    assert again["sequence"].upper() == doc["sequence"].upper()
    assert again["circular"] == doc["circular"]
    key = lambda a: (a["type"], a["name"], a["start"], a["end"], a["strand"])
    assert sorted(map(key, again["features"])) == sorted(map(key, doc["features"]))
    assert sorted(map(key, again["primers"])) == sorted(map(key, doc["primers"]))


def test_feature_across_origin_uses_start_greater_than_end():
    seq = "ACGT" * 25  # 100 bp
    record = SeqRecord(Seq(seq), annotations={"topology": "circular", "molecule_type": "DNA"})
    loc = CompoundLocation([SimpleLocation(90, 100, 1), SimpleLocation(0, 10, 1)])
    record.features.append(SeqFeature(loc, type="misc_feature", qualifiers={"label": ["wrap"]}))
    doc = record_to_ove(record)
    [feat] = doc["features"]
    assert (feat["start"], feat["end"]) == (90, 9)
    assert "locations" not in feat
    back = ove_to_record(doc)
    assert str(back.features[0].location) == str(loc)


def test_minus_strand_join_round_trips():
    seq = "ACGT" * 25
    record = SeqRecord(Seq(seq), annotations={"topology": "linear", "molecule_type": "DNA"})
    loc = CompoundLocation([SimpleLocation(50, 60, -1), SimpleLocation(10, 20, -1)])
    record.features.append(SeqFeature(loc, type="CDS", qualifiers={"label": ["split"]}))
    doc = record_to_ove(record)
    [feat] = doc["features"]
    assert feat["locations"] == [{"start": 10, "end": 19}, {"start": 50, "end": 59}]
    assert (feat["start"], feat["end"], feat["strand"]) == (10, 59, -1)
    back = ove_to_record(doc)
    assert str(back.features[0].location) == str(loc)


def test_ove_dict_keyed_annotations_are_accepted():
    doc = {"name": "t", "sequence": "ACGTACGTAC", "circular": False,
           "features": {"a": {"name": "f", "type": "CDS", "start": 0, "end": 2, "strand": 1}},
           "primers": {}}
    record = ove_to_record(doc)
    assert [(f.type, str(f.location)) for f in record.features] == [("CDS", "[0:3](+)")]


def test_ab1_is_read_as_sequence():
    [doc] = fileio.read_documents("3730.ab1", (DATA / "3730.ab1").read_bytes())
    assert len(doc["sequence"]) > 500


def test_unsupported_extension_is_rejected():
    with pytest.raises(fileio.UnsupportedFile):
        fileio.read_documents("notes.docx", b"hello")


def test_dna_file_colours_primer_tails_and_description_are_read():
    doc = read_doc("sample-hybridization-params.dna")
    colours = {f["name"]: f.get("color") for f in doc["features"]}
    assert colours["hht2-ORF"] == "#993366"
    primers = {p["name"]: p.get("bases") for p in doc["primers"]}
    assert primers["XhoI-hht2(US)-Fwd"] == "AAAAAACTCGAGGAAAAGCTTCAAC"
    assert doc["description"].startswith("SJAG_03543.1")


def test_white_feature_colours_are_ignored():
    doc = read_doc("pFA-KanMX4.dna")
    colours = {f["name"]: f.get("color") for f in doc["features"]}
    assert colours["SP6 promoter"] is None  # stored as #ffffff, would be invisible
    assert colours["ori"] == "#ffff00"


def locus_line(text: str) -> list[str]:
    return text.splitlines()[0].split()


def test_saved_genbank_has_todays_date_and_keeps_the_division():
    from datetime import datetime

    doc = read_doc("MT612431.1.gb")
    assert doc["division"] == "SYN"
    locus = locus_line(fileio.write_genbank({**doc, "division": "VRL"}))
    assert locus[-2:] == ["VRL", fileio.genbank_date(datetime.now())]
    # A construct made here, or a document whose division is unknown, is synthetic.
    plain = {"name": "pNew", "sequence": "ACGT" * 10, "circular": True}
    assert locus_line(fileio.write_genbank(plain))[-2] == "SYN"
    assert locus_line(fileio.write_genbank({**plain, "division": "UNK"}))[-2] == "SYN"


def test_genbank_date_is_the_same_in_every_locale():
    from datetime import datetime

    assert fileio.genbank_date(datetime(2008, 8, 7)) == "07-AUG-2008"
    assert fileio.genbank_date(datetime(2026, 10, 20)) == "20-OCT-2026"


@pytest.mark.parametrize(("name", "locus"), [
    ("pUC19", "pUC19"),
    ("colleague-pFA-KanMX4", "colleague-pFA"),
    ("pAAV CAG GFP WPRE bGH", "pAAV_CAG_GFP"),
    ("Abcdefghijklmnopqrstuvwxyz", "Abcdefghijklmnop"),
    ("tab\tand\nnewline", "tab_and_newline"),
    ("   ", "Untitled"),
])
def test_long_names_are_shortened_at_a_boundary_and_kept_in_the_definition(name, locus):
    text = fileio.write_genbank({"name": name, "sequence": "ACGT" * 10, "circular": False})
    assert locus_line(text)[1] == locus
    definition = text.splitlines()[1]
    if locus != "_".join(name.split()) and name.strip():
        assert definition.startswith(f"DEFINITION  {name}")


def test_primer_sequence_survives_genbank_round_trip():
    doc = read_doc("sample-hybridization-params.dna")
    [again] = fileio.read_documents("x.gb", fileio.write_genbank(doc).encode())
    primers = {p["name"]: p.get("bases") for p in again["primers"]}
    assert primers["XmaI-hht2(Cter)-Rev"] == "AAAAAACCCGGGAGAACGCTCACCACG"


def test_regulatory_features_open_as_their_class_and_save_in_their_own_form():
    # Current NCBI records write promoters as `regulatory` + /regulatory_class="promoter"
    # (PZ629150.1); the editor colours them as promoters, and saving keeps the record's form.
    from pathlib import Path

    from vettore import fileio

    data = Path(__file__).parent / "data" / "PZ629150.1.gb"
    doc = fileio.read_file(data.name, data.read_bytes()).documents[0]
    by_name = {f["name"]: f for f in doc["features"]}
    ctnt = next(f for name, f in by_name.items() if name.startswith("cardiac troponin T"))
    kozak = by_name["Kozak sequence"]
    assert ctnt["type"] == "promoter" and kozak["type"] == "regulatory"
    saved = fileio.write_genbank(doc)
    assert "     regulatory      144..556" in saved
    assert "original_type" not in saved
    again = fileio.read_file("again.gb", saved.encode()).documents[0]
    assert [f["type"] for f in again["features"]] == [f["type"] for f in doc["features"]]
