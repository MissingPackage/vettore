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


def test_primer_sequence_survives_genbank_round_trip():
    doc = read_doc("sample-hybridization-params.dna")
    [again] = fileio.read_documents("x.gb", fileio.write_genbank(doc).encode())
    primers = {p["name"]: p.get("bases") for p in again["primers"]}
    assert primers["XmaI-hht2(Cter)-Rev"] == "AAAAAACCCGGGAGAACGCTCACCACG"
