"""Writing .dna files: what Vettore writes must read back the same with Biopython and with sgffp.

Each fixture is opened, written as .dna and read back three ways: Biopython's reader (sequence,
topology, feature and primer locations), sgffp (names, types, strands, colours, primer
sequences with their 5' tails) and Vettore's own reader (the whole document).
"""

import io
from pathlib import Path

import pytest
from Bio import SeqIO
from sgffp import SgffReader

from vettore import fileio
from vettore.convert import record_to_ove

DATA = Path(__file__).parent / "data"
FIXTURES = ["pFA-KanMX4.dna", "sample-d.dna", "sample-e.dna", "sample-hybridization-params.dna",
            "MT612431.1.gb", "PZ629150.1.gb", "exchange-pBR322.embl", "exchange-sample-b.xdna"]


def location_key(ann: dict) -> tuple:
    return (ann["type"], ann["name"], ann["start"], ann["end"], ann["strand"],
            tuple((loc["start"], loc["end"]) for loc in ann.get("locations") or []))


def primer_key(ann: dict) -> tuple:
    return (ann["name"], ann["start"], ann["end"], ann["strand"])


def open_doc(name: str) -> dict:
    [doc] = fileio.read_documents(name, (DATA / name).read_bytes())
    return doc


@pytest.mark.parametrize("name", FIXTURES)
def test_biopython_reads_back_sequence_topology_and_locations(name):
    doc = open_doc(name)
    blob = fileio.write_dna(doc)
    record = next(SeqIO.parse(io.BytesIO(blob), fileio.DNA_FILE_FORMAT))
    assert str(record.seq) == doc["sequence"]
    assert record.annotations["topology"] == ("circular" if doc["circular"] else "linear")
    again = record_to_ove(record)
    assert sorted(map(location_key, again["features"])) == sorted(map(location_key, doc["features"]))
    assert sorted(map(primer_key, again["primers"])) == sorted(map(primer_key, doc["primers"]))


@pytest.mark.parametrize("name", FIXTURES)
def test_sgffp_reads_back_names_strands_colours_and_primers(name):
    doc = open_doc(name)
    sg = SgffReader.from_bytes(fileio.write_dna(doc))
    assert sg.sequence.value == doc["sequence"]
    assert sg.sequence.topology == ("circular" if doc["circular"] else "linear")

    def ove(ann):
        strand = "-" if ann["strand"] == -1 else "+"
        ranges = [(loc["start"], loc["end"]) for loc in ann.get("locations") or []] or [
            (ann["start"], ann["end"])]
        return (ann["name"], ann["type"], strand, tuple(ranges), ann.get("color"))

    def sgff(feature):
        ranges = tuple((s.start, s.end - 1) for s in feature.segments)
        return (feature.name, feature.type, feature.strand, ranges, feature.segments[0].color)

    assert sorted(map(sgff, sg.features), key=str) == sorted(map(ove, doc["features"]), key=str)
    written = {p.name: p for p in sg.primers}
    assert set(written) == {p["name"] for p in doc["primers"]}
    for ann in doc["primers"]:
        primer = written[ann["name"]]
        [site] = primer.binding_sites
        assert (site.start, site.end - 1) == (ann["start"], ann["end"])
        assert site.bound_strand == ("-" if ann["strand"] == -1 else "+")
        if ann.get("bases"):
            assert primer.sequence == ann["bases"]


@pytest.mark.parametrize("name", FIXTURES)
def test_vettore_reopens_its_own_dna_files_unchanged(name):
    doc = open_doc(name)
    [again] = fileio.read_documents("again.dna", fileio.write_dna(doc))

    def full(ann):
        return (*location_key(ann), ann.get("color"))

    assert again["sequence"] == doc["sequence"]
    assert again["circular"] == doc["circular"]
    assert sorted(map(full, again["features"])) == sorted(map(full, doc["features"]))
    assert (sorted((*primer_key(p), p.get("bases")) for p in again["primers"])
            == sorted((*primer_key(p), p.get("bases")) for p in doc["primers"]))
    assert again["description"] == doc["description"]


def test_tailed_primers_keep_their_tails():
    doc = open_doc("sample-hybridization-params.dna")
    sg = SgffReader.from_bytes(fileio.write_dna(doc))
    tails = {p.name: p.sequence for p in sg.primers}
    assert tails["XhoI-hht2(US)-Fwd"] == "AAAAAACTCGAGGAAAAGCTTCAAC"
    assert tails["XmaI-hht2(Cter)-Rev"] == "AAAAAACCCGGGAGAACGCTCACCACG"


# A 60-bp stretch of pBR322 (NCBI J01749.1, bases 1-60) carrying the hard cases.
EDGE = {
    "name": "edge µ-cases",
    "description": "Ampersand & <angle> brackets",
    "sequence": "TTCTCATGTTTGACAGCTTATCATCGATAAGCTTTAATGCGGTAGTTTATCACAGTTAAA",
    "circular": True,
    "features": [
        # Across the origin, with a colour.
        {"name": "wrap", "type": "misc_feature", "start": 55, "end": 4, "strand": 1,
         "color": "#ff0000"},
        # A minus-strand join, with a note containing a control character.
        {"name": "split", "type": "CDS", "start": 10, "end": 39, "strand": -1,
         "locations": [{"start": 10, "end": 19}, {"start": 30, "end": 39}],
         "notes": {"note": ["tab\there, bell\x07"]}},
    ],
    "primers": [
        # Minus-strand primer with a 5' tail, across the origin.
        {"name": "rev-tail", "start": 50, "end": 5, "strand": -1, "bases": "GGATCCTTTAA"},
        {"name": "fwd", "start": 20, "end": 35, "strand": 1},
    ],
}


def test_hard_cases_round_trip():
    blob = fileio.write_dna(EDGE)
    record = next(SeqIO.parse(io.BytesIO(blob), fileio.DNA_FILE_FORMAT))
    locations = {f.qualifiers["label"][0]: str(f.location) for f in record.features}
    assert locations["wrap"] == "join{[55:60](+), [0:5](+)}"
    assert locations["split"] == "join{[30:40](-), [10:20](-)}"
    assert locations["fwd"] == "[20:36](+)"
    assert locations["rev-tail"] == "join{[50:60](-), [0:6](-)}"
    [again] = fileio.read_documents("edge.dna", blob)
    assert again["description"] == EDGE["description"]
    wrap = next(f for f in again["features"] if f["name"] == "wrap")
    assert (wrap["start"], wrap["end"], wrap.get("color")) == (55, 4, "#ff0000")
    split = next(f for f in again["features"] if f["name"] == "split")
    assert split["locations"] == [{"start": 10, "end": 19}, {"start": 30, "end": 39}]
    rev = next(p for p in again["primers"] if p["name"] == "rev-tail")
    assert (rev["start"], rev["end"], rev["strand"], rev["bases"]) == (50, 5, -1, "GGATCCTTTAA")
    sg = SgffReader.from_bytes(blob)
    # The annealing part of the minus-strand primer is the reverse complement of its site.
    site = EDGE["sequence"][50:] + EDGE["sequence"][:6]
    from Bio.Seq import reverse_complement
    assert sg.primers[0].binding_sites[0].annealed_bases == reverse_complement(site)


def test_empty_document_is_not_written():
    with pytest.raises(fileio.SaveError, match="no sequence"):
        fileio.write_dna({"name": "x", "sequence": "", "circular": False})
