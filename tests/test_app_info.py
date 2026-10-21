"""Application routes: version, licence and the sample constructs on the welcome screen."""

import re
from pathlib import Path

from Bio import SeqIO
from fastapi.testclient import TestClient

from vettore import app_info
from vettore.server import TOKEN_HEADER, create_app

TOKEN = "t0ken"
SAMPLES = Path(app_info.__file__).parent / "data" / "samples"


def client() -> TestClient:
    return TestClient(create_app(TOKEN, None), base_url="http://127.0.0.1")


def get(path: str):
    return client().get(path, headers={TOKEN_HEADER: TOKEN})


def test_about_names_the_licence_and_where_the_log_is():
    data = get("/api/app/about").json()
    assert data["version"]
    assert data["licence"] == "GPL-3.0-or-later"
    assert data["log_path"].endswith(str(Path(".vettore") / "vettore.log"))
    names = {entry["name"] for entry in data["third_party"]}
    assert {"pydna", "Biopython", "primer3-py"} <= names
    # The GPL-2.0-or-later component is the reason the whole build is GPL; it must be listed.
    assert any(e["name"] == "primer3-py" and "GPL-2.0" in e["licence"] for e in data["third_party"])


def test_samples_are_listed_with_a_label_and_one_line_of_description():
    samples = get("/api/app/samples").json()["samples"]
    assert len(samples) >= 2
    for sample in samples:
        assert sample["id"] and sample["label"] and sample["description"]
        assert sample["accession"]
        assert len(sample["description"]) < 160


def test_samples_carry_an_outline_for_their_mini_map():
    samples = {s["id"]: s for s in get("/api/app/samples").json()["samples"]}
    puc = samples["puc19"]
    assert puc["length"] == 2686 and puc["circular"] is True
    assert 0 < len(puc["features"]) <= 12
    assert all(f["type"] != "source" for f in puc["features"])
    assert all(0 <= f["start"] < puc["length"] and 0 <= f["end"] < puc["length"] for f in puc["features"])


def test_opening_a_sample_gives_an_unsaved_document_named_after_the_construct():
    result = get("/api/app/samples/puc19").json()
    assert result["path"] is None          # unsaved: Save must ask where to write
    [doc] = result["documents"]
    assert doc["name"] == "pUC19"          # not the record's locus name (SYNPUC19V)
    assert len(doc["sequence"]) == 2686
    assert doc["circular"]


def test_every_sample_opens_and_matches_its_genbank_record():
    for sample in app_info.SAMPLES:
        record = SeqIO.read(SAMPLES / sample.filename, "genbank")
        assert record.id == sample.accession, sample.id
        [doc] = get(f"/api/app/samples/{sample.id}").json()["documents"]
        assert doc["sequence"].upper() == str(record.seq).upper(), sample.id
        assert doc["name"] == sample.label
        assert doc["circular"] is True, sample.id


def test_the_aav_samples_carry_two_itrs_that_the_packaging_check_can_find():
    """The AAV samples are there to be used with the AAV packaging check, which needs two
    ITRs it can recognise from the annotation alone. "ITR" is matched as a word: the SaCas9
    protein translation contains the letters ITR inside a peptide."""
    itr = re.compile(r"\bITR\b|inverted terminal repeat", re.IGNORECASE)
    for sample_id in ("paav-ctnt-sacas9", "paav-ef1a-dio"):
        [doc] = get(f"/api/app/samples/{sample_id}").json()["documents"]
        itrs = []
        for feature in doc["features"]:
            notes = {k: v for k, v in (feature.get("notes") or {}).items() if k != "translation"}
            if itr.search(feature["name"]) or itr.search(str(notes)):
                itrs.append(feature)
        assert len(itrs) == 2, sample_id
        assert itrs[0]["end"] < itrs[1]["start"], sample_id


def test_unknown_sample_says_which_one_is_missing():
    res = get("/api/app/samples/not-a-sample")
    assert res.status_code == 404
    assert "not-a-sample" in res.json()["detail"]


def test_sample_files_are_documented_with_their_provenance():
    readme = (SAMPLES / "README.md").read_text(encoding="utf-8")
    for sample in app_info.SAMPLES:
        assert sample.accession in readme, sample.accession
        assert sample.filename in readme, sample.filename
