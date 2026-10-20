"""No request makes the engine fail with a stack trace: every answer the window can show is one
plain sentence. The cases come from a fuzzer run against every endpoint before 0.2.0."""

import base64
import logging
import os
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from vettore import fileio
from vettore.server import TOKEN_HEADER, create_app

TOKEN = "t0ken"
H = {TOKEN_HEADER: TOKEN}
S = "ATGC" * 100
TECHNICAL = ("Traceback", "Error", "error(", "NoneType", "object has no attribute", "Errno",
             "pydantic", "Internal Server Error", "KeyError", "invalid literal")


@pytest.fixture(autouse=True)
def private_home(tmp_path, monkeypatch):
    monkeypatch.setenv("VETTORE_HOME", str(tmp_path / "vettore-home"))


def client(**kwargs) -> TestClient:
    return TestClient(create_app(TOKEN, None), base_url="http://127.0.0.1", **kwargs)


def doc(sequence: str, circular: bool = False, **extra) -> dict:
    return {"name": "t", "sequence": sequence, "circular": circular, "features": [],
            "primers": [], **extra}


def plain(res) -> str:
    """The answer is a single string a person can read, not a 500 or a list of objects."""
    assert res.status_code < 500, res.text
    detail = res.json()["detail"]
    assert isinstance(detail, str), detail
    assert not any(marker in detail for marker in TECHNICAL), detail
    return detail


@pytest.mark.parametrize(("path", "body"), [
    ("open-path", {}),
    ("open-path", {"path": None}),
    ("unsaved", {"count": "x"}),
    ("tm", {}),
    ("cloning/gateway", {"reaction": "XX", "inputs": []}),
    ("cloning/mutagenesis", {"doc": doc(S), "kind": "bogus", "start": 10, "end": 10}),
])
def test_a_malformed_request_is_one_plain_sentence(path, body):
    detail = plain(client().post(f"/api/{path}", headers=H, json=body))
    assert "could not read" in detail and "vettore.log" in detail


def test_an_unexpected_problem_is_logged_and_explained(monkeypatch, caplog):
    def broken(*_args, **_kwargs):
        raise RuntimeError("internal detail the user must not see")

    monkeypatch.setattr("vettore.gel.simulate", broken)
    with caplog.at_level(logging.ERROR, logger="vettore.server"):
        res = client(raise_server_exceptions=False).post(
            "/api/gel", headers=H, json={"lanes": []})
    assert res.status_code == 500
    detail = res.json()["detail"]
    assert detail.startswith("Simulating the gel failed because of an unexpected problem "
                             "inside Vettore.")
    assert "vettore.log" in detail and "internal detail" not in detail
    assert "internal detail the user must not see" in caplog.text     # the log has it all


@pytest.mark.skipif(os.name == "nt" or os.geteuid() == 0, reason="needs POSIX permissions")
def test_a_file_vettore_may_not_read_is_explained(tmp_path):
    folder = tmp_path / "private"
    folder.mkdir()
    (folder / "secret.gb").write_text("LOCUS", encoding="utf-8")
    folder.chmod(0)
    try:
        res = client().post("/api/open-path", headers=H,
                            json={"path": str(folder / "secret.gb")})
    finally:
        folder.chmod(0o755)
    assert plain(res).startswith("“secret.gb” could not be opened: Vettore is not allowed")


def test_opening_a_folder_or_nothing_is_explained(tmp_path):
    c = client()
    assert plain(c.post("/api/open-path", headers=H, json={"path": ""})) == (
        "No file was chosen. Choose a sequence file to open.")
    assert plain(c.post("/api/open-path", headers=H, json={"path": str(tmp_path)})).startswith(
        f"“{tmp_path.name}” is a folder, not a file.")


@pytest.mark.parametrize("route", ["save-path", "save-binary"])
def test_saving_without_a_file_name_or_to_a_folder_is_explained(tmp_path, route):
    body = ({"sequenceData": doc("ACGT")} if route == "save-path"
            else {"content_b64": base64.b64encode(b"png").decode()})
    c = client()
    for empty in ("", "/"):
        detail = plain(c.post(f"/api/{route}", headers=H, json={**body, "path": empty}))
        assert detail.startswith("No file name was given.")
    detail = plain(c.post(f"/api/{route}", headers=H, json={**body, "path": str(tmp_path)}))
    assert detail.startswith(f"“{tmp_path.name}” is a folder.")
    detail = plain(c.post(f"/api/{route}", headers=H,
                          json={**body, "path": str(tmp_path / "a\x00b.gb")}))
    assert "cannot be used as a file name" in detail
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("content", ["!!!", "QUJ"])
def test_damaged_base64_is_never_written(tmp_path, content):
    target = tmp_path / "x.png"
    detail = plain(client().post("/api/save-binary", headers=H,
                                 json={"path": str(target), "content_b64": content}))
    assert detail == "“x.png” did not reach the engine intact; save it again."
    assert not target.exists()


def test_features_that_cannot_be_placed_are_left_out_with_a_warning(tmp_path):
    """A save never fails because of one bad annotation: the good ones are written, the others
    are named in a warning."""
    features = [{"name": "good", "start": 0, "end": 9}, {"name": "past the end", "start": 10,
                                                         "end": 500},
                {"name": "no position"}, {"name": "text", "start": "a", "end": 3},
                "not a feature", {"name": "wraps", "start": 30, "end": 5}]
    primers = [{"name": "p", "start": -5, "end": 2}]
    data = doc(S, features=features, primers=primers)
    for name in ("x.gb", "x.dna"):
        res = client().post("/api/save-path", headers=H,
                            json={"path": str(tmp_path / name), "sequenceData": data})
        assert res.status_code == 200, res.text
        warnings = res.json()["warnings"]
        assert warnings[0] == ("5 features were left out of the saved file because their "
                               "positions are not on the sequence: “past the end”, “no "
                               "position”, “text”, “feature”, “wraps”.")
        assert warnings[1].startswith("1 primer was left out of the saved file because its "
                                      "position is")
        [again] = fileio.read_documents(name, (tmp_path / name).read_bytes())
        assert [f["name"] for f in again["features"]] == ["good"]
    res = client().post("/api/export", headers=H, json={"sequenceData": data, "format": "genbank"})
    assert res.status_code == 200 and len(res.json()["warnings"]) == 2
    clean = client().post("/api/save-path", headers=H,
                          json={"path": str(tmp_path / "c.gb"), "sequenceData": doc(S)})
    assert "warnings" not in clean.json()


def test_a_wrapping_feature_is_kept_on_a_circular_sequence(tmp_path):
    data = doc(S, circular=True, features=[{"name": "wraps", "start": 390, "end": 5}])
    res = client().post("/api/save-path", headers=H,
                        json={"path": str(tmp_path / "w.gb"), "sequenceData": data})
    assert "warnings" not in res.json()
    [again] = fileio.read_documents("w.gb", (tmp_path / "w.gb").read_bytes())
    assert [(f["start"], f["end"]) for f in again["features"]] == [(390, 5)]


@pytest.mark.parametrize("sequence", ["ACGTXYZ!", 5])
def test_a_sequence_that_could_not_be_reopened_is_not_saved(tmp_path, sequence):
    for name in ("x.gb", "x.dna", "x.fasta"):
        res = client().post("/api/save-path", headers=H,
                            json={"path": str(tmp_path / name), "sequenceData": doc(sequence)})
        detail = plain(res)
        assert detail.startswith("“t” has “X” at position 5" if sequence != 5
                                 else "“t” could not be saved")
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize(("path", "body", "expected"), [
    ("tm", {"sequence": "A"}, "A Tm needs at least 2 bases"),
    ("alignment/align", {"reference": {"sequence": 123}, "sequences": [doc(S[:200])]},
     "too short to align to"),
    ("alignment/align", {"reference": doc(S), "sequences": [{"sequence": 5},
                                                            {"name": "y", "sequence": None}]},
     "left out: it is shorter than 15 bp"),
    ("cloning/digest", {"doc": doc("HELLO WORLD 123 !!"), "enzymes": ["EcoRI"]},
     "“t” has “E” at position 2"),
    ("cloning/digest", {"doc": doc(S, features=[{"name": "f"}]), "enzymes": ["EcoRI"]},
     "EcoRI does not cut “t”"),
    ("cloning/delete-fragment", {"doc": doc("GAATTC" + S + "GAATTC" + S, True),
                                 "enzymes": ["EcoRI"], "fragment_index": 99},
     "Fragment 100 does not exist (the digest gives 2 fragments)"),
    ("cloning/design-pcr-primers", {"template": doc(S), "start": 100, "end": 200,
                                    "target_tm": 500},
     "The target Tm must be between 40 and 80 °C"),
    ("cloning/design-pcr-primers", {"template": doc(S), "start": -5, "end": 9999},
     "is not on the template “t”, which runs from 1 to 400"),
    ("sanger/align", {"reference": {"sequence": S}, "reads": [{"filename": "a.ab1",
                                                               "content_b64": "QUJ"}]},
     "“a.ab1” did not reach the engine intact"),
    ("sanger/align", {"reference": {"sequence": "Z" * 100}, "reads": []}, "not a DNA letter"),
    ("features/detect", {"sequenceData": {"sequence": 5}}, "too short to annotate"),
    ("export", {"sequenceData": doc("ACGT"), "format": "embl"}, "cannot export"),
])
def test_specific_cases_get_specific_answers(path, body, expected):
    res = client().post(f"/api/{path}", headers=H, json=body)
    if res.status_code == 200:
        assert any(expected in w for w in res.json()["warnings"]), res.json()
    else:
        assert expected in plain(res)


def test_the_log_file_follows_the_vettore_home(tmp_path):
    from vettore import app_info

    assert app_info.log_path() == Path(os.environ["VETTORE_HOME"]) / "vettore.log"
