import base64
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from vettore.server import TOKEN_HEADER, create_app

DATA = Path(__file__).parent / "data"
TOKEN = "t0ken"


@pytest.fixture(autouse=True)
def private_home(tmp_path, monkeypatch):
    """Recent files and the last folder go to a temporary home, not the user's."""
    monkeypatch.setenv("VETTORE_HOME", str(tmp_path / "vettore-home"))


def client() -> TestClient:
    return TestClient(create_app(TOKEN, None), base_url="http://127.0.0.1")


def test_api_requires_token():
    c = client()
    assert c.get("/api/health").status_code == 403
    assert c.get("/api/health", headers={TOKEN_HEADER: "wrong"}).status_code == 403
    assert c.get("/api/health", headers={TOKEN_HEADER: TOKEN}).status_code == 200


def test_foreign_host_is_rejected():
    c = TestClient(create_app(TOKEN, None), base_url="http://evil.example")
    assert c.get("/api/health", headers={TOKEN_HEADER: TOKEN}).status_code == 403


def test_open_content_and_save(tmp_path):
    c = client()
    content = base64.b64encode((DATA / "pFA-KanMX4.dna").read_bytes()).decode()
    res = c.post("/api/open-content", headers={TOKEN_HEADER: TOKEN},
                 json={"filename": "pFA-KanMX4.dna", "content_b64": content})
    assert res.status_code == 200
    [doc] = res.json()["documents"]
    out = tmp_path / "copy.gb"
    res = c.post("/api/save-path", headers={TOKEN_HEADER: TOKEN},
                 json={"path": str(out), "sequenceData": doc})
    assert res.status_code == 200
    assert out.read_text().startswith("LOCUS")


def test_corrupt_file_returns_422():
    c = client()
    res = c.post("/api/open-content", headers={TOKEN_HEADER: TOKEN},
                 json={"filename": "broken.dna", "content_b64": base64.b64encode(b"junk").decode()})
    assert res.status_code == 422
    assert "broken.dna" in res.json()["detail"]


def test_startup_files_are_served_once(tmp_path):
    bad = tmp_path / "missing.dna"
    c = TestClient(create_app(TOKEN, None, [DATA / "sample-d.dna", bad]),
                   base_url="http://127.0.0.1")
    first = c.get("/api/startup-files", headers={TOKEN_HEADER: TOKEN}).json()
    assert [len(r["documents"]) for r in first["results"]] == [1]
    assert len(first["errors"]) == 1
    second = c.get("/api/startup-files", headers={TOKEN_HEADER: TOKEN}).json()
    assert second == {"results": [], "errors": []}


def test_recent_files_track_open_and_save(tmp_path, monkeypatch):
    monkeypatch.setenv("VETTORE_HOME", str(tmp_path / "home"))
    c = client()
    src = tmp_path / "a.dna"
    src.write_bytes((DATA / "sample-d.dna").read_bytes())
    [doc] = c.post("/api/open-path", headers={TOKEN_HEADER: TOKEN},
                   json={"path": str(src)}).json()["documents"]
    out = tmp_path / "a.gb"
    c.post("/api/save-path", headers={TOKEN_HEADER: TOKEN},
           json={"path": str(out), "sequenceData": doc})
    paths = c.get("/api/recent", headers={TOKEN_HEADER: TOKEN}).json()["paths"]
    assert paths == [str(out), str(src)]


def test_tm_endpoint():
    c = client()
    res = c.post("/api/tm", headers={TOKEN_HEADER: TOKEN}, json={"sequence": "AGCTTGCATGCCTGCAGGTCG"})
    assert res.status_code == 200
    assert 60 < res.json()["tm"] < 72
    assert "SantaLucia" in res.json()["method"]
    assert {p["key"] for p in res.json()["presets"]} == {"standard", "neb_q5", "neb_taq"}
    bad = c.post("/api/tm", headers={TOKEN_HEADER: TOKEN}, json={"sequence": "ACGTNNACGT"})
    assert bad.status_code == 422


def test_save_binary(tmp_path):
    c = client()
    out = tmp_path / "map.png"
    res = c.post("/api/save-binary", headers={TOKEN_HEADER: TOKEN},
                 json={"path": str(out), "content_b64": base64.b64encode(b"\x89PNG").decode()})
    assert res.status_code == 200 and out.read_bytes() == b"\x89PNG"


def test_unsaved_count_is_forwarded():
    seen = []
    c = TestClient(create_app(TOKEN, None, [], seen.append), base_url="http://127.0.0.1")
    assert c.post("/api/unsaved", headers={TOKEN_HEADER: TOKEN}, json={"count": 2}).status_code == 200
    assert seen == [2]


def test_desktop_api_exposes_only_methods_to_javascript():
    # pywebview walks every public non-callable attribute of the js_api object.
    from vettore.app import DesktopApi

    api = DesktopApi()
    public = [n for n in vars(api) if not n.startswith("_")]
    assert public == []


H = {TOKEN_HEADER: TOKEN}


def b64(path: Path) -> str:
    return base64.b64encode(path.read_bytes()).decode()


def test_open_content_reports_format_warnings_and_writability():
    c = client()
    res = c.post("/api/open-content", headers=H, json={
        "filename": "exchange-plasmids.fasta", "content_b64": b64(DATA / "exchange-plasmids.fasta")})
    body = res.json()
    assert res.status_code == 200
    assert body["format"] == "fasta" and body["format_label"] == "FASTA"
    assert [len(d["sequence"]) for d in body["documents"]] == [2686, 4361, 4245]
    assert body["writable"] is False and body["warnings"] == [] and body["path"] is None
    res = c.post("/api/open-content", headers=H, json={
        "filename": "exchange-quirks-features.gb",
        "content_b64": b64(DATA / "exchange-quirks-features.gb")})
    assert len(res.json()["warnings"]) == 2


def test_open_text_recognises_pasted_records():
    c = client()
    genbank = (DATA / "MT612431.1.gb").read_text()
    res = c.post("/api/open-text", headers=H, json={"text": genbank})
    assert res.status_code == 200
    [doc] = res.json()["documents"]
    assert (doc["name"], len(doc["sequence"]), doc["circular"]) == ("MT612431", 7204, True)
    res = c.post("/api/open-text", headers=H, json={"text": "  1 acgtacgt\n 9 gg\n", "name": "ins"})
    assert res.json()["documents"][0]["sequence"] == "acgtacgtgg"
    assert res.json()["documents"][0]["name"] == "ins"
    res = c.post("/api/open-text", headers=H, json={"text": "ACGT\nACQT"})
    assert res.status_code == 422
    assert "line 2, column 3 has “Q”" in res.json()["detail"]


def test_open_errors_are_plain_sentences(tmp_path):
    c = client()
    res = c.post("/api/open-path", headers=H, json={"path": str(tmp_path / "gone.gb")})
    assert res.status_code == 404
    assert res.json()["detail"].startswith("“gone.gb” was not found in")
    broken = tmp_path / "broken.gb"
    broken.write_bytes(b"LOCUS       x 10 bp DNA\nFEATURES             Location/Qualifiers\n"
                       b"     CDS             join(\nORIGIN\n")
    res = c.post("/api/open-path", headers=H, json={"path": str(broken)})
    assert res.status_code == 422
    detail = res.json()["detail"]
    assert detail.startswith("“broken.gb”")
    assert "Error" not in detail and "Traceback" not in detail


def test_save_path_picks_the_format_from_the_extension(tmp_path, monkeypatch):
    monkeypatch.setenv("VETTORE_HOME", str(tmp_path / "home"))
    c = client()
    [doc] = c.post("/api/open-content", headers=H, json={
        "filename": "pFA-KanMX4.dna", "content_b64": b64(DATA / "pFA-KanMX4.dna")}
    ).json()["documents"]

    def save(path, fmt=None):
        res = c.post("/api/save-path", headers=H,
                     json={"path": str(tmp_path / path), "sequenceData": doc, "format": fmt})
        assert res.status_code == 200, res.text
        return res.json()

    assert save("a")["path"] == str(tmp_path / "a.gb")  # extension appended
    assert (tmp_path / "a.gb").read_text().startswith("LOCUS")
    assert save("b", "fasta") == {"path": str(tmp_path / "b.fasta"), "format": "fasta",
                                  "format_label": "FASTA"}
    assert (tmp_path / "b.fasta").read_text().startswith(">pFA-KanMX4")
    assert save("c.fa", "genbank")["format"] == "fasta"  # the typed extension wins
    dna = save("d.dna")
    assert dna["format_label"] == ".dna"
    assert (tmp_path / "d.dna").read_bytes()[:5] == b"\x09\x00\x00\x00\x0e"
    assert save("e", "dna")["path"] == str(tmp_path / "e.dna")
    assert save("pAAV.CMV.v2")["path"] == str(tmp_path / "pAAV.CMV.v2.gb")
    # Reopening what was saved gives the same document back.
    for name in ("a.gb", "d.dna"):
        [again] = c.post("/api/open-path", headers=H,
                         json={"path": str(tmp_path / name)}).json()["documents"]
        assert again["sequence"] == doc["sequence"] and len(again["features"]) == 9
    assert not list(tmp_path.glob(".*vettore-tmp"))
    # Dialogs start in the folder used last.
    from vettore import settings
    assert settings.last_folder() == str(tmp_path)


def test_save_failure_is_explained(tmp_path):
    c = client()
    doc = {"name": "x", "sequence": "ACGT", "circular": False, "features": [], "primers": []}
    res = c.post("/api/save-path", headers=H,
                 json={"path": str(tmp_path / "missing-folder" / "x.gb"), "sequenceData": doc})
    assert res.status_code == 422
    assert res.json()["detail"].startswith("“x.gb” could not be saved in")
    assert "does not exist" in res.json()["detail"]
    res = c.post("/api/save-path", headers=H,
                 json={"path": str(tmp_path / "x.dna"),
                       "sequenceData": {**doc, "sequence": ""}})
    assert res.status_code == 422 and "no sequence" in res.json()["detail"]


def test_save_never_truncates_an_existing_file(tmp_path, monkeypatch):
    c = client()
    target = tmp_path / "keep.gb"
    target.write_text("original")
    doc = {"name": "x", "sequence": "ACGT", "circular": False, "features": [], "primers": []}

    def fail(*_args):
        raise OSError(28, "No space left on device")

    monkeypatch.setattr("vettore.server.os.replace", fail)
    res = c.post("/api/save-path", headers=H, json={"path": str(target), "sequenceData": doc})
    assert res.status_code == 422 and "disk is full" in res.json()["detail"]
    assert target.read_text() == "original"
    assert not list(tmp_path.glob(".*vettore-tmp"))


def test_export_dna_for_browser_download():
    c = client()
    doc = {"name": "x", "sequence": "ACGTACGT", "circular": True, "features": [], "primers": []}
    res = c.post("/api/export", headers=H, json={"sequenceData": doc, "format": "dna"})
    assert base64.b64decode(res.json()["content_b64"])[:5] == b"\x09\x00\x00\x00\x0e"
    res = c.post("/api/export", headers=H, json={"sequenceData": doc, "format": "fasta"})
    assert res.json()["content"].startswith(">x")


def test_startup_file_errors_are_plain(tmp_path):
    junk = tmp_path / "junk.docx"
    junk.write_bytes(b"PK\x03\x04\x00\x00binary")
    c = TestClient(create_app(TOKEN, None, [junk, tmp_path / "gone.dna"]),
                   base_url="http://127.0.0.1")
    errors = c.get("/api/startup-files", headers=H).json()["errors"]
    assert errors[0].startswith("“junk.docx” is not a sequence file Vettore can read")
    assert errors[1].startswith("“gone.dna” was not found")


def test_desktop_dialog_filters_are_valid_for_pywebview():
    from webview.util import parse_file_type

    from vettore import app, fileio

    for spec in app.FILE_TYPES + tuple(app.SAVE_TYPES.values()):
        parse_file_type(spec)  # raises on a filter pywebview rejects
    for ext in fileio.FORMATS:
        assert f"*{ext}" in app.FILE_TYPES[0]
    assert app.save_file_types("fasta")[0].startswith("FASTA")
    assert app.save_file_types("png") == ("PNG image (*.png)",)


def test_dropped_files_are_split_into_paths_and_unresolved_names():
    from vettore.app import DROP_EVENT, drop_script, dropped_paths

    event = {"type": "drop", "dataTransfer": {"files": [
        {"name": "a.gb", "pywebviewFullPath": "/data/a.gb"},
        {"name": "b.dna"},
    ]}}
    payload = dropped_paths(event)
    assert payload == {"paths": ["/data/a.gb"], "unresolved": ["b.dna"]}
    assert dropped_paths({"type": "drop"}) == {"paths": [], "unresolved": []}
    script = drop_script(payload)
    assert script.startswith(f'window.dispatchEvent(new CustomEvent("{DROP_EVENT}"')
    assert '"/data/a.gb"' in script
