import base64
from pathlib import Path

from fastapi.testclient import TestClient

from vettore.server import TOKEN_HEADER, create_app

DATA = Path(__file__).parent / "data"
TOKEN = "t0ken"


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
