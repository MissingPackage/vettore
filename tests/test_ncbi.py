"""File → Open from NCBI…: the request sent (the accession only), the record opened, and a
plain sentence for every way it can fail. NCBI itself is replaced by a stand-in: the tests run
offline, like the rest of Vettore."""

import http.client
import io
import socket
import ssl
import urllib.error
import urllib.parse
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from vettore import ncbi
from vettore.server import TOKEN_HEADER, create_app

DATA = Path(__file__).parent / "data"
TOKEN = "t0ken"


@pytest.fixture(autouse=True)
def private_home(tmp_path, monkeypatch):
    monkeypatch.setenv("VETTORE_HOME", str(tmp_path / "vettore-home"))


class Answer:
    """What the stand-in for NCBI answers, and the URLs it was asked for."""

    def __init__(self, monkeypatch, body: bytes = b"", error: Exception | None = None):
        self.body, self.error, self.urls = body, error, []
        monkeypatch.setattr(ncbi.urllib.request, "urlopen", self.urlopen)

    def urlopen(self, url, timeout=None, context=None):
        self.urls.append(url)
        if self.error:
            raise self.error
        return io.BytesIO(self.body)


def fetch(accession: str):
    client = TestClient(create_app(TOKEN, None), base_url="http://127.0.0.1")
    return client.post("/api/ncbi/fetch", headers={TOKEN_HEADER: TOKEN},
                       json={"accession": accession})


def test_a_record_opens_as_an_unsaved_document_named_by_its_accession(monkeypatch):
    Answer(monkeypatch, (DATA / "MT612431.1.gb").read_bytes())
    res = fetch(" mt612431.1 ")
    assert res.status_code == 200
    body = res.json()
    assert body["path"] is None
    assert body["accession"] == "MT612431.1"
    assert body["format"] == "genbank"
    [doc] = body["documents"]
    assert doc["name"] == "MT612431.1"
    assert len(doc["sequence"]) == 7204
    assert doc["circular"] is True
    assert any("GFP" in f["name"] for f in doc["features"])


def test_only_the_accession_is_sent(monkeypatch):
    answer = Answer(monkeypatch, (DATA / "MT612431.1.gb").read_bytes())
    fetch("MT612431.1")
    [url] = answer.urls
    parsed = urllib.parse.urlsplit(url)
    assert parsed.scheme == "https"
    assert parsed.netloc == "eutils.ncbi.nlm.nih.gov"
    assert urllib.parse.parse_qs(parsed.query) == {
        "db": ["nuccore"], "id": ["MT612431.1"], "rettype": ["gbwithparts"], "retmode": ["text"]}


@pytest.mark.parametrize("accession", [
    "MT612431.1", "M77789", "NC_001422.1", "NZ_CP009072.1", "AAAA02000001", "pz629150.1"])
def test_accessions_of_every_kind_are_accepted(accession):
    assert ncbi.clean_accession(accession) == accession.upper()


@pytest.mark.parametrize("typed, message", [
    ("", "Type an accession number"),
    ("   ", "Type an accession number"),
    ("pUC19", "is not a GenBank accession number"),
    ("MT612431.1 or M77789.2", "is not a GenBank accession number"),
    ("ACGTACGTACGT", "is not a GenBank accession number"),
])
def test_what_is_not_an_accession_is_never_sent(monkeypatch, typed, message):
    answer = Answer(monkeypatch, b"LOCUS")
    res = fetch(typed)
    assert res.status_code == 422
    assert message in res.json()["detail"]
    assert answer.urls == []


@pytest.mark.parametrize("error", [
    urllib.error.URLError(socket.gaierror(-3, "Temporary failure in name resolution")),
    urllib.error.URLError(ConnectionRefusedError(111, "Connection refused")),
    TimeoutError("timed out"),
])
def test_offline_says_so_plainly(monkeypatch, error):
    Answer(monkeypatch, error=error)
    res = fetch("MT612431.1")
    assert res.status_code == 503
    detail = res.json()["detail"]
    assert "could not reach NCBI" in detail
    assert "everything else in Vettore works without it" in detail


def test_an_unknown_accession_is_named(monkeypatch):
    Answer(monkeypatch, error=urllib.error.HTTPError(
        "https://eutils.ncbi.nlm.nih.gov", 400, "Bad Request", None, io.BytesIO(b"")))
    res = fetch("ZZ999999.9")
    assert res.status_code == 404
    assert "NCBI has no nucleotide record ZZ999999.9" in res.json()["detail"]


def test_an_error_text_sent_with_status_200_is_not_opened(monkeypatch):
    Answer(monkeypatch, b"Error: F a i l e d  t o  u n d e r s t a n d  i d:  ZZ999999.9\n")
    res = fetch("ZZ999999.9")
    assert res.status_code == 404
    assert "NCBI has no nucleotide record" in res.json()["detail"]


def test_ncbi_busy_asks_to_try_again(monkeypatch):
    Answer(monkeypatch, error=urllib.error.HTTPError(
        "https://eutils.ncbi.nlm.nih.gov", 429, "Too Many Requests", None, io.BytesIO(b"")))
    res = fetch("MT612431.1")
    assert res.status_code == 502
    assert "Try again in a few minutes" in res.json()["detail"]


def test_a_genome_is_refused_with_what_to_do(monkeypatch):
    monkeypatch.setattr(ncbi, "MAX_BYTES", 1000)
    Answer(monkeypatch, b"LOCUS" + b"a" * 2000)
    res = fetch("NC_000001.11")
    assert res.status_code == 422
    assert "too large to open in Vettore" in res.json()["detail"]


def test_an_accession_copied_with_words_around_it_is_picked_out():
    assert ncbi.clean_accession("GenBank: MT612431.1") == "MT612431.1"
    assert ncbi.clean_accession("MT612431.1 (circular)") == "MT612431.1"
    # Only the accession is ever sent, whatever surrounds it.
    assert ncbi.clean_accession("MT612431.1; rm -rf") == "MT612431.1"


@pytest.mark.parametrize("accession", ["NP_000509.1", "XP_011544321", "WP_000123456.1",
                                       "AAB12345.1", "CAA1234567"])
def test_a_protein_accession_is_refused_before_anything_is_sent(monkeypatch, accession):
    answer = Answer(monkeypatch, b"LOCUS")
    res = fetch(accession)
    assert res.status_code == 422
    assert "is a protein accession" in res.json()["detail"]
    assert answer.urls == []


def test_an_accession_without_version_is_named_after_the_version_received(monkeypatch):
    Answer(monkeypatch, (DATA / "MT612431.1.gb").read_bytes())
    body = fetch("MT612431").json()
    assert body["accession"] == "MT612431.1"
    assert body["documents"][0]["name"] == "MT612431.1"


def test_an_outage_page_is_not_reported_as_a_missing_record(monkeypatch):
    Answer(monkeypatch, b"<html><body>Service unavailable</body></html>")
    res = fetch("MT612431.1")
    assert res.status_code == 502
    assert "Try again in a few minutes" in res.json()["detail"]


def test_a_download_cut_half_way_counts_as_offline(monkeypatch):
    Answer(monkeypatch, error=http.client.IncompleteRead(b"LOCUS", 1000))
    res = fetch("MT612431.1")
    assert res.status_code == 503
    assert "could not reach NCBI" in res.json()["detail"]


def test_a_certificate_that_cannot_be_verified_is_named(monkeypatch):
    Answer(monkeypatch, error=urllib.error.URLError(
        ssl.SSLCertVerificationError(1, "certificate verify failed")))
    res = fetch("MT612431.1")
    assert res.status_code == 503
    assert "could not verify NCBI's certificate" in res.json()["detail"]


def test_a_genome_is_refused_from_its_locus_line_before_the_download(monkeypatch):
    class Stream(io.BytesIO):
        read_after_locus = False

        def read(self, *args):
            Stream.read_after_locus = True
            return super().read(*args)

    locus = b"LOCUS       NC_000913            4641652 bp    DNA     circular CON 09-MAR-2022\n"
    monkeypatch.setattr(ncbi.urllib.request, "urlopen",
                        lambda *a, **k: Stream(locus + b"x" * 100))
    res = fetch("NC_000913.3")
    assert res.status_code == 422
    assert "4,641,652 bp" in res.json()["detail"]
    assert not Stream.read_after_locus


def test_a_sequencing_project_master_record_points_to_its_contigs(monkeypatch):
    Answer(monkeypatch, b"LOCUS       AAAA02000000          8 rc    DNA     linear   PLN\n"
                        b"VERSION     AAAA02000000.1\n"
                        b"WGS         AAAA02000001-AAAA02050231\n//\n")
    res = fetch("AAAA02000000")
    assert res.status_code == 422
    assert "master record" in res.json()["detail"]


def test_a_record_that_cannot_be_read_does_not_talk_about_a_file(monkeypatch):
    Answer(monkeypatch, b"LOCUS       MT612431    7204 bp    DNA     circular SYN\n"
                        b"VERSION     MT612431.1\nORIGIN\n        1 zzzz\n")
    res = fetch("MT612431.1")
    assert res.status_code == 422
    detail = res.json()["detail"]
    assert "could not read the record NCBI sent" in detail
    assert "colleague" not in detail
