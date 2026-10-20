"""Smoke-test a packaged Vettore build: start it headless, open and save files through the API.

Usage: python packaging/smoke_test.py <path-to-executable>
"""

import base64
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

PORT = int(os.environ.get("VETTORE_SMOKE_PORT", "8799"))
TOKEN = "smoke-token"
FIXTURES = Path(__file__).resolve().parents[1] / "tests" / "data"
DATA = FIXTURES / "pFA-KanMX4.dna"

# One file per format the engine reads, with the sequence lengths of its records.
EXCHANGE = [
    ("pFA-KanMX4.dna", [3941]),
    ("MT612431.1.gb", [7204]),
    ("exchange-pBR322.embl", [4361]),
    ("exchange-sample-b.xdna", [1000]),
    ("exchange-artificial.gck", [12]),
    ("exchange-plasmids.fasta", [2686, 4361, 4245]),
    ("exchange-pUC19-numbered.seq", [2686]),
    ("exchange-quirks-header.gb", [4361]),
]


def round_trips(folder: Path) -> int:
    """Open each format, save it as GenBank, FASTA and .dna, and open what was saved."""
    checked = 0
    for name, lengths in EXCHANGE:
        opened = call("open-content", {
            "filename": name,
            "content_b64": base64.b64encode((FIXTURES / name).read_bytes()).decode(),
        })
        docs = opened["documents"]
        assert [len(d["sequence"]) for d in docs] == lengths, (name, [len(d["sequence"]) for d in docs])
        doc = docs[0]
        for ext in (".gb", ".fasta", ".dna"):
            saved = call("save-path", {"path": str(folder / f"{Path(name).stem}{ext}"),
                                       "sequenceData": doc})
            [again] = call("open-path", {"path": saved["path"]})["documents"]
            assert again["sequence"].upper() == doc["sequence"].upper(), (name, ext)
            if ext != ".fasta":  # FASTA keeps only the sequence
                assert again["circular"] == doc["circular"], (name, ext)
                assert len(again["features"]) == len(doc["features"]), (name, ext)
            checked += 1
    pasted = call("open-text", {"text": ">pasted\nACGTACGTNN\n"})
    assert pasted["documents"][0]["sequence"] == "ACGTACGTNN"
    return checked


def call(path: str, body: dict | None = None) -> dict:
    req = urllib.request.Request(
        f"http://127.0.0.1:{PORT}/api/{path}",
        data=None if body is None else json.dumps(body).encode(),
        headers={"X-Vettore-Token": TOKEN, "Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=10) as res:
        return json.loads(res.read())


def hand_over(exe: str, env: dict) -> None:
    """A double-click while Vettore runs: ``Vettore <file>`` gives the file to the running copy
    (the server under test, started with --single-instance) and exits. This script plays the
    page that takes it."""
    sample = (FIXTURES / "MT612431.1.gb").resolve()
    second = subprocess.Popen([exe, str(sample)], env=env)
    try:
        taken: list = []
        deadline = time.monotonic() + 120
        while not taken and time.monotonic() < deadline:
            taken = call("startup-files")["results"]
            time.sleep(0.2)
        assert [r["path"] for r in taken] == [str(sample)], taken
        assert second.wait(timeout=60) == 0, second.returncode
    finally:
        if second.poll() is None:   # it opened a window instead of handing the file over
            second.kill()


def main() -> int:
    exe = sys.argv[1]
    work = Path(tempfile.mkdtemp(prefix="vettore-smoke-"))
    # Recent files and the last folder go to a throwaway home, not the runner's.
    env = {**os.environ, "VETTORE_HOME": str(work / "home")}
    proc = subprocess.Popen([exe, "--no-window", "--single-instance", "--port", str(PORT),
                             "--token", TOKEN], env=env)
    try:
        deadline = time.monotonic() + 90
        while True:
            try:
                assert call("health") == {"status": "ok"}
                break
            except OSError:
                if proc.poll() is not None:
                    print(f"executable exited early with code {proc.returncode}")
                    return 1
                if time.monotonic() > deadline:
                    print("server did not start within 90 s")
                    return 1
                time.sleep(1)
        result = call("open-content", {
            "filename": DATA.name,
            "content_b64": base64.b64encode(DATA.read_bytes()).decode(),
        })
        [doc] = result["documents"]
        assert len(doc["sequence"]) == 3941, len(doc["sequence"])
        assert len(doc["features"]) == 9, len(doc["features"])
        with urllib.request.urlopen(f"http://127.0.0.1:{PORT}/index.html", timeout=10) as res:
            assert b'id="root"' in res.read()

        # The cloning engine pulls in pydna/primer3, whose packaging is easy to get wrong.
        vector = "AAAA" + "GAATTC" + "CCCC" * 50 + "GGATCC" + "TTTT" * 50
        insert = "GG" + "GAATTC" + "ACGT" * 25 + "GGATCC" + "GG"
        cloned = call("cloning/restriction-insert", {
            "vector": {"name": "pV", "sequence": vector, "circular": True},
            "vector_enzymes": ["EcoRI", "BamHI"],
            "insert": {"name": "ins", "sequence": insert, "circular": False},
            "insert_enzymes": ["EcoRI", "BamHI"],
        })
        product = cloned["products"][0]["doc"]
        assert product["circular"] and len(product["sequence"]) == 316, len(product["sequence"])
        # The Gateway dialogs pick their inputs by the att sites pydna finds in each document.
        entry = call("open-content", {
            "filename": "cloning-gateway-entry.gb",
            "content_b64": base64.b64encode(
                (FIXTURES / "cloning-gateway-entry.gb").read_bytes()).decode(),
        })["documents"][0]
        sites = call("cloning/gateway-sites", {"docs": [entry]})["sites"]
        assert sites == [["attL1", "attL2"]], sites
        tm = call("tm", {"sequence": "AGCTTGCATGCCTGCAGGTCG"})["tm"]
        assert 55 < tm < 70, tm
        # primer3 reads its thermodynamic tables from files next to its compiled module and
        # swallows the error at import: a build without them only fails when primer design runs.
        status = call("primer-design/status")
        assert status["ready"], status
        designed = call("primer-design/design", {"template": doc, "start": 458, "end": 1267})
        assert len(designed["pairs"]) == 5, designed
        fragment = {"name": "frag", "sequence": doc["sequence"][1000:1600], "circular": False}
        aligned = call("alignment/align", {"reference": doc, "sequences": [fragment]})
        [placed] = aligned["sequences"]
        assert (placed["start"], placed["identity"]) == (1001, 100.0), placed
        # Sanger reads are aligned with Biopython's BLASTN scoring matrix, a data file.
        reference = (FIXTURES.parents[1] / "frontend" / "e2e" / "fixtures"
                     / "sanger-reference.txt").read_text(encoding="utf-8")
        trace = base64.b64encode((FIXTURES / "3730.ab1").read_bytes()).decode()
        sanger = call("sanger/align", {
            "reference": {"name": "ref", "sequence": reference.strip(), "circular": False},
            "reads": [{"filename": "3730.ab1", "content_b64": trace}],
        })
        assert sanger["reads"] and sanger["reads"][0]["identity"] > 99, sanger["reads"]
        # Engine data files (vettore/data) must be in the bundle.
        sample = call("app/samples")["samples"][0]
        sample_doc = call(f"app/samples/{sample['id']}")["documents"][0]
        assert len(sample_doc["sequence"]) > 1000, len(sample_doc["sequence"])
        detected = call("features/detect", {"sequenceData": sample_doc})
        assert detected["catalogue_size"] >= 20 and detected["features"], detected.get("catalogue_size")
        # Opening sequence files on double-click is offered by the packaged app only. Checked
        # without writing anything: a dry run lists the writes, and the status is unchanged.
        assoc = call("desktop/file-associations")
        assert assoc["available"] is True, assoc
        planned = call("desktop/file-associations", {"enabled": True, "dry_run": True})["plan"]
        assert any(Path(exe).name in line for line in planned), planned
        assert all(any(ext in line for line in planned) for ext in assoc["extensions"]), planned
        assert call("desktop/file-associations") == assoc
        hand_over(exe, env)
        saved_files = work / "saved"
        saved_files.mkdir()
        checked = round_trips(saved_files)
        print("smoke test passed: server up, .dna parsed (3941 bp, 9 features), UI served, "
              f"restriction cloning, Tm, primer design, alignment, Sanger reads, samples and feature "
              "detection work, file associations reachable (dry run), a second start hands its "
              f"file to the running copy, {checked} open-and-save "
              "round trips "
              f"({len(EXCHANGE)} formats x GenBank/FASTA/.dna)")
        return 0
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
