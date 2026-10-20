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

PORT, TOKEN = 8799, "smoke-token"
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


def main() -> int:
    exe = sys.argv[1]
    work = Path(tempfile.mkdtemp(prefix="vettore-smoke-"))
    # Recent files and the last folder go to a throwaway home, not the runner's.
    env = {**os.environ, "VETTORE_HOME": str(work / "home")}
    proc = subprocess.Popen([exe, "--no-window", "--port", str(PORT), "--token", TOKEN], env=env)
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
        tm = call("tm", {"sequence": "AGCTTGCATGCCTGCAGGTCG"})["tm"]
        assert 55 < tm < 70, tm
        saved_files = work / "saved"
        saved_files.mkdir()
        checked = round_trips(saved_files)
        print("smoke test passed: server up, .dna parsed (3941 bp, 9 features), UI served, "
              f"restriction cloning and Tm work, {checked} open-and-save round trips "
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
