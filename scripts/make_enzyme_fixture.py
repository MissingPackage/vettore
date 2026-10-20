"""Write the cross-engine fixture: Biopython cut positions for every bundled enzyme.

The frontend test (frontend/src/__tests__/enzymeDatabase.test.ts) checks that Open Vector
Editor's cut-site engine, fed with frontend/src/data/enzymes.json, cuts at exactly these
positions. Run: python -I scripts/make_enzyme_fixture.py > frontend/src/__tests__/fixtures/cuts.json
"""

import json
import random
import sys
from pathlib import Path

from Bio import Restriction, SeqIO
from Bio.Seq import Seq

from vettore.fileio import DNA_FILE_FORMAT

ROOT = Path(__file__).resolve().parents[1]


def sequences() -> list[dict]:
    rng = random.Random(20261019)
    out = []
    for name, circular in (("pFA-KanMX4.dna", True), ("sample-d.dna", False)):
        record = SeqIO.read(ROOT / "tests" / "data" / name, DNA_FILE_FORMAT)
        out.append({"name": name, "sequence": str(record.seq).upper(), "circular": circular})
    for i in range(3):
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        out.append({"name": f"random-{i}", "sequence": seq, "circular": i == 0})
    return out


def main() -> None:
    db = json.loads((ROOT / "frontend" / "src" / "data" / "enzymes.json").read_text())["enzymes"]
    cases = []
    for s in sequences():
        seq = Seq(s["sequence"])
        cuts = {}
        for entry in db.values():
            enzyme = getattr(Restriction, entry["name"])
            # Biopython: 1-based position of the first base after the top-strand cut.
            positions = enzyme.search(seq, linear=not s["circular"])
            if positions:
                cuts[entry["name"]] = sorted((p - 1) % len(seq) for p in positions)
        cases.append({**s, "cuts": cuts})
    json.dump(cases, sys.stdout, indent=0)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
