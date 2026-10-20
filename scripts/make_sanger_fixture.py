"""Write the reference used by the Sanger end-to-end test.

It embeds a good-quality slice of the real 3730.ab1 read between filler sequence, so the
alignment in the browser must come back at 100 % identity with no discrepancies.
Run: python -I scripts/make_sanger_fixture.py > frontend/e2e/fixtures/sanger-reference.txt
"""

import sys
from pathlib import Path

from vettore.sanger import mott_trim, read_trace

ROOT = Path(__file__).resolve().parents[1]
FILLER = "ATTGCCTACGGTACATGCAACGTTGACCATAGCATTACGGACCATTGACCAGTTACGGATCAGCTTAGGCATAC"


def main() -> None:
    path = ROOT / "tests" / "data" / "3730.ab1"
    trace = read_trace(path.name, path.read_bytes())
    start, _ = mott_trim(trace.quality)
    core = trace.sequence[start + 60 : start + 520]
    sys.stdout.write(FILLER * 4 + core + FILLER * 5 + "\n")


if __name__ == "__main__":
    main()
