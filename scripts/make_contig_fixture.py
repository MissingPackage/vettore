"""Write the reference used by the Sanger contig end-to-end test.

The forward read of the contig fixtures (tests/data/contig-*.ab1, a public forward/reverse
pair) between filler, its N calls replaced by what the reverse read calls there: the sequence
both reads agree on, so the contig in the browser starts from no difference.
Run: python -I scripts/make_contig_fixture.py > frontend/e2e/fixtures/contig-reference.txt
"""

import sys
from pathlib import Path

from vettore.sanger import alignment_payload, mott_trim, read_trace

ROOT = Path(__file__).resolve().parents[1]
FILLER = "ATTGCCTACGGTACATGCAACGTTGACCATAGCATTACGGACCATTGACCAGTTACGGATCAGCTTAGGCATAC"


def main() -> None:
    forward = ROOT / "tests" / "data" / "contig-ACHLO006-09_1_F.ab1"
    reverse = ROOT / "tests" / "data" / "contig-ACHLO006-09_2_R.ab1"
    trace = read_trace(forward.name, forward.read_bytes())
    start, end = mott_trim(trace.quality)
    draft = FILLER * 3 + trace.sequence[start:end] + FILLER * 3
    reads = [(path.name, path.read_bytes()) for path in (forward, reverse)]
    payload = alignment_payload("draft", draft, False, reads,
                                document={"name": "draft", "sequence": draft, "circular": False})
    sys.stdout.write(payload["contig"]["document"]["sequence"] + "\n")


if __name__ == "__main__":
    main()
