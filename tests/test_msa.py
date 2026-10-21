"""Progressive multiple alignment (``vettore.msa``).

Accuracy is measured with the sum-of-pairs score of benchmark suites (BAliBASE: Thompson,
Plewniak & Poch 1999, Bioinformatics 15:87): the fraction of the residue pairs aligned in a
reference alignment that the computed alignment also aligns.

- Simulated families: a random root sequence evolves down a tree with substitutions,
  insertions and deletions, every base keeping the identity of the root column it descends
  from, so the true alignment is known exactly.
- ``tests/data/msa-opuntia.aln``: seven *Opuntia* chloroplast sequences (GenBank
  AF191658-AF191665) and the ClustalW 1.81 alignment published with them in Biopython's test
  suite (``Tests/Clustalw/opuntia.aln``).
- AAV capsid genes: AAV2 VP1 (from ``analysis-NC_001401.2.gb``, NCBI) and AAV9 VP1
  (``analysis-AY530579.1.gb``, NCBI), whose start and stop codons and the start of VP3 must
  line up.
"""

import random
import time
from pathlib import Path

import pytest
from Bio import AlignIO, SeqIO
from cloning_helpers import clean_dna, doc, rc
from fastapi.testclient import TestClient

from vettore import msa
from vettore.alignment import AlignmentError, rotate
from vettore.server import TOKEN_HEADER, create_app

DATA = Path(__file__).parent / "data"
TOKEN = "t0ken"


# ---- simulated evolution with a known true alignment ------------------------------------

def evolve(root: list[tuple[int, str]], rate: float, rng: random.Random, next_id: list[int]
           ) -> list[tuple[int, str]]:
    """Substitutions (70 %), deletions and insertions of 1-6 bases (15 % each)."""
    out = []
    i = 0
    while i < len(root):
        r = rng.random()
        if r < rate * 0.7:
            out.append((root[i][0], rng.choice("ACGT".replace(root[i][1], ""))))
        elif r < rate * 0.85:
            i += rng.randint(1, 6)
            continue
        elif r < rate:
            out.append(root[i])
            for _ in range(rng.randint(1, 6)):
                next_id[0] += 1
                out.append((next_id[0], rng.choice("ACGT")))
        else:
            out.append(root[i])
        i += 1
    return out


def family(n: int, length: int, rate: float, seed: int) -> list[list[tuple[int, str]]]:
    """n sequences in clades of three: each clade's ancestor evolves from the root at
    ``rate`` per base, each sequence from its clade's ancestor at half that."""
    rng = random.Random(seed)
    root = list(enumerate(rng.choice("ACGT") for _ in range(length)))
    next_id = [length]
    leaves: list[list[tuple[int, str]]] = []
    while len(leaves) < n:
        clade = evolve(root, rate, rng, next_id)
        leaves += [evolve(clade, rate / 2, rng, next_id) for _ in range(min(3, n - len(leaves)))]
    return leaves


def aligned_pairs(rows: list[list[int | None]]) -> set[tuple[int, int, int, int]]:
    """(row a, residue of a, row b, residue of b) for every pair of residues in one column."""
    pairs = set()
    for column in zip(*rows):
        present = [(r, k) for r, k in enumerate(column) if k is not None]
        for x in range(len(present)):
            for y in range(x + 1, len(present)):
                pairs.add((*present[x], *present[y]))
    return pairs


def residue_rows(gapped: list[str]) -> list[list[int | None]]:
    rows = []
    for row in gapped:
        k = -1
        out: list[int | None] = []
        for c in row:
            if c == "-":
                out.append(None)
            else:
                k += 1
                out.append(k)
        rows.append(out)
    return rows


def true_rows(leaves: list[list[tuple[int, str]]]) -> list[list[int | None]]:
    ids = sorted({i for leaf in leaves for i, _ in leaf})
    rows = []
    for leaf in leaves:
        where = {i: k for k, (i, _) in enumerate(leaf)}
        rows.append([where.get(i) for i in ids])
    return rows


def sum_of_pairs(reference: list[list[int | None]], test: list[list[int | None]]) -> float:
    expected = aligned_pairs(reference)
    return len(expected & aligned_pairs(test)) / len(expected)


def rows_by_request(result: dict) -> list[str]:
    return [s["row"] for s in sorted(result["sequences"], key=lambda s: s["index"])]


@pytest.mark.parametrize("n, length, rate, floor", [
    (6, 2000, 0.03, 0.97),
    (8, 2000, 0.08, 0.90),
    (12, 1500, 0.10, 0.83),
])
def test_simulated_families_recover_the_true_alignment(n, length, rate, floor):
    leaves = family(n, length, rate, seed=n)
    docs = [doc(f"s{i}", "".join(b for _, b in leaf), False) for i, leaf in enumerate(leaves)]
    result = msa.multiple_alignment(docs)
    rows = rows_by_request(result)
    for row, leaf in zip(rows, leaves):
        assert row.replace("-", "") == "".join(b for _, b in leaf)
    assert sum_of_pairs(true_rows(leaves), residue_rows(rows)) >= floor


def test_opuntia_matches_the_published_clustal_alignment():
    reference = AlignIO.read(DATA / "msa-opuntia.aln", "clustal")
    docs = [doc(r.id, str(r.seq).replace("-", ""), False) for r in reference]
    rows = rows_by_request(msa.multiple_alignment(docs))
    score = sum_of_pairs(residue_rows([str(r.seq) for r in reference]), residue_rows(rows))
    assert score >= 0.99


def _vp1(path: str, start: int | None = None, end: int | None = None) -> str:
    record = SeqIO.read(DATA / path, "genbank")
    return str(record.seq[start:end] if start is not None else record.seq)


def test_aav_capsid_genes_line_up_codon_landmarks():
    aav2 = _vp1("analysis-NC_001401.2.gb", 2202, 4410)        # VP1 CDS 2203..4410
    aav9 = _vp1("analysis-AY530579.1.gb")
    vp3 = aav2[2809 - 2203:]                                   # VP3 starts at 2809
    result = msa.multiple_alignment([doc("AAV2 VP1", aav2, False), doc("AAV9 VP1", aav9, False),
                                     doc("AAV2 VP3", vp3, False)])
    a2, a9, v3 = rows_by_request(result)
    column = lambda row, k: [i for i, c in enumerate(row) if c != "-"][k]
    assert column(a2, 0) == column(a9, 0)                      # ATG of VP1
    assert column(a2, len(aav2) - 3) == column(a9, len(aav9) - 3)   # stop codon
    assert column(v3, 0) == column(a2, 2809 - 2203)           # ATG of VP3 inside VP1
    assert aav9.startswith("ATG") and aav2.startswith("ATG")


# ---- strand, origin, unrelated sequences ----------------------------------------------

BASE = clean_dna(3000, 11)


def test_reverse_complemented_and_rotated_circles_are_turned_to_the_first():
    turned = rc(rotate(BASE, 1000))
    result = msa.multiple_alignment([doc("a", BASE, True), doc("b", turned, True),
                                     doc("fragment", BASE[500:1500], False)])
    by_name = {s["name"]: s for s in result["sequences"]}
    assert by_name["b"]["strand"] == -1
    assert by_name["fragment"]["strand"] == 1
    assert by_name["a"]["row"] == by_name["b"]["row"]
    assert by_name["fragment"]["first_column"] == 501
    assert result["conservation"] == [100] * 3000
    assert result["consensus"] == BASE
    # The consensus track comes first, then the sequences in guide-tree order.
    assert result["alignmentTracks"][0]["sequenceData"]["name"] == "Consensus"


def test_an_unrelated_sequence_is_left_out_with_the_reason():
    result = msa.multiple_alignment([doc("a", BASE, False), doc("b", BASE[:2500], False),
                                     doc("stranger", clean_dna(2000, 99), False)])
    assert [s["name"] for s in result["sequences"]] == ["a", "b"] or \
        sorted(s["name"] for s in result["sequences"]) == ["a", "b"]
    assert any("“stranger” was left out" in w for w in result["warnings"])


def test_a_sequence_unrelated_to_the_first_is_placed_through_another():
    left, right = BASE[:1500], BASE[1500:]
    result = msa.multiple_alignment([doc("left", left, False), doc("right", rc(right), False),
                                     doc("whole", BASE, False)])
    assert len(result["sequences"]) == 3
    assert not result["warnings"]
    by_name = {s["name"]: s for s in result["sequences"]}
    assert by_name["right"]["strand"] == by_name["whole"]["strand"] * -1


def test_limits_are_explained():
    with pytest.raises(AlignmentError, match="at least two"):
        msa.multiple_alignment([doc("a", BASE, False)])
    with pytest.raises(AlignmentError, match="at most 50"):
        msa.multiple_alignment([doc(f"s{i}", BASE[:100], False) for i in range(51)])
    with pytest.raises(AlignmentError, match="15,000 bp"):
        msa.multiple_alignment([doc("a", BASE, False), doc("long", BASE * 6, False)])
    with pytest.raises(AlignmentError, match="add up to"):
        msa.multiple_alignment([doc(f"s{i}", BASE * 4, False) for i in range(13)])


# ---- differences between two rows ------------------------------------------------------

def test_compare_reports_in_the_first_sequence_numbering_on_its_strand():
    variant = BASE[:100] + ("A" if BASE[100] != "A" else "C") + BASE[101:2000] + BASE[2010:]
    result = msa.multiple_alignment([doc("ref", rc(rotate(BASE, 700)), True),
                                     doc("variant", variant, True)])
    by_name = {s["name"]: s for s in result["sequences"]}
    # "variant" is turned to the strand of "ref"; compare in variant's own numbering.
    out = msa.compare(by_name["variant"], by_name["ref"])
    kinds = sorted((d["kind"], d["start"]) for d in out["differences"])
    assert kinds == [("insertion", 2000), ("mismatch", 101)]
    assert out["counts"]["inserted_bases"] == 10
    # And the other way round: ref is reverse-complemented and rotated in its file.
    back = msa.compare(by_name["ref"], by_name["variant"])
    assert back["counts"] == {"mismatches": 1, "insertions": 0, "deletions": 1,
                              "inserted_bases": 0, "deleted_bases": 10}
    deletion = next(d for d in back["differences"] if d["kind"] == "deletion")
    ref_seq = rc(rotate(BASE, 700))
    assert ref_seq[deletion["start"] - 1:deletion["end"]] == deletion["reference"]


# ---- speed (the limits promised in the dialog) -----------------------------------------

@pytest.mark.parametrize("n, length", [(10, 14_500), (50, 2_000)])
def test_the_largest_alignments_finish_in_time(n, length):
    leaves = family(n, length, 0.05, seed=3)
    docs = [doc(f"s{i}", "".join(b for _, b in leaf)[:15_000], False)
            for i, leaf in enumerate(leaves)]
    began = time.perf_counter()
    msa.multiple_alignment(docs)
    assert time.perf_counter() - began < 30


# ---- HTTP: a background job with progress ----------------------------------------------

@pytest.fixture
def client() -> TestClient:
    return TestClient(create_app(TOKEN, None), base_url="http://127.0.0.1")


def test_job_runs_in_the_background_and_reports_progress(client):
    headers = {TOKEN_HEADER: TOKEN}
    body = {"sequences": [doc("a", BASE, False), doc("b", BASE[200:], False)]}
    job = client.post("/api/alignment/multiple/start", json=body, headers=headers).json()["job"]
    for _ in range(200):
        state = client.get(f"/api/alignment/multiple/job/{job}", headers=headers).json()
        if state["state"] != "running":
            break
        assert 0 <= state["progress"] <= 1 and state["message"]
        time.sleep(0.05)
    assert state["state"] == "done"
    rows = state["result"]["sequences"]
    compared = client.post("/api/alignment/multiple/compare", headers=headers,
                           json={"a": rows[0], "b": rows[1]}).json()
    assert compared["differences"] == []
    # Collected once: the job is gone afterwards.
    assert client.get(f"/api/alignment/multiple/job/{job}", headers=headers).status_code == 404


def test_over_the_limits_is_refused_at_once(client):
    response = client.post("/api/alignment/multiple/start", headers={TOKEN_HEADER: TOKEN},
                           json={"sequences": [doc("a", BASE, False)]})
    assert response.status_code == 422
    assert "at least two" in response.json()["detail"]


def test_limits_route(client):
    limits = client.get("/api/alignment/multiple/limits", headers={TOKEN_HEADER: TOKEN}).json()
    assert limits == {"max_sequences": 50, "max_length": 15_000, "max_total": 150_000,
                      "min_length": 15}
