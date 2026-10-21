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


# ---- regressions found by the adversarial review -----------------------------------------

def test_strand_comes_from_the_closest_relative_not_an_inverted_shared_element():
    # A and C share only a 150 bp element, inverted in C (a FLEx-like cassette); C matches B
    # end to end. C must be turned by B, forward.
    element, r1, y = clean_dna(150, 41), clean_dna(1500, 42), clean_dna(2000, 43)
    variant = y[:1000] + ("A" if y[1000] != "A" else "C") + y[1001:]
    result = msa.multiple_alignment([doc("A", element + r1, False),
                                     doc("C", variant + rc(element), False),
                                     doc("B", r1 + y, False)])
    by_name = {s["name"]: s for s in result["sequences"]}
    assert by_name["C"]["strand"] == 1
    out = msa.compare(by_name["B"], by_name["C"])
    assert out["counts"]["mismatches"] == 1 and out["counts"]["insertions"] == 0


def test_a_fragment_across_the_origin_of_a_circle_lines_up_in_one_piece():
    fragment = BASE[-600:] + BASE[:600]
    result = msa.multiple_alignment([doc("P", BASE, True), doc("F", fragment, False)])
    by_name = {s["name"]: s for s in result["sequences"]}
    assert result["columns"] == len(BASE)
    out = msa.compare(by_name["F"], by_name["P"])
    assert (out["start"], out["end"], out["overlap"]) == (1, 1200, 1200)
    # The circle was cut elsewhere, and says so.
    assert by_name["P"]["origin_moved_to"] is not None


def test_two_circles_with_an_indel_near_the_origin_start_together():
    target = 2996 - 1005                                    # q[k] is BASE[1005 + k]
    q = list(rotate(BASE, 1000)[5:-5])
    q[target] = "A" if q[target] != "A" else "C"
    result = msa.multiple_alignment([doc("P", BASE, True), doc("Q", "".join(q), True)])
    by_name = {s["name"]: s for s in result["sequences"]}
    out = msa.compare(by_name["P"], by_name["Q"])
    deletion, mismatch = sorted(out["differences"], key=lambda d: d["kind"])
    # BASE[995:1005] deleted (1-based 996-1005); equal flanking bases may slide it along.
    assert deletion["kind"] == "deletion" and len(deletion["reference"]) == 10
    assert 990 <= deletion["start"] <= 1000
    assert (mismatch["kind"], mismatch["start"]) == ("mismatch", 2997)


def test_n_runs_are_unknown_bases_not_indels_or_differences():
    masked = BASE[:1000] + "N" * 200 + BASE[1200:]
    result = msa.multiple_alignment([doc("a", BASE, False), doc("b", masked, False)])
    assert result["columns"] == len(BASE)
    by_name = {s["name"]: s for s in result["sequences"]}
    out = msa.compare(by_name["a"], by_name["b"])
    assert out["differences"] == [] and out["unknown"] == 200


def test_anchors_keep_the_longest_identical_stretch():
    long, blocks = clean_dna(3000, 51), [clean_dna(40, 60 + i) for i in range(4)]
    chain = msa.anchors(long + "".join(blocks), "".join(blocks) + long)
    assert max(n for _, _, n in chain) == 3000


def test_unrelated_inserts_in_one_backbone_are_not_reported_as_mutations():
    backbone = clean_dna(1500, 71)
    gfp, cherry = clean_dna(700, 72), clean_dna(700, 73)
    result = msa.multiple_alignment([doc("pGFP", backbone[:700] + gfp + backbone[700:], False),
                                     doc("pCherry", backbone[:700] + cherry + backbone[700:],
                                         False)])
    by_name = {s["name"]: s for s in result["sequences"]}
    out = msa.compare(by_name["pGFP"], by_name["pCherry"])
    assert out["counts"]["mismatches"] < 5
    [gap] = out["unrelated"]
    assert abs(gap["start"] - 701) < 20 and abs(gap["end"] - 1400) < 20
    assert out["identity"] > 99


def test_a_divergent_homologous_gene_is_compared_all_along():
    aav2 = _vp1("analysis-NC_001401.2.gb", 2202, 4410)
    aav9 = _vp1("analysis-AY530579.1.gb")
    result = msa.multiple_alignment([doc("AAV2", aav2, False), doc("AAV9", aav9, False)])
    by_name = {s["name"]: s for s in result["sequences"]}
    out = msa.compare(by_name["AAV2"], by_name["AAV9"])
    assert out["unrelated"] == [] and (out["start"], out["end"]) == (1, len(aav2))


def test_conservation_never_counts_gaps_as_agreement():
    inserted = BASE[:500] + "ACGTACGTAC" + BASE[500:1500]
    rows = [doc(f"s{i}", BASE[:1500], False) for i in range(9)] + [doc("ins", inserted, False)]
    result = msa.multiple_alignment(rows)
    low = [c for c in result["conservation"] if c is not None and c < 50]
    assert len(low) == 10 and set(low) == {10}


def test_columns_covered_by_one_sequence_have_no_conservation():
    result = msa.multiple_alignment([doc("whole", BASE, False), doc("part", BASE[1000:2000],
                                                                     False)])
    assert result["conservation"][:1000] == [None] * 1000
    assert result["conservation"][1000:2000] == [100] * 1000


def test_a_new_alignment_stops_the_one_still_running():
    first = msa.start([doc("a", BASE, False), doc("b", BASE, False)])
    second = msa.start([doc("a", BASE, False), doc("b", BASE, False)])
    assert msa.JOBS[first].cancelled or msa.JOBS[first].done.is_set()
    msa.JOBS[second].done.wait(30)


# ---- regressions found by the second adversarial review ----------------------------------

def _two(a_seq, b_seq, circular=True):
    result = msa.multiple_alignment([doc("a", a_seq, circular), doc("b", b_seq, circular)])
    by_name = {s["name"]: s for s in result["sequences"]}
    return result, msa.compare(by_name["a"], by_name["b"])


def snp(seq: str, i: int) -> str:
    return seq[:i] + ("A" if seq[i] != "A" else "C") + seq[i + 1:]


def test_plasmids_with_inserts_of_different_length_are_compared_all_round():
    backbone = clean_dna(3000, 81)
    gfp, cherry = clean_dna(720, 82), clean_dna(1400, 83)
    mutated = snp(snp(backbone, 120), 2500)
    result, out = _two(backbone[:1000] + gfp + backbone[1000:],
                       mutated[:1000] + cherry + mutated[1000:])
    assert result["columns"] < 3000 + 720 + 1400
    starts = sorted(d["start"] for d in out["differences"] if d["kind"] == "mismatch")
    assert starts == [121, 2501 + 720]
    [gap] = out["unrelated"]
    assert abs(gap["start"] - 1001) < 5 and gap["a_bases"] > 700


@pytest.mark.parametrize("change", ["snp1", "del5", "del50", "ins50"])
def test_a_difference_next_to_base_1_of_two_circles_is_never_lost(change):
    marker = 2000
    other = snp(BASE, marker)
    other = {"snp1": snp(other, 0), "del5": other[:4] + other[14:],
             "del50": other[:20] + other[70:],
             "ins50": other[:20] + clean_dna(50, 84) + other[20:]}[change]
    _, out = _two(BASE, other)
    kinds = sorted(d["kind"] for d in out["differences"])
    expected = {"snp1": ["mismatch", "mismatch"], "del5": ["deletion", "mismatch"],
                "del50": ["deletion", "mismatch"], "ins50": ["insertion", "mismatch"]}[change]
    assert kinds == expected and out["unrelated"] == []
    assert (out["start"], out["end"]) == (1, len(BASE))


@pytest.mark.parametrize("size", [80, 150, 400])
def test_long_deletions_are_deletions_not_unrelated_stretches(size):
    _, out = _two(BASE, BASE[:1000] + BASE[1000 + size:], circular=False)
    assert out["counts"]["deletions"] == 1 and out["counts"]["deleted_bases"] == size
    assert out["unrelated"] == []


def test_unrelated_ends_of_linear_sequences_are_named():
    shared = BASE[:1200]
    _, out = _two(clean_dna(800, 85) + shared, clean_dna(800, 86) + shared, circular=False)
    [gap] = out["unrelated"]
    assert gap["start"] < 200 and gap["end"] > 700   # before it, b has not begun yet


def test_circles_open_inside_what_they_share_and_fragments_stay_whole():
    # P's base 1 is not in Q (Q holds P as a linear stretch): the circles are opened inside
    # P, which both share, and a fragment of Q stays in one piece.
    p = clean_dna(1500, 87)
    q = clean_dna(2500, 88) + p + clean_dna(2000, 89)
    fragment = q[3700:4300]
    result = msa.multiple_alignment([doc("P", p, True), doc("Q", q, True),
                                     doc("F", fragment, False)])
    by_name = {s["name"]: s for s in result["sequences"]}
    row = by_name["F"]["row"]
    assert "-" not in row[by_name["F"]["first_column"] - 1:by_name["F"]["last_column"]]
    out = msa.compare(by_name["P"], by_name["Q"])
    [extra] = out["differences"]           # Q is P with 4.5 kb more: one insertion
    assert extra["kind"] == "insertion" and len(extra["sequence"]) == 4500
    assert out["overlap"] == len(p)


def test_identical_circles_keep_their_origin():
    result = msa.multiple_alignment([doc("a", BASE, True), doc("b", snp(BASE, 1500), True)])
    assert all(s["origin_moved_to"] is None for s in result["sequences"])


@pytest.mark.parametrize("insert_at", [0, 200, 2800])
def test_a_root_starting_in_or_near_an_insert_the_other_lacks(insert_at):
    insert = clean_dna(700, 90)
    other = snp(snp(BASE, 10), 2900)
    _, out = _two(BASE[:insert_at] + insert + BASE[insert_at:], other)
    mism = sorted(d["start"] for d in out["differences"] if d["kind"] == "mismatch")
    expected = sorted(k + 1 + (700 if k >= insert_at else 0) for k in (10, 2900))
    assert mism == expected
    deletions = [d for d in out["differences"] if d["kind"] == "deletion"]
    assert len(deletions) == 1 and len(deletions[0]["reference"]) == 700
    assert out["unrelated"] == []


# ---- regressions found by the third adversarial review -----------------------------------

@pytest.mark.parametrize("change, size", [("ins", 500), ("del", 290), ("ins", 310),
                                          ("del", 500)])
def test_an_indel_next_to_base_1_of_two_circles_is_one_event(change, size):
    other = (BASE[:20] + clean_dna(size, 91) + BASE[20:] if change == "ins"
             else BASE[:20] + BASE[20 + size:])
    _, out = _two(BASE, other)
    [event] = out["differences"]
    assert event["kind"] == ("insertion" if change == "ins" else "deletion")
    assert len(event["sequence"] or event["reference"]) == size
    assert out["unrelated"] == []


def test_an_adapter_tail_is_not_compared_rather_than_listed_as_mutations():
    _, out = _two(BASE[500:1400] + clean_dna(150, 92), BASE, circular=False)
    assert out["differences"] == []
    [tail] = out["unrelated"]
    assert tail["start"] == 901 and tail["a_bases"] == 150


@pytest.mark.parametrize("stuffer", [30, 95, 120])
def test_an_insert_replaced_by_an_unrelated_stuffer_is_one_named_replacement(stuffer):
    a = BASE[:1500] + clean_dna(1000, 93) + BASE[1500:]
    b = BASE[:1500] + clean_dna(stuffer, 94) + BASE[1500:]
    _, out = _two(a, b, circular=False)
    assert out["differences"] == []
    [swap] = out["unrelated"]
    assert abs(swap["start"] - 1501) <= 2 and swap["a_bases"] >= 995


# ---- regressions found by the fourth adversarial review ----------------------------------

def _substitute(seq: str, rate: float, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT".replace(c, "")) if rng.random() < rate else c for c in seq)


@pytest.mark.parametrize("circular", [True, False])
def test_two_members_never_get_an_indel_the_others_needed(circular):
    # Substitutions only, plus a third member with a real deletion: the gaps the others
    # carry for each other must not read as an insertion and a deletion between two.
    rng = random.Random(41)
    members = [_substitute(BASE[:2500], rate, 400 + i) for i, rate in enumerate(
        (0.01, 0.03, 0.06, 0.03, 0.06))]
    members.append(BASE[:1200] + BASE[1230:2500])
    docs = [doc(f"m{i}", rotate(s, rng.randrange(len(s))) if circular else s, circular)
            for i, s in enumerate(members)]
    by_name = {s["name"]: s for s in msa.multiple_alignment(docs)["sequences"]}
    for i in range(5):
        for j in range(i + 1, 5):
            out = msa.compare(by_name[f"m{i}"], by_name[f"m{j}"])
            assert {d["kind"] for d in out["differences"]} <= {"mismatch"}, (i, j)


@pytest.mark.parametrize("circular", [True, False])
@pytest.mark.parametrize("at", [0, 1, 797, 798, 799])
def test_a_snp_at_the_end_of_an_amplicon_is_reported(circular, at):
    amplicon = snp(BASE[1000:1800], at)
    result = msa.multiple_alignment([doc("a", BASE, circular), doc("b", amplicon, False)])
    by_name = {s["name"]: s for s in result["sequences"]}
    out = msa.compare(by_name["a"], by_name["b"])
    assert [(d["kind"], d["start"]) for d in out["differences"]] == [("mismatch", 1001 + at)]
    assert out["unrelated"] == []


def test_a_divergent_homologous_end_is_compared_not_named_unrelated():
    other = _substitute(BASE[:1920], 0.1, 42) + _substitute(BASE[1920:2000], 0.35, 43)
    _, out = _two(BASE[:2000], other, circular=False)
    assert out["unrelated"] == []
    assert (out["start"], out["end"]) == (1, 2000)
    # At 35 % divergence a pair alignment may trade a few mismatches for an indel pair;
    # what matters is that the end's differences are listed, not hidden.
    expected = sum(1 for x, y in zip(BASE[:2000], other) if x != y)
    assert out["counts"]["mismatches"] >= expected - 6
    assert any(d["start"] > 1950 for d in out["differences"])


@pytest.mark.parametrize("seed", range(8))
def test_unrelated_inserts_leave_no_differences_at_the_junctions(seed):
    a = BASE[:1500] + clean_dna(600, 600 + seed) + BASE[1500:]
    b = BASE[:1500] + clean_dna(400, 700 + seed) + BASE[1500:]
    _, out = _two(a, b, circular=False)
    assert out["differences"] == []
    [swap] = out["unrelated"]
    assert swap["a_bases"] >= 590


# ---- regressions found by the fifth adversarial review -----------------------------------


@pytest.mark.parametrize("seed", range(12))
def test_divergent_linear_pairs_are_compared_without_error(seed):
    # Seed 0 is the reviewer's case: realigning the first stretch moved gaps to the front.
    rng = random.Random(7)
    a = "".join(rng.choice("ACGT") for _ in range(2100)) if seed == 0 else clean_dna(1500, 800 + seed)
    b = _substitute(a, 0.15, 31) if seed == 0 else _substitute(a, 0.12 + 0.005 * seed, 900 + seed)
    for x, y in ((a, b), (b, a)):
        _, out = _two(x, y, circular=False)
        assert out["start"] is not None and out["unrelated"] == []
        assert out["counts"]["mismatches"] > 100


@pytest.mark.parametrize("side", ["5'", "3'"])
def test_a_short_barcode_tail_is_not_compared_rather_than_listed(side):
    barcode = clean_dna(10, 95)
    amplicon = barcode + BASE[1000:1800] if side == "5'" else BASE[1000:1800] + barcode
    _, out = _two(amplicon, BASE, circular=False)
    assert out["differences"] == []
    [tail] = out["unrelated"]
    assert tail["a_bases"] == 10
