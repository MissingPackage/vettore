"""S5.2 Sanger contig: several reads of one clone, one consensus over the reference.

Real traces: a forward and a reverse read of the same template, published with
sangeranalyseR (provenance in tests/data/README.md), and Biopython's 3730.ab1. Known
differences are planted in the reference, so the clone (what the reads show) stays real.
"""

from functools import cache
from pathlib import Path

from cloning_helpers import clean_dna, doc, feature, rc

from vettore.contig import CONFIRMING_QUALITY, _read_calls
from vettore.sanger import Trace, align_trace, alignment_payload, mott_trim, read_trace

DATA = Path(__file__).parent / "data"
FORWARD = DATA / "contig-ACHLO006-09_1_F.ab1"
REVERSE = DATA / "contig-ACHLO006-09_2_R.ab1"
UNRELATED = DATA / "3730.ab1"
FILLER = "ATTGCCTACGGTACATGCAACGTTGACCATAGCATTACGGACCATTGACCAGTTACGGATCAGCTTAGGCATAC"
LEFT = 3 * len(FILLER)                     # the clone starts at reference index 222


def pair() -> list[tuple[str, bytes]]:
    return [(FORWARD.name, FORWARD.read_bytes()), (REVERSE.name, REVERSE.read_bytes())]


def contig(reference: str, circular: bool = False, reads=None, features=()) -> dict:
    payload = alignment_payload("pCOI", reference, circular, reads or pair(),
                                document=doc("pCOI", reference, circular, features))
    return payload["contig"]


@cache
def clone() -> str:
    """The forward read's trimmed calls between filler, its N calls replaced by what the
    reverse read calls there: the sequence both reads agree on."""
    forward = read_trace(FORWARD.name, FORWARD.read_bytes())
    start, end = mott_trim(forward.quality)
    draft = FILLER * 3 + forward.sequence[start:end] + FILLER * 3
    return contig(draft)["document"]["sequence"]


@cache
def forward_calls() -> dict:
    forward = read_trace(FORWARD.name, FORWARD.read_bytes())
    aligned = align_trace(clone(), False, forward)
    return _read_calls(aligned, len(clone()))[0]


def strand_runs(result: dict, strands: int) -> list[tuple[int, int]]:
    return [(r["start"], r["end"]) for r in result["coverage"] if r["strands"] == strands]


def mutate(sequence: str, index: int) -> tuple[str, str]:
    replacement = "A" if sequence[index] != "A" else "G"
    return sequence[:index] + replacement + sequence[index + 1:], replacement


# ---- real forward and reverse reads ------------------------------------------------------

def test_the_pair_covers_the_insert_on_both_strands_and_one_end_on_one():
    result = contig(clone())
    assert result["reads"] == 2
    [(both_start, both_end)] = strand_runs(result, 2)
    [(one_start, one_end)] = strand_runs(result, 1)
    assert both_start == LEFT + 1 and one_start == both_end + 1
    assert strand_runs(result, 0) == [(1, LEFT), (one_end + 1, len(clone()))]
    assert result["covered"] == [{"start": LEFT + 1, "end": one_end}]


def test_reads_that_agree_with_the_reference_give_it_back_unchanged():
    result = contig(clone())
    assert result["differences"] == []
    assert result["document"]["sequence"] == clone()
    # The reads only differ where one of them makes no call (N): nothing to settle.
    assert {d["kind"] for d in result["disagreements"]} <= {"no call"}


def test_a_difference_both_reads_show_is_confirmed():
    index = LEFT + 300                                   # inside the two-strand stretch
    reference, replacement = mutate(clone(), index)
    [difference] = contig(reference)["differences"]
    assert difference == {"kind": "mismatch", "position": index + 1, "reference": replacement,
                          "consensus": clone()[index], "reads": 2,
                          "quality": difference["quality"], "confirmed": True}


def test_one_read_confirms_a_difference_only_at_q30_or_better():
    result = contig(clone())
    [(start, end)] = strand_runs(result, 1)
    single = [p for p in range(start - 1, end - 10) if p in forward_calls()
              and forward_calls()[p].base != "N"]
    good = next(p for p in single if forward_calls()[p].quality >= CONFIRMING_QUALITY + 5)
    weak = next(p for p in single if forward_calls()[p].quality < CONFIRMING_QUALITY)
    for index, confirmed in ((good, True), (weak, False)):
        reference, _ = mutate(clone(), index)
        [difference] = contig(reference)["differences"]
        assert difference["position"] == index + 1
        assert difference["reads"] == 1
        assert difference["quality"] == forward_calls()[index].quality
        assert difference["confirmed"] is confirmed


def test_an_extra_base_in_the_clone_moves_the_features_after_it():
    index = LEFT + 400
    reference = clone()[:index] + clone()[index + 1:]   # the clone carries one more base
    before = feature("before", 10, 100)
    after = feature("after", index + 50, index + 120)
    across = feature("across", index - 20, index + 20)
    result = contig(reference, features=[before, after, across])
    [difference] = result["differences"]
    assert difference["kind"] == "insertion" and difference["confirmed"]
    document = result["document"]
    assert document["sequence"] == clone()
    by_name = {f["name"]: f for f in document["features"]}
    assert (by_name["before"]["start"], by_name["before"]["end"]) == (10, 100)
    assert (by_name["after"]["start"], by_name["after"]["end"]) == (index + 51, index + 121)
    assert (by_name["across"]["start"], by_name["across"]["end"]) == (index - 20, index + 21)
    [variation] = [f for f in document["features"] if f["type"] == "variation"]
    inserted = difference["consensus"]
    assert document["sequence"][variation["start"]:variation["end"] + 1] == inserted
    assert variation["name"] == f"+{inserted} (Sanger)"


def test_bases_missing_from_the_clone_are_one_deletion():
    index = LEFT + 250
    reference = clone()[:index] + "GTC" + clone()[index:]
    result = contig(reference)
    [difference] = result["differences"]
    assert difference["kind"] == "deletion" and difference["confirmed"]
    assert len(difference["reference"]) == 3
    assert result["document"]["sequence"] == clone()
    [variation] = [f for f in result["document"]["features"] if f["type"] == "variation"]
    assert variation["end"] - variation["start"] == 1         # the two bases around the gap


def test_verified_stretches_leave_out_the_positions_no_read_called():
    """The forward read calls N at two places only it covers: the stretch is verified around
    them, and they are marked instead."""
    result = contig(clone())
    assert result["document"]["name"] == "pCOI consensus"
    features = result["document"]["features"]
    verified = sorted((f["start"], f["end"], f["name"]) for f in features
                      if f["name"].startswith("Sanger verified"))
    unclear = sorted(f["start"] for f in features if f["name"] == "Sanger: no clear call")
    assert result["unclear"] == 2 and len(unclear) == 2
    [covered] = result["covered"]
    assert verified[0][0] == covered["start"] - 1 and verified[-1][1] == covered["end"] - 1
    assert verified[0][2] == "Sanger verified (1 read)"  # one stretch is read by one read
    # Bases backed only by a call under Q20 are not verified either.
    weak = [p for p, c in forward_calls().items() if p >= both_end(result)
            and c.base != "N" and c.quality < 20]
    assert weak and result["doubtful"] == len(weak)
    for position in unclear + weak:
        assert not any(start <= position <= end for start, end, _ in verified)


def both_end(result: dict) -> int:
    [(_, end)] = strand_runs(result, 2)
    return end


def test_a_difference_to_check_is_annotated_as_such_and_not_verified():
    result = contig(clone())
    [(start, end)] = strand_runs(result, 1)
    weak = next(p for p in range(start - 1, end - 10) if p in forward_calls()
                and forward_calls()[p].base != "N"
                and forward_calls()[p].quality < CONFIRMING_QUALITY)
    reference, _ = mutate(clone(), weak)
    features = contig(reference)["document"]["features"]
    [variation] = [f for f in features if f["type"] == "variation"]
    assert variation["name"].endswith("(Sanger, check the trace)")
    assert variation["start"] == weak
    assert not any(f["start"] <= weak <= f["end"] for f in features
                   if f["name"].startswith("Sanger verified"))


def test_the_same_trace_chosen_twice_is_counted_once():
    reads = [*pair(), ("copy.ab1", FORWARD.read_bytes())]
    payload = alignment_payload("pCOI", clone(), False, reads)
    [warning] = payload["warnings"]
    assert "copy.ab1" in warning and "counted once" in warning
    assert payload["contig"]["reads"] == 2


def test_a_read_from_another_plasmid_is_listed_apart_and_left_out():
    reads = [*pair(), (UNRELATED.name, UNRELATED.read_bytes())]
    payload = alignment_payload("pCOI", clone(), False, reads,
                                document=doc("pCOI", clone(), False))
    [warning] = payload["warnings"]
    assert "probably comes from another plasmid" in warning
    assert payload["contig"]["reads"] == 2
    assert payload["contig"]["differences"] == []


def test_a_circular_reference_wraps_the_coverage_across_the_origin():
    shift = LEFT + 200                                   # the origin now falls inside the reads
    reference = clone()[shift:] + clone()[:shift]
    result = contig(reference, circular=True)
    assert result["differences"] == []
    [(start, end)] = [(r["start"], r["end"]) for r in result["covered"]]
    assert start > end                                   # one stretch, across the origin
    assert result["document"]["sequence"] == reference


# ---- the voting rules, on reads whose qualities are set by hand ----------------------------

REFERENCE = clean_dna(200, 301) + clean_dna(300, 302) + clean_dna(200, 303)


def made_read(name: str, sequence: str, quality: int | list[int], strand: int = 1) -> Trace:
    qualities = quality if isinstance(quality, list) else [quality] * len(sequence)
    peaks = [{"aTrace": [1], "tTrace": [1], "gTrace": [1], "cTrace": [1]}] * len(sequence)
    trace = Trace(name, sequence, qualities, peaks)
    return trace if strand == 1 else trace.reverse_complement()


def made_contig(*traces: Trace) -> dict:
    from vettore.contig import contig_payload

    reads = [align_trace(REFERENCE, False, t) for t in traces]
    return contig_payload(REFERENCE, False, reads, doc("pMade", REFERENCE, False))


def with_base(sequence: str, index: int, base: str) -> str:
    return sequence[:index] + base + sequence[index + 1:]


def test_reads_that_disagree_are_settled_by_the_better_call_and_flagged():
    insert = REFERENCE[200:500]
    other = "A" if insert[150] != "A" else "C"
    good = made_read("good", insert, 45)
    bad = made_read("bad", with_base(insert, 150, other), [45] * 150 + [12] + [45] * 149, -1)
    result = made_contig(good, bad)
    assert result["differences"] == []                   # the Q45 call is the reference base
    [flag] = result["disagreements"]
    assert flag["kind"] == "base" and flag["position"] == 351
    assert flag["consensus"] == insert[150]
    assert {(c["read"], c["call"], c["quality"]) for c in flag["calls"]} == {
        ("good", insert[150], 45), ("bad", other, 12)}


def test_a_lone_low_quality_call_is_to_be_checked():
    insert = REFERENCE[200:500]
    other = "A" if insert[150] != "A" else "C"
    quality = [45] * 150 + [CONFIRMING_QUALITY - 1] + [45] * 149
    [difference] = made_contig(made_read("only", with_base(insert, 150, other), quality))[
        "differences"]
    assert difference["confirmed"] is False and difference["reads"] == 1


def test_two_reads_confirm_a_difference_together_only_at_q20_or_better():
    insert = REFERENCE[200:500]
    other = "A" if insert[150] != "A" else "C"
    mutated = with_base(insert, 150, other)
    for low, confirmed in ((22, True), (15, False)):
        quality = [45] * 150 + [low] + [45] * 149
        [difference] = made_contig(made_read("a", mutated, quality),
                                   made_read("b", mutated, quality, -1))["differences"]
        assert difference["reads"] == 2
        assert difference["confirmed"] is confirmed


def test_many_good_reads_outvote_one_better_read():
    """Summed quality, as in phrap: three Q35 reads of the reference base beat one Q40 read."""
    insert = REFERENCE[200:500]
    other = "A" if insert[150] != "A" else "C"
    result = made_contig(made_read("a", insert, 35), made_read("b", insert, 35, -1),
                         made_read("c", insert[:290], 35),
                         made_read("odd", with_base(insert, 150, other), 40, -1))
    assert result["differences"] == []
    [flag] = result["disagreements"]
    assert flag["consensus"] == insert[150]


def test_a_difference_a_good_read_contradicts_is_not_confirmed():
    insert = REFERENCE[200:500]
    other = "A" if insert[150] != "A" else "C"
    result = made_contig(made_read("alt", with_base(insert, 150, other), 50),
                         made_read("ref", insert, [45] * 150 + [25] + [45] * 149, -1))
    [difference] = result["differences"]
    assert difference["consensus"] == other and difference["confirmed"] is False


def test_a_read_never_wins_half_of_its_own_deletion():
    """Read A (Q40) lacks 2 bases; read B (Q45, then Q35) has them. Voted base by base the
    consensus kept one and lost the other: a 1 bp frameshift no read shows."""
    insert = REFERENCE[200:500]
    deleted = insert[:150] + insert[152:]
    b_quality = [40] * 150 + [45, 35] + [40] * 148
    result = made_contig(made_read("A", deleted, 40), made_read("B", insert, b_quality, -1))
    # Read A's whole deletion (Q40 at its flanks) against read B's weaker of the two bases
    # (Q35): the deletion wins as one event, and B, a good read, keeps it unconfirmed.
    [difference] = result["differences"]
    assert difference["kind"] == "deletion" and difference["reference"] == insert[150:152]
    assert difference["confirmed"] is False
    assert len(result["document"]["sequence"]) == len(REFERENCE) - 2


def test_a_mixed_iupac_call_is_no_call():
    insert = REFERENCE[200:500]
    result = made_contig(made_read("a", with_base(insert, 150, "R"), 40),
                         made_read("b", with_base(insert, 150, "R"), 40, -1))
    assert result["differences"] == []
    assert [d["kind"] for d in result["disagreements"]] == ["no call"]


def test_a_deletion_across_the_origin_is_one_event():
    from vettore.contig import contig_payload

    reference = clean_dna(600, 304)
    clone_ = reference[2:-2]                     # 4 bases missing around the origin
    read = clone_[-250:] + clone_[:250]
    reads = [align_trace(reference, True, t) for t in (made_read("a", read, 40),
                                                         made_read("b", read, 40, -1))]
    result = contig_payload(reference, True, reads, doc("pRing", reference, True))
    [difference] = result["differences"]
    assert difference["kind"] == "deletion"
    assert difference["reference"] == reference[-2:] + reference[:2]
    assert len(result["document"]["sequence"]) == len(reference) - 4


def test_an_n_is_no_call_and_never_enters_the_consensus():
    insert = REFERENCE[200:500]
    # Base callers give an N a very low quality; a good-quality N is doubt (tested below).
    result = made_contig(made_read("a", with_base(insert, 150, "N"), [40] * 150 + [4] + [40] * 149),
                         made_read("b", insert, 20, -1))
    assert result["differences"] == []
    assert result["disagreements"] == []
    alone = made_contig(made_read("a", with_base(insert, 150, "N"), [40] * 150 + [4] + [40] * 149))
    assert alone["differences"] == []
    assert [d["kind"] for d in alone["disagreements"]] == ["no call"]
    assert alone["document"]["sequence"] == REFERENCE


def test_coverage_counts_strands_not_reads():
    insert = REFERENCE[200:500]
    result = made_contig(made_read("a", insert[:200], 40), made_read("b", insert[:200], 40),
                         made_read("c", rc(insert[100:]), 40))
    assert [(r["start"], r["end"], r["strands"]) for r in result["coverage"]] == [
        (1, 200, 0), (201, 300, 1), (301, 400, 2), (401, 500, 1), (501, 700, 0)]


def test_e2e_contig_fixture_is_up_to_date():
    import subprocess
    import sys

    root = Path(__file__).resolve().parents[1]
    fixture = root / "frontend" / "e2e" / "fixtures" / "contig-reference.txt"
    out = subprocess.run([sys.executable, "-I", str(root / "scripts" / "make_contig_fixture.py")],
                         capture_output=True, text=True, encoding="utf-8", check=True).stdout
    assert out == fixture.read_text(encoding="utf-8"), (
        "Regenerate: python -I scripts/make_contig_fixture.py > "
        "frontend/e2e/fixtures/contig-reference.txt")
    assert out.strip() == clone()


# ---- contrary evidence that does not get a full vote (second review) -----------------------

def verified_at(result: dict, position: int) -> bool:
    """Whether reference position (0-based) lies in a verified stretch of the consensus."""
    return any(r["start"] - 1 <= position <= r["end"] - 1 for r in result["verified"])


def test_a_read_ending_inside_a_deletion_still_speaks_against_it():
    insert = REFERENCE[200:500]
    deleted = insert[:150] + insert[154:]                 # A lacks reference 350–353
    partial = insert[:152]                                # B ends at 351, showing 350 and 351
    result = made_contig(made_read("A", deleted, 40), made_read("B", partial, 40, -1))
    [difference] = result["differences"]
    assert difference["kind"] == "deletion" and difference["confirmed"] is False
    assert [d["kind"] for d in result["disagreements"]] == ["base"]
    assert not verified_at(result, 350)


def test_a_read_with_an_n_inside_a_deletion_still_speaks_against_it():
    insert = REFERENCE[200:500]
    deleted = insert[:150] + insert[152:]                 # A lacks reference 350–351
    with_n = with_base(insert, 151, "N")                  # B shows 350, and N at 351
    result = made_contig(made_read("A", deleted, 40), made_read("B", with_n, 40, -1))
    [difference] = result["differences"]
    assert difference["confirmed"] is False
    assert result["disagreements"]


def test_a_tie_between_good_reads_is_not_verified():
    insert = REFERENCE[200:500]
    other = "A" if insert[150] != "A" else "C"
    result = made_contig(made_read("ref", insert, 30),
                         made_read("alt", with_base(insert, 150, other), 30, -1))
    assert result["differences"] == []                    # the plasmid's base is shown
    assert result["doubtful"] >= 1
    assert not verified_at(result, 350)


def test_the_reference_base_is_not_verified_against_a_good_contrary_read():
    insert = REFERENCE[200:500]
    other = "A" if insert[150] != "A" else "C"
    result = made_contig(made_read("ref1", insert, 25), made_read("ref2", insert[:280], 25),
                         made_read("alt", with_base(insert, 150, other), 40, -1))
    assert result["differences"] == []                    # 25 + 25 outweighs 40
    assert not verified_at(result, 350)
    features = result["document"]["features"]
    assert any(f["name"] == "Sanger: check the trace" and f["start"] <= 350 <= f["end"]
               for f in features)


def test_a_tied_insertion_slot_is_not_verified():
    insert = REFERENCE[200:500]
    extra = insert[:150] + "G" + insert[150:]
    result = made_contig(made_read("plain", insert, 40), made_read("extra", extra, 40, -1))
    assert not verified_at(result, 350)
    assert all(not d["confirmed"] for d in result["differences"])


def test_verified_counts_the_reads_that_agree_and_skips_weak_bases():
    insert = REFERENCE[200:500]
    result = made_contig(made_read("a", with_base(insert, 150, "N"), [40] * 150 + [4] + [40] * 149),
                         made_read("b", insert, [40] * 160 + [12] + [40] * 139, -1))
    names = {(r["start"], r["end"]): r["reads"] for r in result["verified"]}
    assert any(start <= 351 <= end and n == 1 for (start, end), n in names.items())
    assert verified_at(result, 360)                      # Q12 in b, but a has Q40 there
    weak = made_contig(made_read("only", insert, [40] * 160 + [12] + [40] * 139))
    assert not verified_at(weak, 360)


# ---- third review ----------------------------------------------------------------------------

def test_extra_bases_called_n_are_never_verified():
    insert = REFERENCE[200:500]
    extra = insert[:150] + "N" + insert[150:]
    both = made_contig(made_read("a", extra, 40), made_read("b", extra, 40, -1))
    assert not verified_at(both, 350)
    assert both["doubtful"] >= 1 and both["disagreements"]
    against = made_contig(made_read("n", extra, 40), made_read("plain", insert, 40, -1))
    assert not verified_at(against, 350)


def test_a_good_mixed_call_on_a_reference_base_is_doubt():
    insert = REFERENCE[200:500]
    result = made_contig(made_read("plain", insert, 40),
                         made_read("mixed", with_base(insert, 150, "R"), 40, -1))
    assert result["differences"] == []
    assert not verified_at(result, 350)


def test_verified_counts_only_good_calls():
    insert = REFERENCE[200:500]
    result = made_contig(made_read("good", insert, 40),
                         made_read("poor", insert, [40] * 140 + [8] * 20 + [40] * 140, -1))
    assert any(r["start"] <= 351 <= r["end"] and r["reads"] == 1 for r in result["verified"])


def test_a_flag_on_deleted_bases_stays_in_the_document():
    insert = REFERENCE[200:500]
    deleted = insert[:150] + insert[152:]
    result = made_contig(made_read("a", deleted, 40), made_read("b", deleted, 40, -1),
                         made_read("c", insert, [40] * 150 + [25, 25] + [40] * 148))
    [difference] = result["differences"]
    assert difference["kind"] == "deletion" and difference["confirmed"] is False
    checks = [f for f in result["document"]["features"] if f["name"] == "Sanger: check the trace"]
    assert checks and checks[0]["end"] - checks[0]["start"] == 1


# ---- fourth review: alignments set by hand -------------------------------------------------

def hand_read(name: str, reference_gapped: str, read_gapped: str, offset: int,
              quality: list[int] | int = 40, strand: int = 1):
    """A read whose alignment to the reference is given column by column."""
    from vettore.sanger import AlignedRead

    bases = read_gapped.replace("-", "")
    qualities = quality if isinstance(quality, list) else [quality] * len(bases)
    trace = Trace(name, bases, qualities, [{}] * len(bases), strand)
    span = len(reference_gapped.replace("-", ""))
    return AlignedRead(trace, reference_gapped, read_gapped, offset, offset + span - 1, 99.0, [],
                       (0, len(bases)), len(bases), False, offset)


def test_a_weak_extra_base_never_hides_a_good_contrary_call_beside_it():
    from vettore.contig import contig_payload

    reference = REFERENCE[:300]
    other = {"A": "C", "C": "G", "G": "T", "T": "A"}[reference[101]]
    # a: Q40, a mismatch at 101. b: Q40 at 101 with the reference base, but 100 deleted and a
    # Q15 extra base before 101. c: a Q5 deletion of 101, which makes 100–101 one unit.
    a = hand_read("a", reference[50:151], reference[50:101] + other + reference[102:151], 50)
    b_read = reference[50:100] + "-" + "T" + reference[101:151]
    b_quality = [40] * len(b_read.replace("-", ""))
    b_quality[50] = 15
    b = hand_read("b", reference[50:100] + reference[100] + "-" + reference[101:151], b_read, 50,
                  b_quality, -1)
    c_read = reference[60:101] + "-" + reference[102:141]
    c_quality = [40] * len(c_read.replace("-", ""))
    c_quality[40] = c_quality[41] = 5
    c = hand_read("c", reference[60:141], c_read, 60, c_quality)
    result = contig_payload(reference, False, [a, b, c], doc("p", reference, False))
    mismatches = [d for d in result["differences"] if d["kind"] == "mismatch"]
    assert mismatches and all(d["confirmed"] is False for d in mismatches)
    assert not verified_at(result, 101)


def test_a_deleted_stretch_across_the_origin_maps_to_nothing():
    from vettore.contig import _new_range

    kept = [True] * 10
    kept[8] = kept[9] = kept[0] = False
    start_of, cursor = [], 0
    for k in kept:
        start_of.append(cursor)
        cursor += k
    assert _new_range(8, 0, start_of, kept, cursor, True) is None
    assert _new_range(7, 1, start_of, kept, cursor, True) == (6, 0)



# ---- fifth review ------------------------------------------------------------------------------

def test_a_join_no_read_runs_across_is_not_verified():
    """A forward read ending exactly where a reverse read begins: an insertion there would be
    clipped off both read ends by the local alignment, so the join is not verified."""
    from vettore.contig import contig_payload

    reference = REFERENCE[:400]
    forward = hand_read("F", reference[20:200], reference[20:200], 20)
    reverse = hand_read("R", reference[200:372], reference[200:372], 200, strand=-1)
    result = contig_payload(reference, False, [forward, reverse], doc("p", reference, False))
    assert result["differences"] == []
    assert [d["kind"] for d in result["disagreements"]] == ["not read across"]
    assert result["disagreements"][0]["position"] == 201
    assert result["doubtful"] == 1
    assert not any(r["start"] <= 200 and r["end"] >= 201 for r in result["verified"])


def test_a_read_round_a_whole_small_ring_does_not_read_its_own_join():
    from vettore.contig import contig_payload

    ring = clean_dna(150, 306)
    read = ring[40:] + ring[:40]
    whole = hand_read("whole", read, read, 40)
    result = contig_payload(ring, True, [whole], doc("pSmall", ring, True))
    assert [d["kind"] for d in result["disagreements"]] == ["not read across"]
    assert result["disagreements"][0]["position"] == 41
