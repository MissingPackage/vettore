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
    assert result["covered"] == [{"start": LEFT + 1, "end": one_end, "reads": 2}]


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
              and forward_calls()[p].value != "N"]
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


def test_the_covered_stretch_is_annotated_as_verified():
    result = contig(clone())
    [verified] = [f for f in result["document"]["features"]
                  if f["name"] == "Sanger verified (2 reads)"]
    span = result["covered"][0]
    assert (verified["start"], verified["end"]) == (span["start"] - 1, span["end"] - 1)
    assert result["document"]["name"] == "pCOI consensus"


def test_unconfirmed_differences_are_not_annotated():
    result = contig(clone())
    [(start, end)] = strand_runs(result, 1)
    weak = next(p for p in range(start - 1, end - 10) if p in forward_calls()
                and forward_calls()[p].value != "N"
                and forward_calls()[p].quality < CONFIRMING_QUALITY)
    reference, _ = mutate(clone(), weak)
    features = contig(reference)["document"]["features"]
    assert not [f for f in features if f["type"] == "variation"]


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


def test_two_low_quality_reads_together_confirm_a_difference():
    insert = REFERENCE[200:500]
    other = "A" if insert[150] != "A" else "C"
    mutated = with_base(insert, 150, other)
    quality = [45] * 150 + [15] + [45] * 149
    [difference] = made_contig(made_read("a", mutated, quality),
                               made_read("b", mutated, quality, -1))["differences"]
    assert difference["confirmed"] is True and difference["reads"] == 2


def test_an_n_is_no_call_and_never_enters_the_consensus():
    insert = REFERENCE[200:500]
    result = made_contig(made_read("a", with_base(insert, 150, "N"), 40),
                         made_read("b", insert, 20, -1))
    assert result["differences"] == []
    assert result["disagreements"] == []
    alone = made_contig(made_read("a", with_base(insert, 150, "N"), 40))
    assert alone["differences"] == []
    assert [d["kind"] for d in alone["disagreements"]] == ["no call"]
    assert alone["document"]["sequence"] == REFERENCE


def test_coverage_counts_strands_not_reads():
    insert = REFERENCE[200:500]
    result = made_contig(made_read("a", insert[:200], 40), made_read("b", insert[:200], 40),
                         made_read("c", rc(insert[100:]), 40))
    assert [(r["start"], r["end"], r["strands"]) for r in result["coverage"]] == [
        (1, 200, 0), (201, 300, 1), (301, 400, 2), (401, 500, 1), (501, 700, 0)]
