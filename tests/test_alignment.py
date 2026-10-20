"""Star alignment of several sequences to a reference (``vettore.alignment``).

Hand-worked cases plant known substitutions, deletions and insertions into site-free random
DNA (``cloning_helpers.clean_dna``) and check that they come back at the right reference
coordinates, whatever the strand or origin of the other sequence. Indels are planted where
their placement is unambiguous (the bases on either side differ from the indel's ends), so
that only one optimal alignment exists.

Real-data fixtures, downloaded from NCBI nuccore with E-utilities (efetch, rettype=gbwithparts)
on 2026-10-20 and stored unchanged:
- ``tests/data/analysis-NC_001401.2.gb``: adeno-associated virus 2, complete genome
  (VP1 CDS annotated at 2203..4410);
- ``tests/data/analysis-AY530579.1.gb``: adeno-associated virus 9 (isolate hu.14) capsid
  protein VP1 gene, complete cds (2211 bp).
"""

from pathlib import Path

import pytest
from Bio import SeqIO
from cloning_helpers import clean_dna, doc, rc
from fastapi.testclient import TestClient

from vettore import alignment as al
from vettore.server import TOKEN_HEADER, create_app

DATA = Path(__file__).parent / "data"
TOKEN = "t0ken"
REF = clean_dna(1500, 301)


def unambiguous_deletion(seq: str, near: int, size: int) -> int:
    """A start position near ``near`` where deleting ``size`` bases has one placement."""
    for p in range(near, near + 200):
        if seq[p - 1] != seq[p + size - 1] and seq[p] != seq[p + size]:
            return p
    raise AssertionError("no unambiguous deletion site")


def unambiguous_insertion(seq: str, near: int, inserted: str) -> int:
    """A position q near ``near`` where inserting after base q has one placement."""
    for q in range(near, near + 200):
        if inserted[-1] != seq[q] and inserted[0] != seq[q + 1]:
            return q
    raise AssertionError("no unambiguous insertion site")


SUB = 150
SUB_BASE = "A" if REF[SUB] != "A" else "C"
DEL = unambiguous_deletion(REF, 600, 4)
INSERTED = "GATTACA"
INS = unambiguous_insertion(REF, 1000, INSERTED)
# The variant, built from the reference by hand: one substitution, a 4-base deletion and a
# 7-base insertion, far enough apart not to interact.
VARIANT = (REF[:SUB] + SUB_BASE + REF[SUB + 1:DEL] + REF[DEL + 4:INS + 1] + INSERTED
           + REF[INS + 1:])
EXPECTED = [
    {"kind": "mismatch", "start": SUB + 1, "end": SUB + 1, "reference": REF[SUB],
     "sequence": SUB_BASE},
    {"kind": "deletion", "start": DEL + 1, "end": DEL + 4, "reference": REF[DEL:DEL + 4],
     "sequence": ""},
    {"kind": "insertion", "start": INS + 1, "end": INS + 2, "reference": "",
     "sequence": INSERTED},
]


def payload(reference: dict, *others: dict) -> dict:
    return al.alignment_payload(reference, list(others))


def test_planted_differences_come_back_at_their_reference_positions():
    out = payload(doc("ref", REF, False), doc("variant", VARIANT, False))
    [s] = out["sequences"]
    assert s["differences"] == EXPECTED
    assert s["strand"] == 1 and s["mode"] == "semi-global"
    assert (s["start"], s["end"], s["coverage"]) == (1, 1500, 100.0)
    # BLAST-style identity: identical columns over all columns, gaps included.
    assert s["identity"] == round(100 * (1500 - 1 - 4) / (1500 + len(INSERTED)), 2)
    assert s["counts"] == {"mismatches": 1, "insertions": 1, "deletions": 1,
                           "inserted_bases": 7, "deleted_bases": 4}


def test_the_reverse_strand_is_found_and_differences_stay_on_the_reference_strand():
    out = payload(doc("ref", REF, False), doc("variant", rc(VARIANT), False))
    [s] = out["sequences"]
    assert s["strand"] == -1
    assert s["differences"] == EXPECTED
    assert out["alignmentTracks"][1]["sequenceData"]["name"] == "variant (reverse complement)"


def test_a_circular_sequence_is_rotated_to_the_reference_origin():
    rotated = VARIANT[900:] + VARIANT[:900]
    out = payload(doc("ref", REF, True), doc("variant", rotated, True))
    [s] = out["sequences"]
    assert s["differences"] == EXPECTED
    assert s["mode"] == "global"
    # VARIANT[0] (which faces reference base 1) sits at index len - 900 of the rotated copy.
    assert s["origin_moved_to"] == len(VARIANT) - 900 + 1
    assert (s["start"], s["end"], s["spans_origin"]) == (1, 1500, False)


def test_a_rotated_reverse_complemented_circle_gives_the_same_answer():
    rotated = rc(VARIANT[400:] + VARIANT[:400])
    [s] = payload(doc("ref", REF, True), doc("variant", rotated, True))["sequences"]
    assert s["strand"] == -1
    assert s["differences"] == EXPECTED
    # Reference base 1 faces the base that was VARIANT[0] (complemented on this strand).
    assert rotated[s["origin_moved_to"] - 1] == rc(VARIANT[0])


def test_a_deletion_across_the_origin_of_two_circles_is_one_event():
    length = len(REF)
    # Delete the last `tail` and the first `head` bases of the circular reference, choosing
    # the split so that the deletion has a single placement.
    tail, head = next((t, 5 - t) for t in range(1, 5)
                      if REF[length - t - 1] != REF[5 - t - 1] and REF[length - t] != REF[5 - t])
    out = payload(doc("ref", REF, True), doc("short", REF[head:length - tail], True))
    [s] = out["sequences"]
    assert s["differences"] == [{"kind": "deletion", "start": length - tail + 1, "end": head,
                                 "reference": REF[-tail:] + REF[:head], "sequence": ""}]


def test_an_insertion_at_the_origin_of_two_circles_is_reported_between_last_and_first_base():
    inserted = next(s for s in ("TTGCA", "ACGGT", "GACTT", "CAGTG")
                    if s[-1] != REF[-1] and s[0] != REF[0])
    [s] = payload(doc("ref", REF, True), doc("longer", REF + inserted, True))["sequences"]
    assert s["differences"] == [{"kind": "insertion", "start": len(REF), "end": 1,
                                 "reference": "", "sequence": inserted}]


def test_a_fragment_across_the_origin_of_a_circular_reference():
    length = len(REF)
    fragment = (REF[length - 300:] + REF)[:800]
    [s] = payload(doc("ref", REF, True), doc("frag", fragment, False))["sequences"]
    assert (s["start"], s["end"]) == (length - 300 + 1, 500)
    assert s["spans_origin"] is True
    assert s["coverage"] == round(100 * 800 / length, 1)
    assert s["differences"] == [] and s["identity"] == 100.0


def test_a_fragment_inside_a_longer_reference_lands_where_it_belongs():
    [s] = payload(doc("ref", REF, False), doc("frag", REF[700:1100], False))["sequences"]
    assert (s["start"], s["end"], s["mode"]) == (701, 1100, "semi-global")


def test_a_short_linear_reference_inside_a_circular_plasmid():
    # Reference: an insert; the other sequence: the plasmid carrying it, origin elsewhere.
    plasmid = REF[1100:] + REF[:1100]
    [s] = payload(doc("insert", REF[300:900], False), doc("plasmid", plasmid, True))["sequences"]
    assert (s["start"], s["end"], s["coverage"]) == (1, 600, 100.0)
    assert s["identity"] == 100.0 and s["aligned_length"] == 600


def test_a_circular_reference_contained_in_a_longer_linear_sequence():
    flank_left, flank_right = clean_dna(200, 77), clean_dna(150, 78)
    linear = flank_left + REF[500:] + REF[:500] + flank_right
    [s] = payload(doc("ref", REF, True), doc("linear", linear, False))["sequences"]
    assert (s["start"], s["end"], s["coverage"]) == (1, 1500, 100.0)
    assert s["differences"] == [] and s["aligned_length"] == 1500


def test_a_partly_related_sequence_is_aligned_locally():
    other = REF[200:500] + clean_dna(1200, 999)
    [s] = payload(doc("ref", REF, False), doc("other", other, False))["sequences"]
    assert s["mode"] == "local"
    assert (s["start"], s["end"], s["identity"]) == (201, 500, 100.0)
    assert s["aligned_length"] == 300


def test_an_unrelated_sequence_is_left_out_with_a_reason():
    out = payload(doc("ref", REF, False), doc("variant", VARIANT, False),
                  doc("stranger", clean_dna(800, 4242), False))
    assert [s["name"] for s in out["sequences"]] == ["variant"]
    assert "stranger" in out["warnings"][0] and "left out" in out["warnings"][0]


def test_a_chance_match_is_left_out_without_running_the_alignment(monkeypatch):
    """Two unrelated plasmids that share one 16 bp stretch: a row with '100 % identity' over
    16 bp would read as 'similar'. They are left out, and the alignment itself (seconds on
    plasmid-sized sequences) is never run."""
    shared = REF[700:716]
    stranger = clean_dna(1000, 4243) + shared + clean_dna(1000, 4244)
    aligner = al._aligner
    sizes = []

    def small_alignments_only(mode):
        real = aligner(mode)

        class Guarded:
            def align(self, target, query):
                sizes.append(len(target) * len(query))
                assert len(target) < 500 and len(query) < 500, "whole alignment was run"
                return real.align(target, query)

        return Guarded()

    monkeypatch.setattr(al, "_aligner", small_alignments_only)
    out = payload(doc("ref", REF, False), doc("stranger", stranger, True))
    assert out["sequences"] == [] and out["alignmentTracks"] == []
    [warning] = out["warnings"]
    assert warning.startswith("“stranger” was left out: it has nothing recognisably in common "
                              "with “ref”")
    assert "their best match pairs only 16 bp" in warning
    assert sizes      # the shared stretch was looked at, in a small window


def test_a_short_oligo_and_a_shared_cassette_are_still_aligned():
    """Short is not the same as meaningless: a 20 bp oligo found whole in the reference, and a
    120 bp cassette shared by two otherwise unrelated plasmids, are real matches."""
    oligo = REF[300:320]
    cassette = clean_dna(900, 4245) + REF[1000:1120] + clean_dna(900, 4246)
    out = payload(doc("ref", REF, False), doc("oligo", oligo, False),
                  doc("cassette", cassette, False))
    rows = {s["name"]: s for s in out["sequences"]}
    assert (rows["oligo"]["start"], rows["oligo"]["end"]) == (301, 320)
    assert rows["cassette"]["mode"] == "local"
    assert (rows["cassette"]["start"], rows["cassette"]["end"]) == (1001, 1120)
    assert out["warnings"] == []


def test_tracks_share_one_coordinate_system():
    out = payload(doc("ref", REF, False), doc("variant", VARIANT, False),
                  doc("frag", REF[700:1100], False), doc("rc", rc(VARIANT), True))
    tracks = out["alignmentTracks"]
    assert len(tracks) == 4
    assert len({len(t["alignmentData"]["sequence"]) for t in tracks}) == 1
    assert tracks[0]["alignmentData"]["sequence"].replace("-", "") == REF
    for track in tracks[1:]:
        assert (track["alignmentData"]["sequence"].replace("-", "")
                == track["sequenceData"]["sequence"])
    # The insertion of the variant opened columns in every track.
    assert "-" * len(INSERTED) in tracks[0]["alignmentData"]["sequence"]
    assert [s["index"] for s in out["sequences"]] == [0, 1, 2]


def test_reference_features_are_shown_except_those_across_the_origin():
    features = [{"name": "inside", "start": 10, "end": 99, "strand": 1, "forward": True},
                {"name": "wraps", "start": 1400, "end": 20, "strand": 1, "forward": True}]
    reference = doc("ref", REF, True, features=features)
    out = payload(reference, doc("variant", VARIANT, True))
    names = [f["name"] for f in out["alignmentTracks"][0]["sequenceData"]["features"]]
    assert names == ["inside"]


def record(name: str) -> str:
    return str(SeqIO.read(DATA / name, "genbank").seq)


def test_aav2_vp1_lands_on_its_annotated_cds_in_the_genome():
    genome = record("analysis-NC_001401.2.gb")
    vp1 = genome[2202:4410]                       # CDS 2203..4410 in the record
    for sequence, strand in ((vp1, 1), (rc(vp1), -1)):
        [s] = payload(doc("AAV2", genome, False), doc("VP1", sequence, False))["sequences"]
        assert (s["start"], s["end"], s["strand"]) == (2203, 4410, strand)
        assert s["identity"] == 100.0 and s["differences"] == []


def test_aav9_capsid_aligns_end_to_end_with_aav2_capsid():
    genome = record("analysis-NC_001401.2.gb")
    aav9 = record("analysis-AY530579.1.gb")
    out = payload(doc("AAV2", genome, False), doc("AAV9 VP1", aav9, False),
                  doc("AAV9 VP1 reversed", rc(aav9), False))
    forward, reverse = out["sequences"]
    # The AAV9 VP1 ORF (2211 bp) is homologous to AAV2 VP1 (2208 bp) along its whole length.
    assert (forward["start"], forward["end"]) == (2203, 4410)
    assert forward["mode"] == "semi-global" and forward["aligned_length"] == 2211
    assert 70 < forward["identity"] < 90          # divergent serotypes, not the same gene
    assert reverse["strand"] == -1 and reverse["differences"] == forward["differences"]
    # Every difference falls inside the capsid region.
    assert all(2203 <= d["start"] <= 4410 for d in forward["differences"])


def test_unrelated_plasmids_are_left_out_against_a_viral_genome():
    """The reported case: the AAV2 genome as reference and plasmids that share nothing with it
    gave rows like '93.75 % identity' over 16 bp. They are left out with the reason."""
    genome = record("analysis-NC_001401.2.gb")
    plasmids = list(SeqIO.parse(DATA / "exchange-plasmids.fasta", "fasta"))
    pbr322 = SeqIO.read(DATA / "exchange-pBR322.embl", "embl")
    others = [doc(p.id, str(p.seq).upper(), True) for p in plasmids[:2]]
    others.append(doc("pBR322", str(pbr322.seq).upper(), True))
    out = payload(doc("AAV2", genome, False), *others)
    assert out["sequences"] == []
    assert len(out["warnings"]) == len(others)
    assert all("nothing recognisably in common with “AAV2”" in w for w in out["warnings"])


# ---- limits and HTTP ----------------------------------------------------------------------

def client() -> TestClient:
    return TestClient(create_app(TOKEN, None), base_url="http://127.0.0.1")


def post(body: dict):
    return client().post("/api/alignment/align", headers={TOKEN_HEADER: TOKEN}, json=body)


def test_limits_route():
    res = client().get("/api/alignment/limits", headers={TOKEN_HEADER: TOKEN})
    assert res.json() == {"max_sequences": al.MAX_SEQUENCES, "max_length": al.MAX_LENGTH,
                          "min_length": al.MIN_LENGTH}


def test_http_route_aligns():
    res = post({"reference": doc("ref", REF, False), "sequences": [doc("v", VARIANT, False)]})
    assert res.status_code == 200
    assert res.json()["sequences"][0]["differences"] == EXPECTED


@pytest.mark.parametrize("sequences, message", [
    ([], "at least one sequence"),
    ([doc(f"s{i}", REF, False) for i in range(al.MAX_SEQUENCES)], "at most"),
    ([doc("huge", "ACGT" * (al.MAX_LENGTH // 4 + 1), False)], "longer than"),
    ([doc("odd", REF[:100] + "XZ" + REF[100:], False)], "not DNA"),
])
def test_requests_beyond_the_limits_are_refused_with_a_reason(sequences, message):
    res = post({"reference": doc("ref", REF, False), "sequences": sequences})
    assert res.status_code == 422
    assert message in res.json()["detail"]


def test_a_too_short_reference_is_refused():
    res = post({"reference": doc("ref", "ACGT", False), "sequences": [doc("v", REF, False)]})
    assert res.status_code == 422 and "too short" in res.json()["detail"]


def test_rna_letters_are_read_as_dna():
    [s] = payload(doc("ref", REF, False),
                  doc("rna", REF[:500].replace("T", "U").lower(), False))["sequences"]
    assert s["identity"] == 100.0
