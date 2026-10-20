from pathlib import Path

import pytest
from Bio import SeqIO
from cloning_helpers import clean_dna, doc, rc
from fastapi.testclient import TestClient

from vettore.sanger import (
    DEFAULT_CUTOFF,
    TraceError,
    align_trace,
    alignment_payload,
    mott_trim,
    read_trace,
)
from vettore.server import TOKEN_HEADER, create_app

DATA = Path(__file__).parent / "data"
AB1 = DATA / "3730.ab1"
TOKEN = "t0ken"


def trace():
    return read_trace(AB1.name, AB1.read_bytes())


def test_reading_a_trace_exposes_bases_quality_and_one_peak_bundle_per_base():
    t = trace()
    assert len(t.sequence) == 1165
    assert len(t.quality) == len(t.sequence) == len(t.base_traces)
    assert set(t.base_traces[0]) == {"aTrace", "tTrace", "gTrace", "cTrace"}
    assert all(len(b["aTrace"]) == len(b["cTrace"]) > 0 for b in t.base_traces[:50])
    assert t.annotations["machine_model"]


def test_trimming_agrees_with_biopythons_abi_trim():
    t = trace()
    start, end = mott_trim(t.quality)
    biopython = SeqIO.read(AB1, "abi-trim")
    # Biopython anchors the start at the first positive cumulative score rather than at the
    # start of the best segment, so the two can differ by a base at each edge.
    assert abs((end - start) - len(biopython)) <= 2
    assert t.sequence[start:end].startswith(str(biopython.seq)[:40])


def test_trimming_removes_the_bad_quality_head():
    t = trace()
    start, _ = mott_trim(t.quality)
    assert start > 0
    assert min(t.quality[:start]) < 20


def test_reverse_complement_also_flips_the_peaks():
    t = trace().slice(0, 30)
    r = t.reverse_complement()
    assert r.sequence == rc(t.sequence)
    assert r.strand == -1
    assert r.base_traces[0]["aTrace"] == list(reversed(t.base_traces[-1]["tTrace"]))
    assert r.quality == list(reversed(t.quality))


def reference_around(read_core: str, left: int = 300, right: int = 400) -> str:
    return clean_dna(left, 201) + read_core + clean_dna(right, 202)


def test_a_perfect_read_aligns_with_no_discrepancies():
    t = trace()
    start, _ = mott_trim(t.quality)
    core = t.sequence[start + 40 : start + 500]
    aligned = align_trace(reference_around(core), False, t)
    assert aligned.identity == 100.0
    assert aligned.discrepancies == []
    assert aligned.start == 300
    assert aligned.trace.strand == 1


def test_a_planted_mismatch_is_reported_at_the_right_position():
    t = trace()
    start, _ = mott_trim(t.quality)
    core = t.sequence[start + 40 : start + 500]
    reference = reference_around(core)
    position = 400                      # 0-based, inside the read
    replacement = "A" if reference[position] != "A" else "G"
    mutated = reference[:position] + replacement + reference[position + 1 :]
    aligned = align_trace(mutated, False, t)
    mismatches = [d for d in aligned.discrepancies if d["kind"] == "mismatch"]
    assert [(m["position"], m["reference"]) for m in mismatches] == [(position + 1, replacement)]
    assert mismatches[0]["read"] == reference[position]


def test_a_difference_carries_the_read_quality_of_its_base():
    t = trace()
    start, _ = mott_trim(t.quality)
    core = t.sequence[start + 40 : start + 500]
    reference = reference_around(core)
    position = 400                      # reference index; read index start + 40 + 100
    replacement = "A" if reference[position] != "A" else "G"
    mutated = reference[:position] + replacement + reference[position + 1 :]
    [mismatch] = [d for d in align_trace(mutated, False, t).discrepancies
                  if d["kind"] == "mismatch"]
    assert mismatch["quality"] == t.quality[start + 140]


def test_column_qualities_give_a_gap_the_weaker_neighbour():
    from vettore.sanger import _column_qualities

    assert _column_qualities("AC-GT", [30, 12, 25, 40]) == [30, 12, 12, 25, 40]
    assert _column_qualities("--AC", [20, 35]) == [20, 20, 20, 35]
    assert _column_qualities("AC--", [20, 35]) == [20, 35, 35, 35]


def test_a_merged_indel_carries_its_lowest_quality():
    from vettore.sanger import _discrepancies

    diffs = _discrepancies("AC---GT", "ACTTAGT", 0, 100, [30, 30, 25, 9, 31, 30, 30])
    assert diffs == [{"kind": "insertion", "position": 3, "reference": "", "read": "TTA",
                      "quality": 9}]


def test_a_planted_deletion_is_reported_as_one_event():
    t = trace()
    start, _ = mott_trim(t.quality)
    core = t.sequence[start + 40 : start + 500]
    reference = reference_around(core)
    mutated = reference[:400] + reference[403:]        # the read now has 3 extra bases
    aligned = align_trace(mutated, False, t)
    events = [d for d in aligned.discrepancies if d["kind"] == "insertion"]
    assert len(events) == 1 and len(events[0]["read"]) == 3


def test_a_read_on_the_reverse_strand_is_detected():
    t = trace()
    start, _ = mott_trim(t.quality)
    core = t.sequence[start + 40 : start + 500]
    reference = clean_dna(200, 203) + rc(core) + clean_dna(200, 204)
    aligned = align_trace(reference, True, t)
    assert aligned.trace.strand == -1
    assert aligned.identity > 95


def test_payload_satisfies_the_editor_invariants():
    reference = reference_around(trace().sequence[60:520])
    payload = alignment_payload("pTest", reference, False, [(AB1.name, AB1.read_bytes())])
    ref_track, read_track = payload["alignmentTracks"]
    assert len(ref_track["alignmentData"]["sequence"]) == len(read_track["alignmentData"]["sequence"])
    assert ref_track["alignmentData"]["sequence"].replace("-", "") == reference
    assert (read_track["alignmentData"]["sequence"].replace("-", "")
            == read_track["sequenceData"]["sequence"])
    assert read_track["chromatogramData"] is read_track["sequenceData"]["chromatogramData"]
    chromatogram = read_track["sequenceData"]["chromatogramData"]
    assert len(chromatogram["baseCalls"]) == len(read_track["sequenceData"]["sequence"])
    assert len(chromatogram["basePos"]) == len(chromatogram["baseCalls"])
    assert 0 < payload["coverage_percent"] <= 100
    assert payload["reads"][0]["identity"] == 100.0
    aligned_quality = read_track["sequenceData"]["chromatogramData"]["qualNums"]
    assert payload["reads"][0]["mean_quality"] == round(sum(aligned_quality) / len(aligned_quality), 1)


def test_several_reads_share_one_coordinate_system():
    reference = reference_around(trace().sequence[60:520])
    payload = alignment_payload("pTest", reference, False,
                                [("a.ab1", AB1.read_bytes()), ("b.ab1", AB1.read_bytes())])
    tracks = payload["alignmentTracks"]
    assert len(tracks) == 3
    assert len({len(t["alignmentData"]["sequence"]) for t in tracks}) == 1
    assert payload["first_covered"] == payload["reads"][0]["start"]


def test_an_insertion_in_one_read_opens_a_column_in_every_track():
    t = trace()
    start, _ = mott_trim(t.quality)
    core = t.sequence[start + 40 : start + 500]
    reference = reference_around(core)
    mutated = reference[:400] + reference[403:]       # the read carries 3 extra bases
    payload = alignment_payload("pTest", mutated, False, [(AB1.name, AB1.read_bytes())])
    reference_track, read_track = payload["alignmentTracks"]
    gapped = reference_track["alignmentData"]["sequence"]
    assert gapped.replace("-", "") == mutated
    assert "---" in gapped                            # the insertion opened columns
    assert len(gapped) == len(read_track["alignmentData"]["sequence"])


def test_an_unreadable_file_is_reported_as_a_warning_not_a_crash():
    payload = alignment_payload("pTest", clean_dna(500, 205), False, [("broken.ab1", b"junk")])
    assert payload["alignmentTracks"] == []
    assert "broken.ab1" in payload["warnings"][0]


def test_a_read_from_another_plasmid_is_refused_not_given_a_verdict():
    """3730.ab1 against an unrelated 10 kb plasmid: the best local alignment is a chance match
    of a few dozen bases, which must not be reported as a read with differences."""
    from vettore.fileio import read_documents

    path = DATA / "cloning-gateway-destination.gb"
    [plasmid] = read_documents(path.name, path.read_bytes())
    payload = alignment_payload("pAGRIKOLA", plasmid["sequence"].upper(), True,
                                [(AB1.name, AB1.read_bytes())])
    assert payload["reads"] == [] and payload["alignmentTracks"] == []
    [warning] = payload["warnings"]
    assert warning.startswith("Only ")
    assert "good-quality bases of" in warning and "match “pAGRIKOLA”" in warning
    assert "probably comes from another plasmid" in warning
    assert "active tab" in warning


def test_a_read_running_past_the_end_of_a_short_linear_reference_is_kept():
    """A linear reference shorter than the read (an insert on its own): the read bases past its
    ends are not evidence against it."""
    t = trace()
    start, _ = mott_trim(t.quality)
    core = t.sequence[start + 300: start + 380]          # 80 bases, inside the read
    aligned = align_trace(core, False, t)
    assert aligned.identity == 100.0
    assert (aligned.start, aligned.end) == (0, 79)


def test_reading_a_non_trace_file_raises():
    with pytest.raises(TraceError):
        read_trace("x.ab1", b"not a trace")


def test_http_route():
    import base64

    reference = reference_around(trace().sequence[60:520])
    client = TestClient(create_app(TOKEN, None), base_url="http://127.0.0.1")
    body = {"reference": doc("pTest", reference, False),
            "reads": [{"filename": AB1.name,
                       "content_b64": base64.b64encode(AB1.read_bytes()).decode()}],
            "cutoff": DEFAULT_CUTOFF}
    res = client.post("/api/sanger/align", headers={TOKEN_HEADER: TOKEN}, json=body)
    assert res.status_code == 200
    assert res.json()["reads"][0]["strand"] == 1
    short = client.post("/api/sanger/align", headers={TOKEN_HEADER: TOKEN},
                        json={"reference": doc("x", "ACGT", False), "reads": []})
    assert short.status_code == 422


def test_e2e_reference_fixture_is_up_to_date():
    import subprocess
    import sys

    root = Path(__file__).resolve().parents[1]
    fixture = root / "frontend" / "e2e" / "fixtures" / "sanger-reference.txt"
    out = subprocess.run([sys.executable, "-I", str(root / "scripts" / "make_sanger_fixture.py")],
                         capture_output=True, text=True, encoding="utf-8", check=True).stdout
    assert out == fixture.read_text(encoding="utf-8"), (
        "Regenerate: python -I scripts/make_sanger_fixture.py > "
        "frontend/e2e/fixtures/sanger-reference.txt")
