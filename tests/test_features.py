import collections
import functools
import json
from pathlib import Path
from typing import ClassVar
from unittest import mock

import pytest
from Bio import Align, SeqIO
from Bio.Seq import Seq
from cloning_helpers import clean_dna, rc

from vettore import features as features_module
from vettore.features import (
    DEFAULT_MIN_COVERAGE,
    DEFAULT_MIN_IDENTITY,
    SHORT_FEATURE,
    catalogue,
    detect,
    detect_payload,
    shown_identity,
)

DATA = Path(__file__).parent / "data"


def plasmid(accession: str):
    record = SeqIO.read(DATA / f"{accession}.gb", "genbank")
    return str(record.seq).upper(), record.annotations.get("topology") == "circular"


def names(matches) -> list[str]:
    return [m.name for m in matches]


def cited(source: str) -> tuple[str, int, int, int]:
    """(accession, start, end, strand) of a ``source`` such as ``X64410.1:429..598`` or, for
    an element a record annotates on its minus strand, ``LT727488.1:complement(3997..4178)``.
    Coordinates are 1-based and inclusive, as the record's FEATURES table prints them."""
    accession, _, span = source.partition(":")
    strand = 1
    if span.startswith("complement(") and span.endswith(")"):
        span, strand = span[len("complement("):-1], -1
    start, _, end = span.partition("..")
    assert accession and start.isdigit() and end.isdigit(), source
    return accession, int(start), int(end), strand


def test_every_catalogue_entry_cites_a_record_and_a_range():
    for feature in catalogue():
        accession, start, end, _ = cited(feature.source)
        assert accession, feature.name
        assert end - start + 1 == feature.length, feature.name
        assert set(feature.sequence) <= set("ACGTRYSWKMBDHVN"), feature.name


def test_loxp_matches_the_canonical_site():
    [loxp] = [f for f in catalogue() if f.name == "loxP"]
    assert loxp.sequence == "ATAACTTCGTATAGCATACATTATACGAAGTTAT"


def test_features_are_found_across_records_not_only_in_their_own_source():
    """AmpR, the origin and WPRE come from other accessions, so finding them here is a real
    test of the matcher rather than of the database."""
    sequence, circular = plasmid("MT612431.1")
    matches = detect(sequence, circular)
    found = names(matches)
    assert "AmpR" in found and "ori" in found and "WPRE" in found and "EGFP" in found
    for match in matches:
        if match.name in {"AmpR", "WPRE", "EGFP", "ori"}:
            assert not match.source.startswith("MT612431"), match.name
            assert match.identity >= 98.5


def test_a_dio_cassette_shows_both_lox_pairs_in_both_orientations():
    sequence, circular = plasmid("MT612431.1")
    matches = detect(sequence, circular)
    lox = [m for m in matches if m.name in {"loxP", "lox2272"}]
    assert len(lox) == 4
    assert {m.strand for m in lox} == {1, -1}
    assert sorted({m.name for m in lox}) == ["lox2272", "loxP"]


def test_a_gene_on_the_minus_strand_is_reported_as_reverse():
    sequence, circular = plasmid("MT612431.1")
    [egfp] = [m for m in detect(sequence, circular) if m.name == "EGFP"]
    assert egfp.strand == -1
    # The EGFP here is fused at both ends (Gephyrin.FingR-GFP-CCR5TC): linker codons stand
    # where the catalogue has ATG and AAG TAA. The annotation is the EGFP that is there, codons
    # 2 to 238 in its own frame (VSKGEE…MDELY in the record's translation), not padded with
    # the linkers to the catalogue's 720 bp; the record's own "GFP" (2147..2858) is 712.
    cds = rc(sequence[egfp.start : egfp.end + 1])
    assert cds.startswith("GTGAGCAAGGGCGAGGAGCTGTTCACC") and cds.endswith("GACGAGCTGTAC")
    assert egfp.end - egfp.start + 1 == 711
    assert egfp.part == (4, 714) and egfp.entry_length == 720
    assert egfp.identity == 100.0


# Both records carry an SV40 polyA cloned between sites that differ from the ones around the
# catalogue entry; one base inside it differs too. (record, start, end), 0-based inclusive.
SV40_POLYA_IN_RECORDS = [("HG530138.1", 6964, 7185), ("MT891341.1", 2008, 2229)]


@pytest.mark.parametrize(("accession", "start", "end"), SV40_POLYA_IN_RECORDS)
def test_an_element_is_annotated_without_the_junction_it_was_cloned_into(accession, start,
                                                                        end):
    """The bases around the element are the vector's, not the catalogue entry's (6 and 4 of
    them, nearly all different), so they are not annotated as SV40 polyA, and the identity
    shown is that of the 222 bp that are: 221 of them match. With the threshold at 98 % the
    element is still found, and at the identity the threshold was tested on."""
    sequence = record_sequence(accession)
    for threshold in (95.0, 98.0):
        [polya] = [m for m in detect(sequence, False, min_identity=threshold)
                   if m.name == "SV40 polyA"]
        assert (polya.start, polya.end) == (start, end), threshold
        assert polya.identity == pytest.approx(100.0 * 221 / 222), threshold
        assert polya.part == (7, 228) and polya.entry_length == 232, threshold


@pytest.mark.parametrize("threshold", [90.0, 95.0, 98.0, 99.0, 100.0])
def test_no_match_is_shown_below_the_identity_threshold(threshold):
    """The threshold filters the number the table shows, on every fixture record."""
    sequences = [(record_sequence(path.name.removeprefix("features-").removesuffix(".fasta")),
                  False) for path in sorted(DATA.glob("features-*.fasta"))]
    sequences += [plasmid(accession) for accession in ("MT612431.1", "PZ629150.1")]
    for sequence, circular in sequences:
        for match in detect(sequence, circular, min_identity=threshold):
            assert match.identity >= threshold, (match.name, match.identity)
            row = match.as_dict(len(sequence))
            assert row["identity"] >= threshold, (match.name, row["identity"])


def test_a_coding_sequence_runs_from_its_own_start_codon_in_frame():
    """NeoR/KanR here starts ATG GGG ATT (an NcoI site put a codon in after the start), so
    the alignment begins on the last G of that codon. The annotation is the ORF the record
    itself annotates, ATG to TGA, with the codon put in counted as an insertion: it
    translates as the record's MGIEQDG… and covers the whole 795 bp entry."""
    sequence, circular = plasmid("cloning-multisite-destination-r4r2")
    [neo] = [m for m in detect(sequence, circular) if m.name == "NeoR/KanR"]
    cds = rc(sequence[neo.start : neo.end + 1])
    assert len(cds) == 798
    protein = str(Seq(cds).translate())
    assert protein.startswith("MGIEQDGLHAGSPA") and protein.endswith("LDEFF*")
    assert neo.part == (1, 795) and not neo.frameshift
    # 798 columns: the ATG matches, the 3 bases put in count as differences, and 4 bases of
    # the rest differ from the catalogue as they did before.
    assert neo.identity == pytest.approx(100.0 * 791 / 798)


def test_a_match_tells_which_part_of_the_entry_it_covers():
    """mCherry in MT891341.1 has GTC where the catalogue starts ATG, and stops with TAG where
    the catalogue has TAA. The stop is the vector's own and is annotated; the first codon is
    not mCherry's start and is left out, and the row says the match is bases 4 to 711 of the
    711 bp entry rather than letting the identity pass for the whole of it."""
    sequence = record_sequence("MT891341.1")
    [mcherry] = [m for m in detect(sequence, False) if m.name == "mCherry"]
    assert sequence[mcherry.end - 2 : mcherry.end + 1] == "TAG"
    row = mcherry.as_dict(len(sequence))
    assert (row["length"], row["identity"]) == (708, shown_identity(100.0 * 707 / 708))
    assert (row["catalogueFrom"], row["catalogueTo"], row["catalogueLength"]) == (4, 711, 711)


@pytest.mark.parametrize("strand", [1, -1])
def test_a_different_stop_codon_of_the_vector_is_part_of_the_cds(strand):
    [egfp] = [f for f in catalogue() if f.name == "EGFP"]
    assert egfp.sequence.endswith("TAA")
    insert = egfp.sequence[:-3] + "TAG"
    filler = clean_dna(400, 311)
    construct = filler + (insert if strand == 1 else rc(insert)) + filler
    [match] = [m for m in detect(construct, False) if m.name == "EGFP"]
    assert (match.start, match.end, match.strand) == (400, 400 + 719, strand)
    assert match.part == (1, 720) and match.identity == pytest.approx(100.0 * 719 / 720)


@pytest.mark.parametrize(("linker", "start"), [("ATGGGC", 1), ("ATGCCC", 4), ("ATGTAA", 4),
                                               ("ATG", 1)])
def test_an_atg_before_the_match_is_the_start_only_in_an_ncoi_context(linker, start):
    """EGFP without its ATG, after a linker. ATG GGC is the NcoI start (CCATGG) with a codon
    put in; ATG CCC is a Met in a linker, and ATG TAA would put a stop in the protein: neither
    is EGFP's start."""
    [egfp] = [f for f in catalogue() if f.name == "EGFP"]
    filler = clean_dna(400, 313)
    for strand in (1, -1):
        insert = linker + egfp.sequence[3:]
        construct = filler + (insert if strand == 1 else rc(insert)) + filler
        [match] = [m for m in detect(construct, False) if m.name == "EGFP"]
        assert match.part[0] == start, (linker, strand)


@pytest.mark.parametrize(("edit", "frameshift"), [(1, True), (2, True), (3, False)])
def test_an_indel_that_breaks_the_reading_frame_is_flagged(edit, frameshift):
    """One or two bases missing in the middle of EGFP leave a 99 % match whose protein is
    garbage after that point; a whole codon missing does not."""
    [egfp] = [f for f in catalogue() if f.name == "EGFP"]
    filler = clean_dna(400, 312)
    construct = filler + egfp.sequence[:360] + egfp.sequence[360 + edit :] + filler
    [match] = [m for m in detect(construct, False, min_identity=95.0) if m.name == "EGFP"]
    assert match.frameshift is frameshift
    assert match.as_dict(len(construct))["frameshift"] is frameshift


def test_an_insertion_and_a_deletion_that_make_up_for_each_other_are_flagged():
    """One base in at 300, one out at 330: the frame is right again after, but the ten codons
    between are another protein."""
    [egfp] = [f for f in catalogue() if f.name == "EGFP"]
    seq = egfp.sequence
    filler = clean_dna(400, 314)
    edited = seq[:300] + ("A" if seq[300] != "A" else "C") + seq[300:330] + seq[331:]
    [match] = [m for m in detect(filler + edited + filler, False) if m.name == "EGFP"]
    assert len(edited) == len(seq) and match.frameshift


def test_one_mismatch_in_a_long_element_is_not_shown_as_100_percent():
    assert shown_identity(100.0 * 3344 / 3345) == 99.9
    assert shown_identity(98.0) == 98.0
    assert shown_identity(100.0 * 49 / 50) == 98.0
    assert shown_identity(100.0) == 100.0


def test_a_truncated_itr_is_not_padded_to_full_length():
    sequence, circular = plasmid("MT612431.1")
    itrs = [m for m in detect(sequence, circular) if m.type == "repeat_region"]
    assert all(100 <= len(m.span(len(sequence))) <= 135 for m in itrs)
    assert all(m.start <= m.end for m in itrs)       # no false wrap around the origin


def test_itrs_are_found_at_both_ends_of_an_aav_vector():
    for accession in ("PZ629150.1", "MT612431.1"):
        sequence, circular = plasmid(accession)
        itrs = [m for m in detect(sequence, circular) if m.type == "repeat_region"]
        assert len(itrs) == 2, accession
        assert itrs[0].start < 200, accession
        assert itrs[1].start > len(sequence) / 2, accession


def test_the_same_region_is_not_annotated_twice():
    for accession in ("PZ629150.1", "MT612431.1"):
        sequence, circular = plasmid(accession)
        matches = detect(sequence, circular)
        for i, first in enumerate(matches):
            for second in matches[i + 1 :]:
                if first.type != second.type:
                    continue
                overlap = len(first.span(len(sequence)) & second.span(len(sequence)))
                smaller = min(len(first.span(len(sequence))), len(second.span(len(sequence))))
                assert overlap < 0.9 * smaller, f"{first.name} / {second.name} in {accession}"


def test_a_tag_inside_a_coding_sequence_is_still_reported():
    sequence, circular = plasmid("PZ629150.1")
    found = names(detect(sequence, circular))
    assert "SaCas9" in found and "3×HA" in found


def test_random_sequence_gets_no_annotations():
    assert detect(clean_dna(3000, 301), False) == []


def test_a_short_site_must_match_exactly():
    [loxp] = [f for f in catalogue() if f.name == "loxP"]
    assert loxp.length <= SHORT_FEATURE
    filler = clean_dna(500, 302)
    assert names(detect(filler + loxp.sequence + filler, False)) == ["loxP"]
    mutated = loxp.sequence[:10] + ("A" if loxp.sequence[10] != "A" else "C") + loxp.sequence[11:]
    assert detect(filler + mutated + filler, False) == []


def test_a_feature_crossing_the_origin_of_a_plasmid_is_found():
    [wpre] = [f for f in catalogue() if f.name == "WPRE"]
    tail, head = wpre.sequence[:200], wpre.sequence[200:]
    sequence = head + clean_dna(1000, 303) + tail
    [match] = [m for m in detect(sequence, circular=True) if m.name == "WPRE"]
    assert match.start > match.end                      # it wraps
    assert match.start == len(sequence) - 200


def test_payload_flags_which_features_are_new():
    sequence, circular = plasmid("PZ629150.1")
    existing = [{"name": "SaCas9", "start": 580, "end": 3924}]
    payload = detect_payload(sequence, circular, existing=existing)
    rows = {row["name"]: row for row in payload["features"]}
    assert rows["SaCas9"]["isNew"] is False
    assert rows["U6 promoter (human)"]["isNew"] is True
    assert payload["new_count"] == len(payload["features"]) - 1
    assert payload["catalogue_size"] == len(catalogue())


def test_an_element_annotated_under_another_name_is_not_offered_again():
    """The published destination vector carries CDS “nptI” at 820..1635 (complement) and an
    undirected rep_origin “colE1” at 18..672: the catalogue's KanR and ori fall on them and
    must not be offered as new features, which would add overlapping duplicates."""
    from vettore.fileio import read_documents

    path = Path(__file__).parent / "data" / "cloning-gateway-destination.gb"
    [document] = read_documents(path.name, path.read_bytes())
    payload = detect_payload(document["sequence"], True, existing=document["features"])
    rows = {row["name"]: row for row in payload["features"]}
    assert rows["KanR"]["isNew"] is False
    assert rows["KanR"]["annotatedAs"] in {"nptI", "aminoglycoside 3'-phosphotransferase"}
    assert rows["ori"]["isNew"] is False and rows["ori"]["annotatedAs"] == "colE1"
    assert payload["new_count"] == sum(1 for row in payload["features"] if row["isNew"])


def test_the_unannotated_fixture_is_the_aav_plasmid_without_its_features():
    """The end-to-end test adds detected features to a plasmid that has none; the published
    record annotates every element, so it would rightly offer nothing to add."""
    from vettore.fileio import read_documents

    data = Path(__file__).parent / "data"
    [bare] = read_documents("x.gb", (data / "PZ629150.1-unannotated.gb").read_bytes())
    [full] = read_documents("y.gb", (data / "PZ629150.1.gb").read_bytes())
    assert bare["sequence"] == full["sequence"] and bare["circular"] is True
    assert bare["features"] == [] and len(full["features"]) == 16
    payload = detect_payload(bare["sequence"], True, existing=bare["features"])
    assert payload["new_count"] == len(payload["features"]) == 9
    annotated = detect_payload(full["sequence"], True, existing=full["features"])
    assert annotated["new_count"] == 0
    rows = {row["name"]: row["annotatedAs"] for row in annotated["features"]}
    assert rows["SaCas9"] == "SaCas9 fusion protein"


def test_a_large_feature_does_not_hide_the_elements_inside_it():
    sequence, circular = plasmid("PZ629150.1")
    everything = [{"name": "backbone", "start": 0, "end": len(sequence) - 1}]
    payload = detect_payload(sequence, circular, existing=everything)
    assert all(row["isNew"] and row["annotatedAs"] is None for row in payload["features"])


def test_malformed_existing_features_are_ignored():
    sequence, circular = plasmid("PZ629150.1")
    junk = [None, "x", {"start": "a", "end": 3}, {"start": 5}, {"locations": "bad"},
            {"name": "SaCas9", "start": 580, "end": 3924}, {"start": -4, "end": 10**9}]
    rows = {row["name"]: row for row in detect_payload(sequence, circular,
                                                       existing=junk)["features"]}
    assert rows["SaCas9"]["annotatedAs"] == "SaCas9"
    assert detect_payload(sequence, circular, existing="bad")["new_count"] > 0


def test_the_database_file_is_the_one_the_builder_produces():
    """The builder needs the network, so this only checks the file's own consistency."""
    text = (Path(__file__).resolve().parents[1] / "src" / "vettore" / "data"
            / "common_features.json").read_text(encoding="utf-8")
    payload = json.loads(text)
    assert payload["source"].startswith("NCBI")
    assert len(payload["features"]) == len(catalogue())
    for entry in payload["features"]:
        assert entry["length"] == len(entry["sequence"])


def test_http_route():
    from fastapi.testclient import TestClient

    from vettore.server import TOKEN_HEADER, create_app

    sequence, circular = plasmid("PZ629150.1")
    client = TestClient(create_app("t0ken", None), base_url="http://127.0.0.1")
    res = client.post("/api/features/detect", headers={TOKEN_HEADER: "t0ken"},
                      json={"sequenceData": {"name": "p", "sequence": sequence,
                                             "circular": circular, "features": []}})
    assert res.status_code == 200
    body = res.json()
    assert {f["name"] for f in body["features"]} >= {"SaCas9", "U6 promoter (human)", "AmpR"}
    assert all(f["isNew"] for f in body["features"])
    short = client.post("/api/features/detect", headers={TOKEN_HEADER: "t0ken"},
                        json={"sequenceData": {"sequence": "ACGT"}})
    assert short.status_code == 422


def test_the_lightly_annotated_sample_vector_gets_its_backbone_named():
    """pUC19 is offered on the welcome screen as the construct to try feature detection on, so
    detection has to find something useful in it."""
    from vettore import app_info

    [sample] = [s for s in app_info.SAMPLES if s.id == "puc19"]
    [doc] = app_info.sample_documents(sample)
    found = names(detect(doc["sequence"].upper(), doc["circular"]))
    assert "AmpR" in found and "ori" in found


def test_an_unannotated_aav_cassette_gives_the_packaging_check_two_itrs():
    """The AAV packaging check falls back on this search when a document has no ITR feature."""
    [itr] = [f for f in catalogue() if f.name == "AAV2 ITR"]
    filler = clean_dna(2000, 11)
    matches = [m for m in detect(itr.sequence + filler + itr.sequence, circular=True)
               if m.type == "repeat_region"]
    assert len(matches) == 2
    assert matches[0].start == 0
    assert matches[1].start == len(itr.sequence) + len(filler)
    # ITR to ITR, ITRs included — the number the check reports.
    assert matches[1].end - matches[0].start + 1 == 2 * len(itr.sequence) + len(filler)


# ---- names and labels -------------------------------------------------------------------
# fidelity-block-cmv.gb and fidelity-mislabelled-cmv.gb are a CMV-EGFP plasmid assembled from
# catalogue entries, in this order: CMV enhancer (PV031350.1:96..475) and CMV promoter
# (PV031350.1:476..679), contiguous in that record, EGFP (U55762.1:679..1398), SV40 polyA
# (MH282430.1:2999..3230), ori (PZ629150.1:6431..7019), AmpR promoter and AmpR
# (PV031350.1:2514..3479, contiguous). One block annotation covers the CMV enhancer and
# promoter (1..584), as vector vendors write it: "CMV promoter" in the first file, "PGK
# promoter" in the second. A real case: a vendor's AAV vector labels its CMV enhancer and
# promoter "PGK promoter" while its own /official_designation says CMV.


def fixture_payload(name: str, existing=None) -> dict:
    from vettore.fileio import read_documents

    [document] = read_documents(name, (DATA / name).read_bytes())
    features = document["features"] if existing is None else existing
    return detect_payload(document["sequence"], document["circular"], existing=features)


def rows_by_name(payload: dict) -> dict[str, dict]:
    return {row["name"]: row for row in payload["features"]}


def test_the_cmv_fixtures_are_built_from_the_catalogue():
    from vettore.fileio import read_documents

    parts = {f.name: f.sequence for f in catalogue()}
    order = ["CMV enhancer", "CMV promoter", "EGFP", "SV40 polyA", "ori", "AmpR promoter",
             "AmpR"]
    for name in ("fidelity-block-cmv.gb", "fidelity-mislabelled-cmv.gb"):
        [document] = read_documents(name, (DATA / name).read_bytes())
        assert document["sequence"].upper() == "".join(parts[p] for p in order), name
        [block] = [f for f in document["features"] if (f["start"], f["end"]) == (0, 583)]
        assert block["type"] == "promoter"


def test_an_element_inside_a_block_named_after_it_is_already_annotated():
    """The CMV promoter is 204 of the block's 584 bp: by size alone it would be offered again,
    on top of the block that already annotates it."""
    payload = fixture_payload("fidelity-block-cmv.gb")
    rows = rows_by_name(payload)
    assert payload["new_count"] == 0
    for element in ("CMV enhancer", "CMV promoter"):
        assert rows[element]["isNew"] is False, element
        assert rows[element]["annotatedAs"] == "CMV promoter", element
        assert rows[element]["labelWarning"] is None, element
    assert payload["label_warnings"] == 0


def test_a_block_named_after_an_element_annotates_it_however_long_it_is():
    # One block from the CMV enhancer to the end of EGFP: 1304 bp, named after both.
    block = [{"name": "CMV-EGFP", "start": 0, "end": 1303}]
    rows = rows_by_name(fixture_payload("fidelity-block-cmv.gb", existing=block))
    assert rows["CMV promoter"]["annotatedAs"] == "CMV-EGFP"
    assert rows["EGFP"]["annotatedAs"] == "CMV-EGFP"
    assert rows["SV40 polyA"]["isNew"] is True


def test_a_cmv_promoter_labelled_pgk_is_flagged_and_not_ticked():
    payload = fixture_payload("fidelity-mislabelled-cmv.gb")
    rows = rows_by_name(payload)
    assert rows["CMV enhancer"]["labelWarning"] == (
        "Annotated as “PGK promoter”, but the sequence is the CMV enhancer (100 % identity). "
        "Check the label.")
    assert rows["CMV promoter"]["labelWarning"] == (
        "Annotated as “PGK promoter”, but the sequence is the CMV promoter (100 % identity). "
        "Check the label.")
    for element in ("CMV enhancer", "CMV promoter"):
        assert rows[element]["isNew"] is False and rows[element]["annotatedAs"] == "PGK promoter"
    assert [r["name"] for r in payload["features"] if r["labelWarning"]] == [
        "CMV enhancer", "CMV promoter"]
    assert payload["label_warnings"] == 1            # one label, two elements under it


def test_a_name_that_is_not_an_element_raises_no_label_warning():
    mine = [{"name": "my insert", "start": 0, "end": 583}]
    payload = fixture_payload("fidelity-mislabelled-cmv.gb", existing=mine)
    rows = rows_by_name(payload)
    assert payload["label_warnings"] == 0
    assert rows["CMV enhancer"]["annotatedAs"] == "my insert"       # 380 of 584 bp
    assert rows["CMV promoter"]["isNew"] is True                     # 204 of 584 bp


def test_a_short_element_at_the_edge_of_a_named_feature_does_not_cast_doubt():
    """An "EGFP" feature that runs on over the SV40 polyA (232 of its 952 bp) is still EGFP:
    the polyA is offered as a new feature, and the label is not questioned."""
    block = [{"name": "EGFP", "start": 584, "end": 1535}]
    payload = fixture_payload("fidelity-block-cmv.gb", existing=block)
    rows = rows_by_name(payload)
    assert payload["label_warnings"] == 0
    assert rows["SV40 polyA"]["labelWarning"] is None and rows["SV40 polyA"]["isNew"] is True
    assert rows["EGFP"]["annotatedAs"] == "EGFP"


def test_a_cag_block_is_not_a_wrong_label_for_its_cmv_enhancer():
    """CAG is the CMV enhancer, the chicken β-actin promoter and the rabbit β-globin intron
    (Niwa, Yamamura & Miyazaki 1991, Gene 108:193-199)."""
    [enhancer] = [f for f in catalogue() if f.name == "CMV enhancer"]
    sequence = enhancer.sequence + clean_dna(1500, 41)
    cag = [{"name": "CAG", "start": 0, "end": 1599}]
    [row] = detect_payload(sequence, True, existing=cag)["features"]
    assert row["name"] == "CMV enhancer" and row["labelWarning"] is None


def test_a_label_is_only_doubted_for_a_nearly_exact_match():
    """At 97.9 % identity a variant is likelier than a wrong label."""
    [enhancer] = [f for f in catalogue() if f.name == "CMV enhancer"]
    bases = list(enhancer.sequence)
    for position in (30, 50, 100, 120, 180, 200, 250, 270):     # away from the search seeds
        bases[position] = "A" if bases[position] != "A" else "C"
    sequence = "".join(bases) + clean_dna(1500, 42)
    pgk = [{"name": "PGK promoter", "start": 0, "end": 379}]
    [row] = detect_payload(sequence, True, existing=pgk)["features"]
    assert row["name"] == "CMV enhancer" and 95 <= row["identity"] < 98
    assert row["labelWarning"] is None and row["annotatedAs"] == "PGK promoter"


@pytest.mark.parametrize("name", ["PZ629150.1.gb", "MT612431.1.gb", "cloning-gateway-entry.gb",
                                  "cloning-gateway-destination.gb",
                                  "cloning-gateway-expression.gb", "fidelity-notes.gb",
                                  "fidelity-block-cmv.gb"])
def test_published_and_correct_labels_raise_no_warning(name):
    payload = fixture_payload(name)
    assert payload["label_warnings"] == 0
    assert all(row["labelWarning"] is None for row in payload["features"])


@pytest.mark.parametrize(("label", "element", "expected"), [
    ("CMV", "CMV enhancer", True), ("CMV promoter", "CMV enhancer", True),
    ("EF1A", "EF-1α promoter (human)", True),
    ("EF-1a promoter", "EF-1α promoter (human)", True),
    ("hEF1a", "EF-1α promoter (human)", True), ("SV40 late pA", "SV40 polyA", True),
    ("BGH pA", "bGH polyA", True), ("5' ITR", "AAV2 ITR", True),
    ("hSyn1 promoter", "hSyn promoter", True), ("CMV-promoter", "CMV promoter", True),
    ("HA tag", "3×HA", True), ("EGFP:T2A:Puro", "EGFP", True),
    ("CMV promoter", "SV40 polyA", False), ("pUC ori", "f1 ori", False),
    ("Exon 1a", "EF-1α promoter (human)", False), ("PGK promoter", "CMV promoter", False),
    ("promoter", "CMV promoter", False), ("my insert", "CMV enhancer", False),
])
def test_names_are_matched_on_the_words_that_identify_an_element(label, element, expected):
    from vettore.features import _names_element

    assert _names_element(label, element) is expected


def test_name_words_ignore_generic_words_and_numbers():
    from vettore.features import name_words

    assert name_words("SV40 late polyA signal") == {"sv40"}
    assert name_words("EF-1α promoter (human)", pieces=False) == {"ef1a"}
    assert name_words("5' ITR") == {"itr"}
    assert name_words("promoter") == set()


def test_holds_lists_every_catalogue_element_found_inside_another():
    """A catalogue element found inside another catalogue entry must be listed in HOLDS, so a
    block named after the larger one is never reported as a wrong label for the smaller."""
    from vettore.features import HOLDS

    entries = catalogue()
    found = set()
    for outer in entries:
        others = tuple(f for f in entries if f is not outer)
        found |= {(outer.name, m.name) for m in detect(outer.sequence, False, features=others)}
    listed = {(outer, inner) for outer, inners in HOLDS.items() for inner in inners}
    assert found <= listed
    names = {f.name for f in entries}
    assert all(outer in names and inner <= names for outer, inner in HOLDS.items())


def test_the_full_length_cag_runs_from_the_cmv_enhancer_to_the_splice_acceptor():
    """The full-length CAG entry (ON012604.1:4799..6471) starts with the CMV enhancer, which
    is why a search finds the enhancer inside it and HOLDS lists the pair, and ends at the
    acceptor of the chimeric intron — about 1.7 kb, the whole promoter."""
    parts = {f.name: f.sequence for f in catalogue()}
    enhancer = parts["CMV enhancer"]
    cag = parts["CAG promoter"]
    assert cag[: len(enhancer)] == enhancer
    assert cag.endswith("CCTACAG")                  # the 3′ splice site of the intron
    assert 1600 <= len(cag) <= 1700


def test_the_truncated_cag_promoter_starts_inside_the_cmv_enhancer():
    """The truncated entry (ON785707.1:189..1125) begins with the last 249 bp of the CMV
    enhancer, too little for a search to find: that containment is listed in HOLDS by hand."""
    parts = {f.name: f.sequence for f in catalogue()}
    head, tail = parts["CAG promoter (truncated)"][:249], parts["CMV enhancer"][-249:]
    assert sum(a != b for a, b in zip(head, tail, strict=True)) <= 1


# ---- the promoters the catalogue gained -------------------------------------------------
# Every entry is a subsequence of the record its ``source`` names, and tests/data carries two
# records for each: the entry's own, so the cited coordinates can be checked offline, and an
# independent one — another laboratory, another vector — that holds the same element, so
# finding it there tests the sequence rather than the record it came from. The provenance of
# each fixture is in tests/data/README.md.

CROSS_CHECKS = {
    "PGK promoter (human)": "MT921016.1",
    "PGK promoter (mouse)": "MT891341.1",
    "UbC promoter (human)": "MN389527.1",
    "H1 promoter (human)": "FJ687158.1",
    "7SK promoter (human)": "MQ296444.1",
    "SFFV promoter": "MP040844.1",
    "CAG promoter": "HG530138.1",
    # The ubiquitous promoters a construct carries while it is checked in a cell line, before
    # it goes into neurons, and the parts those cassettes are built from.
    "EFS promoter (human)": "MP933339.1",
    "CBh promoter": "QQ167606.1",
    "SV40 early promoter": "M83237.1",
    "β-actin promoter (human)": "AY952326.1",
    "HSV TK promoter": "AJ277960.1",
    "RSV LTR promoter": "M83237.1",
    "TRE promoter (tight)": "OQ427855.1",
    "Cre": "AB542060.1",
    "FLPe": "KM521212.1",
    "FLPo": "MH493814.1",
    "loxN": "KP686439.1",
    "F3": "HQ333206.1",
    "EMCV IRES": "OQ427855.1",
    "P2A": "MG880225.1",
    "T2A": "PZ285977.1",
    "WPRE3": "MQ017553.1",
    "hGH polyA": "MP470633.1",
    "rabbit β-globin polyA": "OP697988.1",
    "PuroR (pac)": "PZ267690.1",
    "BlastR (bsd)": "PZ267689.1",
    "HygR (hph)": "U40398.1",
    "ZeoR (ble)": "MW471125.1",
    "Myc tag": "MF966512.1",
    "V5 tag": "AY682091.1",
}


def entry(name: str):
    [feature] = [f for f in catalogue() if f.name == name]
    return feature


def record_sequence(accession: str) -> str:
    record = SeqIO.read(DATA / f"features-{accession}.fasta", "fasta")
    return str(record.seq).upper()


@pytest.mark.parametrize("element", sorted(CROSS_CHECKS))
def test_a_catalogue_entry_is_its_record_at_the_coordinates_it_cites(element):
    """The entry's own record is here too, so the claim its ``source`` makes — these are the
    record's bases at these coordinates — is checked without the network, and the element is
    found there whole and exact."""
    feature = entry(element)
    accession, start, end, strand = cited(feature.source)
    record = record_sequence(accession)
    span = record[start - 1 : end]
    assert (rc(span) if strand == -1 else span) == feature.sequence, element
    matches = [m for m in detect(record, False) if m.name == element]
    # A vector can carry the same element twice (two loxN sites in one cassette), so the copy
    # the source cites must be among the matches rather than the only one; each copy found is
    # whole and exact.
    assert matches, element
    assert all(m.identity == 100.0 for m in matches), element
    assert all(m.end - m.start + 1 == feature.length for m in matches), element
    assert (start - 1, end - 1, strand) in {(m.start, m.end, m.strand) for m in matches}


def shared_span(feature, found: str) -> tuple[float, float]:
    """(identity, coverage) of a catalogue ``feature`` against a copy of it found in another
    record, measured here rather than taken from the match: identity over the bases the two
    share, coverage over the entry's own length. The two are reported apart because a second
    record often carries the element with a few bases less at one end — the last 19 bp of the
    EFS promoter, 27 of FLPe and FLPo, 22 of hph — and averaging a missing flank into an
    identity would hide how well the bases that are there agree.

    Computed with a plain local alignment, not with ``Match.identity``, so that this check
    tests the sequences and not the engine that found them.
    """
    aligner = Align.PairwiseAligner(scoring="blastn")
    aligner.mode = "local"
    alignment = aligner.align(feature.sequence, found)[0]
    reference, query = alignment[0], alignment[1]
    same = sum(1 for a, b in zip(reference, query) if a == b and a != "-")
    identity = 100.0 * same / len(reference)
    aligned = sum(1 for a in reference if a != "-")
    return identity, 100.0 * aligned / feature.length


@pytest.mark.parametrize(("element", "accession"), sorted(CROSS_CHECKS.items()))
def test_a_catalogue_entry_is_found_in_an_independent_record(element, accession):
    feature = entry(element)
    assert not feature.source.startswith(accession), element
    sequence = record_sequence(accession)
    matches = [m for m in detect(sequence, False) if m.name == element]
    assert matches, element
    assert all(m.identity >= 98.0 for m in matches), element
    found = max(m.end - m.start + 1 for m in matches)
    assert found >= 0.9 * feature.length, element
    # The same claim again, measured outside the engine: the bases the two records share
    # agree at 98 % or better, over at least nine tenths of the entry.
    best = max(matches, key=lambda m: m.identity)
    copy = sequence[best.start : best.end + 1]
    identity, coverage = shared_span(feature, rc(copy) if best.strand == -1 else copy)
    assert identity >= 98.0, f"{element}: {identity:.1f} % over the shared bases"
    assert coverage >= 90.0, f"{element}: {coverage:.1f} % of the entry found"


def test_every_cross_checked_element_is_found_side_by_side_in_one_construct():
    """Every cross-checked element in one sequence, promoters, recombinases, markers and tags
    together: each is reported once, whole, at its own position and at full identity, so none
    of them swallows or hides another."""
    filler = [clean_dna(300, 60 + i) for i in range(len(CROSS_CHECKS) + 1)]
    pieces: list[str] = []
    starts: dict[str, int] = {}
    for index, name in enumerate(CROSS_CHECKS):
        pieces.append(filler[index])
        starts[name] = len("".join(pieces))
        pieces.append(entry(name).sequence)
    sequence = "".join(pieces) + filler[-1]
    found = {m.name: m for m in detect(sequence, False)}
    for name in CROSS_CHECKS:
        assert found[name].identity == 100.0, name
        assert found[name].start == starts[name], name
        assert found[name].end - found[name].start + 1 == entry(name).length, name


def test_a_plasmid_built_from_the_catalogue_shows_one_full_length_cag():
    """The full-length CAG opens with the whole CMV enhancer, so a search finds both. Only the
    CAG is reported: the enhancer is part of it, not a second element next to it."""
    parts = {f.name: f.sequence for f in catalogue()}
    order = ["AAV2 ITR (truncated)", "CAG promoter", "EGFP", "bGH polyA",
             "AAV2 ITR (truncated)", "ori"]
    sequence = "".join(parts[name] for name in order)
    found = {m.name: m for m in detect(sequence, circular=True)}
    cag = found["CAG promoter"]
    assert cag.identity == 100.0
    assert cag.start == len(parts["AAV2 ITR (truncated)"])
    assert cag.end - cag.start + 1 == len(parts["CAG promoter"])
    assert "CMV enhancer" not in found and "CAG promoter (truncated)" not in found


def test_a_block_named_cag_already_annotates_the_promoter_under_it():
    """A vector annotated with one “CAG” block, as vendors write it: the CAG is ticked off as
    annotated, nothing inside it is offered again and the label is not questioned."""
    parts = {f.name: f.sequence for f in catalogue()}
    head = parts["AAV2 ITR (truncated)"]
    sequence = head + parts["CAG promoter"] + parts["EGFP"] + parts["bGH polyA"]
    block = [{"name": "CAG", "start": len(head),
              "end": len(head) + len(parts["CAG promoter"]) - 1}]
    payload = detect_payload(sequence, True, existing=block)
    rows = rows_by_name(payload)
    assert rows["CAG promoter"]["isNew"] is False
    assert rows["CAG promoter"]["annotatedAs"] == "CAG"
    assert payload["label_warnings"] == 0
    assert all(row["labelWarning"] is None for row in payload["features"])
    assert "CMV enhancer" not in rows


def test_a_pgk_promoter_labelled_cmv_is_flagged():
    """The other way round from the CMV block labelled “PGK promoter”: now the sequence is the
    PGK promoter and the label names the CMV one."""
    pgk = entry("PGK promoter (human)")
    sequence = pgk.sequence + clean_dna(1500, 44)
    mislabelled = [{"name": "CMV promoter", "start": 0, "end": pgk.length - 1}]
    payload = detect_payload(sequence, True, existing=mislabelled)
    [row] = payload["features"]
    assert row["name"] == "PGK promoter (human)"
    assert row["labelWarning"] == (
        "Annotated as “CMV promoter”, but the sequence is the PGK promoter (human) "
        "(100 % identity). Check the label.")
    assert row["isNew"] is False and row["annotatedAs"] == "CMV promoter"
    assert payload["label_warnings"] == 1


@pytest.mark.parametrize(("label", "element", "expected"), [
    ("PGK promoter", "PGK promoter (human)", True),
    ("PGK promoter", "PGK promoter (mouse)", True),
    ("hPGK promoter", "PGK promoter (human)", True),
    ("UBC", "UbC promoter (human)", True), ("UbC-EGFP", "UbC promoter (human)", True),
    ("H1", "H1 promoter (human)", True), ("H1-shRNA", "H1 promoter (human)", True),
    ("7SK promoter", "7SK promoter (human)", True), ("SFFV", "SFFV promoter", True),
    ("CAG", "CAG promoter", True), ("CAG", "CAG promoter (truncated)", True),
    ("U6 promoter", "H1 promoter (human)", False),
    ("PGK promoter", "UbC promoter (human)", False),
    ("SFFV promoter", "CMV promoter", False), ("CMV promoter", "PGK promoter (human)", False),
])
def test_the_new_promoters_are_named_by_the_labels_vectors_use(label, element, expected):
    from vettore.features import _names_element

    assert _names_element(label, element) is expected


def test_the_species_in_a_name_does_not_tell_two_promoters_apart():
    """The catalogue names the species, so the table says which promoter a hit is. The word
    is not one that identifies an element: a label that says only "PGK promoter" names both
    entries, and is not read as naming some element the catalogue does not have."""
    from vettore.features import name_words, named_elements

    assert name_words("PGK promoter (human)", pieces=False) == {"pgk"}
    assert name_words("PGK promoter (mouse)", pieces=False) == {"pgk"}
    assert name_words("U6 promoter (human)", pieces=False) == {"u6"}
    assert named_elements("PGK promoter") == {"PGK promoter (human)", "PGK promoter (mouse)"}


def test_the_human_and_mouse_pgk_promoters_are_two_sequences_not_one():
    """Two entries earn their place only if a search tells them apart: the human promoter is
    not found in the mouse one, nor the other way round."""
    human, mouse = entry("PGK promoter (human)"), entry("PGK promoter (mouse)")
    assert names(detect(human.sequence, False)) == ["PGK promoter (human)"]
    assert names(detect(mouse.sequence, False)) == ["PGK promoter (mouse)"]


def _plasmid_with(element: str, label: str) -> tuple[str, list[dict]]:
    """A linear construct: random flanks around a catalogue element annotated as ``label``."""
    import random

    sequence = next(f.sequence for f in catalogue() if f.name == element)
    rng = random.Random(7)
    left = "".join(rng.choice("ACGT") for _ in range(300))
    right = "".join(rng.choice("ACGT") for _ in range(300))
    annotation = {"name": label, "type": "promoter", "start": 300,
                  "end": 300 + len(sequence) - 1, "strand": 1}
    return left + sequence + right, [annotation]


@pytest.mark.parametrize("element, label, wrong", [
    ("PGK promoter (mouse)", "PGK promoter (human)", True),
    ("PGK promoter (mouse)", "hPGK", True),
    ("PGK promoter (mouse)", "mPGK promoter", False),
    ("PGK promoter (mouse)", "PGK promoter", False),
    ("PGK promoter (human)", "murine PGK promoter", True),
    ("U6 promoter (human)", "mU6", True),
    ("U6 promoter (human)", "hU6 promoter", False),
])
def test_a_label_naming_the_other_species_is_flagged(element, label, wrong):
    sequence, existing = _plasmid_with(element, label)
    out = detect_payload(sequence, False, DEFAULT_MIN_IDENTITY, existing)
    [row] = [r for r in out["features"] if r["name"] == element]
    assert row["isNew"] is False and row["annotatedAs"] == label
    assert bool(row["labelWarning"]) is wrong
    assert out["label_warnings"] == (1 if wrong else 0)
    if wrong:
        assert row["labelWarning"].startswith(f"Annotated as “{label}”, but the sequence is the {element}")


# ---- the elements of the cell-line pretests ---------------------------------------------
# The constructs are not used in neurons only: they are checked in a cell line first, and
# those pretests drive the payload from a ubiquitous promoter (S5.4, group 1b). Each of those
# promoters is in the catalogue, named like the rest.

PRETEST_PROMOTERS = ("CMV promoter", "SV40 early promoter", "β-actin promoter (human)",
                     "EFS promoter (human)", "HSV TK promoter", "RSV LTR promoter")


def test_the_catalogue_carries_the_promoters_of_the_cell_line_pretests():
    held = {f.name for f in catalogue()}
    assert set(PRETEST_PROMOTERS) <= held
    # The CMV/chicken β-actin hybrids the pretests also use, and the inducible promoter.
    assert {"CAG promoter", "CAG promoter (truncated)", "CBh promoter",
            "TRE promoter (tight)"} <= held


@pytest.mark.parametrize("element", PRETEST_PROMOTERS)
def test_each_pretest_promoter_is_found_whole_on_either_strand(element):
    feature = entry(element)
    ahead = clean_dna(200, 31) + feature.sequence + clean_dna(200, 32)
    behind = clean_dna(200, 33) + rc(feature.sequence) + clean_dna(200, 34)
    [forward] = [m for m in detect(ahead, False) if m.name == element]
    [reverse] = [m for m in detect(behind, False) if m.name == element]
    assert (forward.strand, forward.start, forward.identity) == (1, 200, 100.0)
    assert (reverse.strand, reverse.start, reverse.identity) == (-1, 200, 100.0)
    assert forward.end - forward.start + 1 == feature.length


def test_the_rsv_promoter_is_stored_as_a_vector_reads_it_though_its_record_is_reversed():
    """LT727488.1 annotates the RSV U3 promoter on its minus strand, so the entry is the
    reverse complement of the cited span and its source says ``complement(...)``. What the
    catalogue holds is the promoter the way a vector carries it, pointing at its payload."""
    feature = entry("RSV LTR promoter")
    accession, start, end, strand = cited(feature.source)
    assert (accession, start, end, strand) == ("LT727488.1", 3997, 4178, -1)
    record = record_sequence(accession)
    assert feature.sequence == rc(record[start - 1 : end])
    # In the independent record the promoter reads forward, which is the orientation stored.
    [match] = [m for m in detect(record_sequence("M83237.1"), False)
               if m.name == "RSV LTR promoter"]
    assert match.strand == 1 and match.identity >= 98.0


def test_the_human_actb_promoter_is_not_the_chicken_one_inside_cag_and_cbh():
    """CAG and CBh are built on the *chicken* β-actin promoter; the new entry is the human
    ACTB promoter. Neither is found in the other, which is what makes them separate entries
    rather than one element named twice."""
    human = entry("β-actin promoter (human)")
    for hybrid in ("CAG promoter", "CAG promoter (truncated)", "CBh promoter"):
        assert not [m for m in detect(entry(hybrid).sequence, False) if m.name == human.name]
        assert not [m for m in detect(human.sequence, False) if m.name == hybrid]


def test_the_efs_core_is_not_reported_beside_the_whole_ef1a_promoter():
    """EFS is the intron-less core of the EF-1α promoter and sits inside the longer entry, so
    a cassette carrying the whole promoter shows one element; the core on its own is EFS."""
    full, core = entry("EF-1α promoter (human)"), entry("EFS promoter (human)")
    assert names(detect(clean_dna(300, 41) + full.sequence + clean_dna(300, 42), False)) \
        == ["EF-1α promoter (human)"]
    assert names(detect(clean_dna(300, 43) + core.sequence + clean_dna(300, 44), False)) \
        == ["EFS promoter (human)"]


def test_a_block_named_ef1a_is_not_a_wrong_label_for_the_efs_core():
    """HOLDS lists the pair, so the label that names the whole promoter is read as naming its
    core too."""
    sequence, existing = _plasmid_with("EFS promoter (human)", "EF-1α promoter")
    out = detect_payload(sequence, False, DEFAULT_MIN_IDENTITY, existing)
    [row] = [r for r in out["features"] if r["name"] == "EFS promoter (human)"]
    assert row["isNew"] is False and row["labelWarning"] is None


def test_wpre3_and_the_full_length_wpre_are_told_apart():
    """WPRE3 is the shortened element AAV cassettes carry, not a piece of the long one: each
    is reported where it is and neither stands for the other."""
    for element in ("WPRE", "WPRE3"):
        sequence = clean_dna(200, 45) + entry(element).sequence + clean_dna(200, 46)
        assert names(detect(sequence, False)) == [element]


def test_flpe_and_flpo_are_two_entries_a_search_tells_apart():
    """Vectors carry both the thermostable yeast gene and the mammalian codon-optimised one,
    and they share too little for one entry to find the other."""
    for element in ("FLPe", "FLPo"):
        assert names(detect(entry(element).sequence, False)) == [element]


def test_the_new_recombination_sites_are_exact_and_not_each_other():
    """A single base in a lox or Flp site changes what the recombinase does with it, so the
    new sites are matched exactly and no site is reported for another."""
    for element in ("loxN", "F3"):
        feature = entry(element)
        assert feature.length <= SHORT_FEATURE, element
        found = names(detect(clean_dna(200, 47) + feature.sequence + clean_dna(200, 48), False))
        assert found == [element], found
    for site in ("loxP", "lox2272", "loxN"):
        others = {"loxP", "lox2272", "loxN"} - {site}
        assert not others & set(names(detect(entry(site).sequence, False)))


@pytest.mark.parametrize(("label", "element", "expected"), [
    ("pac", "PuroR (pac)", True), ("PuroR", "PuroR (pac)", True),
    ("hph", "HygR (hph)", True), ("HygR", "HygR (hph)", True),
    ("bsd", "BlastR (bsd)", True), ("BlastR", "BlastR (bsd)", True),
    ("Sh ble", "ZeoR (ble)", True), ("ZeoR", "ZeoR (ble)", True),
    ("ACTB promoter", "β-actin promoter (human)", True),
    ("hACTB", "β-actin promoter (human)", True),
    ("beta-actin promoter", "β-actin promoter (human)", True),
    ("EFS", "EFS promoter (human)", True), ("CBh", "CBh promoter", True),
    ("myc tag", "Myc tag", True), ("V5", "V5 tag", True),
    ("P2A", "P2A", True), ("T2A", "T2A", True),
    ("pac", "NeoR/KanR", False), ("hph", "PuroR (pac)", False),
    ("P2A", "T2A", False), ("Myc tag", "V5 tag", False),
    ("CBh", "CAG promoter", False), ("WPRE", "WPRE3", False),
])
def test_the_markers_sites_and_tags_are_named_by_the_labels_vectors_use(label, element,
                                                                       expected):
    from vettore.features import _names_element

    assert _names_element(label, element) is expected


@pytest.mark.parametrize(("element", "label", "wrong"), [
    ("β-actin promoter (human)", "hACTB promoter", False),
    ("β-actin promoter (human)", "ACTB promoter", False),
    ("β-actin promoter (human)", "mouse Actb promoter", True),
    ("EFS promoter (human)", "mEFS", True),
    ("EFS promoter (human)", "hEFS promoter", False),
])
def test_the_species_check_covers_the_new_entries(element, label, wrong):
    sequence, existing = _plasmid_with(element, label)
    out = detect_payload(sequence, False, DEFAULT_MIN_IDENTITY, existing)
    [row] = [r for r in out["features"] if r["name"] == element]
    assert row["isNew"] is False and row["annotatedAs"] == label
    assert bool(row["labelWarning"]) is wrong


@pytest.mark.parametrize(("element", "label"), [
    # The CAG and CBh promoters are hybrids, and the labels their users write name the pieces
    # they are made of. Each of those names belongs to a different catalogue entry — the human
    # β-actin promoter, the rabbit β-globin polyA signal — and before the kind and species
    # checks every one of these raised a wrong-label warning on a correctly labelled construct.
    ("CAG promoter", "chicken beta-actin promoter"),
    ("CAG promoter", "chicken β-actin promoter"),
    ("CAG promoter", "rabbit beta-globin intron"),
    ("CBh promoter", "chicken ACTB promoter"),
    # WPRE3 is 246 bp of the WPRE, base for base, so its author's loose "WPRE" is not wrong.
    ("WPRE3", "WPRE"),
])
def test_a_label_naming_what_a_hybrid_is_made_of_is_not_a_wrong_label(element, label):
    sequence, existing = _plasmid_with(element, label)
    out = detect_payload(sequence, False, DEFAULT_MIN_IDENTITY, existing)
    [row] = [r for r in out["features"] if r["name"] == element]
    assert not row["labelWarning"], row["labelWarning"]


def test_a_label_of_the_wrong_kind_or_species_is_still_caught():
    """The two new filters narrow the check; they do not switch it off. A block that really is
    the human β-actin promoter and calls itself the chicken one is still a wrong label."""
    from vettore.features import _kind_conflict, _names_another_element, _species_conflict

    assert _species_conflict("chicken β-actin promoter", "β-actin promoter (human)")
    assert _kind_conflict("rabbit β-globin intron", "rabbit β-globin polyA")
    assert not _kind_conflict("rabbit β-globin polyA", "rabbit β-globin polyA")
    # A promoter label over an enhancer is left alone on purpose: a CMV enhancer-and-promoter
    # block is labelled either way and neither label is wrong.
    assert not _kind_conflict("CMV promoter", "CMV enhancer")
    assert _names_another_element("hSyn1 promoter", "CamKIIα promoter")

    sequence, existing = _plasmid_with("β-actin promoter (human)", "chicken ACTB promoter")
    out = detect_payload(sequence, False, DEFAULT_MIN_IDENTITY, existing)
    [row] = [r for r in out["features"] if r["name"] == "β-actin promoter (human)"]
    assert row["labelWarning"]


# ---- detection speed --------------------------------------------------------------------
# S5.4 asks that detection stay under two seconds for a 15 kb plasmid with the grown
# catalogue. (`specs/vettore.md` lists S5.4 as the catalogue's growth to human and mouse
# elements; the two-second criterion attached to it is in the S5-S6 roadmap, which lives in the
# notes repository, so the seconds are quoted in this file and nowhere under `specs/`.) The cost
# is not in the catalogue's size — a feature absent from the sequence is ruled out by an exact
# seed search and costs almost nothing — but in the elements that *are* there, each of which is
# aligned where its seeds land.
#
# Three mistakes are easy here, and earlier versions of these tests made all three.
#
# The first is to measure a plasmid of elements copied *exactly* from the catalogue. That is
# the best case, not the worst: an exact copy is taken as it stands and costs no alignment at
# all. Promoters, ITRs and reporters differ between vectors by a few bases — the reason this
# module aligns rather than searching — so the case to measure is a plasmid whose elements are
# nearly the catalogue's.
#
# The second is to assume the costly divergence is a changed base. It is not. The cost is the
# number of candidate offsets, not the number of mismatches: a substitution at a position no
# seed covers leaves all six seeds in one register, so the element is aligned once, while an
# **indel** splits the seeds into two registers and the element is aligned twice. So the
# fixtures below diverge by deleting.
#
# What the cost tracks is the number of *places the engine aligns an element*, which used to be
# the number of surviving seed registers and is now the number of copies: registers a few bases
# apart are merged and aligned once (`MERGE_GAP`). That change is what brought the all-in-one
# plasmid inside S5.4's two seconds, and it flattened the curve with it. On the AAV fixture
# (36 entries, one container):
#
#   sites     0      1      2      3      4      5      6
#   cells   0.34 M 15.49M 15.44M 15.55M 29.57M 15.57M 29.56M
#   seconds  0.03   0.39   0.39   0.41   0.75   0.41   0.73
#   same?    yes    yes    yes    yes    yes    yes    NO
#
# A single changed base is 15.3 M, so a deleted base and a changed one now cost about the
# same: that factor of two was the registers, and it is gone. The curve is flat at one pass per
# copy, give or take the AAV2 ITR — a hairpin, so self-similar enough that `_merge_limit` will
# not merge its registers, which is the slow creep from 15.44 to 15.57 M across the row — and
# except where the fixture's even spacing drops a cut inside a seed window. The engine places
# its seeds at `round(i * (length - 18) / 5)` and the fixture its cuts at
# `length // (sites + 1) * k`, which for four cuts coincide, and the registers that survive
# are then further apart than the merge is willing to join, so SaCas9 is aligned twice and the
# plasmid costs 29.6 M. At six sites an element is no longer found at all, which is the
# boundary its own test asserts. The bounds below pin the points where the plasmid still holds
# together.
#
# The third is to write a divergence the engine itself reports as a different document. A
# **single** deleted base inside a coding sequence is a frameshift, and detection says so:
# `Match.frameshift` is set, the dialog leaves the match unticked and shows an amber warning. One
# deleted base per element flags SpCas9, SaCas9, tdTomato, mCherry, EGFP and AmpR, so a
# fixture built that way is not the same plasmid to the user, and a nuclease with a frameshift
# is not a vector anyone builds. Deleting a whole **codon** costs the same (15.48 M on the
# AAV fixture against 15.49 M) and is a real variant, so that is what `_plasmid` does to a CDS.
#
# What the plasmid holds matters as much as how far it diverges, and a fixture that packs
# every awkward entry together measures a vector nobody builds. The three below are a single
# AAV nuclease vector, a dense all-in-one expression plasmid, and that packed plasmid, kept as
# a stress case and named as one.

#: S5.4's criterion, in the seconds it is written in. **Nothing asserts it**, and saying so is
#: the point of this constant: four attempts to assert it, and four adversarial rounds that
#: measured each one, established that no quantity a clock can produce here separates the code
#: from the machine. What the suite asserts instead is the *work*, in cells of dynamic
#: programming, by `ALIGNMENT_WORK`, which is the same number on every machine — three of whose
#: bounds are at points recorded above two seconds, so it is a guard against the work changing
#: and not a restatement of the criterion. The seconds are recorded below.
DETECTION_BUDGET = 2.0

# ---- what detection costs, measured and not asserted --------------------------------------
#
# Every figure here is a measurement, kept because the criterion is written in seconds and
# somebody will want to know whether it is met. None of it is asserted, and the four rounds
# that got it to that state are worth recording, because the obvious fixes all look sound until
# they are measured.
#
# Twelve warm runs per point, 36 entries, on one container idle and with four busy loops on its
# four cores, and the same points on a second container (its figures are the seconds column of
# `ALIGNMENT_WORK`):
#
#   plasmid                      idle           contended      the other container
#   AAV, 1 codon per element     0.46-0.59 s    0.51-0.82 s    0.39 s
#   AAV, 5 codons                0.44-0.62 s    0.54-0.72 s    0.41 s
#   all-in-one, 1 codon          1.49-2.10 s    1.84-2.26 s    1.19 s
#   all-in-one, 5 codons         1.70-2.10 s    2.11-2.61 s    1.53 s
#   all-in-one, 4 codons         2.53-3.10 s    3.31-3.77 s    1.95-2.07 s
#   packed stress, 1 codon       -              -              2.75 s
#
# Read those figures carefully, because the favourable reading is wrong. PR #13 brought the
# all-in-one plasmid from 2.3 s with a codon lost per element and 5.5 s with five down to the
# ranges above, which is what it set out to do. But the criterion is met by the **median** run
# and not by the range: on the slower container the 1-codon and 5-codon points both reach
# 2.10 s idle and the 5-codon point is over budget across its whole contended range, so the
# margin is 1.1x to 1.3x at the fastest runs and nothing at the slowest. The four-site point —
# where the fixture's even cuts land on the engine's seed positions, so the merge cannot help —
# straddles the budget on the faster container, at 1.95 to 2.07 s, and is over it on the
# slower. The same engine is 7 to 23 % quicker on the faster container, point by point against
# this one's fastest idle runs.
#
# Two things keep that from being a failure of S5.4. The fixtures are deliberately the worst
# case, every element diverged so that none takes the exact-copy shortcut, which is not a
# plasmid anyone builds; and the real records in `tests/data` are far cheaper. The dearest of
# them all is `cloning-multisite-entry-l4r1.gb`, 8.6 kb and circular, at 0.62 to 0.71 s
# across warm runs on the slower container, and every other record is at or under 0.50 s. But
# the honest statement is the one below: a test asserting two seconds would report which
# container ran it.
#
# **Why none of it is asserted.** An assertion at two seconds fails on the slower container
# with nothing regressed, which is what sent this down four rounds:
#
#   - the fastest of three timings hides a cost paid once per process, and still failed about
#     one run in ten on a container a third slower, because the jitter here is slow drift that
#     three consecutive runs share rather than spikes a minimum filters out;
#   - a loose ceiling (3 s, then 8 s) cannot see a regression that costs a factor: matching the
#     seeds base by base in Python instead of with `str.find` is 0.47 -> 2.12 s on the AAV
#     fixture at three deleted sites, with identical annotations and identical cells, and
#     2.12 s is inside the slowest honest run of the all-in-one points. A ceiling placed where
#     the fixtures need it is a ceiling that regression fits under;
#   - the share of the time spent outside `PairwiseAligner.align` looked like the
#     machine-independent answer and is not. It is a property of how much of the document needs
#     aligning, not of the engine, so it ranges as widely as the documents do: measured here,
#     healthy, 28.2-30.9 % on the all-in-one fixture, 33.8-36.4 % on the AAV one, 30-38 % on the
#     dearest real record but 55-60 % on two cheaper ones — and **78.7-79.5 % on the verbatim
#     AAV fixture**, where almost nothing is aligned. Against that, the base-by-base scan drops
#     to 53.2-55.4 % on the packed stress plasmid. So a healthy run reaches 79 % while a
#     regressed one sits at 53 %: the two ranges overlap, and no threshold fits between them.
#
# What survives is deterministic. `ALIGNMENT_WORK` pins the alignment work; the annotation tests
# pin what is found; `test_the_catalogue_is_parsed_once_however_many_documents_are_opened` pins
# the one cost too small for any clock; and
# `test_the_document_is_not_walked_in_python_once_per_element` pins the scan the clock could not
# see, by charging the document for what is asked of it rather than timing it. That last one is
# what the fourth refutation above was reaching for, and it needs no threshold at all: the
# engine asks the document for 267 things and a base-by-base scan for 3,962,927 of its bases.

#: A single-AAV nuclease vector, with the elements of pX601 the catalogue carries, in that
#: construct's own order: the AAV ends, the CMV–SaCas9–bGH polyA cassette, the U6 promoter of
#: the sgRNA cassette downstream of it, and the plasmid backbone. (The catalogue has no sgRNA
#: scaffold, so the U6 promoter is all there is of that cassette.) SaCas9 rather than SpCas9
#: because that is the point of the design: ITR to ITR is 4291 bp **including both ITRs**,
#: inside AAV's ~4.7 kb packaging limit, where the same vector built on CAG and SpCas9 would
#: be 6687 bp on the same basis. SaCas9's entry carries its 3×HA tag internally, in frame at
#: base 3262 and immediately before its stop codon, so the tag is not listed separately:
#: appended after the stop it would never be translated. Detection reports it as a tenth
#: match inside the ninth element, which is why the counts below are 10 and not 9.
AAV_ELEMENTS = (
    "AAV2 ITR", "CMV promoter", "SaCas9", "bGH polyA", "U6 promoter (human)",
    "AAV2 ITR (truncated)", "AmpR promoter", "AmpR", "ori",
)

#: The second ITR of a real AAV plasmid is inverted. It is also the only element any of these
#: fixtures puts on the minus strand.
AAV_INVERTED = ("AAV2 ITR (truncated)",)

#: A dense all-in-one expression plasmid, not a viral genome: the U6 promoter of an sgRNA
#: cassette at the head, as lentiCRISPR vectors carry it (the catalogue has no scaffold, so
#: the promoter is all there is of that cassette), then two complete cassettes — each a
#: promoter, a coding sequence and a terminator, with the WPRE where it belongs, in the 3' UTR
#: between the nuclease and its polyA signal — and a backbone of one marker and two origins.
#: 10.4 kb of catalogue elements, one nuclease, and nothing that has to fit a packaging limit.
#: This is the plausible heavy case: the one that used to miss S5.4's budget, at 2.3 s with a
#: codon lost per element. It now meets it on one container, at 1.2 s, and sits on the line on a
#: slower one, 1.5 to 2.1 s across runs with the budget at 2.0 — which is why nothing asserts
#: the seconds. The block under `DETECTION_BUDGET` records them in full.
ALL_IN_ONE_ELEMENTS = (
    "U6 promoter (human)", "EF-1α promoter (human)", "SpCas9", "WPRE", "bGH polyA",
    "CMV promoter", "tdTomato", "SV40 polyA", "AmpR promoter", "AmpR", "ori", "f1 ori",
)

#: The stress case, not a layout: both nucleases at once, and with them the pairs close enough
#: to seed on each other (SpCas9/SaCas9, tdTomato/mCherry, the two ITRs), an element that
#: contains another (the CAG promoter opens with the whole CMV enhancer, SaCas9 ends with
#: 3×HA), and the short exact sites. 14 kb of elements in a 15 kb plasmid. Named rather than
#: taken by length, so that adding a catalogue entry cannot quietly change what is measured.
PACKED_ELEMENTS = (
    "SpCas9", "SaCas9", "CAG promoter", "tdTomato", "mCherry", "EGFP", "AmpR", "ori",
    "AAV2 ITR", "AAV2 ITR (truncated)", "3×HA", "loxP", "lox2272", "FRT",
)

#: The entries S5.4 added, packed the same way: the long recombinases and markers first, then
#: the promoters of the cell-line pretests, then the short sites and tags. 10.3 kb of elements
#: none of which the other fixtures carry, so this is the one that says what the new half of
#: the catalogue costs on its own. Named rather than taken by length, for the same reason as
#: PACKED_ELEMENTS.
GROWN_ELEMENTS = (
    "FLPe", "FLPo", "Cre", "HygR (hph)", "CBh promoter", "PuroR (pac)", "EMCV IRES",
    "β-actin promoter (human)", "hGH polyA", "BlastR (bsd)", "ZeoR (ble)",
    "SV40 early promoter", "TRE promoter (tight)", "EFS promoter (human)", "WPRE3",
    "RSV LTR promoter", "HSV TK promoter", "rabbit β-globin polyA", "P2A", "T2A",
    "loxN", "F3", "Myc tag", "V5 tag",
)

#: Each fixture with the points on the cost curve the tests pin for it: none, one deleted site
#: per element, three, and five where five is still annotated as the verbatim plasmid is. The
#: stress plasmid is not pinned at five because it loses an element there — see
#: `LOSES_AN_ELEMENT`, which is where that boundary is asserted.
FIXTURES = (
    ("aav", AAV_ELEMENTS, AAV_INVERTED, (0, 1, 3, 5)),
    ("all-in-one", ALL_IN_ONE_ELEMENTS, (), (0, 1, 3, 4, 5)),
    ("stress", PACKED_ELEMENTS, (), (0, 1, 3)),
    # Pinned only to one site: the entries S5.4 added include T2A at 54 bp, barely above
    # SHORT_FEATURE, and three deletions in 54 bases leave it no seed at all.
    ("grown", GROWN_ELEMENTS, (), (0, 1)),
)

#: Where each fixture loses an element, and which: the deleted-site count at which the
#: shortest entry the fixture aligns has no seed left anywhere near its copy, measured, and
#: explained in the test that asserts it — the identity floor is never reached. The
#: stress plasmid gets there first because it carries 3×HA at 81 bp, which is also why it is
#: pinned only to three sites.
LOSES_AN_ELEMENT = {
    "aav": (6, "AmpR promoter"),         # 105 bp
    "all-in-one": (6, "AmpR promoter"),  # 105 bp
    "stress": (5, "3×HA"),               #  81 bp, the copy that stands on its own
    "grown": (3, "T2A"),                 #  54 bp, the shortest entry S5.4 aligns
}
#: Both halves of the "grown" boundary are measured, not reasoned: at three sites exactly
#: the T2A goes and nothing else, and nothing is reported that the verbatim plasmid does not
#: report. At two sites nothing is lost at all, and the T2A — cut at 18 and 36 of its 54
#: bases, which leaves the first seed window free — sits at 96.30 %, 1.30 points above the
#: 95 % floor. That is the thinnest margin in this fixture: P2A has 1.97 and the rabbit
#: β-globin polyA 2.98. Worth writing down because adding another short entry to
#: GROWN_ELEMENTS is what would break it, and the failure would read as a budget problem.

#: One case per (fixture, deleted sites), named so a failure says which.
POINTS = tuple(
    pytest.param(label, element_names, inverted, sites, id=f"{label}-{sites}")
    for label, element_names, inverted, counts in FIXTURES for sites in counts
)
DIVERGED_POINTS = tuple(point for point in POINTS if point.values[3])


def _cuts(feature, sites: int) -> tuple[list[int], int]:
    """Where ``_diverged`` cuts ``feature`` and how wide each cut is.

    A coding sequence loses a whole codon at each site, so its reading frame survives and
    detection does not flag it; anything else loses a single base. An entry the engine
    searches rather than aligns is left alone, because deleting a base of a 34 bp loxP does
    not make it a near miss, it makes it absent — the engine's own predicate is
    ``feature.length <= SHORT_FEATURE``.
    """
    if not sites or feature.length <= SHORT_FEATURE:
        return [], 0
    width = 3 if feature.type == "CDS" else 1
    step = feature.length // (sites + 1)
    cuts = [step * index for index in range(1, sites + 1)]
    return [at - at % 3 if width == 3 else at for at in cuts], width


def _diverged(feature, sites: int) -> str:
    """``feature``'s sequence with ``sites`` stretches deleted, evenly spaced."""
    sequence = feature.sequence
    cuts, width = _cuts(feature, sites)
    for at in reversed(cuts):
        sequence = sequence[:at] + sequence[at + width:]
    return sequence


def _losses(element_names: tuple[str, ...],
            sites: int) -> dict[tuple[str, str | None], tuple[int, ...]]:
    """How many bases each copy of each element lost, keyed by (name, the element it sits
    inside or ``None``), largest loss first.

    Keyed by where the copy is and not by name alone, because an entry can appear twice over:
    SaCas9 carries 3×HA at base 3262, so a plasmid that lists both the nuclease and the tag
    has the tag reported twice, and a cut in SaCas9 shortens the copy inside it only if it
    lands inside it. Keyed by name alone, the copy inside SaCas9 could inherit the standalone
    copy's loss and degrade by three points unnoticed; keyed by name with the losses matched
    to the copies in whatever order made them pass, it could still do so by trading values
    with the standalone copy, which is why the test pairs them by position.
    """
    catalogued = {feature.name: feature for feature in catalogue()}
    lost: dict[tuple[str, str | None], list[int]] = {}
    for name in element_names:
        cuts, width = _cuts(catalogued[name], sites)
        lost.setdefault((name, None), []).append(width * len(cuts))
    for name in element_names:
        host = catalogued[name]
        cuts, width = _cuts(host, sites)
        deleted = {cut + offset for cut in cuts for offset in range(width)}
        for inner in catalogue():
            if inner.name == host.name:
                continue
            at = host.sequence.find(inner.sequence)
            # Every copy, and the bases gone from each rather than the cuts that start
            # inside it, so that a cut straddling either end counts for what it removes.
            # No fixture reaches that case: SaCas9's 3×HA sits at 3262..3342 of 3345, and a
            # cut lands inside it only past about forty sites.
            while at >= 0:
                lost.setdefault((inner.name, host.name), []).append(
                    len(deleted & set(range(at, at + inner.length))))
                at = host.sequence.find(inner.sequence, at + 1)
    return {key: tuple(sorted(counts, reverse=True)) for key, counts in lost.items()}


@functools.cache
def _plasmid(element_names: tuple[str, ...], sites: int = 0,
             inverted: tuple[str, ...] = ()) -> str:
    """``element_names`` end to end, diverged by ``sites`` deletions each, padded to 15 kb.

    Whatever ``sites`` is, the plasmid carries the same elements as the verbatim one: all of
    them found, in the same order, as many of them, none flagged, each covering its whole
    entry, and each at an identity no lower than the bases that element lost. It is a
    genuinely shorter plasmid, so every element's reported size, coordinates and identity move
    with the bases it lost, and the dialog says "N bp missing" against the catalogue entry —
    that is the engine being right, not a different annotation.
    ``test_a_plasmid_is_annotated_the_same_however_its_elements_diverge`` holds all of it.
    """
    catalogued = {feature.name: feature for feature in catalogue()}
    unknown = set(element_names) - set(catalogued)
    assert not unknown, f"not in the catalogue: {sorted(unknown)}"
    pieces = []
    for name in element_names:
        sequence = _diverged(catalogued[name], sites)
        pieces.append(rc(sequence) if name in inverted else sequence)
    assert sum(map(len, pieces)) <= 15_000, f"{element_names} no longer fit the budget"
    return ("".join(pieces) + clean_dna(15_000, 11))[:15_000]


#: The budget in the work itself: cells of dynamic programming, summed over every alignment
#: pass. Cells rather than passes because a pass is not a fixed cost — it is the window's
#: length times the query's — so a count of passes can stay flat while the work behind it
#: grows. Figures for the catalogue as it stands and for the 60 entries PR #10 grows it to,
#: keyed by (fixture, deleted sites). The bounds are sized for the grown one, since the
#: criterion is written for it; the seconds are one container's on the catalogue as it stands,
#: for scale only — they move by a third between runs on one container and by up to a quarter
#: between the two containers this project has run on, point by point, which is why the bounds
#: are in cells and not in them. The ranges are recorded beside `DETECTION_BUDGET`.
#:
#: Each bound sits a sixth to a fifth above the grown catalogue's figure — 17 to 19 %, the
#: spread being the rounding to a tidy number — except ("aav", 0): a
#: verbatim AAV vector costs almost nothing precisely because every element is exact, so 2 M
#: is nearly six times its cost and still an eighth of the 15.3 M it costs without the
#: shortcut that makes it cheap.
#:
#: Growing the catalogue from 36 entries to 60 moves only the stress and all-in-one figures,
#: and by a few per cent: on the stress plasmid 3.8 M cells, and all of that is the CBh
#: promoter, whose seeds land on the CAG promoter there and whose alignment is then thrown
#: away. Thrown away for want of coverage, not identity: CBh and CAG agree on 96.2 % of the
#: bases they align, which clears the 95 % the engine asks, but the alignment reaches only
#: 87.8 % of the CBh entry, short of the 90 % a promoter must cover.
#: `test_the_cbh_promoter_is_not_reported_on_a_cag_promoter` pins both figures, because an
#: engine that let this pair through would call every CAG promoter a CBh as well. The cost is
#: in the elements a sequence carries, not in how many entries the catalogue holds.
#:
#: Each shortcut the engine takes was removed on a scratch copy to see which bounds notice.
#: The one that reports an exact copy without aligning it is reached only where an element is
#: exact, so removing it fails the three `sites == 0` bounds and leaves the other nine alone;
#: removing the circular-copy skip fails every bound but `("aav", 0)`, which is too cheap for
#: it to show; and removing the merge of neighbouring guesses fails seven of the nine bounds
#: with a deleted site in them. Not `("all-in-one", 4)`, where the merge costs more rather
#: than less, so its bound is above its own per-offset cost by construction (84.1 M against
#: 103 M); and not `("stress", 1)`, at 130.0 M against 136 M. So these bounds are a
#: regression guard and not the guard for the shortcuts themselves: that is
#: `test_the_copy_of_an_element_a_circular_search_appends_is_not_aligned`,
#: `test_an_element_present_verbatim_is_found_without_aligning_it` and
#: `test_the_registers_of_one_diverged_copy_are_aligned_once`, one each.
ALIGNMENT_WORK = {
    #  fixture,    sites         36 entries         60 entries      seconds
    ("aav",        0):   2_000_000,  #    0.34 M / 15 p     0.34 M / 15 p     0.03
    ("aav",        1):  18_500_000,  #   15.49 M / 29 p    15.49 M / 29 p     0.39
    ("aav",        3):  18_500_000,  #   15.55 M / 31 p    15.55 M / 31 p     0.40
    ("aav",        5):  18_500_000,  #   15.57 M / 32 p    15.57 M / 32 p     0.41
    ("all-in-one", 0):  26_000_000,  #   22.01 M /  9 p    22.24 M / 12 p     0.48
    ("all-in-one", 1):  60_000_000,  #   50.38 M / 22 p    50.61 M / 25 p     1.19
    ("all-in-one", 3):  74_000_000,  #   62.69 M / 30 p    62.92 M / 33 p     1.46
    ("all-in-one", 4): 103_000_000,  #   87.27 M / 32 p    87.80 M / 39 p     2.07
    ("all-in-one", 5):  74_000_000,  #   62.15 M / 30 p    62.53 M / 35 p     1.53
    ("stress",     0):  52_000_000,  #   40.42 M / 39 p    44.22 M / 44 p     0.93
    ("stress",     1): 136_000_000,  #  111.75 M / 58 p   115.55 M / 63 p     2.75
    ("stress",     3): 147_000_000,  #  121.05 M / 65 p   124.85 M / 70 p     3.04
    # The S5.4 fixture exists only in the grown catalogue, so it has one column.
    ("grown",      0):   8_000_000,  #        —           6.65 M / 15 p     0.20
    ("grown",      1):  24_000_000,  #        —          20.64 M / 41 p     0.57
}


def _alignment_work(sequence: str, circular: bool = True, features=None) -> tuple[int, int]:
    """(passes, cells) of dynamic programming detection spends on ``sequence``. This is the
    measure a clock cannot give: the same numbers on any machine, where wall-clock time on a
    shared CI runner says as much about the runner as about the code."""
    aligner = features_module._aligner()
    passes = cells = 0
    original = aligner.align

    def counted(target, query, *args, **kwargs):
        nonlocal passes, cells
        passes += 1
        cells += len(target) * len(query)
        return original(target, query, *args, **kwargs)

    with mock.patch.object(aligner, "align", counted):
        detect(sequence, circular=circular, features=features)
    return passes, cells


@functools.cache
def _annotated(element_names: tuple[str, ...], sites: int,
               inverted: tuple[str, ...]) -> tuple:
    """``detect`` on a fixture, remembered: eleven pinned points and eight annotation cases ask
    for the same plasmids, and each detection is seconds. The matches come back shared, since
    ``Match`` is not frozen, so read them and do not write to them."""
    return tuple(detect(_plasmid(element_names, sites, inverted), circular=True))


@functools.cache
def _fixture_work(element_names: tuple[str, ...], sites: int,
                  inverted: tuple[str, ...]) -> tuple[int, int]:
    """``_alignment_work`` on a fixture, remembered for the same reason."""
    return _alignment_work(_plasmid(element_names, sites, inverted))


@pytest.mark.parametrize("label, element_names, inverted, sites", POINTS)
def test_the_alignment_work_of_a_15kb_plasmid_is_pinned(label, element_names, inverted, sites):
    """The work each point on the curve costs, pinned so that it cannot grow unnoticed."""
    bound = ALIGNMENT_WORK[label, sites]       # before the measurement, so a gap fails fast
    passes, cells = _fixture_work(element_names, sites, inverted)
    # An upper bound alone is satisfied by work that vanished rather than grew, and that is not
    # hypothetical: `_alignment_work` counts the aligner `_aligner()` hands out, so an engine
    # that built a fresh one per call would be counted as doing nothing at all and would pass
    # every bound in the table. Every point here aligns something.
    assert passes, f"{label}, {sites} deleted sites: no alignment was counted"
    assert cells <= bound, f"{label}, {sites} deleted sites: {cells / 1e6:.2f} M cells"


@pytest.mark.parametrize("label, element_names, inverted, sites", DIVERGED_POINTS)
def test_a_plasmid_is_annotated_the_same_however_its_elements_diverge(label, element_names,
                                                                     inverted, sites):
    """A diverged fixture is only worth measuring if it carries the same elements as the
    verbatim one, so this is the correctness half of every number above, and it is what the
    comment's claim rests on. Absolute, not relative: an engine that stopped finding things
    would make the cells go *down* and pass every bound."""
    exact = _annotated(element_names, 0, inverted)
    assert set(element_names) <= set(names(exact)), (
        f"missing from the verbatim plasmid: {sorted(set(element_names) - set(names(exact)))}"
    )
    diverged = _annotated(element_names, sites, inverted)
    # The same elements, in the same order, and as many of them: a tag inside a coding
    # sequence (SaCas9 carries 3×HA) is reported too, and must not quietly disappear.
    assert names(diverged) == names(exact)
    # None of them flagged, where a single deleted base rather than a codon would flag every
    # coding sequence and the dialog would leave them unticked.
    assert [match.frameshift for match in diverged] == [match.frameshift for match in exact]
    # Not one of them flagged — asserted outright, and not only as parity with the verbatim
    # plasmid, which an engine that flagged both would satisfy while the words above were
    # false.
    assert not any(match.frameshift for match in diverged)
    assert not any(match.frameshift for match in exact)
    losses = _losses(element_names, sites)
    for match in diverged:
        # The whole entry, not part of it: an annotation that dropped bases off an element's
        # head would otherwise satisfy the identity below by having less to be wrong about.
        assert match.part == (1, match.entry_length), f"{match.name} is partial: {match.part}"
    # Then, per copy of each element, how far below 100 % it may read and how long it may be
    # — both derived from the bases that copy lost rather than from a constant, so a shorter
    # entry in some later fixture cannot fail them for no reason. The floors are exact for
    # every element these fixtures cut — measured slack zero everywhere — because the cuts
    # are spaced far apart, so each deletion is one gap column and costs nothing else; cuts
    # clustered together would cost a mismatch as well and would need a column of slack.
    #
    # A copy is matched to its floor by *where it is*, not by sorting the two lists together:
    # the tag inside SaCas9 loses nothing while the standalone copy loses a codon, and sorted
    # pairing is satisfied by those two values in either order, so the nested copy could
    # still have degraded by three points as long as the standalone one read 100 %.
    spans = {id(match): match.span(15_000) for match in diverged}

    def host_of(match):
        """The shortest other element whose span covers this one, or None for a copy that
        stands on its own."""
        covering = [other for other in diverged
                    if other.name != match.name and spans[id(match)] <= spans[id(other)]]
        if not covering:
            return None
        return min(covering, key=lambda other: len(spans[id(other)])).name

    grouped: dict[tuple[str, str | None], list] = {}
    for match in diverged:
        grouped.setdefault((match.name, host_of(match)), []).append(match)
    for (name, host), copies in sorted(grouped.items(), key=lambda item: item[0][0]):
        floors = losses.get((name, host), (0,) * len(copies))
        assert len(floors) >= len(copies), (
            f"{name} in {host}: {len(copies)} copies, floors {floors}")
        # Within one host the copies are indistinguishable, so there pair worst with worst:
        # fewer bases lost means a higher identity and a longer span.
        copies.sort(key=lambda match: match.identity, reverse=True)
        for match, lost in zip(copies, sorted(floors)):
            floor = 100.0 * (match.entry_length - lost) / match.entry_length
            assert match.identity >= floor - 1e-9, (
                f"{name} at {match.start} reads {match.identity:.3f} %, {lost} bases deleted "
                f"of {match.entry_length}, so no lower than {floor:.3f} %"
            )
            # And the reported size moves with those bases, which is what the dialog turns
            # into "N bp missing" against the catalogue entry. No match in these fixtures
            # crosses the origin, so this does not exercise a wrapping span.
            assert len(match.span(15_000)) == match.entry_length - lost, (
                f"{name} at {match.start} spans {len(match.span(15_000))} bp, expected "
                f"{match.entry_length - lost}"
            )


@pytest.mark.parametrize("label, element_names, inverted, counts", FIXTURES)
def test_enough_deleted_sites_stop_an_element_being_found_at_all(label, element_names,
                                                                 inverted, counts):
    """The other end of the curve, and the reason the pinned points stop where they do.

    What ends an element is **seed exhaustion**, not the identity floor. The engine aligns an
    element only where one of its six seeds occurs in the plasmid, and at each of these three
    boundaries nothing is left to align: the AmpR promoter has no candidate offset at all,
    and 3×HA still has five, but every one of them comes from a seed hit inside the copy
    nested in SaCas9, 6 kb from the standalone copy and far outside the alignment window
    around an offset. Dropping
    `min_identity` to 0 brings neither of them back, so the floor is never reached and is not
    the reason. (Were they anchored they would fail it: 94.29 % for the promoter, 93.83 % for
    the tag. That is arithmetic about the deletions, not the engine's cause.)

    A seed window free of a cut is therefore sufficient but not necessary. With evenly spaced
    cuts the end windows last longest — the AmpR promoter still has its last window intact at
    five sites, with the spacing already below `SEED`, which is why "spacing below `SEED`" is
    not the rule either — and 3×HA is found at four sites with no window of its own left,
    because it is three exact 27 bp repeats and the junction in front of it supplies a hit
    from two bases outside the copy.

    The casualties are the shortest entries each fixture aligns, because a cut spacing that
    leaves no window of a 105 bp entry still leaves both end windows of a nuclease, which
    also loses only 18 of its 3345 bases to six sites (99.46 %, since a CDS is cut by whole
    codons). "Shortest" only orders them: at seven sites the AAV fixture also loses its
    truncated ITR at 131 bp.

    This is what makes the cost a sawtooth rather than a ramp, and why a bound on one point
    is not a ceiling.
    """
    sites, casualty = LOSES_AN_ELEMENT[label]
    exact = collections.Counter(names(_annotated(element_names, 0, inverted)))
    lost = collections.Counter(names(_annotated(element_names, sites, inverted)))
    # One copy of one named element, and nothing gained: a bare inequality of the two lists
    # would be satisfied by an element appearing, or by the order changing.
    assert exact - lost == collections.Counter({casualty: 1}), f"{label}: {exact} -> {lost}"
    assert lost - exact == collections.Counter()
    # And one site fewer is still the same document, so the line is exactly here, and the
    # last point the fixture pins is on the right side of it. This is the fragile assertion
    # of the file: one site before it vanishes the casualty is still anchored but sits just
    # above the identity threshold — 95.24 % for the AmpR promoter, 95.06 % for 3×HA — so
    # either a change to the seed spacing or a change to the scoring moves this line.
    before = collections.Counter(names(_annotated(element_names, sites - 1, inverted)))
    assert before == exact
    assert max(counts) < sites


def test_the_catalogue_is_parsed_once_however_many_documents_are_opened():
    """The regression no clock in this file can see, counted instead of timed. Parsing
    `common_features.json` costs 0.06 to 0.17 ms across 200 runs, and the whole cache miss —
    reading the file, parsing it and building the 36 dataclasses — 0.15 to 0.66 ms. That is
    three to five orders of magnitude under the seconds recorded beside `DETECTION_BUDGET`, and
    so far under the jitter of any timed assertion that a version parsing it per document, or
    per element, would be invisible while the user paid for it on every file they open.

    So this counts the parses rather than their cost, which is the same assertion on every
    machine: two documents, one parse. Counting `json.loads` rather than the `lru_cache`'s own
    miss counter on purpose — a second loader, or a merge with a user catalogue, would leave
    that counter at one and still reparse. Counting it on the `json` module rather than on this
    module's reference to it, so that a loader written as `json.load(fh)` is counted too, since
    that calls `loads` itself. What slips through is a loader that goes round `loads`
    altogether: a review planted a user catalogue merged in on every detection through
    `json.JSONDecoder().decode` and this test passed while reparsing on every document. That is
    the limit of the method and worth knowing before trusting it; what it would hide is another
    0.06 to 0.17 ms per file opened.

    One more thing to expect rather than debug: parsing the catalogue eagerly at import is a
    *correct* change, and it fails this test on `cache_clear()` with an `AttributeError` rather
    than on the count, because there is then no `lru_cache` to clear.

    The second document is what makes the count mean anything, since `detect` asks for the
    catalogue once per call. Both carry an element, so the per-element and the aligning paths
    are walked and not only the seed search: a review found that two random documents reach no
    alignment at all, which would have left a reparse on the aligning path uncounted.

    The aligner's cache needs nothing here. A review built a second aligner two ways — the
    cache removed, and `_match_long` constructing its own — and five existing tests fail on
    either, because they patch `align` on the object `_aligner()` returns; and
    `test_the_alignment_work_of_a_15kb_plasmid_is_pinned` now fails too rather than counting
    zero passes and passing."""
    entry_in_it = entry("EGFP").sequence
    first = clean_dna(400, 12) + entry_in_it + clean_dna(400, 13)
    changed = list(entry_in_it)
    changed[200] = "A" if changed[200] != "A" else "C"
    second = clean_dna(400, 14) + "".join(changed) + clean_dna(400, 15)
    features_module.catalogue.cache_clear()
    parses = 0
    original = features_module.json.loads

    def counted(*args, **kwargs):
        nonlocal parses
        parses += 1
        return original(*args, **kwargs)

    with mock.patch.object(features_module.json, "loads", counted):
        found = [detect(first, circular=False), detect(second, circular=False)]
    assert parses == 1, f"the catalogue was parsed {parses} times for two documents"
    # And both documents really were annotated, so that a `detect` which found nothing — and
    # therefore walked none of the paths a reparse could sit on — cannot pass this test.
    assert [names(matches) for matches in found] == [["EGFP"], ["EGFP"]], found


def _bases_in(result: object) -> int:
    """How many of the document's bases a call on it handed back."""
    if isinstance(result, (str, bytes, bytearray)):
        return len(result)
    if isinstance(result, (list, tuple)):
        return sum(_bases_in(part) for part in result)
    return 0


class _WatchesHowItIsRead(str):
    """A document that charges for what Python asks of it *through its own methods*.

    A call that goes through the instance costs at least one, and one that hands back the
    document's bases costs however many it hands back. So `document.find(seed)`, which scans the
    whole thing in C and returns a number, costs one — that is the idiom the engine is built on,
    six times per catalogue entry — while `document[i]` costs one base, `document.upper()` or
    `document.encode()` cost the document's whole length, and `document.startswith(seed, i)` at
    every position costs one per position.

    Three ways to the bases go round the instance entirely and are charged **nothing**, even
    once per element, and no `str` subclass can close them:

      * a C-level constructor or codec, which reads the characters with no method call at all —
        `bytes(document, "ascii")`, `bytearray(...)`, and anything built on them, which includes
        Biopython, since `Seq.__init__` encodes a `str` that way: `str(Seq(document))` is
        charged 0, and this module's own `_reverse_complement` is written like that;
      * an unbound call, `str.__getitem__(document, i)` or `str.upper(document)`, which finds
        the method on the type;
      * a call on some *other* object that takes the document as an argument, such as
        `"".join([document])`.

    And once a copy exists by any route, walking the copy is not charged, because the document
    is no longer involved. The test's own docstring records what all of that leaves open."""

    __slots__ = ()
    bases: ClassVar[list[int]] = [0]

    @classmethod
    def _charge(cls, amount: int) -> None:
        cls.bases[0] += amount

    def __getattribute__(self, name):
        attribute = str.__getattribute__(self, name)
        if name.startswith("_") or not callable(attribute):
            return attribute

        def charged(*args, **kwargs):
            result = attribute(*args, **kwargs)
            _WatchesHowItIsRead._charge(max(1, _bases_in(result)))
            return result

        return charged

    # The dunders syntax reaches through the type, so `__getattribute__` never sees them.
    def __getitem__(self, item):
        taken = str.__getitem__(self, item)
        _WatchesHowItIsRead._charge(len(taken))
        return taken

    def __iter__(self):
        _WatchesHowItIsRead._charge(len(self))
        return str.__iter__(self)

    def __str__(self):
        _WatchesHowItIsRead._charge(len(self))
        return str.__str__(self)

    def __add__(self, other):
        _WatchesHowItIsRead._charge(len(self))
        return str.__add__(self, other)

    def __mul__(self, other):
        _WatchesHowItIsRead._charge(len(self) * max(1, other))
        return str.__mul__(self, other)

    def __mod__(self, other):
        _WatchesHowItIsRead._charge(len(self))
        return str.__mod__(self, other)

    def __radd__(self, other):
        _WatchesHowItIsRead._charge(len(self))
        return str.__add__(other, str.__str__(self))

    def __rmul__(self, other):
        _WatchesHowItIsRead._charge(len(self) * max(1, other))
        return str.__mul__(self, other)

    def __contains__(self, other):
        _WatchesHowItIsRead._charge(1)
        return str.__contains__(self, other)


def _a_base_at_a_time(target: str, feature: str) -> set[int]:
    """`_candidate_offsets` rewritten with the per-base loop in Python, exactness preserved.

    The reference regression, kept so that the counting below can be checked against something
    rather than trusted. It returns the same offsets as the engine does."""
    offsets: set[int] = set()
    for position, seed in features_module._seeds(feature):
        width = len(seed)
        for start in range(len(target) - width + 1):
            for index in range(width):
                if target[start + index] != seed[index]:
                    break
            else:
                offsets.add(max(0, start - position))
    return offsets


def _bases_taken_from_the_document(document: str, elements, scan=None) -> int:
    """What `scan` takes out of `document` in Python, looking for `elements`.

    `scan` defaults to the engine's own `_candidate_offsets`."""
    scan = scan or features_module._candidate_offsets
    _WatchesHowItIsRead.bases[0] = 0
    watched = _WatchesHowItIsRead(document)
    for element in elements:
        scan(watched, element.sequence)
    return _WatchesHowItIsRead.bases[0]


def test_the_document_is_not_walked_in_python_once_per_element():
    """The other regression no clock can see, and this one asserted rather than recorded.

    `_candidate_offsets` scans the whole document for six seeds of every catalogue entry longer
    than `SHORT_FEATURE`, which is the largest single part of what detection does on a plasmid
    whose elements are not diverged — the common case, since an element present verbatim is
    taken whole without being aligned at all. `str.find` does that scan in C. Pulling it into
    Python costs 40x there (0.040 s to 1.60 s on the verbatim AAV fixture, medians), while
    leaving every annotation and, where the rewrite preserves exactness, every cell of
    `ALIGNMENT_WORK` untouched. On the real records the factor is widest where there is least
    alignment to dilute it, so the cheap ones suffer worst and the dearest of them, which pays
    most of its time inside the aligner, costs under three times as much. Those figures, and
    how much they move between containers, are in `_candidate_offsets`. Four rounds of timed
    assertions could not separate any of it from a slower machine.

    It can be asserted without a clock, because the property that matters is not a duration. It
    is that the document is not read out of the C layer **once per element**: that is what
    separates the regression from a correct implementation that also reads the document in
    Python, and it is why neither bound below is a constant. Measured here, over the 33 entries
    and a 15 kb document, in 10 ms (the test as a whole costs 0.21 s, most of it the
    instrument's own check below):

        as shipped, `str.find`                            267 things asked of the document
        a k-mer index built once per document         269,694 bases    correct, and quicker
        a compiled-regex IUPAC reader                       0          correct
        base by base, exactness preserved           3,962,927          the regression
        base by base over a `.upper()` copy           495,000          the regression
        base by base over `document + ""`             495,000          the regression
        base by base over `document.encode()`         495,000          the regression
        base by base over `iter(document)`            495,000          the regression
        `document.startswith(seed, i)` per position  2,966,634          the regression

    All nine find exactly the same offsets. The index is not a straw man: on this container it
    is 1.4x quicker than `str.find` on the verbatim fixture (0.029 s against 0.042 s) and level
    with it at one deleted site, so a test that forbade it would be sending a developer to C to
    undo a speed-up. Note that it reads *more* bases per pass than any of the regressions — a
    window of `SEED` bases at every position — and still costs the least, because it reads them
    once where each of the others reads the whole document again for every entry. So the bound
    is one windowed pass over the document, and what it forbids is the per-element repetition:
    at 33 entries the cheapest of the regressions is already 1.7x above it, and every extra
    entry takes it further over while the index stays put. (Which is also why
    `assert len(long_enough) > 20` above is load-bearing and not decoration: at nineteen entries
    a full read per entry would fit inside the bound.) The count is then taken again on the same
    document doubled, against a bound that doubles with it, so that a scan cannot slip under by
    the document being short.

    **What it does not catch**, measured rather than assumed, every route below verified to
    pass while costing 46x to 75x on the verbatim fixture:

      * a copy taken by any of the three routes that go round the document's own methods, even
        *per element*: a C-level constructor or codec, including `str(Seq(document))`, which is
        the idiom this module's own `_reverse_complement` uses; an unbound
        `str.__getitem__(document, i)`;
        and `"".join([document])`. `_WatchesHowItIsRead`'s docstring has the list, and no `str`
        subclass can close it. A developer adding a protein-level search writes
        `str(Seq(target).translate())` and this test will not see the walk that follows;
      * **one** full-length copy by *any* route, including the charged ones, followed by a walk
        over the copy. The copy is charged once and the walking never is, because by then the
        document is not involved.

    So the honest summary is narrower than it first looks: this forbids *asking the document
    itself* for its bases once per element. It is a guard on the shape the regression has taken
    every time it has been written here, not a proof that the per-base loop is in C. The layer
    underneath is `ALIGNMENT_WORK`, which catches the tolerant scans whatever route they take,
    and the one tolerant scan that slips past this guard is caught there. This test also calls
    `_candidate_offsets` directly, because `detect` hands on a plain `str` from
    `sequence.upper()`, so a per-base loop introduced in `_match_long` or in a new helper is
    not seen here either.

    **False positives to expect**, and they are the mirror of the second hole: a scan that
    copies the document once **per element** through a charged route and then does all its
    per-base work in C is charged 495,000 and fails, although it is correct. That covers
    `encode()` with `bytes.find`, and `translate()` with `str.find`, which is the obvious way
    to canonicalise ambiguity codes and costs 1.05x the engine. Hoisting the copy out of the
    loop into a per-document cache — which is what either wants anyway — brings it to 15,000
    and passes.

    It forbids the per-element loop and not the tolerance. A reader of ambiguity codes built on
    a compiled regex — which is what `_candidate_offsets`'s own docstring recommends — is
    charged nothing at all, because `re` takes the document into C, and it finds the same 25
    offsets, so it passes. A mismatch-tolerant scan is another matter: it finds *more* starting
    points, so `ALIGNMENT_WORK` fails on it, which is the alarm working and not a stale
    bound."""
    long_enough = [e for e in catalogue() if e.length > SHORT_FEATURE]
    assert len(long_enough) > 20, len(long_enough)
    document = _plasmid(AAV_ELEMENTS, 1, AAV_INVERTED)
    # A counter that charges nothing looks exactly like an engine that reads nothing, which is
    # the hole `assert passes` closes in `test_the_alignment_work_of_a_15kb_plasmid_is_pinned`.
    # So before a low count below is trusted, check the instrument against the regression it
    # exists to catch — and against *both* of the ways it charges, since a review found that a
    # control going only through `__getitem__` passes with the rest of the instrument blind.
    # Two entries are enough, and keep this to 0.2 s rather than the 3.7 s all 33 would cost.
    # Entries the fixture actually carries, because the equality below is vacuous on an element
    # that is absent: the two longest are SpCas9 and SaCas9, and this is the SaCas9 vector.
    a_couple = [element for element in long_enough
                if features_module._candidate_offsets(document, element.sequence)][:2]
    subscripted = _bases_taken_from_the_document(document, a_couple, _a_base_at_a_time)
    assert subscripted > 100_000, (
        f"the instrument charged a base-by-base scan over {len(a_couple)} entries only "
        f"{subscripted} bases, so a low count below would mean nothing")

    def _over_a_copy(target: str, feature: str) -> set[int]:
        """The same scan over a full-length copy, which reaches the document by a method."""
        return _a_base_at_a_time(target.upper(), feature)

    copied = _bases_taken_from_the_document(document, a_couple, _over_a_copy)
    assert copied >= 2 * len(document), (
        f"the instrument charged {copied} bases for {len(a_couple)} full-length copies of a "
        f"{len(document)} bp document, so it is not watching how the document is copied")
    # And the reference regression really is the engine's scan, so that the control above
    # cannot be satisfied by a scan that is merely expensive.
    found = [features_module._candidate_offsets(document, element.sequence)
             for element in a_couple]
    assert len(a_couple) == 2 and all(found), (a_couple, found)
    assert [_a_base_at_a_time(document, element.sequence) for element in a_couple] == found, (
        "the reference scan is not the engine's scan")
    for target in (document, document * 2):
        #: One windowed pass over it — `SEED` bases at each position, which is what building an
        #: index costs — plus one base per seed per entry.
        budget = (features_module.SEED + 1) * len(target) + 6 * len(long_enough)
        taken = _bases_taken_from_the_document(target, long_enough)
        assert taken <= budget, (
            f"{taken} bases of a {len(target)} bp document were taken out in Python for "
            f"{len(long_enough)} entries, against {budget}: the per-base work belongs in C, or "
            f"at least in a single pass over the document rather than one pass per element")


def test_a_second_deleted_site_per_element_is_nearly_free():
    """What the merge bought, asserted rather than described: the cost barely depends on how
    many times an element has been cut, because the registers of one copy are one alignment
    however many there are. Five deleted sites per element cost a quarter more than one does,
    where they used to cost three times as much (90.1 M against 30.4 M on the AAV fixture,
    197.1 M against 78.6 M on the all-in-one one).

    A quarter rather than nothing, and the quarter is what checking the merge costs instead
    of trusting it. Taking the all-in-one fixture at **three** deleted sites: 49.8 M cells if
    a run's one alignment is trusted, 62.7 M as the engine actually does it. Every one of
    those re-queues is a run whose alignment came back *below the thresholds*, which answers
    no guess in it and so sends all of them back — tdTomato, a tandem dimer, for two runs and
    four extra passes of a 1431 bp entry (9.8 M cells), and mCherry for two runs and five of a
    711 bp one (3.0 M). The other branch, an alignment that drifts onto a placement none of
    its guesses is at, fires nowhere here; the one place it fires on any fixture is a 3×HA run
    on the packed stress plasmid at three deleted sites.

    The four-site trough is left out on purpose: there the fixture's cuts land in the seeds
    themselves, the surviving registers end up further apart than `MERGE_GAP`, and SpCas9 is
    aligned twice — its two guesses are 12 bases apart, because all four cuts fall ahead of
    the last surviving seed and each of them takes a codon out. (SaCas9 is aligned once here;
    on the AAV fixture it is SaCas9 the same arithmetic splits.) That is a property of evenly
    spaced cuts, not of divergence — it is one of the two points the tests pin where the
    engine does more work than aligning every guess separately, 87.3 M cells against 84.1 M
    here and 40.4 M against 39.4 M on the packed stress plasmid with nothing deleted at all,
    and the point furthest over S5.4's two seconds on a plausible vector: 1.95 to 2.07 s across
    runs on one container, against 1.86 to 1.96 s per guess, and 2.53 to 3.10 s on a slower
    one, where the 1-codon and 5-codon points only touch the budget rather than clear it. So
    this is the one point where the merge neither saves work nor leaves room against the
    budget, which is one reason no test asserts the seconds. `ALIGNMENT_WORK` records the point
    so it cannot drift unseen, and `MERGE_GAP`'s own comment records what widening the gap
    would cost and buy there.
    """
    for element_names, inverted, slack in ((AAV_ELEMENTS, AAV_INVERTED, 1.01),
                                           (ALL_IN_ONE_ELEMENTS, (), 1.30)):
        _, one = _fixture_work(element_names, 1, inverted)
        for sites in (3, 5):
            _, several = _fixture_work(element_names, sites, inverted)
            assert several <= slack * one, (
                f"{sites} deleted sites cost {several / 1e6:.2f} M cells against "
                f"{one / 1e6:.2f} M for one")


def test_a_plasmid_with_nothing_in_the_catalogue_is_ruled_out_without_aligning():
    """The guard for the catalogue's growth: a sequence holding none of its elements is
    dismissed on the seed search alone, so entries can be added without slowing every
    document down."""
    sequence = clean_dna(15_000, 12)
    assert detect(sequence, circular=True) == []
    assert _alignment_work(sequence) == (0, 0)


def test_the_copy_of_an_element_a_circular_search_appends_is_not_aligned():
    """A circular search runs against the sequence doubled, so every element is found twice
    and the copy past the end is dropped as a duplicate. It must not be aligned first: that
    was half the work detection used to do."""
    feature = max(catalogue(), key=lambda f: f.length)
    # A copy with a few changed bases, so that it has to be aligned: an exact copy is taken as
    # it stands and would hide the difference this test is looking for.
    bases = list(feature.sequence)
    for position in range(40, len(bases), 500):
        bases[position] = "A" if bases[position] != "A" else "C"
    sequence = clean_dna(300, 15) + "".join(bases) + clean_dna(300, 16)
    linear, _ = _alignment_work(sequence, circular=False, features=(feature,))
    circular, _ = _alignment_work(sequence, circular=True, features=(feature,))
    assert linear > 0, "the mutated copy should need an alignment"
    # Passes, not cells: a circular search looks at a doubled target, so the window around the
    # genuine copy reaches further and costs a little more. What must not appear is a second
    # pass for the appended copy.
    assert circular == linear


def test_an_element_present_verbatim_is_found_without_aligning_it():
    """An exact copy is its own best local alignment, so it is reported without a dynamic
    programming pass. The test pins the shortcut by failing if the aligner is used at all."""
    feature = max(catalogue(), key=lambda f: f.length)
    sequence = clean_dna(400, 13) + feature.sequence + clean_dna(400, 14)
    [found] = detect(sequence, circular=False, features=(feature,))
    assert found.identity == 100.0
    assert found.start == 400
    assert found.end == 400 + feature.length - 1
    # Only this element is searched for, so any pass at all would be this one being aligned.
    assert _alignment_work(sequence, circular=False, features=(feature,)) == (0, 0)


def test_the_cbh_promoter_is_not_reported_on_a_cag_promoter():
    """Why the growth above costs what it costs, and why the cost is worth paying. CBh is a
    CAG made shorter, so its seeds hit every CAG and the engine aligns the two in full. What
    keeps it from reporting a CBh on every CAG is coverage, not identity: the two agree on
    almost every base they align, and the alignment simply does not reach enough of the CBh
    entry. Both figures are pinned, because a change that loosened either one would start
    reporting an element the construct does not carry.
    """
    from vettore.features import COVERAGE_BY_TYPE, DEFAULT_MIN_COVERAGE

    cbh, cag = entry("CBh promoter"), entry("CAG promoter")
    identity, coverage = shared_span(cbh, cag.sequence)
    needed = COVERAGE_BY_TYPE.get("promoter", DEFAULT_MIN_COVERAGE)
    assert identity == pytest.approx(96.2, abs=0.1), identity
    assert identity >= DEFAULT_MIN_IDENTITY, "identity alone would let the pair through"
    assert coverage == pytest.approx(87.8, abs=0.1), coverage
    assert coverage < needed, f"{coverage:.1f} % covered, {needed} % asked"
    assert "CBh promoter" not in names(detect(cag.sequence, circular=False))


def test_both_copies_of_a_tandem_tag_are_reported():
    """A 6×HA tag is two 3×HA in a row, and 3×HA is itself three identical 27 bp repeats, so
    the seeds hit in every register and the guesses at where the tag starts land closer
    together than the search margin. Both copies must still be reported: an element that is
    there verbatim and goes unreported is the one failure the engine cannot afford, and
    grouping those guesses as one place to look caused exactly that.
    """
    tag = entry("3×HA")
    head = clean_dna(414, 17)
    sequence = head + tag.sequence + tag.sequence + clean_dna(414, 18)
    starts = {m.start for m in detect(sequence, circular=False) if m.name == tag.name}
    assert len(head) in starts, "the first tag is missing"
    assert len(head) + tag.length in starts, "the second tag is missing"
    for match in (m for m in detect(sequence, circular=False) if m.name == tag.name):
        assert match.identity == 100.0
        assert match.end - match.start + 1 == tag.length


def test_the_registers_of_one_diverged_copy_are_aligned_once():
    """The merge, pinned where it is cheapest to see: one copy of one element, alone in a
    flank, with a base taken out of it once, twice and three times. Each deletion shifts every
    seed after it, so the engine guesses at two, three and two starts for the same copy, one
    to three bases apart (offsets 59-60, 58-60 and 57/60 as measured: the third cut falls
    inside a seed window and destroys that seed, which is why the third case has two guesses
    rather than four and why they are three apart). All of them are one alignment, and the
    match that comes out of it still covers the whole entry at the copy's own start."""
    feature = entry("AmpR promoter")
    for sites in (1, 2, 3):
        sequence = clean_dna(60, 21) + _diverged(feature, sites) + clean_dna(60, 22)
        passes, _ = _alignment_work(sequence, circular=False, features=(feature,))
        assert passes == 1, f"{sites} deleted sites took {passes} alignments"
        [found] = detect(sequence, circular=False, features=(feature,))
        assert found.start == 60
        assert found.part == (1, feature.length)
        assert len(found.span(len(sequence))) == feature.length - sites


def test_two_diverged_copies_of_one_element_stay_two_matches():
    """And the merge must not join two copies, which is how a wider version of this idea
    failed a review on PR #5: it lost the second half of a 6×HA tag. Here both 3×HA copies
    have lost a base, in different places, so neither is reported through the exact-copy
    shortcut and both go through the alignment the merge decides the shape of."""
    tag = entry("3×HA")
    first = _diverged(tag, 1)
    second = tag.sequence[:20] + tag.sequence[21:]
    sequence = clean_dna(60, 23) + first + second + clean_dna(60, 24)
    starts = {match.start for match in detect(sequence, circular=False, features=(tag,))}
    assert 60 in starts, f"the first copy is missing: {sorted(starts)}"
    assert 60 + len(first) in starts, f"the second copy is missing: {sorted(starts)}"


def test_only_guesses_that_can_be_one_copy_are_merged():
    """The rule itself, on numbers rather than sequence, because what it must never do is join
    two copies. Neighbours within the gap are one run; a wider step starts a new one; and a run
    stops growing at the span, so no window stretches over more than a seed's worth of
    uncertainty however many registers line up. Whether those two figures are safe for a given
    entry is the next test's business, not this one's."""
    assert features_module._offset_runs([]) == []
    assert features_module._offset_runs([10]) == [(10, 10)]
    assert features_module._offset_runs([10, 11, 13]) == [(10, 13)]
    assert features_module._offset_runs([10, 15]) == [(10, 10), (15, 15)]
    assert features_module._offset_runs([10, 91, 92]) == [(10, 10), (91, 92)]
    # Steps of three, a lost codon each: the run stops at MERGE_SPAN and the next one begins.
    runs = features_module._offset_runs(list(range(10, 40, 3)))
    assert runs == [(10, 28), (31, 37)]
    assert all(last - first <= features_module.MERGE_SPAN for first, last in runs)


def test_how_far_apart_two_guesses_may_be_comes_from_the_entry_itself():
    """What makes the merge safe is not a constant but the entry. A reported copy of an entry
    `L` bases long agrees with the entry at all but `slack` bases — it may leave
    `100 - min_coverage` per cent of them unaligned and differ at `100 - min_identity` per
    cent of the rest — so two copies `shift` apart force the entry to be
    `(L - shift - 2 * slack) / (L - shift)` the same as itself at that shift. Below that, two
    guesses that close cannot be two copies."""
    assert features_module._self_similarity("CA" * 30, 1) == 0.0
    assert features_module._self_similarity("CA" * 30, 2) == 100.0
    # A (CA)30 entry repeats itself two bases over, so nothing beyond one base may be merged.
    assert features_module._merge_limit("CA" * 30, DEFAULT_MIN_IDENTITY, 100.0) == 1
    # A (GGGGS)6 linker repeats itself fifteen bases over, which is beyond what a run spans.
    assert features_module._merge_limit(
        "GGAGGAGGAGGATCA" * 6, DEFAULT_MIN_IDENTITY, 100.0) == 14
    real = entry("EGFP").sequence
    assert features_module._merge_limit(
        real, DEFAULT_MIN_IDENTITY, DEFAULT_MIN_COVERAGE) == features_module.MERGE_SPAN


def test_the_coverage_floor_is_part_of_how_far_apart_two_guesses_may_be():
    """The term that review round 9 found missing, and the one that bites. A `repeat_region`
    entry only has to be 70 % covered, so a copy may ignore nearly a third of the entry and
    still be reported at 100 % identity — which an identity-only bound does not see at all. A
    (CAG)16 polyQ tract with unique flanks is exactly that shape: it is reported on its tract
    alone, the tract repeats itself every three bases, and the identity-only bound merged six
    of its nine copies away."""
    tract = "CAGATTTTCA" + "CAG" * 16 + "AAAGCGGCAC"
    # It is not self-similar enough for an identity-only bound to narrow anything...
    assert features_module._merge_limit(tract, DEFAULT_MIN_IDENTITY, 100.0) \
        == features_module.MERGE_SPAN
    # ...but at the coverage a repeat_region is held to, almost nothing may be merged.
    coverage = features_module.COVERAGE_BY_TYPE["repeat_region"]
    assert features_module._merge_limit(tract, DEFAULT_MIN_IDENTITY, coverage) <= 2


def _records(matches) -> list[tuple]:
    return [(m.name, m.start, m.end, m.strand, round(m.identity, 6), m.part, m.frameshift)
            for m in matches]


def _with_and_without_the_merge(sequence: str, feature, **kwargs) -> tuple[list, list]:
    """What the engine reports, and what it reports with every guess aligned on its own. The
    merge has to reproduce the second exactly, in order as well as in content: `_dedupe` keeps
    whichever of two overlapping copies it meets first, so a reordering is an annotation
    change."""
    merged = _records(detect(sequence, circular=False, features=(feature,), **kwargs))
    with mock.patch.object(features_module, "_offset_runs",
                           lambda offsets, *_, **__: [(o, o) for o in sorted(offsets)]):
        alone = _records(detect(sequence, circular=False, features=(feature,), **kwargs))
    return merged, alone


def test_a_partly_covered_repeat_entry_keeps_every_copy():
    """The same tract end to end: nine overlapping copies in a longer tract, every one of them
    reported, on bases 11-58 of 68 — the tract alone, just over the 70 % floor — and every one
    at 100 % identity. An identity-only limit reported three; the coverage term and the
    diagonal check each recover them, which is why this asserts against the engine's own
    per-guess output rather than against either guard."""
    tract = features_module.CatalogueFeature(
        name="(CAG)16 polyQ tract", type="repeat_region",
        sequence="CAGATTTTCA" + "CAG" * 16 + "AAAGCGGCAC", note="", source="test")
    sequence = clean_dna(60, 32) + "CAG" * 36 + clean_dna(60, 33)
    found = detect(sequence, circular=False, features=(tract,))
    starts = sorted(m.start for m in found)
    assert starts == list(range(60, 60 + 6 * 9, 6)), starts
    merged, alone = _with_and_without_the_merge(sequence, tract)
    assert merged == alone
    [first] = [m for m in found if m.start == 60]
    assert first.identity == 100.0 and first.part == (11, 58)
    assert 100.0 * 48 / tract.length < features_module.DEFAULT_MIN_COVERAGE
    assert 100.0 * 48 / tract.length >= features_module.COVERAGE_BY_TYPE["repeat_region"]


def test_which_catalogue_entries_the_limit_narrows():
    """The catalogue's own self-similarity, pinned at every threshold the dialog offers, so
    that an entry added later cannot quietly change what the merge does. Three shipped entries
    are narrowed and they are the right three: both AAV2 ITRs, which are hairpins and so are
    genuinely self-similar, and 3×FLAG, which is three copies of one 24 bp tag. Everything
    else merges over the full span, which is why no measured figure moved when the coverage
    floor went into the bound.

    What this table is not is the thing that keeps a repeating entry safe. An entry whose
    internal period is longer than `MERGE_SPAN` — 3×HA's is 27 — is never narrowed and so is
    absent here, and that is the shape review round 10 lost a copy of. What answers that is
    `_match_long` checking which guesses its alignment passes through; see
    `test_a_guess_the_merged_alignment_does_not_pass_through_is_aligned_on_its_own`."""
    narrowed = {}
    for min_identity in (90.0, 95.0, 98.0, 100.0):
        for feature in catalogue():
            if feature.length <= SHORT_FEATURE:
                continue
            coverage = features_module.COVERAGE_BY_TYPE.get(
                feature.type, features_module.DEFAULT_MIN_COVERAGE)
            limit = features_module._merge_limit(feature.sequence, min_identity, coverage)
            if limit != features_module.MERGE_SPAN:
                narrowed[(min_identity, feature.name)] = limit
    assert narrowed == {
        (90.0, "AAV2 ITR"): 0, (90.0, "AAV2 ITR (truncated)"): 1, (90.0, "3×FLAG"): 8,
        (95.0, "AAV2 ITR"): 1, (95.0, "AAV2 ITR (truncated)"): 1, (95.0, "3×FLAG"): 14,
        (98.0, "AAV2 ITR"): 1, (98.0, "AAV2 ITR (truncated)"): 3,
        (100.0, "AAV2 ITR"): 3, (100.0, "AAV2 ITR (truncated)"): 3,
    }, narrowed


def test_the_limit_narrows_a_run_however_close_its_steps_are():
    """The derived limit is handed to `_offset_runs` as the span, and the span alone is what
    it narrows: capping the gap by it as well would change nothing, because a step wider than
    the limit already puts the offset further than the limit from the run's first guess.
    (Checked exhaustively: over every subset of 0..25 of size up to six and every limit up to
    the span, `min(MERGE_GAP, limit)` and `MERGE_GAP` give the same runs, without exception.)"""
    close = [10, 11, 12, 13, 14]
    apart = [10, 14, 18, 22, 26]
    assert features_module._offset_runs(close, 4, 18) == [(10, 14)]
    assert features_module._offset_runs(close, 4, 2) == [(10, 12), (13, 14)]
    assert features_module._offset_runs(close, 4, 0) == [(o, o) for o in close]
    assert features_module._offset_runs(apart, 4, 18) == [(10, 26)]
    assert features_module._offset_runs(apart, 4, 2) == [(o, o) for o in apart]
    # Sorted rather than trusted: an unsorted list must not build a backwards run.
    assert features_module._offset_runs([20, 10, 11]) == [(10, 11), (20, 20)]


def test_the_derived_limit_is_what_the_engine_merges_by():
    """`_merge_limit` is handed to `_offset_runs` as the span, and this is the only test that
    would notice if that stopped happening: the diagonal check recovers the copies either way,
    so correctness alone does not pin it, but the work does. A (CA)30 entry has a limit of one,
    so the guesses two bases apart through a long CA array are never joined; merge them and the
    engine pays for the run as well as for the guesses it then puts back."""
    repeat = features_module.CatalogueFeature(
        name="(CA)30", type="misc_feature", sequence="CA" * 30, note="", source="test")
    assert features_module._merge_limit(
        repeat.sequence, DEFAULT_MIN_IDENTITY, DEFAULT_MIN_COVERAGE) == 1
    array = "CA" * 90
    middle = len(array) // 2
    sequence = clean_dna(60, 32) + array[:middle] + "G" + array[middle + 1:] + clean_dna(60, 33)
    passes, _ = _alignment_work(sequence, circular=False, features=(repeat,))
    with mock.patch.object(features_module, "_merge_limit",
                           lambda *_, **__: features_module.MERGE_SPAN):
        ignored, _ = _alignment_work(sequence, circular=False, features=(repeat,))
    assert passes < ignored, (passes, ignored)


def _passes(sequence: str, feature, **kwargs) -> tuple[int, int]:
    """Alignments the engine does on ``sequence``, merged and with every guess on its own."""
    aligner = features_module._aligner()
    original, counted = aligner.align, [0]

    def count(target, query, *args, **kwargs_):
        counted[0] += 1
        return original(target, query, *args, **kwargs_)

    out = []
    for runs in (None, lambda offsets, *_, **__: [(o, o) for o in sorted(offsets)]):
        counted[0] = 0
        with mock.patch.object(aligner, "align", count):
            if runs is None:
                detect(sequence, circular=False, features=(feature,), **kwargs)
            else:
                with mock.patch.object(features_module, "_offset_runs", runs):
                    detect(sequence, circular=False, features=(feature,), **kwargs)
        out.append(counted[0])
    return out[0], out[1]


def test_a_run_whose_alignment_is_rejected_answers_none_of_its_guesses():
    """The other half of the same idea, and the one no construct distinguishes. When a run's
    one alignment comes back *below* the identity or coverage threshold it has answered
    nothing, however many of the run's guesses lie on it, so they all go back. Identity and
    coverage are not monotone in alignment score, so a window that can only score higher can
    still come back with nothing where a narrower one would have reported a match.

    A review looked for a construct where that costs an annotation and found none in 31,680
    checks, so what this test can pin is the work: at a threshold the merged alignment cannot
    meet, the engine does one pass per guess *plus* the run's own wasted one. Treat the run as
    having answered the guesses on it and that becomes one pass fewer than aligning them
    separately, not one more.

    What it does not pin is this half of the condition on its own. On this construct the run's
    alignment also drifts clear of both its guesses, so the diagonal check alone would send
    them back, and deleting `hit is None` from the condition leaves the whole suite green. A
    review then hunted for a construct where the rejection is the only thing holding — a run
    rejected, its diagonals covering a guess, that guess accepted when aligned alone — and
    found none in 5262 runs of more than one guess. So the clause is belt and braces, and what
    this test asserts is the work it costs, not a loss it is known to prevent."""
    tag = entry("3×HA")
    body = tag.sequence * 2
    sequence = clean_dna(60, 99) + body[:65] + "AAAA" + body[65:] + clean_dna(60, 99)
    merged, alone = _passes(sequence, tag, min_identity=100.0)
    assert merged == alone + 1, (merged, alone)
    # And nothing is lost by it: the same matches, in the same order.
    with_merge, without = _with_and_without_the_merge(sequence, tag, min_identity=100.0)
    assert with_merge == without


def test_a_guess_the_merged_alignment_does_not_pass_through_is_aligned_on_its_own():
    """The merge's one real hazard, and the check that answers it. One window yields one
    alignment, and the best alignment over a run's window is not always the best alignment
    over the window of each guess in it: a wider window can hold a better-scoring placement
    elsewhere and drift to it. A 6×HA tag — the two shipped 3×HA entries end to end — with
    four bases inserted in its second repeat does exactly that. The insertion puts two guesses
    four bases apart, the run's window reaches back over the first block and reports that
    instead, and the block at 87 is gone.

    No bound on the entry can see this, because it is a property of the target rather than of
    the entry: 3×HA is at most 44.9 % the same as itself over the shifts the merge uses — that
    at twelve bases, and under 31 % at fourteen of the eighteen — so its limit is the full
    span. So the engine checks instead — the guesses the accepted alignment passes
    through, one per register it reads, are the ones it stands in for, and any other guess in
    the run is aligned on its own."""
    tag = entry("3×HA")
    body = tag.sequence * 2
    sequence = clean_dna(60, 99) + body[:65] + "AAAA" + body[65:] + clean_dna(60, 99)
    limit = features_module._merge_limit(
        tag.sequence, DEFAULT_MIN_IDENTITY, DEFAULT_MIN_COVERAGE)
    assert limit == features_module.MERGE_SPAN, "the bound sees nothing here, by design"
    offsets = sorted(features_module._candidate_offsets(sequence, tag.sequence))
    assert (87, 91) in features_module._offset_runs(offsets, features_module.MERGE_GAP, limit)
    # That run's own window drifts: its best alignment reads the block at 60, passing through
    # neither 87 nor 91, so neither guess is accounted for and both are aligned again.
    margin = max(30, tag.length // 10)
    columns = features_module._columns(
        features_module._aligner().align(
            sequence[87 - margin : 91 + tag.length + margin], tag.sequence)[0], 87 - margin)
    diagonals = {t - q for t, q, _ in columns if t is not None and q is not None}
    assert diagonals.isdisjoint({87, 91}), sorted(diagonals)
    # All three blocks are reported even so, which is the point of the check.
    found = detect(sequence, circular=False, features=(tag,))
    assert sorted(m.start for m in found) == [60, 87, 145]
    merged, alone = _with_and_without_the_merge(sequence, tag)
    assert merged == alone
    # And the drifted answer is dropped rather than reported from that run, so the engine
    # hands `_dedupe` exactly the hits the guesses give one at a time, no extra copy of the
    # block the run drifted onto.
    def raw():
        return sorted(features_module._match_long(
            sequence, tag.sequence, len(sequence), DEFAULT_MIN_IDENTITY, DEFAULT_MIN_COVERAGE))
    with mock.patch.object(features_module, "_offset_runs",
                           lambda offsets, *_, **__: [(o, o) for o in sorted(offsets)]):
        one_at_a_time = raw()
    assert raw() == one_at_a_time


def test_a_repeating_entry_keeps_every_copy():
    """A (CA)30 entry inside a longer, imperfect CA array: every register of the array is a
    copy of the entry, and each must be reported. One base in the middle of the array differs,
    so the copies that cover it are aligned rather than found verbatim, and a fixed gap of four
    joins two of their registers into one run. What keeps the copies is the diagonal check
    rather than the limit — the limit stops them being merged in the first place, and the test
    below pins that separately — and either way the engine must give exactly what aligning
    every guess on its own gives."""
    repeat = features_module.CatalogueFeature(
        name="(CA)30", type="misc_feature", sequence="CA" * 30, note="", source="test")
    array = "CA" * 90
    at = len(array) // 2
    sequence = clean_dna(60, 31) + array[:at] + "G" + array[at + 1:] + clean_dna(60, 32)
    found = detect(sequence, circular=False, features=(repeat,))
    starts = sorted(match.start for match in found)
    assert starts == [60, 68, 76, 84, 152, 160, 168, 176, 188], starts
    merged, alone = _with_and_without_the_merge(sequence, repeat)
    assert merged == alone
    # The copy that straddles the changed base is the one a fixed gap loses.
    [straddling] = [match for match in found if match.start == 188]
    assert straddling.part == (1, repeat.length)
    assert straddling.identity == pytest.approx(96.67, abs=0.01)


def test_guesses_further_apart_than_the_gap_are_aligned_separately():
    """The narrow half of the rule, measured on the case that shows it: a nuclease cut seven
    times loses the seeds that would have bridged its registers, the ones left are six bases
    apart, and it is aligned more than once. Two sites keep every register a codon apart and
    cost one alignment. This is what the four-site trough on the AAV fixture is."""
    feature = entry("SaCas9")
    for sites, alignments in ((2, 1), (5, 1), (7, 3)):
        sequence = clean_dna(60, 25) + _diverged(feature, sites) + clean_dna(60, 26)
        passes, _ = _alignment_work(sequence, circular=False, features=(feature,))
        assert passes == alignments, f"{sites} deleted sites took {passes} alignments"
        assert detect(sequence, circular=False, features=(feature,))


def test_a_tandem_array_reports_a_copy_for_every_repeat():
    """The same thing as the array grows: however many copies are there, each is found."""
    tag = entry("3×HA")
    for copies in (2, 4, 8):
        sequence = clean_dna(414, 19) + tag.sequence * copies + clean_dna(414, 20)
        starts = {m.start for m in detect(sequence, circular=False) if m.name == tag.name}
        assert {414 + i * tag.length for i in range(copies)} <= starts, copies
