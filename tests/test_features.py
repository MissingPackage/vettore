import collections
import functools
import json
import time
from pathlib import Path
from unittest import mock

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from cloning_helpers import clean_dna, rc

from vettore import features as features_module
from vettore.features import (
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


def test_every_catalogue_entry_cites_a_record_and_a_range():
    for feature in catalogue():
        accession, _, span = feature.source.partition(":")
        start, _, end = span.partition("..")
        assert accession and start.isdigit() and end.isdigit(), feature.name
        assert int(end) - int(start) + 1 == feature.length, feature.name
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
}


def entry(name: str):
    [feature] = [f for f in catalogue() if f.name == name]
    return feature


def record_sequence(accession: str) -> str:
    record = SeqIO.read(DATA / f"features-{accession}.fasta", "fasta")
    return str(record.seq).upper()


@pytest.mark.parametrize("element", sorted(CROSS_CHECKS))
def test_a_catalogue_promoter_is_its_record_at_the_coordinates_it_cites(element):
    """The entry's own record is here too, so the claim its ``source`` makes — these are the
    record's bases at these coordinates — is checked without the network, and the element is
    found there whole and exact."""
    feature = entry(element)
    accession, _, span = feature.source.partition(":")
    start, _, end = span.partition("..")
    record = record_sequence(accession)
    assert record[int(start) - 1 : int(end)] == feature.sequence, element
    matches = [m for m in detect(record, False) if m.name == element]
    assert len(matches) == 1, element
    assert matches[0].identity == 100.0, element
    assert matches[0].start == int(start) - 1, element
    assert matches[0].end == int(end) - 1, element


@pytest.mark.parametrize(("element", "accession"), sorted(CROSS_CHECKS.items()))
def test_a_catalogue_promoter_is_found_in_an_independent_record(element, accession):
    feature = entry(element)
    assert not feature.source.startswith(accession), element
    sequence = record_sequence(accession)
    matches = [m for m in detect(sequence, False) if m.name == element]
    assert len(matches) == 1, element
    assert matches[0].identity >= 98.0, element
    found = matches[0].end - matches[0].start + 1
    assert found >= 0.9 * feature.length, element


def test_the_new_promoters_are_all_found_side_by_side_in_one_construct():
    """Seven promoters in one sequence: each is reported once, whole, at its own position and
    at full identity, so none of them swallows or hides another."""
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


# ---- detection speed --------------------------------------------------------------------
# S5.4 asks that detection stay under two seconds for a 15 kb plasmid with the grown
# catalogue. The cost is not in the catalogue's size — a feature absent from the sequence is
# ruled out by an exact seed search and costs almost nothing — but in the elements that *are*
# there, each of which is aligned where its seeds land.
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
# What the cost tracks is the number of *surviving seed registers*, and with six seeds per
# entry the ceiling is six. That is not the same as the number of deleted sites, and the two
# come apart in a way worth knowing, because the fixtures space their cuts evenly and the
# engine spaces its seeds evenly too. On the AAV fixture (36 entries, this machine):
#
#   sites     0      1      2      3      4      5      6
#   cells   0.34 M 30.4 M 45.3 M 60.3 M 30.2 M 90.1 M 89.4 M
#   seconds  0.04   0.85   1.24   1.67   0.83   2.56   2.37
#   same?    yes    yes    yes    yes    yes    yes    NO
#
# A single changed base, for comparison, is 15.3 M. So it is a sawtooth, not a ramp. At four
# sites the even spacing drops every cut inside a seed window — the engine places its seeds at
# `round(i * (length - 18) / 5)` and the fixture its cuts at `length // (sites + 1) * k`, which
# for four cuts coincide — so four of the six seeds are destroyed and the cost halves. Five
# evenly spaced sites is the ceiling: 90.1 M here, 197.1 M on the all-in-one plasmid, both
# still annotated exactly as their verbatim selves. At six an element is no longer found at
# all, and the cost does *not* collapse with it (89.4 M): only the four-site trough is cheap,
# and "divergence heavy enough to destroy a seed is cheaper" is true of a destroyed seed, not
# of heavier divergence in general. The bounds below pin that ceiling where the plasmid still
# holds together, which is every fixture but the stress one, and the six-site edge has a test
# of its own.
#
# The third is to write a divergence the engine itself reports as a different document. A
# **single** deleted base inside a coding sequence is a frameshift, and detection says so:
# `Match.frameshift` is set, the dialog leaves the match unticked and shows an amber warning. One
# deleted base per element flags SpCas9, SaCas9, tdTomato, mCherry, EGFP and AmpR, so a
# fixture built that way is not the same plasmid to the user, and a nuclease with a frameshift
# is not a vector anyone builds. Deleting a whole **codon** costs exactly the same (30.4 M on
# the AAV fixture either way) and is a real variant, so that is what `_plasmid` does to a CDS.
#
# What the plasmid holds matters as much as how far it diverges, and a fixture that packs
# every awkward entry together measures a vector nobody builds. The three below are a single
# AAV nuclease vector, a dense all-in-one expression plasmid, and that packed plasmid, kept as
# a stress case and named as one.

DETECTION_BUDGET = 2.0

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
#: This is the plausible heavy case, and it is the one that misses S5.4's budget.
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

#: Each fixture with the points on the cost curve the tests pin for it: none, one deleted site
#: per element, three, and five where five is still annotated as the verbatim plasmid is. The
#: stress plasmid is not pinned at five because it loses an element there — see
#: `LOSES_AN_ELEMENT`, which is where that boundary is asserted.
FIXTURES = (
    ("aav", AAV_ELEMENTS, AAV_INVERTED, (0, 1, 3, 5)),
    ("all-in-one", ALL_IN_ONE_ELEMENTS, (), (0, 1, 3, 5)),
    ("stress", PACKED_ELEMENTS, (), (0, 1, 3)),
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
}

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
#: criterion is written for it; the seconds are this machine's on the catalogue as it stands,
#: for scale only, and are the same measurements the table above records.
#:
#: Each bound sits about 17 % above the grown catalogue's figure, except ("aav", 0): a
#: verbatim AAV vector costs almost nothing precisely because every element is exact, so 2 M
#: is nearly six times its cost and still only a seventh of the 15.3 M it costs without the
#: shortcut that makes it cheap.
#:
#: Each shortcut the engine takes was removed on a scratch copy to see which bounds notice.
#: The one that reports an exact copy without aligning it is reached only where an element is
#: exact, so removing it fails the three `sites == 0` bounds and leaves the other eight alone;
#: removing the circular-copy skip fails every bound but `("aav", 0)`, which is too cheap for
#: it to show. So these bounds are a regression guard and not the guard for the shortcuts
#: themselves: that is `test_the_copy_of_an_element_a_circular_search_appends_is_not_aligned`
#: and `test_an_element_present_verbatim_is_found_without_aligning_it`, one each.
ALIGNMENT_WORK = {
    #  fixture,    sites         36 entries         60 entries      seconds
    ("aav",        0):   2_000_000,  #    0.34 M / 15 p     0.34 M / 15 p     0.04
    ("aav",        1):  36_000_000,  #   30.38 M / 35 p    30.38 M / 35 p     0.85
    ("aav",        3):  70_000_000,  #   60.28 M / 50 p    60.28 M / 50 p     1.67
    ("aav",        5): 105_000_000,  #   90.09 M / 63 p    90.09 M / 63 p     2.56
    ("all-in-one", 0):  26_000_000,  #   22.01 M /  9 p    22.24 M / 12 p     0.7
    ("all-in-one", 1):  92_000_000,  #   78.56 M / 33 p    78.79 M / 36 p     2.3
    ("all-in-one", 3): 165_000_000,  #  141.22 M / 59 p   141.45 M / 62 p     4.0
    ("all-in-one", 5): 230_000_000,  #  197.14 M / 80 p   197.44 M / 84 p     5.5
    ("stress",     0):  50_000_000,  #   39.36 M / 38 p    42.40 M / 42 p     1.2
    ("stress",     1): 155_000_000,  #  129.98 M / 64 p   133.02 M / 68 p     3.8
    ("stress",     3): 265_000_000,  #  222.89 M / 84 p   225.93 M / 88 p     6.3
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
    """The work each point on the curve costs, pinned so that it cannot grow unnoticed — and,
    for the all-in-one plasmid, so that whoever closes the budget gap sees the number move."""
    bound = ALIGNMENT_WORK[label, sites]       # before the measurement, so a gap fails fast
    _, cells = _fixture_work(element_names, sites, inverted)
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
    and 3×HA still has five, but every one of them is inside the copy nested in SaCas9, 6 kb
    from the standalone copy and far outside the alignment window around an offset. Dropping
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


def test_the_budget_is_met_on_a_single_aav_vector_with_a_codon_lost_per_element():
    """S5.4's two seconds as they are written, in seconds, on the vector the criterion is for.
    The cells tests are what fail on a regression; this one catches a cost that is not
    alignment, so it also times the heaviest element set that has a point inside the budget.
    The margins against the two seconds are 2.4x for the AAV point (0.85 s) and 1.7x for the
    stress one (1.2 s), so a shared CI runner would have to be about two thirds slower than
    this machine for either to be a false failure."""
    for element_names, sites, inverted in ((AAV_ELEMENTS, 1, AAV_INVERTED),
                                           (PACKED_ELEMENTS, 0, ())):
        sequence = _plasmid(element_names, sites, inverted)
        detect(sequence, circular=True)      # warm the catalogue and the aligner
        start = time.perf_counter()
        matches = detect(sequence, circular=True)
        elapsed = time.perf_counter() - start
        assert matches, "the plasmid should carry recognisable elements"
        assert elapsed < DETECTION_BUDGET, f"detection took {elapsed:.2f} s"


def test_the_all_in_one_plasmid_is_still_far_enough_over_the_budget_to_record_it():
    """The gap S5.4 still has, asserted rather than described.

    A dense all-in-one plasmid that has lost one codon per element takes about 2.3 s on this
    machine, against 0.7 s for the same plasmid with its elements exact, and 4.5 s if it has
    lost three. Nothing about it is exotic: 10.4 kb of ordinary elements, one nuclease, and a
    codon variant in each — annotated identically to the verbatim plasmid, as the test above
    holds.

    Closing the gap means changing how an element's alignment is chosen, which can change
    which variant is annotated and with what coordinates, so it is not a constant to tune.
    The floor below is the point of this test: it fails on the day the engine gets fast
    enough, and whoever sees it fail should record that S5.4's budget is met and delete it.
    """
    _, cells = _fixture_work(ALL_IN_ONE_ELEMENTS, 1, ())
    assert cells >= 60_000_000, (
        f"{cells / 1e6:.2f} M cells: the gap has closed. Re-time the fixture, record that "
        f"S5.4's budget is met, and delete this test."
    )


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


def test_a_tandem_array_reports_a_copy_for_every_repeat():
    """The same thing as the array grows: however many copies are there, each is found."""
    tag = entry("3×HA")
    for copies in (2, 4, 8):
        sequence = clean_dna(414, 19) + tag.sequence * copies + clean_dna(414, 20)
        starts = {m.start for m in detect(sequence, circular=False) if m.name == tag.name}
        assert {414 + i * tag.length for i in range(copies)} <= starts, copies
