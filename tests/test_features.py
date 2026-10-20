import json
from pathlib import Path

import pytest
from Bio import SeqIO
from cloning_helpers import clean_dna, rc

from vettore.features import SHORT_FEATURE, catalogue, detect, detect_payload

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
    # Read on the other strand, the match carries EGFP's own first residues.
    assert "GGTGAGCAAGGGCGAGGAGCTGTTCACC" in rc(sequence[egfp.start : egfp.end + 1])
    assert egfp.end - egfp.start + 1 == 720          # the whole CDS, not the aligned core only
    assert 98.0 <= egfp.identity < 100.0             # two bases of the flank differ


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
    assert rows["U6 promoter"]["isNew"] is True
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
    assert {f["name"] for f in body["features"]} >= {"SaCas9", "U6 promoter", "AmpR"}
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
    ("EF1A", "EF-1α promoter", True), ("EF-1a promoter", "EF-1α promoter", True),
    ("hEF1a", "EF-1α promoter", True), ("SV40 late pA", "SV40 polyA", True),
    ("BGH pA", "bGH polyA", True), ("5' ITR", "AAV2 ITR", True),
    ("hSyn1 promoter", "hSyn promoter", True), ("CMV-promoter", "CMV promoter", True),
    ("HA tag", "3×HA", True), ("EGFP:T2A:Puro", "EGFP", True),
    ("CMV promoter", "SV40 polyA", False), ("pUC ori", "f1 ori", False),
    ("Exon 1a", "EF-1α promoter", False), ("PGK promoter", "CMV promoter", False),
    ("promoter", "CMV promoter", False), ("my insert", "CMV enhancer", False),
])
def test_names_are_matched_on_the_words_that_identify_an_element(label, element, expected):
    from vettore.features import _names_element

    assert _names_element(label, element) is expected


def test_name_words_ignore_generic_words_and_numbers():
    from vettore.features import name_words

    assert name_words("SV40 late polyA signal") == {"sv40"}
    assert name_words("EF-1α promoter", pieces=False) == {"ef1a"}
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
# Every entry is a subsequence of the record its ``source`` names. The fixtures below are
# other records — another laboratory, another vector — that carry the same element, so
# finding it there tests the sequence rather than the record it was taken from. The
# provenance of each fixture is in tests/data/README.md.

CROSS_CHECKS = {
    "human PGK promoter": "MT921016.1",
    "mouse PGK promoter": "MT891341.1",
    "UbC promoter": "MN389527.1",
    "H1 promoter": "FJ687158.1",
    "7SK promoter": "MQ296444.1",
    "SFFV promoter": "MP040844.1",
    "CAG promoter": "HG530138.1",
}


def entry(name: str):
    [feature] = [f for f in catalogue() if f.name == name]
    return feature


def cross_check_sequence(accession: str) -> str:
    record = SeqIO.read(DATA / f"features-{accession}.fasta", "fasta")
    return str(record.seq).upper()


@pytest.mark.parametrize(("element", "accession"), sorted(CROSS_CHECKS.items()))
def test_a_catalogue_promoter_is_found_in_an_independent_record(element, accession):
    feature = entry(element)
    assert not feature.source.startswith(accession), element
    sequence = cross_check_sequence(accession)
    matches = [m for m in detect(sequence, False) if m.name == element]
    assert len(matches) == 1, element
    assert matches[0].identity >= 98.0, element
    found = matches[0].end - matches[0].start + 1
    assert found >= 0.9 * feature.length, element


def test_the_promoters_come_back_whole_from_their_own_sequences():
    """Each entry is its record's subsequence at the cited coordinates, so a construct that
    carries the catalogue sequence carries the record's bases: the match must be exact."""
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
    pgk = entry("human PGK promoter")
    sequence = pgk.sequence + clean_dna(1500, 44)
    mislabelled = [{"name": "CMV promoter", "start": 0, "end": pgk.length - 1}]
    payload = detect_payload(sequence, True, existing=mislabelled)
    [row] = payload["features"]
    assert row["name"] == "human PGK promoter"
    assert row["labelWarning"] == (
        "Annotated as “CMV promoter”, but the sequence is the human PGK promoter "
        "(100 % identity). Check the label.")
    assert row["isNew"] is False and row["annotatedAs"] == "CMV promoter"
    assert payload["label_warnings"] == 1


@pytest.mark.parametrize(("label", "element", "expected"), [
    ("PGK promoter", "human PGK promoter", True),
    ("PGK promoter", "mouse PGK promoter", True),
    ("hPGK promoter", "human PGK promoter", True),
    ("UBC", "UbC promoter", True), ("UbC-EGFP", "UbC promoter", True),
    ("H1", "H1 promoter", True), ("H1-shRNA", "H1 promoter", True),
    ("7SK promoter", "7SK promoter", True), ("SFFV", "SFFV promoter", True),
    ("CAG", "CAG promoter", True), ("CAG", "CAG promoter (truncated)", True),
    ("U6 promoter", "H1 promoter", False), ("PGK promoter", "UbC promoter", False),
    ("SFFV promoter", "CMV promoter", False), ("CMV promoter", "human PGK promoter", False),
])
def test_the_new_promoters_are_named_by_the_labels_vectors_use(label, element, expected):
    from vettore.features import _names_element

    assert _names_element(label, element) is expected
