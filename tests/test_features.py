import json
from pathlib import Path

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
