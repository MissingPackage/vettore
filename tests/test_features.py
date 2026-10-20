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


def test_the_database_file_is_the_one_the_builder_produces():
    """The builder needs the network, so this only checks the file's own consistency."""
    text = (Path(__file__).resolve().parents[1] / "src" / "vettore" / "data"
            / "common_features.json").read_text()
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
